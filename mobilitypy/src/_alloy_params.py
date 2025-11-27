from .database import database
import numpy as np

## ==============================================================================
class _AlloyParams:
    '''
    The functions in this class calculates the parameters for alloy from their
    binary components.
    '''
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN'):
        """
        Initiation function of the class _AlloyParams.

        Parameters
        ----------
        compositions : 1D array of float, optional
            The alloy mole fractions. E.g. x values in Si_xGe_1-x. The default is None.
            If None, a composition array is generated using `np.linspace(start=0.01, end=0.99, num=101)`.
        binaries : list of strings (case sensitive), optional
            Name of the corresponding binaries of requested alloy. They should
            match the names in database. All implemented materials name list 
            can be found in the README. For ternary alloy 'compositions' correspond 
            to the 1st binary in the list; for quaternaries 1st binary is 1st composition
            and so on (from left to right). The default is ['AlN', 'GaN'].
        alloy : string (case sensitive), optional
            The alloy name. The name should match the name in database. All   
            implemented materials name list can be found in the README. Case sensitive.
            The default is 'AlGaN'.

        Returns
        -------
        None.

        """
        self.comps_ = compositions
        self.bins_ = binaries
        self.alloy_ = alloy

    def _get_ternary_params(self):
        """
        This function calculates the parameters for a ternary alloy from its
        binary component parameters using quadratic interpolation.
        E.g. for any parameter, P:
            P_SixGe1-x = x*P_Si + (1-x)*P_Ge - x*(1-x)*P_bowing 
            P_bowing is the quadratic bowing parameter for the parameter P.
        Returns
        -------
        Parameters for alloy.

        """
        assert len(self.bins_) == 2, 'Provide two binary compounds'
        bin_1_params_db = database.get(self.bins_[0])
        bin_2_params_db = database.get(self.bins_[1])
        alloy_params_db = database.get(self.alloy_)
        self.alloy_params_ = {}
        for key, bowing in alloy_params_db.items():
            self.alloy_params_[key] = self.comps_ * bin_1_params_db.get(key) +\
            (1-self.comps_) * bin_2_params_db.get(key) - bowing*self.comps_*(1-self.comps_)
        #print(self.alloy_params_)
            
    def _get_alloy_params(self, system='ternary'):
        """
        This function calculates the parameters for a ternary alloy from its
        binary component parameters using quadratic interpolation.

        Parameters
        ----------
        system : string (case sensitive), optional
            Type of the alloy. E.g. 'ternary'. 
            The default is 'ternary'.

        Returns
        -------
        Parameters for alloy.

        """
        if self.comps_ is None:
            self.comps_ = np.linspace(0.01, 0.99, 101)
        elif isinstance(self.comps_, float) or isinstance(self.comps_, int):
            self.comps_ = np.array([self.comps_])
        if system == 'ternary':
            self._get_ternary_params()
            
    @staticmethod        
    def _get_substrate_properties(substrate_name):
        """
        Generate the substrate properties for phsedomorphic strain.

        Parameters
        ----------
        substrate_name : str
            The name of the substrate. The name should be in the database. If 
            the name does not exists in the database return None.

        Returns
        -------
        Dictionary
            The parameters of the substrate. Get from database. If substrate
            name does not exists in the database return None.

        """
        return database.get(substrate_name)
    
    def _get_Poisson_ratio(self):
        """
        Poisson ratio.
        
        epsilon_yy = Poisson_ratio * epsilon_xx
        
        Poisson_ratio here includes 'negative' sign.
        For WZ: Poisson_ratio = -2*C_13/C_33

        Parameters
        ----------
        alloy_params_ : dictionary
            The parameters dictionary for alloy.
        alloy_type :  str, optional (case insensitive)
            The crystal type of alloy. This will be considered when calculating
            parameters like Poisson ratio etc.
            Use following abbreviation name:
                for wurtzite use 'WZ' or 'wz'.
                for zincblende use 'ZB' or 'zb'.
                for diamond use 'DM' or 'dm'.
            The default is 'WZ'. 

        Raises
        ------
        ValueError
            If alloy type is not implemented yet.

        Returns
        -------
        numpy array
            Poisson ratio. 

        """
        if self.alloy_type_.lower() == 'wz':
            # epsilon_zz = -2*C_13/C_33 * epsilon_xx
            return -2*(self.alloy_params_.get('C_13')/self.alloy_params_.get('C_33'))
        else:
            raise ValueError(f'{self.alloy_type_} is not implemented yet. Contact developer.')