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
            P_SixGe1-x = x*P_Si + (1-x)*P_Ge + x*(1-x)*P_bowing 
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