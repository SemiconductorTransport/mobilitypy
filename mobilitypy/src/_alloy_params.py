from .database import database
import numpy as np

## ==============================================================================
class _AlloyParams:
    '''
    The functions in this class calculates the parameters for alloy from their
    binary components.
    '''
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN',
                 alloy_type:str='wz'):
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
        alloy_type :  str, optional 
            The crystal type of the materials. This will be considered when calculating
            parameters like Poisson ratio etc.
            Use following abbreviation name:
                for wurtzite use 'WZ' or 'wz'.
                for zincblende use 'ZB' or 'zb'.
                for diamond use 'DM' or 'dm'.
            The default is 'wz'. 

        Returns
        -------
        None.

        """
        self.comps_ = compositions
        self.bins_ = binaries
        self.alloy_ = alloy
        self.alloy_type_ = alloy_type
 
    def _get_ternary_params(self, use_bin_params:dict=None):
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
        bin_1_params_db = database.get(self.bins_[0]).copy()
        bin_2_params_db = database.get(self.bins_[1]).copy()
        alloy_params_db = database.get(self.alloy_).copy()
        
        #####_update_database_data_locally
        if use_bin_params is not None:
            for bin_name, parms_ in use_bin_params.items():
                for pms_n in parms_:
                    if bin_name == self.bins_[0]:
                        bin_1_params_db[pms_n] = parms_[pms_n]
                    elif bin_name == self.bins_[1]:
                        bin_2_params_db[pms_n] = parms_[pms_n]
                    elif bin_name == self.alloy_:
                        alloy_params_db[pms_n] = parms_[pms_n]
                    else:
                        raise ValueError('Material name should be one of binary or alloy name')
        ######
        self.alloy_params_ = {}
        for key, bowing in alloy_params_db.items():
            self.alloy_params_[key] = self.comps_ * bin_1_params_db.get(key) +\
            (1-self.comps_) * bin_2_params_db.get(key) - bowing*self.comps_*(1-self.comps_)  
        self._cal_square_electromechanical_coupling_const()
        self._get_strain_realted_properties()
        #print(self.alloy_params_)
            
    def _get_alloy_params(self, system='ternary', use_bin_params:dict=None):
        """
        This function calculates the parameters for a ternary alloy from its
        binary component parameters using quadratic interpolation.

        Parameters
        ----------
        system : string (case sensitive), optional
            Type of the alloy. E.g. 'ternary'. 
            The default is 'ternary'.
        use_bin_params : dict, optional
            To use different materials parameters from that given in the database.
            Units should be same as in the database.
            e.g. use_bin_params = {'AlN': {'mass_density': 3000}}

        Returns
        -------
        Parameters for alloy.

        """
        if self.comps_ is None:
            self.comps_ = np.linspace(0., 1.0, 101)
        elif isinstance(self.comps_, float) or isinstance(self.comps_, int):
            self.comps_ = np.array([self.comps_])
        if system == 'ternary':
            self._get_ternary_params(use_bin_params=use_bin_params)
            
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
        return database.get(substrate_name).copy()
    
    def _get_strain_realted_properties(self):
        """
        This function calculated material properties related to strain-stress relation.
        
        Note-1: For Biaxial strain; epsilon_zz = biaxial_distortion_coefficient * epsilon_xx
        Note-2: For WZ structure the Poisson ratio is anisotropic. However, the strain
        field due to edge dislocation uses 'isotropic/average Poisson ratio'. We use
        the definition isotropic_Poisson_ratio = C_12/(C_11+C_12) following reference
        Shi et al., APL 74, 573 (1999). This was what also used in DJ's thesis.

        For WZ: 
            biaxial_distortion_coefficient = -2*C_13/C_33 (including 'negative' sign.)
            isotropic_Poisson_ratio = C_12/(C_11+C_12)

        Raises
        ------
        ValueError
            If alloy type is not implemented yet.

        Returns
        -------
        numpy array
            Distortion coefficient for biaxial strain. 

        """
        if self.alloy_type_.lower() == 'wz':
            # epsilon_zz = -2*C_13/C_33 * epsilon_xx
            self.alloy_params_['biaxial_distortion_coefficient'] =\
                -2*(self.alloy_params_.get('C_13')/self.alloy_params_.get('C_33'))
            # isotropic Poisson ratio = C_12/(C_11 + C_12)
            self.alloy_params_['isotropic_Poisson_ratio'] =\
                self.alloy_params_.get('C_12')/(self.alloy_params_.get('C_11')+self.alloy_params_.get('C_12'))
        else:
            raise ValueError(f'{self.alloy_type_} is not implemented yet. Contact developer.')
            
    def _cal_omega_0_ad(self, a_lp, b_lp, c_lp):
        """
        NOTE: As described in the DJ's book Omega_0 in alloy disorder scattering
        mobility is not 'literally' the unit cell volume, 
        rather it is the effective volume that atom A (or B) occupies
        in the crystal or the average volume per atom.  

        """
        if self.alloy_type_.lower() == 'wz':
            # the average volume per atom in the wurtzite structure is given by
            # np.sqrt(3)/8 * a^2 c. 4 atoms per unit cell. Unit cell vol = sqrt(3)/2 * a^2 c
            return 0.21650635094610965 * a_lp * b_lp * c_lp # sqrt(3)/8 * a^2 c
        elif self.alloy_type_.lower() == 'zb':
            # Here the primitive cell volume is a^3/4. But we have 2 atoms per
            # unit cell. 
            return 0.125 * a_lp * b_lp * c_lp # a^3/8
        else:
            raise ValueError(f'{self.alloy_type_} is not implemented yet. Contact developer.')
                 
    def _cal_square_electromechanical_coupling_const(self):
        e_33, e_31, e_15 = self.alloy_params_['e_33'], self.alloy_params_['e_31'], self.alloy_params_['e_15']
        u_LA, u_TA = self.alloy_params_['LA_phonon_velocity'], self.alloy_params_['TA_phonon_velocity']
        eps_s, rho = self.alloy_params_['static_dielectric_constant'], self.alloy_params_['mass_density']
        eps_00 = 8.8541878188e-12 #  C.V^-1.m^-1
        
        e_31_p_2_e_15 = e_31 + 2.0*e_15
        e_33_e_31_e_15 = e_33 - e_31 - e_15
        
        e_l_2_av = (15.0*e_33*e_33 + 12.0*e_33*e_31_p_2_e_15 + 8.0*e_31_p_2_e_15*e_31_p_2_e_15) #/105
        e_T_2_av = (48.0*e_15*e_15 + 16.0*e_33_e_31_e_15*e_15 + 6.0*e_33_e_31_e_15*e_33_e_31_e_15) #/105
        
        self.alloy_params_['electromechanical_coupling_const'] = \
            ((e_l_2_av/(u_LA*u_LA*eps_00)) + (e_T_2_av/(u_TA*u_TA*eps_00)))/(eps_s*rho*105.0)
        return 