from .database import material_database
import numpy as np
#import warnings

## ==============================================================================
class _AlloyParams:
    '''
    The functions in this class calculates the parameters for alloy from their
    binary components.
    '''
    _alloy_name_map = {'AlNGaN': 'AlGaN','GaNAlN': 'AlGaN', 'AlGaN':'AlGaN'}
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], 
                 alloy_crystal_structure:str='wz', alloy_type:str=None):
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
            can be found in the README. Alloy 'compositions' will be calculated  
            considering compositions of binaries correspond to x, y, 1-x-y etc. 
            from left to right in the list. For alloy_type='CatAni', it will be counted as
            components it is counted as x, 1-x, y, 1-y etc. 
            The default is ['AlN', 'GaN'].
        alloy_crystal_structure :  str, optional [options: 'WZ', 'ZB', 'DM']
            The crystal type of the materials. This will be considered when calculating
            parameters like Poisson ratio etc.
            Use following abbreviation name:
                for wurtzite use 'WZ' or 'wz'.
                for zincblende use 'ZB' or 'zb'.
                for diamond use 'DM' or 'dm'.
            The default is 'wz'. 
        alloy_type : string (case sensitive), optional [options: 'CatAni']
            The alloy type name. Case sensitive. Only needed if alloy is of AxB1-xCxD1-y
            kind. Will be ignored for alloy of type AxB1-x, AxByC1-x-y, AxByCzD1-x-y-z etc.
            The default is None.

        Returns
        -------
        None.

        """
        self.comps_ = compositions
        self.bins_ = list(binaries)
        self.alloy_type_ = alloy_type
        self.alloy_crys_type_ = alloy_crystal_structure.lower()
 
    def _get_two_comp_component_alloy_params(self, params_db):
        """
        This function calculates the parameters for a 2-composition component alloy 
        from itsbinary component parameters using quadratic interpolation.
        E.g. for any parameter, P:
            P_SixGe1-x = x*P_Si + (1-x)*P_Ge - x*(1-x)*P_bowing 
            P_bowing is the quadratic bowing parameter for the parameter P.
        Returns
        -------
        Parameters for alloy.

        """        
        bin_1_params_db = params_db.get(self.bins_[0])
        bin_2_params_db = params_db.get(self.bins_[1])
        alloy_params_db = params_db.get(self.alloy_name)

        self.alloy_params_ = {}
        for key, bowing in alloy_params_db.items():
            self.alloy_params_[key] = self.comps_ * bin_1_params_db.get(key) +\
                                        (1-self.comps_) * bin_2_params_db.get(key)\
                                            - bowing*self.comps_*(1-self.comps_)  
        self._cal_square_electromechanical_coupling_const()
        self._get_strain_realted_properties()
        return 
    
    @staticmethod
    def _update_material_params_locally(use_mat_params, params_db):
        """
        This function allows to update the materials parameters locally before calculation.
        It will not change the original database that comes with the package.

        Parameters
        ----------
        use_mat_params : dict
            To use different materials parameters from that given in the database.
            Simply join the binary names to construct the alloy name. e.g.,
            for binaries=['AlN', 'GaN'] the alloy name is 'AlNGaN' or 'GaNAlN'.
            Material parameters units should be same as in the database.
            e.g. use_mat_params = {'AlN': {'mass_density': 3000}}
        params_db : dict
            The original material parameters as read from database.

        Returns
        -------
        params_db_tmp : dict
            The updated material parameters after update.

        """
        params_db_tmp = params_db.copy()
        # Avoids appending unnecessary parameters in the return database
        allowed_params = material_database['GaN'].keys() # We know one 'default' material in the database
        for mat_name, parms_ in use_mat_params.items():
            if mat_name in params_db_tmp:
                for pms_n, pms_vals in parms_.items():
                    if pms_n in allowed_params:
                        params_db_tmp[mat_name][pms_n] = pms_vals
                    else:
                        print(f"Warning: Ignoring update for {mat_name}:{pms_n}. Does not exist in database.")
            else:
                print(f"Warning: Ignoring update for {mat_name}. Does not exist in database.")
                # warnings.warn(f"Ignoring update for {mat_name} material, as it does not exits in the database.", 
                #               RuntimeWarning, stacklevel=1) 
        return params_db_tmp
        
    def _get_params_from_database(self, use_mat_params):
        """
        This function reads requested materials parameters from database and
        updates them if needed.

        Parameters
        ----------
        use_mat_params : dict, optional
            To use different materials parameters from that given in the database.
            Simply join the binary names to construct the alloy name. e.g.,
            for binaries=['AlN', 'GaN'] the alloy name is 'AlNGaN' or 'GaNAlN'.
            Material parameters units should be same as in the database.
            e.g. use_mat_params = {'AlN': {'mass_density': 3000}}

        Returns
        -------
        dict
            Material parameters.

        """
        alloy_name = ''.join(self.bins_)
        self.alloy_name = _AlloyParams._alloy_name_map.get(alloy_name)
        if self.alloy_name is None: 
            raise ValueError('Requested alloyes are not implemented yet. Contact developer.')
            
        params_db = {self.alloy_name: material_database.get(self.alloy_name).copy()} # These are bowing parameters
        for mat_name in self.bins_:
            params_db[mat_name] = material_database.get(mat_name).copy()

        if use_mat_params is None:
            return params_db  
        else:
            #####_update_database_data_locally
            if alloy_name in use_mat_params:
                use_mat_params[self.alloy_name] = use_mat_params.pop(alloy_name)
            return self._update_material_params_locally(use_mat_params, params_db)
            
    def _get_alloy_params(self, use_mat_params:dict=None):
        """
        This function calculates material parameters for alloy from its
        binary component parameters using interpolation.

        Parameters
        ----------
        use_mat_params : dict, optional
            To use different materials parameters from that given in the database.
            Simply join the binary names to construct the alloy name. e.g.,
            for binaries=['AlN', 'GaN'] the alloy name is 'AlNGaN' or 'GaNAlN'.
            Material parameters units should be same as in the database.
            e.g. use_mat_params = {'AlN': {'mass_density': 3000}}

        Returns
        -------
        Parameters for alloy.

        """
        # Define compositions
        if self.comps_ is None:
            self.comps_ = np.linspace(0., 1.0, 101)
        elif isinstance(self.comps_, float) or isinstance(self.comps_, int):
            self.comps_ = np.array([self.comps_])
        # Get parameters from database and update it if use_mat_params is not none.
        params_dbs = self._get_params_from_database(use_mat_params)

        if len(self.bins_) == 2:
            self._get_two_comp_component_alloy_params(params_dbs)
        else:
            raise ValueError('Requested multiple composition component alloys are not implemented yet. Contact developer.')
            
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
        return material_database.get(substrate_name).copy()
    
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
        if self.alloy_crys_type_ == 'wz':
            # epsilon_zz = -2*C_13/C_33 * epsilon_xx
            self.alloy_params_['biaxial_distortion_coefficient'] =\
                -2*(self.alloy_params_.get('C_13')/self.alloy_params_.get('C_33'))
            # isotropic Poisson ratio = C_12/(C_11 + C_12)
            self.alloy_params_['isotropic_Poisson_ratio'] =\
                self.alloy_params_.get('C_12')/(self.alloy_params_.get('C_11')+self.alloy_params_.get('C_12'))
        else:
            raise ValueError(f'{self.alloy_crys_type_} is not implemented yet. Contact developer.')
            
    def _cal_pseudomorphic_strain(self, substrate:str|float):
        if isinstance(substrate, str):
            substrate_params_dic = self._get_substrate_properties(substrate)
            substrate_lp = substrate_params_dic.get('lattice_a0') # substrate in-plane lattice parameter
        else:    
            substrate_lp = float(substrate)
            
        lattice_a = self.alloy_params_.get('lattice_a0') 
        lattice_c = self.alloy_params_.get('lattice_c0') 
        epsilon_zz = self.alloy_params_.get('biaxial_distortion_coefficient')\
                              *(substrate_lp/lattice_a - 1.0)
        # Re-populate the lattice parameters
        self.alloy_params_['lattice_a0']  = np.array([substrate_lp]*len(lattice_a)) 
        self.alloy_params_['lattice_c0']= lattice_c * (1.0 + epsilon_zz)
        return
            
    def _cal_omega_0_ad(self, a_lp, b_lp, c_lp):
        """
        NOTE: As described in the DJ's book Omega_0 in alloy disorder scattering
        mobility is not 'literally' the unit cell volume, 
        rather it is the effective volume that atom A (or B) occupies
        in the crystal or the average volume per atom.  

        """
        if self.alloy_crys_type_ == 'wz':
            # the average volume per atom in the wurtzite structure is given by
            # np.sqrt(3)/8 * a^2 c. 4 atoms per unit cell. Unit cell vol = sqrt(3)/2 * a^2 c
            return 0.21650635094610965 * a_lp * b_lp * c_lp # sqrt(3)/8 * a^2 c
        elif self.alloy_crys_type_ == 'zb':
            # Here the primitive cell volume is a^3/4. But we have 2 atoms per
            # unit cell. 
            return 0.125 * a_lp * b_lp * c_lp # a^3/8
        else:
            raise ValueError(f'{self.alloy_crys_type_} is not implemented yet. Contact developer.')
                 
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