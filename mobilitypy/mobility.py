from .src import _DataBase, _AlloyParams, _FermiDiracInt, _MobilityCarrier
from .src import _Mobility2DCarrier, _Mobility3DCarrier
from .utilities import _plot_mobilities
import numpy as np

## ==============================================================================
class DataBase(_DataBase):
    """
    The functions in this class use to print/manipulate the material parameters
    in the database.
    """
    def __init__(self):
        self.dtbase = _DataBase()
        
    def print_database(self, for_material:str=None):
        """
        This function prints the information of material parameters from the database.

        Parameters
        ----------
        for_material : string (case sensitive), optional
            The material name for which the parameters will be printed. 
            The name should match the name in database. If None, prints general
            information about the database.
            The default is None.

        Returns
        -------
        None.

        """
        self.dtbase._print_database(for_material=for_material)
        return
        
    def update_database(self, for_material=None, with_new_database=None):
        """
        To update the material parameters in the database. 
        Not implemented yet. Contact developer.

        Parameters
        ----------
        for_material : string (case sensitive), optional
            The material name for which the parameters will be printed. 
            The name should match the name in database. If None, prints general
            information about the database.
            The default is None.
        with_new_database : dictionary, optional
            The material parameters. The parameter names should match in the database.
            Use 'print_database()' function from DataBase class for material parameter
            names in the database.
            The default is None.

        Returns
        -------
        None.

        """
        self.dtbase._update_database(for_material=for_material, 
                                     with_new_database=with_new_database)
        return
        
class AlloyParams(_AlloyParams):
    '''
    The functions in this class calculates the parameters for alloy from their
    binary components.
    '''
    def __init__(self):
        pass
            
    def get_alloy_params(self, system='ternary', compositions=None, binaries=['AlN', 'GaN'], 
                         alloy='AlGaN', alloy_type:str='wz'):
        """
        This function calculates the parameters for a ternary alloy from its
        binary component parameters using quadratic interpolation.
        E.g. for any parameter, P:
            P_SixGe1-x = x*P_Si + (1-x)*P_Ge - x*(1-x)*P_bowing 
            P_bowing is the quadratic bowing parameter for the parameter P.

        Parameters
        ----------
        system : string (case sensitive), optional
            Type of the alloy. E.g. 'ternary'. 
            The default is 'ternary'.
        compositions : 1D array of float, optional
            The alloy mole fractions. E.g. x values in Si_xGe_1-x. The default is None.
            If None, a composition array is generated using `np.linspace(start=0.01, end=0.99, num=101)`.
        binaries : list of strings (case sensitive), optional
            Name of the corresponding binaries of requested alloy. They should
            match the names in database. All implemented materials name list 
            can be found in the README. 
            The default is ['AlN', 'GaN'].
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
        1D float array
            Parameters for alloy.

        """
        _AlloyParams.__init__(self, compositions=compositions, binaries=binaries, 
                              alloy=alloy, alloy_type=alloy_type)
        return self._get_alloy_params(system=system)

class Mobility2DCarrier(_MobilityCarrier, _Mobility2DCarrier):
    """
    The functions in this class calculates the mobility of 2D carrier gas.  
    The mobility models are implemented based on the following references.
    
    Ref-1: J. Bassaler, J. Mehta, I. Abid, L. Konczewicz, S. Juillaguet, S. Contreras, S. Rennesson, 
    S. Tamariz, M. Nemoz, F. Semond, J. Pernot, F. Medjdoub, Y. Cordier, P. Ferrandis, 
    Al-Rich AlGaN Channel High Electron Mobility Transistors on Silicon: A Relevant Approach for High 
    Temperature Stability of Electron Mobility. Adv. Electron. Mater. 2024, 2400069. 
    https://doi.org/10.1002/aelm.202400069

    Ref-2: Zhang, J., Hao, Y., Zhang, J. et al. The mobility of two-dimensional electron gas in AlGaN/GaN 
    heterostructures with varied Al content. Sci. China Ser. F-Inf. Sci. 51, 780–789 (2008). 
    https://doi.org/10.1007/s11432-008-0056-7
    
    Ref-3: Mondal et. al., Interplay of carrier density and mobility in Al-rich (Al,Ga)N-channel HEMTs: 
    Impact on high-power device performance potential. APL Electronic Devices 1, 026117 (2025)
    https://doi.org/10.1063/5.0277051
    """
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN', 
                 system='ternary', pseudomorphic_strain:bool=False, substrate=None,
                 alloy_type='WZ', eps_n_2d=1e-8, print_log=None):
        """
        Initiation function of the class Mobility2DCarrier.
        
        Parameters
        ----------
        compositions : 1D array of float, optional
            The alloy mole fractions. E.g. x values in Si_xGe_1-x. The default is None.
            If None, a composition array is generated using `np.linspace(start=0.01, end=0.99, num=101)`.
        binaries : list of strings (case sensitive), optional
            Name of the corresponding binaries of requested alloy. They should
            match the names in database. All implemented materials name list 
            can be found in the README. 
            The default is ['AlN', 'GaN'].
        alloy : string (case sensitive), optional
            The alloy name. The name should match the name in database. All   
            implemented materials name list can be found in the README. Case sensitive.
            The default is 'AlGaN'.
        system : string (case sensitive), optional
            Type of the alloy. E.g. 'ternary'. 
            The default is 'ternary'.
        pseudomorphic_strain : bool, optional
            Whether to consider pseudomorphic strain.
            The default is False.
        substrate : string or float, optional (unit: Angstrom)
            The substrate name (if string, warning: the name should be in the database) 
            or the substrate in-plane lattice parameter (if float, Angstrom unit).
            The default is None. Error will be raised if substrate=None and 
            pseudomorphic_strain=True.
        alloy_type :  str, optional 
            The crystal type of alloy. This will be considered when calculating
            parameters like Poisson ratio etc.
            Use following abbreviation name:
                for wurtzite use 'WZ' or 'wz'.
                for zincblende use 'ZB' or 'zb'.
                for diamond use 'DM' or 'dm'.
            The default is 'WZ'. 
        eps_n_2d : float, optional (unit: 10^12 cm^-2)
            Carrier density below eps_n_2d will be considered as zero. 
            The default is 1e-8 == 1e4 cm^-2.
        print_log : string, optional => ['high','medium','low', None]
            Determines the level of log to be printed. The default is None.

        Returns
        -------
        None.

        """
        _MobilityCarrier.__init__(self, compositions=compositions, binaries=binaries, 
                                  alloy=alloy, system=system, pseudomorphic_strain=pseudomorphic_strain, 
                                  substrate=substrate,alloy_type=alloy_type,
                                  print_log=print_log, eps_n=eps_n_2d)
        _Mobility2DCarrier.__init__(self)
        
    def calculate_sheet_mobility(self, n_2d=10, rms_roughness=0.1, corr_len=1,  
                                 n_dis=1, f_dis=0.1, T=300, 
                                 alloy_disordered_effect:bool=False,
                                 interface_roughness_effect:bool=False,
                                 dislocation_effect:bool=False,
                                 deformation_potential_effect:bool=False, 
                                 piezoelectric_effect:bool=False,
                                 acoustic_phonon_effect:bool=False,
                                 polar_optical_phonon_effect:bool=False,
                                 total_mobility:bool=True,
                                 calculate_total_mobility_only:bool=False,
                                 return_sc_rates:bool=False,
                                 mobility_model='v2'):
        """
        This function calculates the sheet mobility from different scattering contributions.
        The mobility models are implemented based on the following references.
        
        Ref-1: J. Bassaler, J. Mehta, I. Abid, L. Konczewicz, S. Juillaguet, S. Contreras, S. Rennesson, 
        S. Tamariz, M. Nemoz, F. Semond, J. Pernot, F. Medjdoub, Y. Cordier, P. Ferrandis, 
        Al-Rich AlGaN Channel High Electron Mobility Transistors on Silicon: A Relevant Approach for High 
        Temperature Stability of Electron Mobility. Adv. Electron. Mater. 2024, 2400069. 
        https://doi.org/10.1002/aelm.202400069

        Ref-2: Zhang, J., Hao, Y., Zhang, J. et al. The mobility of two-dimensional electron gas in AlGaN/GaN 
        heterostructures with varied Al content. Sci. China Ser. F-Inf. Sci. 51, 780–789 (2008). 
        https://doi.org/10.1007/s11432-008-0056-7
        
        Ref-3: Mondal et. al., Interplay of carrier density and mobility in Al-rich (Al,Ga)N-channel HEMTs: 
        Impact on high-power device performance potential. APL Electronic Devices 1, 026117 (2025)
        https://doi.org/10.1063/5.0277051
            
        The considered scattering mechanism are:
            Interface roughness mediated (IRF)
            Threading dislocation mediated (DIS)
            Alloy disorder limited (AD)
            Deformation potential mediated (DP)
            Piezoelectric effect (PE)
            Acoustic phonon (AP)
            Polar optical phonon (POP)

        Parameters
        ----------
        n_2d : 1D float array or float, optional (unit: 10^12 cm^-2)
            Array containing carrier density data for compositions. This can be
            a single number as well. Then all compositions will have same carrier
            density. The default is 10 == 1e13 cm^-2.
        rms_roughness : float, optional (unit: nm)
            Interface root-mean-squared roughness for interface-roughness scattering
            contribution. The default is 0.1.
        corr_len : float, optional (unit: nm)
            Correlation length of interface roughness. The default is 1.
        n_dis : float, optional (unit: 10^8 cm^-2)
            Threading dislocation density. The default is 1.
        f_dis : float, optional (unit: unitless)
            Fraction of dislocation that contributes in scattering. 
            The default is 0.1.
        T : float, optional (unit: K)
            Temperature at which mobility calculations will be done. 
            The default is 300K.
        alloy_disordered_effect : bool, optional
            Whether to calculate alloy disordered mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        interface_roughness_effect : bool, optional
            Whether to calculate interface roughness effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        dislocation_effect : bool, optional
            Whether to calculate interface roughness effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        deformation_potential_effect : bool, optional
            Whether to calculate deformation potential effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        piezoelectric_effect : bool, optional
            Whether to calculate piezoelectric effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        acoustic_phonon_effect : bool, optional
            Whether to calculate acoustic phonon effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
            acoustic_phonon_effect combines deformation_potential and piezoelectric_effect.
            NOTE: One can use deformation_potential_effect=True, piezoelectric_effect=True, and
            acoustic_phonon_effect=True, to return all three contrbutions separately in
            the output. In the implementation it has been made sure that this setup
            does NOT double count the  deformation_potential and piezoelectric_effect;
            So, NO WORRIES!
        polar_optical_phonon_effect : bool, optional
            Whether to calculate polar optical phonon effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        total_mobility : bool, optional
           Whether to calculate total mobility. The default is True.
        calculate_total_mobility_only : 
            Calculate only the total mobility. If False the return data also contains individual 
            specified contributions.
        return_sc_rates : float, optional 
            Return the scattering rates values.The default is False.
        mobility_model : str, optional [options: 'v1', 'v2']
            Which mobility model to use. The default is 'v2'.
            The mobility is implemented based on following publications:
            'v1':
                J. Bassaler, J. Mehta, I. Abid, L. Konczewicz, S. Juillaguet, S. Contreras, S. Rennesson, 
                S. Tamariz, M. Nemoz, F. Semond, J. Pernot, F. Medjdoub, Y. Cordier, P. Ferrandis, 
                Al-Rich AlGaN Channel High Electron Mobility Transistors on Silicon: A Relevant Approach for High 
                Temperature Stability of Electron Mobility. Adv. Electron. Mater. 2024, 2400069. 
                https://doi.org/10.1002/aelm.202400069
                NB 1: Here, the dislocation scattering includes scattering from threading edge dislocation
                charge line only.
            'v2':
                Here, the dislocation scattering includes scattering from threading edge dislocation
                charge line plus scattering from strain field from threading edge dislocations.

        Returns
        -------
        pandas dataframe with compositions and mobility (unit: cm^2 V^-1 S^-1) columns.
            Total (or individual contributions) sheet mobility. If return_sc_rates=True,
            then scattering rates (10^12 s^-1) and m_star_by_e (10^-12 V.m^-2.s^2) are also returned.

        """

        self.alloy_disordered_effect_=alloy_disordered_effect
        self.interface_roughness_effect_=interface_roughness_effect
        self.dislocation_effect_=dislocation_effect
        self.deformation_potential_effect_=deformation_potential_effect
        self.piezoelectric_effect_=piezoelectric_effect
        self.acoustic_phonon_effect_=acoustic_phonon_effect
        self.polar_optical_phonon_effect_=polar_optical_phonon_effect
        self.only_total_mobility = calculate_total_mobility_only
        self.total_mobility_=total_mobility
        self.mobility_model_=mobility_model
        return self._calculate_sheet_mobility(n_2d=n_2d, rms_roughness=rms_roughness, 
                                              corr_len=corr_len, n_dis=n_dis, f_dis=f_dis, 
                                              T=T, return_sc_rates=return_sc_rates)
    @staticmethod
    def sc_rate_2_mobility(mstar_by_e, scattering_rate):
        # Scattering rate to mobility calculation 
        """
        This function calculates sheet mobility from scattering rate.
        
        Parameters
        ----------
        mstar_by_e : float/array (unit: 10^-12 V.m^-2.s^2)
            Carrier effective mass in m0 unit multiplied by m_star_by_e. 
            m_star_by_e = m* X m0 / e
            Note: The unit here is in meter. The conversion to cm is taken care 
            of inside the function operation.
        scattering_rate : float/array (unit: 10^12 s^-1)
            Scattering rate. 

        Returns
        -------
        float/array (unit: cm^2 V^-1 S^-1)
            Total (or individual contributions) sheet mobility. 

        """
        # Scattering rate to mobility calculation
        # 1e4 is unit conversion from m^2 to cm^2
        tau = mstar_by_e * scattering_rate
        tau[tau<1e-10] = np.nan
        return 1e4/tau # cm^2 V^-1 S^-1

    def calculate_sheet_resitance(self, n_2d, mobility):
        """
        This function calculates the sheet resistance.

        Parameters
        ----------
        n_2d : float/ndarray (unit: 10^12 cm^-2)
            Array containing carrier density data. 
        mobility : float/ndarray (unit: cm^2 V^-1 s^-1)
            Array containing mobility data.

        Returns
        -------
        float/ndarray (unit: ohm/square)
            Sheet resistance.

        """
        return self._calculate_sheet_resitance(n_2d, mobility)

    def calculate_figure_of_merit(self, n_2d, mobility, temp:float=300,
                                   mode:str='LFOM', T_corect_bandgap:bool=False, 
                                   direct_bandgap:bool=True, indirect_bandgap:bool=False):
        """
        This function calculates the figure-of-merit (FOM). Available FOMs are:
            LFOM: Lateral figure-of-merit (LFOM = e*n_2d*mu*E_cr^2 ). Here, 
                critical_electric_field is assumed related to bandgap following
                Ref: J. L. Hudgins, G. S. Simin, E. Santi and M. A. Khan, 
                "An assessment of wide bandgap semiconductors for power devices," 
                in IEEE Transactions on Power Electronics, vol. 18, no. 3, pp. 907-914, 
                May 2003, doi: 10.1109/TPEL.2003.810840.
        
                direct_bandgap_critical_electric_field = 1.73e5*(bandgap_**2.5) # V/cm
                indirect_bandgap_critical_electric_field = 2.38e5*(bandgap_**2) # V/cm

        Parameters
        ----------
        n_2d : 1D float array (unit: 10^12 cm^-2)
            Array containing carrier density data.
        mobility : 1D float array (unit: cm^2 V^-1 s^-1)
            Array containing mobility data.
        temp : float, optional (unit: K)
            Temperature for band gap correction. The default is 300K.
        mode : str, optional (['LFOM'])
            The figure-of-merit name. The default is 'LFOM'.
        T_corect_bandgap : bool, optional
            Apply temperature correction to bandgap or not. The default is False.
        direct_bandgap : bool, optional
            If the bandgap is direct bandgap or not. The default is True.
        indirect_bandgap : bool, optional
            If the bandgap is indirect bandgap or not. The default is False.

        Returns
        -------
        1D float array (unit: MW/cm^2)
            Figure-of-merit.

        """
        return self._calculate_figure_of_merit(n_2d, mobility, temp=temp, mode=mode,
                                               T_corect_bandgap=T_corect_bandgap,
                                               direct_bandgap=direct_bandgap, 
                                               indirect_bandgap=indirect_bandgap)
    
class Mobility3DCarrier(_MobilityCarrier, _Mobility3DCarrier):
    """
    The functions in this class calculates the mobility of 3D carrier gas.  
    The mobility models are implemented based on the following references.
    
    Note: Some of the equations in the references has prining mistakes. The mistakes
    are corrected in our implementation. 
    
    """
    
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN', 
                 system='ternary', pseudomorphic_strain:bool=True, substrate=None,
                 alloy_type='WZ', eps_n_3d=1e-14, print_log=None):
        """
        Initialization function of the class Mobility3DCarrier.
        
        Parameters
        ----------
        compositions : 1D array of float, optional
            The alloy mole fractions. E.g. x values in Si_xGe_1-x. The default is None.
            If None, a composition array is generated using `np.linspace(start=0.01, end=0.99, num=101)`.
        binaries : list of strings (case sensitive), optional
            Name of the corresponding binaries of requested alloy. They should
            match the names in database. All implemented materials name list 
            can be found in the README. 
            The default is ['AlN', 'GaN'].
        alloy : string (case sensitive), optional
            The alloy name. The name should match the name in database. All   
            implemented materials name list can be found in the README. Case sensitive.
            The default is 'AlGaN'.
        system : string (case sensitive), optional
            Type of the alloy. E.g. 'ternary'. 
            The default is 'ternary'.
        pseudomorphic_strain : bool, optional
            Whether to consider pseudomorphic strain.
            The default is True.
        substrate : string or float, optional (unit: Angstrom)
            The substrate name (if string, warning: the name should be in the database) 
            or the substrate in-plane lattice parameter (if float, Angstrom unit).
            The default is None. Error will be raised if substrate=None and 
            pseudomorphic_strain=True.
        alloy_type :  str, optional (case insensitive)
            The crystal type of alloy. This will be considered when calculating
            parameters like Poisson ratio etc.
            Use following abbreviation name:
                for wurtzite use 'WZ' or 'wz'.
                for zincblende use 'ZB' or 'zb'.
                for diamond use 'DM' or 'dm'.
            The default is 'WZ'. 
        eps_n_3d : float, optional (unit: 1e18 cm^-2)
            Carrier density below eps_n_3d will be considered as zero. 
            The default is 1e-14 1e18 cm^-2 == 1e4 cm^-2.
        print_log : string, optional => ['high','medium','low', None]
            Determines the level of log to be printed. The default is None.

        Returns
        -------
        None.

        """
        if (pseudomorphic_strain == True) and (substrate is None):
            raise ValueError('substrate tag can not be None when pseudomorphic_strain=True.')
        _MobilityCarrier.__init__(self, compositions=compositions, binaries=binaries, 
                                  alloy=alloy, system=system, pseudomorphic_strain=pseudomorphic_strain, 
                                  substrate=substrate,alloy_type=alloy_type,
                                  print_log=print_log, eps_n=eps_n_3d)
        _Mobility3DCarrier.__init__(self)
        
    def calculate_3DEC_props(self, n_d, T:float=300, inverse_half_FD_method:str='minimax_piecewise'):
        """
        This function calculates some general properties of 3D carrier. 

        Parameters
        ----------
        n_d : float or 1D float array (unit: 1e18 cm^-3)
            Volumetric carrier density data. 
            If array, array size should be same as composition arrary. 
        T : float, optional (unit: K)
            Temperature at which Fermi-Dirac integral calculations will be done. 
            The default is 300K.
        inverse_half_FD_method : str, optional [available: 'JD_approx', 'minimax_piecewise']
            The approximate method to calculate the scaled Fermi energy (E_f/k_BT) 
            using inverse Fermi-Dirac integral of order-1/2. The default is JD_approx.
            JD_approx : Joyce-Dixon approximation (APL 31, 354 (1977)).
            minimax_piecewise : minimax approximation (Applied Mathematics and 
                                                       Computation 259, 698 (2015))
            
        Returns : tuple of lists/scalar 
        -------        
        Scaled_Fermi_energy : float or 1d array of float (unit: unitless)
            Fermi energy w.r.t conduction band w.r.t k_BT.
            general case: inverse Fermi-Dirac integral approach.
            degenerate case: assumes metallic ('degenerate') carriers.
            return : [general case, degenerate case]
        Fermi_energy : float or 1d array of float (unit: eV)
            Fermi energy w.r.t conduction band. 
            general case: inverse Fermi-Dirac integral approach.
            degenerate case: assumes metallic ('degenrate') carriers.
            return : [general case, degenerate case]
        Screening_wave_vector : float or 1d array of float (unit: 10^6 cm^-1)
            Screening wave vector. 
            return : [general,Thomas-Fermi,  Debye]
        tau_c_by_tau_q_dis : list of float or 1d array of float list (unit: unitless)
            The tau_c/Tau_q ratio for charged dislocation (classical by quantum scattering time).
            return : [degenerate case, non degenerate case]
        Fermi_wave_vector : float or 1d array of float (unit: 10^6 cm^-1)
         Fermi wave vector.
        pop_wave_vector : float or 1d array of float (unit: 10^6 cm^-1)
            Polar optical phonon wave vector.

        """
        if inverse_half_FD_method not in ['JD_approx', 'minimax_piecewise']:
            raise ValueError(f'Requested {inverse_half_FD_method} method is not implemeted yet. Contact developer.')
        # Remove small values for the n_3d to avoid 0-division
        n_d_ = np.nan if (np.isscalar(n_d) and n_d < self.eps_n_3d) else\
            np.where(n_d < self.eps_n_3d, np.nan, n_d)   
            
        Scaled_Fermi_energy, Fermi_energy, Screening_wave_vector,\
        tau_c_by_tau_q_dis, Fermi_wave_vector, pop_wave_vector = self._3deg_properties(n_d_, 
                                                           self.alloy_params_.get('static_dielectric_constant'),
                                                           self.alloy_params_.get('e_effective_mass'), 
                                                           self.alloy_params_.get('PO_phonon_energy'),
                                                           T,inv_half_FD_method=inverse_half_FD_method
                                                           )
        return (Scaled_Fermi_energy, Fermi_energy, Screening_wave_vector, 
                tau_c_by_tau_q_dis, Fermi_wave_vector, pop_wave_vector)
        
    def calculate_3D_mobility(self, n_3d=1, n_dis:float=1, f_dis:float=0.5, T:float=300,
                              alloy_disordered_effect:bool=False,
                              dislocation_effect:bool=False,
                              piezoelectric_effect:bool=False,
                              acoustic_phonon_effect:bool=False,
                              polar_optical_phonon_effect:bool=False,
                              total_mobility:bool=True,
                              calculate_total_mobility_only:bool=False, 
                              mobility_model_version:str='v1',
                              inverse_half_FD_method:str='minimax_piecewise',
                              FermiDirac_integration_approach:str='minimax_piecewise',
                              carrier_degeneracy_limit:str='general'
                              ):
        """
        This function calculates the sheet mobility from different scattering contributions.
        The mobility models are implemented based on the following references.
            
        The considered scattering mechanism are:
            Alloy disorder limited (AD)
            Threading dislocation mediated (DIS)
            Piezoelectric effect (PE)
            Acoustic deformation potential phonon (ADP)
            Polar optical phonon (POP)

        Parameters
        ----------
        n_3d : 1D float array, optional (unit: 1e18 cm^-3)
            Array containing carrier density data for compositions. Array size
            should be same as composition arrary. The default is 1.
        n_dis : float, optional (unit: 1e8 cm^-2)
            Threading dislocation density. The default is 1.
        f_dis : float, optional (unit: unitless)
            Fraction of dislocation that contributes in scattering. 
            The default is 0.5.
        T : float, optional (unit: K)
            Temperature at which mobility calculations will be done. 
            The default is 300K.
        alloy_disordered_effect : bool, optional
            Whether to calculate alloy disordered mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        dislocation_effect : bool, optional
            Whether to calculate interface roughness effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        piezoelectric_effect : bool, optional
            Whether to calculate piezoelectric effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        acoustic_phonon_effect : bool, optional
            Whether to calculate acoustic phonon effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. 
            This includes only deformation potential mediated scattering.
            The default is False.
        polar_optical_phonon_effect : bool, optional
            Whether to calculate polar optical phonon effect mediated mobility. Or, whether to include 
            this contribution in total mobility calculation. The default is False.
        total_mobility : bool, optional
           Whether to calculate total mobility. The default is True.
        calculate_total_mobility_only : 
            Calculate only the total mobility. If False the return data also contains individual 
            specified contributions.
        mobility_model_version : str, optional [options: 'v1']
            Which mobility model to use. The default is 'v1'.
        inv_half_FD_method : str, optional [options: 'JD_approx', 'minimax_piecewise']
            The approximate method to calculate the scalled Fermi energy (E_f/k_BT) 
            using inverse Fermi-Dirac integral of order-1/2. The default is minimax_piecewise.
            If JD_approx : Joyce-Dixon approximation (APL 31, 354 (1977)).
            If minimax_piecewise : use Fukishima's minimax_piecewise approximations.
        FermiDirac_integration_approach : str, optional [options: 'num', 'minimax_piecewise', 'polylog']
            Compute the Fermi-Dirac integral. The default is minimax_piecewise.
            If num: calculated numerically, using scipy.quad. 
            If polylog: polylogaritm approach is used for FD_order > 1. For FD_order = 1 
            dilogarithm formulation is used.For FD_order=0, analytical solution 
            is used always.
            If minimax_piecewise: use Fukishima's minimax_piecewise approximation.
        carrier_degeneracy_limit : str, optional [options: 'nondegenerate', 'degenerate', 'general']
            Calculate mobilities at different carrier degenracy limit. The default
            is 'general'.
            NB: Degenerate and non-degenerate limits are only implemented for dislocation scattering.
            Contact developer to request for other scattering mechanisms.

        Returns
        -------
        pandas dataframe of compositions and mobilities (unit: cm^2 V^-1 S^-1).
            Total (or individual contributions) sheet mobility.
            
        """
        if carrier_degeneracy_limit != 'general' and self.print_info is not None:
            print('NB: Degenerate and non-degenerate limits are only implemented for dislocation scattering.')
            print('Contact developer to request for other scattering mechanisms.')
        
        self.alloy_disordered_effect_=alloy_disordered_effect
        self.dislocation_effect_=dislocation_effect
        self.piezoelectric_effect_=piezoelectric_effect
        self.acoustic_phonon_effect_=acoustic_phonon_effect
        self.polar_optical_phonon_effect_=polar_optical_phonon_effect
        self.only_total_mobility = calculate_total_mobility_only
        self.total_mobility_=total_mobility
        self.mobility_model_=mobility_model_version
        self.inverse_half_FD_method_ = inverse_half_FD_method 
        self.FD_int_approach_ = FermiDirac_integration_approach
        self.carrier_degenracy_limit_ = carrier_degeneracy_limit
        return self._calculate_3d_mobility(n_3d=n_3d, n_dis=n_dis, f_dis=f_dis, T=T)

class Plottings(_plot_mobilities):  
    """
    Plotting class for mobilitypy.
    """
    def __init__(self, save_figure_dir='.'):
        """
        Intializing mobilitypy Plotting class.

        Parameters
        ----------
        save_figure_dir : str, optional
            Directory where to save the figure. The default is current directory.

        """
        self.save_figure_directory = save_figure_dir
        _plot_mobilities.__init__(self, save_figure_dir=self.save_figure_directory)

    def plot_2d(self, data2plot, fig=None, ax=None, save_file_name=None, CountFig=None, 
                ymin=None, ymax=None, xmax=None, xmin=None, y_scale_log:bool=True, 
                show_right_ticks:bool=False, title_text:str=None, yaxis_label:str='', 
                xaxis_label:str='', color=None, color_map='viridis', ls_2d='-', 
                show_legend:bool=False, show_colorbar:bool=False, colorbar_label:str=None, 
                savefig:bool=True, vmin=None, vmax=None, show_plot:bool=True, **kwargs_savefig):  
        """
        This function plots 2d plot when providing corresponding x and y values as data2plot.

        Parameters
        ----------
        data2plot : 2D numpy array
            2D numpy array with first column as x and 2nd column as y.
        fig : matplotlib.pyplot figure instance, optional
            Figure instance to plot on. The default is None.
        ax : matplotlib.pyplot axis, optional
            Figure axis to plot on. If None, new figure will be created.
            The default is None.
        save_file_name : str, optional
            Name of the figure file. If None, figure will be not saved. 
            The default is None.
        CountFig: int, optional
            Figure count. The default is None.
        ymin : float, optional
            Minimum in y. The default is None.
        ymax : float, optional
            Maximum in y. The default is None.
        xmin : float, optional
            Minimum in x. The default is None.
        xmax : float, optional
            Maximum in x. The default is None.
        y_scale_log : bool, optional
            Use log scale for y-axis. The default is True.
        show_right_ticks : bool, optional
            Show ticks in the right axis of the figure. the default is False.
        title_text : str, optional
            Title of the figure. The default is None.
        yaxis_label : str, optional
            Y-axis label text. The default is ''.
        xaxis_label : str, optional
            x-axis label text. The default is ''.
        color : str/color, optional
            Color of plot. The default is 'gray'.
        color_map: str/ matplotlib colormap
            Colormap for plot. The default is viridis.
        ls_2d : matplotlib line style, optional
            Matplotlib line style. The default is '-'.
        show_legend : bool, optional
            If show legend or not. The default is True.
        show_colorbar : bool, optional
            Plot the colorbar in the figure or not. If fig=None, this is ignored.
            The default is False.
        colorbar_label : str, optional
            Colorbar label. The default is None. If None, ignored.
        vmin, vmax : float, optional
            vmin and vmax define the data range that the colormap covers. 
            By default, the colormap covers the complete value range of the supplied data.
        show_plot : bool, optional
            To show the plot when not saved. The default is True.
        savefig : bool, optional
            Save the plot or not. The default is True.
        **kwargs_savefig : dict
            The matplotlib keywords for savefig function.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            Figure instance. If ax is not None previously generated/passed fig instance
            will be returned. Return None, if no fig instance is inputed along with ax.
        ax : Axis instance
            Figure axis instance.
        CountFig: int or None
            Figure count.

        """
        return self._plot(data2plot, fig=fig, ax=ax, save_file_name=save_file_name, 
                          CountFig=CountFig, ymin=ymin, ymax=ymax, xmax=xmax, xmin=xmin, 
                          y_scale_log=y_scale_log, mode='plane_2d', yaxis_label=yaxis_label, 
                          title_text=title_text, xaxis_label=xaxis_label, color=color, 
                          show_right_ticks=show_right_ticks, show_legend=show_legend, 
                          ls_2d=ls_2d, color_map=color_map, show_colorbar=show_colorbar, 
                          colorbar_label=colorbar_label, savefig=savefig,
                          vmin=vmin, vmax=vmax, show_plot=show_plot, **kwargs_savefig)
    
    def plot_2d_carrier_mobilities(self, mobility_dataframe, fig=None, ax=None, save_file_name=None, CountFig=None, ymin=None, 
                                   ymax=None, xmax=None, xmin=None, y_scale_log:bool=True, mode:str= '2d_carrier_mobility',
                                   title_text:str=None, mobility_model:str='v2', annotate_pos=(0,0), annotatetextoffset=(0,-20),
                                   yaxis_label:str=r'$\mu$ ($\mathrm{cm}^2\mathrm{V}^{-1}\mathrm{s}^{-1}$)',
                                   xaxis_label:str='Composition', color=None, color_map='viridis', show_legend:bool=False, 
                                   show_right_ticks:bool=False, show_colorbar:bool=False, colorbar_label:str=None, 
                                   savefig:bool=True, vmin=None, vmax=None, show_plot:bool=True, **kwargs_savefig):
        """
        This function plots different mobility values with compositions.

        Parameters
        ----------
        mobility_dataframe : pandas dataframe or 2d array
            Pandas dataframe retured from mobility calculations when mode is '2d_carrier_mobility'.
        fig : matplotlib.pyplot figure instance, optional
            Figure instance to plot on. The default is None.
        ax : matplotlib.pyplot axis, optional
            Figure axis to plot on. If None, new figure will be created.
            The default is None.
        save_file_name : str, optional
            Name of the figure file. If None, figure will be not saved. 
            The default is None.
        CountFig: int, optional
            Figure count. The default is None.
        ymin : float, optional
            Minimum in y. The default is None.
        ymax : float, optional
            Maximum in y. The default is None.
        xmin : float, optional
            Minimum in x. The default is None.
        xmax : float, optional
            Maximum in x. The default is None.
        y_scale_log : bool, optional
            Use log scale for y-axis. The default is True.
        mode : str, optional
            Which plotting mode to use. The options are 
            '2d_carrier_mobility': To plot 2d mobility plots
            'plane_2d': general 2d plots.
        mobility_model :  str, optional [options: 'v1', 'v2']
            Which mobility model used to generate results. The data structure is 
            different for different mobility models. The default is 'v2'.
        annotate_pos : tuple, optional
            To add annotation at position on the plot. The default is (0,0).
        annotatetextoffset : tuple, optional
            To offset the annotated text from the annotate position. The default is (0, -20).
        show_right_ticks : bool, optional
            Show ticks in the right axis of the figure. the default is False.
        title_text : str, optional
            Title of the figure. The default is None.
        yaxis_label : str, optional
            Y-axis label text. The default is 'mu (cm^2V^-1s^-1)'.
        xaxis_label : str, optional
            x-axis label text. The default is 'Composition'.
        color : str/color, optional
            Color of plot. The default is 'gray'.
        color_map: str/ matplotlib colormap
            Colormap for plot. The default is viridis.
        show_legend : bool, optional
            If show legend or not. The default is True.
        show_colorbar : bool, optional
            Plot the colorbar in the figure or not. If fig=None, this is ignored.
            The default is False.
        colorbar_label : str, optional
            Colorbar label. The default is None. If None, ignored.
        vmin, vmax : float, optional
            vmin and vmax define the data range that the colormap covers. 
            By default, the colormap covers the complete value range of the supplied data.
        show_plot : bool, optional
            To show the plot when not saved. The default is True.
        savefig : bool, optional
            Save the plot or not. The default is True.
        **kwargs_savefig : dict
            The matplotlib keywords for savefig function.
        
        Raises
        ------
        ValueError
            If plot mode is unknown.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            Figure instance. If ax is not None previously generated/passed fig instance
            will be returned. Return None, if no fig instance is inputed along with ax.
        ax : Axis instance
            Figure axis instance.
        CountFig: int or None
            Figure count.

        """
        return self._plot(mobility_dataframe, fig=fig, ax=ax, save_file_name=save_file_name, 
                          CountFig=CountFig, ymin=ymin, ymax=ymax, xmax=xmax, xmin=xmin, 
                          annotate_pos=annotate_pos, annotatetextoffset=annotatetextoffset,
                          show_right_ticks=show_right_ticks,
                          y_scale_log=y_scale_log, mode= mode, yaxis_label=yaxis_label, 
                          title_text=title_text, xaxis_label=xaxis_label, color=color, 
                          mobility_model=mobility_model, color_map=color_map, 
                          show_legend=show_legend, show_colorbar=show_colorbar, 
                          colorbar_label=colorbar_label, savefig=savefig,
                          vmin=vmin, vmax=vmax, show_plot=show_plot, **kwargs_savefig)
