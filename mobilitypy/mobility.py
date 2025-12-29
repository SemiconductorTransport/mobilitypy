from .src import _AlloyParams, _MobilityCarrier, _Mobility2DCarrier, _Mobility3DCarrier 
from .utilities import _plot_mobilities
import numpy as np

## ==============================================================================
class AlloyParams(_AlloyParams):
    '''
    The functions in this class calculates the parameters for alloy from their
    binary components.
    '''
    def __init__(self):
        pass
            
    def get_alloy_params(self, system='ternary', compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN'):
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

        Returns
        -------
        1D float array
            Parameters for alloy.

        """
        _AlloyParams.__init__(self, compositions=compositions, binaries=binaries, alloy=alloy)
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
                 system='ternary', psedomorphic_strain=False, substrate=None,
                 alloy_type='WZ', eps_n_2d=1e-10, print_log=None):
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
        psedomorphic_strain : bool, optional
            Whether to consider pseudomorphic strain.
            The default is False.
        substrate : string or float, optional (unit: Angstrom)
            The substrate name (if string, warning: the name should be in the database) 
            or the substrate in-plane lattice parameter (if float, Angstrom unit).
            The default is None. Error will be raised if substrate=None and 
            psedomorphic_strain=True.
        alloy_type :  str, optional (case insensitive)
            The crystal type of alloy. This will be considered when calculating
            parameters like Poisson ratio etc.
            Use following abbreviation name:
                for wurtzite use 'WZ' or 'wz'.
                for zincblende use 'ZB' or 'zb'.
                for diamond use 'DM' or 'dm'.
            The default is 'WZ'. 
        eps_n_2d : float, optional (unit: nm^-2)
            Carrier density below eps_n_2d will be considered as zero. 
            The default is 1e-10 nm^-2 == 1e4 cm^-2.
        print_log : string, optional => ['high','medium','low', None]
            Determines the level of log to be printed. The default is None.

        Returns
        -------
        None.

        """
        if (psedomorphic_strain == True) and (substrate is None):
            raise ValueError('substrate tag can not be None when psedomorphic_strain=True.')
        _MobilityCarrier.__init__(self, compositions=compositions, binaries=binaries, 
                                  alloy=alloy, system=system, psedomorphic_strain=psedomorphic_strain, 
                                  substrate=substrate,alloy_type=alloy_type,
                                  print_log=print_log, eps_n=eps_n_2d)
        _Mobility2DCarrier.__init__(self)
        
    def calculate_sheet_mobility(self, n_2d=0.1, rms_roughness=0.1, corr_len=1, n_dis=1, f_dis=0.1, 
                                 T=300, alloy_disordered_effect:bool=False,
                                 interface_roughness_effect:bool=False,
                                 dislocation_effect:bool=False,
                                 deformation_potential_effect:bool=False, 
                                 piezoelectric_effect:bool=False,
                                 acoustic_phonon_effect:bool=False,
                                 polar_optical_phonon_effect:bool=False,
                                 total_mobility:bool=True,
                                 calculate_total_mobility_only:bool=False,
                                 return_sc_rates:bool=False,
                                 mobility_model='Bassaler'):
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
        
        Units:
            c_lattice => in nm
            a_lattice => in nm
            sc_potential => in eV
            n_2d => in nm^-2
            rms_roughness => nm
            corr_len => nm
            n_dis => nm^-2
            f_dis => unit less
            E_pop => eV

        Parameters
        ----------
        n_2d : 1D float array or float, optional (unit: nm^-2)
            Array containing carrier density data for compositions. This can be
            a single number as well. Then all compositions will have same carrier
            density. The default is 0.1.
        rms_roughness : float, optional (unit: nm)
            Interface root-mean-squared roughness for interface-roughness scattering
            contribution. The default is 0.1.
        corr_len : float, optional (unit: nm)
            Correlation length of interface roughness. The default is 1.
        n_dis : float, optional (unit: nm^-2)
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
        mobility_model : str, optional
            Which mobility model to use. The default is 'Bassaler'.
            The mobility is implemented based on following publications:
            'Bassaler':
                J. Bassaler, J. Mehta, I. Abid, L. Konczewicz, S. Juillaguet, S. Contreras, S. Rennesson, 
                S. Tamariz, M. Nemoz, F. Semond, J. Pernot, F. Medjdoub, Y. Cordier, P. Ferrandis, 
                Al-Rich AlGaN Channel High Electron Mobility Transistors on Silicon: A Relevant Approach for High 
                Temperature Stability of Electron Mobility. Adv. Electron. Mater. 2024, 2400069. 
                https://doi.org/10.1002/aelm.202400069

        Returns
        -------
        pandas dataframe with compositions and mobility (unit: cm^2 V^-1 S^-1) columns.
            Total (or individual contributions) sheet mobility. If return_sc_rates=True,
            then scattering rates are also returned.

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
    def sc_rate_2_mobility(mstar0_by_e, inverse_scattering):
        # Scattering rate to mobility calculation 
        """
        This function calculates sheet mobility from scattering rate.
        
        Parameters
        ----------
        mstar0_by_e : float/array
            Carrier effective mass in m0 unit multiplied by m0_by_e. 
            mstar0_by_e = m* X m0 / e
        inverse_scattering : float/array
            Scattering rate. 

        Returns
        -------
        float/array (unit: cm^2 V^-1 S^-1)
            Total (or individual contributions) sheet mobility. 

        """
        # Scattering rate to mobility calculation
        # 1e4 is unit conversion from m^2 to cm^2
        # unit: cm^2 V^-1 S^-1
        tau = mstar0_by_e * inverse_scattering
        tau[tau<1e-30] = np.nan
        return 1e4/tau

    def calculate_sheet_resitance(self, n_2d, mobility):
        """
        This function calculates the sheet resistance.
        
        Units:
        n_2d => in nm^-2
        e => 1.602176634e-19 C
        mobility (mu) => cm^2 V^-1 s^-1
        
        1 coulomb/volt = 1 second/ohm
        1 ohm = 1 C^-1.V.s

        R = 1/(e * n_2d * mu) ohm/square
          = 1/(1.602176634e-19*1e14 *n_2d * mu C.cm^-2.cm^2.V^-1.S^-1) 
          = 62415.09074/(n_2d * mu) ohm/square

        Parameters
        ----------
        n_2d : 1D float array (unit: nm^-2)
            Array containing carrier density data for compositions. This can be
            a single number as well. Then all compositions will have same carrier
            density.
        mobility : 1D float array (unit: cm^2 V^-1 s^-1)
            Array containing mobility data for compositions.

        Returns
        -------
        1D float array (unit: ohm/square)
            Sheet resistance for compositions.

        """
        return self._calculate_sheet_resitance(n_2d, mobility)

    def calculate_figure_of_merit(self, n_2d, mobility, temp:float=300,
                                   mode:str='LFOM', T_corect_bandgap:bool=False, 
                                   direct_bandgap:bool=True, indirect_bandgap:bool=False):
        """
        This function calculates the figure-of-merit (FOM). Available FOMs are
        LFOM: Lateral figure-of-merit
        
        Ref: J. L. Hudgins, G. S. Simin, E. Santi and M. A. Khan, 
        "An assessment of wide bandgap semiconductors for power devices," 
        in IEEE Transactions on Power Electronics, vol. 18, no. 3, pp. 907-914, 
        May 2003, doi: 10.1109/TPEL.2003.810840.
        
        direct_bandgap_critical_electric_field = 1.73e5*(bandgap_**2.5) # V/cm
        indirect_bandgap_critical_electric_field = 2.38e5*(bandgap_**2.5) # V/cm
        
        Units:
        bandgap_ => in eV.
        temp => in K
        n_2d => in nm^-2
        E_cr => in V/cm
        e => 1.602176634e-19 C
        mobility (mu) => cm^2 V^-1 s^-1
        
        
        LFOM = e*n_2d*mu*E_cr^2 = 1.602e-19 C * 1e14 cm^-2 * cm^2 V^-1 s^-1 * V^2cm^-2
                                = 1.602e-5 CVs^-1cm^-2
                                = 1.602e-5 Wcm^-2    #1 watts = 1 coulombs*volt/second
                                = 1.602e-11 MW/cm^2
        
        Parameters
        ----------
        n_2d : 1D float array (unit: nm^-2)
            Array containing carrier density data for compositions. This can be
            a single number as well. Then all compositions will have same carrier
            density.
        mobility : 1D float array (unit: cm^2 V^-1 s^-1)
            Array containing mobility data for compositions.
        temp : float, optional (unit: K)
            Temperature for band gap correction. The default is 300K.
        mode : str, optional (['LFOM'])
            The figure-of-merit name. The default is 'LFOM'.
        T_corect_bandgap : bool, optional
            Apply temperature correction to bandgap or not. The default is False.
        direct_bandgap : bool, optional
            If the bandgap is direct bandgap or not. The default is True.
        indirect_bandgap : bool, optional
            If the bandgap is indirect bandgap or not.. The default is False.

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
    
    Ref-1: Rajan et al., Appl. Phys. Lett. 88, 042103 (2006) => alloy disorder, polar optical phonon
    Ref-2: DJ. and UKM., PRB 66, 241307(R) (2002) and DJ. et al., PRB 67, 153306 (2003)  => Dislocation
    Ref-3: Debdeep Jena's thesis, Chapter-6 APPENDIX, Sec. Three-dimensional carriers => Acoustic phonon 
    
    """
    
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN', 
                 system='ternary', psedomorphic_strain=False, substrate=None,
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
        psedomorphic_strain : bool, optional
            Whether to consider pseudomorphic strain.
            The default is False.
        substrate : string or float, optional (unit: Angstrom)
            The substrate name (if string, warning: the name should be in the database) 
            or the substrate in-plane lattice parameter (if float, Angstrom unit).
            The default is None. Error will be raised if substrate=None and 
            psedomorphic_strain=True.
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
        if (psedomorphic_strain == True) and (substrate is None):
            raise ValueError('substrate tag can not be None when psedomorphic_strain=True.')
        _MobilityCarrier.__init__(self, compositions=compositions, binaries=binaries, 
                                  alloy=alloy, system=system, psedomorphic_strain=psedomorphic_strain, 
                                  substrate=substrate,alloy_type=alloy_type,
                                  print_log=print_log, eps_n=eps_n_3d)
        _Mobility3DCarrier.__init__(self)
        
    def calculate_3D_mobility(self, n_3d=1, n_dis=1, f_dis=0.5, T=300,
                              alloy_disordered_effect:bool=False,
                              dislocation_effect:bool=False,
                              piezoelectric_effect:bool=False,
                              acoustic_phonon_effect:bool=False,
                              polar_optical_phonon_effect:bool=False,
                              total_mobility:bool=True,
                              calculate_total_mobility_only:bool=False
                              ):
        """
        This function calculates the sheet mobility from different scattering contributions.
        The mobility models are implemented based on the following references.
        
        Ref-1: Rajan et al., Appl. Phys. Lett. 88, 042103 (2006) => alloy disorder, polar optical phonon
        Ref-2: DJ. and UKM., PRB 66, 241307(R) (2002) and DJ. et al., PRB 67, 153306 (2003)  => Dislocation
        Ref-3: Debdeep Jena's thesis, Chapter-6 APPENDIX, Sec. Three-dimensional carriers => Acoustic phonon 
            
        The considered scattering mechanism are:
            Alloy disorder limited (AD)
            Threading dislocation mediated (DIS)
            Piezoelectric effect (PE)
            Acoustic phonon (AP) : Deformation potential mediated
            Polar optical phonon (POP)
        
        Units:
            n_3d => in 1e18 cm^-3
            c_lattice => in A
            a_lattice => in A
            sc_potential => in eV
            n_dis => 1e10 cm^-2
            f_dis => unit less
            E_pop => eV

        Parameters
        ----------
        n_3d : 1D float array, optional (unit: 1e18 cm^-3)
            Array containing carrier density data for compositions. Array size
            should be same as composition arrary. The default is 1.
        n_dis : float, optional (unit: 1e10 cm^-2)
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

        Returns
        -------
        pandas dataframe of compositions and mobilities (unit: cm^2 V^-1 S^-1).
            Total (or individual contributions) sheet mobility.
            
        """
        
        self.alloy_disordered_effect_=alloy_disordered_effect
        self.dislocation_effect_=dislocation_effect
        self.piezoelectric_effect_=piezoelectric_effect
        self.acoustic_phonon_effect_=acoustic_phonon_effect
        self.polar_optical_phonon_effect_=polar_optical_phonon_effect
        self.only_total_mobility = calculate_total_mobility_only
        self.total_mobility_=total_mobility
        return self._calculate_3d_mobility(n_3d=n_3d, n_dis=n_dis, f_dis=f_dis, T=T)

    def Dislocation_3D_Tc_Tq_ratio(self, eps_s, n_3d, m_star):
        """
        Calculate the ratio of classical (or momentum) tp quantum scattering times
        due to charged dislocation scattering.
        
        Ref: DJ. and UKM., PRB 66, 241307(R) (2002) and DJ. et al., PRB 67, 153306 (2003) 

        Parameters
        ----------
        eps_s : float or 1d array of float (unit: epsilon_0)
            Static dielectic constants of the material. In the unit of vacumm permitivity.
        n_3d : float or 1d array of float (unit: 1E18 cm^-3 )
            3DEG density.
        m_star : float or 1d array of float (unit: m0)
            Carrier effective mass. In the unit of m0.

        Returns
        -------
        float or 1d array of float (unitless)
            The tau_c/Tau_q ratio (classical by quantum scattering time).

        """
        return self._ratio_dis_tc_tq(eps_s, n_3d, m_star)

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
                                   title_text:str=None, mobility_model:str='Bassaler', annotate_pos=(0,0), annotatetextoffset=(0,-20),
                                   yaxis_label:str=r'$\mu$ ($\mathrm{cm}^2\mathrm{V}^{-1}\mathrm{s}^{-1}$)',
                                   xaxis_label:str='Composition', color=None, color_map='viridis', show_legend:bool=False, 
                                   show_right_ticks:bool=False, show_colorbar:bool=False, colorbar_label:str=None, 
                                   savefig:bool=True, vmin=None, vmax=None, show_plot:bool=True, **kwargs_savefig):
        """
        This function plots the results.

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
        mobility_model :  str, optional
            Which mobility model used to generate results. The data structure is 
            different for different mobility models. The default is 'Bassaler'.
        annotate_pos : tuple, optional
            To add annotation at position on the plot. The default is (0,0).
        annotatetextoffset : tuple, optional
            To offset the annotated text from the annotate position. The default is (0, -20).
        show_right_ticks : bool, optional
            Show ticks in the right axis of the figure. the default is False.
        title_text : str, optional
            Title of the figure. The default is None.
        yaxis_label : str, optional
            Y-axis label text. The default is 'mu (cm^2V^-1s-1$)'.
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
