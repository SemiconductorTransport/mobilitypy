import numpy as np
import pandas as pd
import scipy as sc
from ._alloy_params import _AlloyParams
from ._constants import *

## ==============================================================================
class _Mobility2DCarrier(_AlloyParams):
    '''
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
    
    Ref-3: Mondal et. al., TBA
    '''
    
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN', 
                 system='ternary', print_log=None, eps_n_2d=1e-10):
        """
        Initiation function of the class _Mobility2DCarrier.
        
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
        system : string (case sensitive), optional
            Type of the alloy. E.g. 'ternary'. 
            The default is 'ternary'.
        print_log : string, optional => ['high','medium','low', None]
            Determines the level of log to be printed. The default is None.
        eps_n_2d : float, optional
            Carrier density below eps_n_2d will be considered as zero. 
            The default is 1e-10.

        Returns
        -------
        None.

        """
        self.print_info = print_log
        if self.print_info is not None: self.print_info = self.print_info.lower()

        self.eps_n_2d = eps_n_2d
            
        _AlloyParams.__init__(self, compositions=compositions, binaries=binaries, alloy=alloy)
        self._get_alloy_params(system=system)

    @staticmethod
    def _calculate_sheet_resitance(n_2d, mobility):
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
        n_2d : 1D float array (nm^-2)
            Array containing carrier density data for compositions. This can be
            a single number as well. Then all compositions will have same carrier
            density.
        mobility : 1D float array (cm^2 V^-1 s^-1)
            Array containing mobility data for compositions.

        Returns
        -------
        1D float array (ohm/square)
            Sheet resistance for compositions.

        """
        return 62415.09074/(n_2d * mobility)

    @staticmethod
    def _apply_Varshni_T_correction_2_bandgap(bandgap_0, temp:float=300, 
                                              bandgap_alpha:float=0, bandgap_beta:float=0):
        """
        This functions applies Varshni's formula for temperature correction to band gap.
        Eg(T) = Eg(T=0) - [aT^2/(T+b)]

        Parameters
        ----------
        bandgap_0 : 1D float array (eV)
            Band gap values at 0K temperature.
        temp : float, optional (K)
            Temperature in K. The default is 300K.
        bandgap_alpha : float, optional (eV/K)
            Temperature correction coefficient alpha. The default is 0.
        bandgap_beta : float, optional (K)
            Temperature correction coefficient beta. The default is 0.

        Returns
        -------
        1D float array (eV)
            The temperature corrected band gap values.

        """
        return bandgap_0 - (bandgap_alpha*temp*temp/(temp+bandgap_beta))

    def _calculate_figure_of_merit(self, n_2d, mobility, temp:float=300, mode:str='LFOM', 
                                   T_corect_bandgap:bool=False,
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
        n_2d : 1D float array (nm^-2)
            Array containing carrier density data for compositions. This can be
            a single number as well. Then all compositions will have same carrier
            density.
        mobility : 1D float array (cm^2 V^-1 s^-1)
            Array containing mobility data for compositions.
        temp : float, optional (K)
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
        1D float array (MW/cm^2)
            Figure-of-merit.

        """
        assert mode in ['LFOM'], 'Requested mode is not implemented yet' 
        bandgap_ = self.alloy_params_.get('bandgap')
        bandgap_alpha_ = self.alloy_params_.get('bandgap_alpha')
        bandgap_beta_ = self.alloy_params_.get('bandgap_beta')
        if T_corect_bandgap:
            bandgap_ = self._apply_Varshni_T_correction_2_bandgap(bandgap_, temp=temp,
                                                                  bandgap_alpha=bandgap_alpha_,
                                                                  bandgap_beta=bandgap_beta_)
        
        if direct_bandgap:
            critical_electric_field = 1.73e5*(bandgap_**2.5) # V/cm
        elif indirect_bandgap:
            critical_electric_field = 2.38e5*(bandgap_**2) # V/cm
        #print(bandgap_, n_2d)
        if mode == 'LFOM': #unit: MW/cm^2
            return 1.602176634e-11 * n_2d * mobility * critical_electric_field * critical_electric_field
        
    def _calculate_sheet_mobility(self, n_2d=0.1, rms_roughness=0.1, corr_len=1, 
                                  n_dis=1, f_dis=0.1, T=300, return_sc_rates:bool=False):
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
        
        Ref-3: Mondal et. al., TBA
            
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
        rms_roughness => nm^-1
        corr_len => nm^-1
        n_dis => nm^-2
        f_dis => unit less
        E_pop => eV

        Parameters
        ----------
        n_2d : 1D float array or float (nm^-2)
            Array containing carrier density data for compositions. This can be
            a single number as well. Then all compositions will have same carrier
            density.
        rms_roughness : float, optional (nm^-1)
            Interface root-mean-squared roughness for interface-roughness scattering
            contribution. The default is 0.1.
        corr_len : float, optional (nm^-1)
            Correlation length of interface roughness. The default is 1.
        n_dis : float, optional (nm^-2)
            Threading dislocation density. The default is 1.
        f_dis : float, optional
            Fraction of dislocation that contributes in scattering. 
            The default is 0.1.
        T : float, optional (K)
            Temperature at which mobility calculations will be done. 
            The default is 300K.
        return_sc_rates : float, optional 
            Return the scattering rates values.The default is False.

        Returns
        -------
        pandas dataframe of compositions and mobilities
            Total (or individual contributions) sheet mobility. If return_sc_rates=True,
            then scattering rates are also returned.

        """
        e_effective_mass = self.alloy_params_.get('e_effective_mass') 
        static_dielectric_constant = self.alloy_params_.get('static_dielectric_constant') 
        high_frequency_dielectric_constant = self.alloy_params_.get('high_frequency_dielectric_constant')
        lattice_a = self.alloy_params_.get('lattice_a0') * 0.1 # angstrom to nm
        lattice_c = self.alloy_params_.get('lattice_c0') * 0.1 # angstrom to nm
        sc_potential = self.alloy_params_.get('alloy_scattering_potential') 
        LA_velocity = self.alloy_params_.get('LA_phonon_velocity')
        mass_densitty = self.alloy_params_.get('mass_density')
        deformation_pot = self.alloy_params_.get('deformation_potential')
        electromech_coupling_sqr = self.alloy_params_.get('electromechanical_coupling_const')
        POP_energy = self.alloy_params_.get('PO_phonon_energy')

        if isinstance(n_2d, int) or isinstance(n_2d, float):
            n_2d = [n_2d] * len(self.comps_)
        
        mobility = {}
        for ii in range(len(self.comps_)):
            #print(n_2d[ii])
            mobility[ii] = {'comp': f'{self.comps_[ii]:.3f}'}
            self._set_params(e_effective_mass[ii], static_dielectric_constant[ii], 
                             high_frequency_dielectric_constant[ii],
                             lattice_c[ii], lattice_a[ii], sc_potential[ii], self.comps_[ii],
                             n_2d[ii], rms_roughness, corr_len, n_dis, f_dis, 
                             T, electromech_coupling_sqr[ii], deformation_pot[ii], mass_densitty[ii], LA_velocity[ii], POP_energy[ii])
            self._print_database_params()
            # mobility unit: cm^2 V^-1 S^-1
            if self.print_info is not None: print(f'- Composition: {self.comps_[ii]:.5f}')
            
            if return_sc_rates: mobility[ii]['m0_by_e'] = self.m0_by_e_
                
            total_inv_sc = 0
            if self.only_total_mobility:
                if self.print_info is not None: print('\t-- Calculating only total mobility')
                
                if self.alloy_disordered_effect_: total_inv_sc += self._inv_tau_ado()                   
                if self.interface_roughness_effect_: total_inv_sc += self._inv_tau_ifr()                   
                if self.dislocation_effect_: total_inv_sc += self._inv_tau_dis()
                if self.polar_optical_phonon_effect_: total_inv_sc += self._inv_tau_pop()
                if self.acoustic_phonon_effect_: # 1/tau_AP = 1/tau_DP + 1/tau_PE
                    total_inv_sc = total_inv_sc + self._inv_tau_pe()+self._inv_tau_dp()
                else:
                    if self.deformation_potential_effect_: total_inv_sc += self._inv_tau_dp()
                    if self.piezoelectric_effect_: total_inv_sc += self._inv_tau_pe()
                mobility[ii]['TOT'] = self._mobility_calculator(total_inv_sc) 
                if return_sc_rates: mobility[ii]['TOT_sc'] = total_inv_sc
            else:
                if self.alloy_disordered_effect_:
                    if self.print_info is not None: print('\t-- Calculating alloy-disordered mobility')
                    inv_sc = self._inv_tau_ado()
                    total_inv_sc += inv_sc
                    mobility[ii]['AD'] = self._mobility_calculator(inv_sc)
                    if return_sc_rates: mobility[ii]['AD_sc'] = inv_sc
                    
                if self.interface_roughness_effect_:
                    if self.print_info is not None: print('\t-- Calculating interface roughness effect mobility')
                    inv_sc = self._inv_tau_ifr()
                    total_inv_sc += inv_sc
                    mobility[ii]['IFR'] = self._mobility_calculator(inv_sc)
                    if return_sc_rates: mobility[ii]['IFR_sc'] = inv_sc
                    
                if self.dislocation_effect_:
                    if self.print_info is not None: print('\t-- Calculating dislocation effect mobility')
                    inv_sc = self._inv_tau_dis()
                    total_inv_sc += inv_sc
                    mobility[ii]['DIS'] = self._mobility_calculator(inv_sc)
                    if return_sc_rates: mobility[ii]['DIS_sc'] = inv_sc
                    
                if self.polar_optical_phonon_effect_:
                    if self.print_info is not None: print('\t-- Calculating polar optical phonon effect mobility')
                    inv_sc = self._inv_tau_pop()
                    total_inv_sc += inv_sc
                    mobility[ii]['POP'] = self._mobility_calculator(inv_sc)
                    if return_sc_rates: mobility[ii]['POP_sc'] = inv_sc
 
                if self.acoustic_phonon_effect_:
                    inv_sc_dp = self._inv_tau_dp()
                    inv_sc_pe = self._inv_tau_pe()
                    if self.print_info is not None: print('\t-- Calculating acoustic effect mobility')
                    inv_sc = inv_sc_dp + inv_sc_pe
                    total_inv_sc += inv_sc # 1/tau_AP = 1/tau_DP + 1/tau_PE
                    mobility[ii]['AP'] = self._mobility_calculator(inv_sc)
                    if return_sc_rates: mobility[ii]['AP_sc'] = inv_sc
                else:
                    if self.deformation_potential_effect_:
                        inv_sc_dp = self._inv_tau_dp()
                        total_inv_sc += inv_sc_dp
                    if self.piezoelectric_effect_:
                        total_inv_sc += inv_sc_pe
                        inv_sc_pe = self._inv_tau_pe()
                    
                if self.deformation_potential_effect_:
                    if self.print_info is not None: print('\t--- Calculating deformation potential effect mobility')
                    mobility[ii]['DP'] = self._mobility_calculator(inv_sc_dp)
                    if return_sc_rates: mobility[ii]['DP_sc'] = inv_sc_dp
                    
                if self.piezoelectric_effect_:
                    if self.print_info is not None: print('\t--- Calculating piezoelectric effect mobility')
                    mobility[ii]['PE'] = self._mobility_calculator(inv_sc_pe)
                    if return_sc_rates: mobility[ii]['PE_sc'] = inv_sc_pe
 
                if self.total_mobility_:
                    if self.print_info is not None: print('\t-- Calculating total mobility')
                    mobility[ii]['TOT'] = self._mobility_calculator(total_inv_sc)
                    if return_sc_rates: mobility[ii]['TOT_sc'] = total_inv_sc
                    
            if self.print_info is not None: print(f'{"="*72}')
        return pd.DataFrame.from_dict(mobility, orient='index')
        
    def _set_params(self, m_star, eps_s, eps_h, c_lattice, a_lattice, sc_potential, 
                    alloy_composition, n_2d, rms_roughness, corr_len, n_dis, f_dis, 
                    T, K_square, E_D, mass_density, v_LA, E_pop):
        """
        This function sets the parameters for mobility calculations.
        """
        self.m_star_ = m_star
        self.eps_s_ = eps_s
        self.eps_h_ = eps_h
        self.c_lp = c_lattice
        self.a_lp = a_lattice
        self.sc_potential_ = sc_potential
        self.comp_ = alloy_composition 
        self.n_2d_ = n_2d
        self.corr_len_ = corr_len
        self.rms_roughness_ = rms_roughness
        self.n_dislocation_ = n_dis
        self.f_dislocation_ = f_dis
        self.temp_ = T
        self.K_sqr = K_square
        self.E_d = E_D
        self.mass_density_ = mass_density
        self.v_LA = v_LA
        self.E_pop = E_pop
        self._get_derived_params()

    def _get_derived_params(self):
        #-------------- derived parameters --------------
        tmp_ = self.m_star_ / self.eps_s_ # unit-less
        self.k_F = np.sqrt(2*pi_*self.n_2d_) # nm^-1
        self.q_TF = fact_q_TF * tmp_ # nm^-1
        self.b_ = fact_b*(self.n_2d_*tmp_)**(1/3) # nm^-1
        self.fact_1 =  fact_irf_dis * tmp_ /self.eps_s_ # s^-1
        self.omega = 0.8660254037844386 * self.a_lp**2 * self.c_lp # sqrt(3)/2 * a^2 c => nm^3
        self.m0_by_e_ = self.m_star_ * m0_by_e
        self.k_0 = fact_pop_k0*np.sqrt(self.m_star_*self.E_pop) # nm^-1

    def _print_database_params(self):
        """
        This function prints the log of model descriptions.

        Returns
        -------
        None.

        """
        if self.print_info == 'high':
            print(f'- Composition={self.comp_:.5f}')
            print(f'\t-- a={self.a_lp:.5f} nm | c={self.c_lp:.5f} nm | m*={self.m_star_:.5f} m0 | eps_s={self.eps_s_:.5f} eps0 | eps_h={self.eps_h_:.5f} eps0')
            print(f'\t-- Mass density={self.mass_density_:.2f} | scattering potential={self.sc_potential_:.2f} eV | 2deg density={self.n_2d_:.5f} nm^-2 | T={self.temp_:.1f} K')
            print(f'\t-- Interface rms roughness={self.rms_roughness_:.3f} nm | correlation length={self.corr_len_:.3f} nm')
            print(f'\t-- Dislocation density={self.n_dislocation_:.4f} nm^-2 | dislocation occupancy={self.f_dislocation_:.1f}')
            print(f'\t-- Electromechanical coupling coefficient={self.K_sqr:.5f} | deformation potential={self.E_d:.5f}')
            print(f'\t-- Longitudinal acoustic phonon velocity={self.v_LA:.2f} m/s | polar optical phonon energy={self.E_pop:.5f} eV')
            print(f'\t-- Fermi wave vector={self.k_F} | b={self.b_}')
            print('')

    def _form_factor(self, x, mode=None):
        """
        This function calculates Fang-Howard form-factors.

        Parameters
        ----------
        x : float
            Scattering states. x=sing(theta/2), theta=scattering angle.
        mode : string, optional ['IRF', 'DIS', 'DP', 'PE', 'POP']
            FW form factor which scattering mechanism. 
            The default is None. If None, None is returned.

        Returns
        -------
        float/None
            FW form factor.

        """
        # Fang-Howard form-factor
        if mode == 'IRF':
            # eta(u) = b/(b+2*k_f*u)
            eta = self.b_/(self.b_ + 2*self.k_F*x)
            # G(eta) = (2*eta^3 + 3*eta^2 + 3*eta) / 8
            return eta*(eta*(2*eta+3)+3)/8
        elif mode == 'DIS':
            return 1
        elif mode in ['DP', 'PE']:
            # eta(u) = b/(b+2*k_f*u)
            eta = self.b_/(self.b_ + 2*self.k_F*x)
            return eta*eta*eta
        elif mode == 'POP':
            # eta(u) = b/(b+k_0)
            eta = self.b_/(self.b_ + self.k_0)
            # G(eta) = (2*eta^3 + 3*eta^2 + 3*eta) / 8
            return eta*(eta*(2*eta+3)+3)/8
        else:
            return None

    def _int_f_str_denomenator(self, x, mode=None):
        return (x + self.q_TF*self._form_factor(x, mode=mode)/2/self.k_F)**2 * np.sqrt(1 - x**2)

    def _int_f_phon_(self, x, mode=None):
        return x**3/((2*self.k_F*x + self.q_TF*self._form_factor(x, mode=mode))**2 * np.sqrt(1 - x**2))

    # ----- interface roughness ----------------
    def _inv_tau_ifr_f(self, x):
        return (x**4 * np.exp(-(self.corr_len_ * self.k_F * x)**2) / 
                self._int_f_str_denomenator(x, mode='IRF'))

    def _inv_tau_ifr_int(self):
        return sc.integrate.quad(self._inv_tau_ifr_f, 0, 1)[0]

    def _inv_tau_ifr(self):
        if self.n_2d_ < self.eps_n_2d: return 0
        fact_2 = (self.rms_roughness_ * self.corr_len_ * self.n_2d_)**2 / 8 
        return self.fact_1 * fact_2 * self._inv_tau_ifr_int() 
        
    # ----------- dislocation -------------------
    def _inv_tau_dis_f(self, x):
        return 1/self._int_f_str_denomenator(x, mode='DIS')

    def _inv_tau_dis_int(self):
        return sc.integrate.quad(self._inv_tau_dis_f, 0, 1)[0]

    def _inv_tau_dis(self):
        if self.n_2d_ < self.eps_n_2d: return 0
        fact_2 = self.n_dislocation_ * self.f_dislocation_**2 / (4*pi_* self.k_F**4 * self.c_lp**2)
        return self.fact_1 * fact_2 * self._inv_tau_dis_int()

    # ----------- alloy disordered -------------------
    def _inv_tau_ado(self):
        if (self.comp_ < 1e-8) or (self.n_2d_ < self.eps_n_2d) or ((1-self.comp_)<1e-8): return 0
        fact_2 = self.m_star_ * self.omega * self.sc_potential_**2 * self.comp_ * (1-self.comp_) * self.b_
        #print(fact_alloy, self.m_star_ ,self.omega ,self.sc_potential_, self.comp_ ,self.b_)
        return fact_alloy * fact_2

    # ----------- deformation potential -------------------
    def _inv_tau_dp_f(self, x):
        return x*self._int_f_phon_(x, mode='DP')

    def _inv_tau_dp_int(self):
        return sc.integrate.quad(self._inv_tau_dp_f, 0, 1)[0]
        
    def _inv_tau_dp(self):
        if self.n_2d_ < self.eps_n_2d: return 0
        fact_2 = (3*self.E_d*self.E_d*self.m_star_*self.temp_*self.k_F*self.k_F*self.b_)/(self.mass_density_*self.v_LA*self.v_LA) 
        return fact_phonon * 1e9*fact_2 * self._inv_tau_dp_int() 

    # ----------- piezoelectric -------------------
    def _inv_tau_pe_f(self, x):
        return self._form_factor(x, mode='PE')*self._int_f_phon_(x, mode='PE')

    def _inv_tau_pe_int(self):
        return sc.integrate.quad(self._inv_tau_pe_f, 0, 1)[0]
        
    def _inv_tau_pe(self):
        if self.n_2d_ < self.eps_n_2d: return 0
        fact_2 = (4*self.k_F*self.K_sqr*self.m_star_*self.temp_)/(eps_0*self.eps_s_)
        return fact_phonon * 1e-9*fact_2 * self._inv_tau_pe_int() 

    # ----------- polar optical phonon -------------------
    def _inv_tau_pop(self):
        if self.n_2d_ < self.eps_n_2d: return 0
        eps_star = 1/(1/self.eps_h_ - 1/self.eps_s_)
        fact_2 = self.m_star_ * self.E_pop * self._form_factor(None, mode='POP')/eps_star/self.k_0 
        yy = fact_pop_y*self.n_2d_/self.m_star_/self.temp_
        fact_3 = yy / ((np.exp(self.E_pop*e_charge/k_B/self.temp_) - 1) * (1+yy-np.exp(-yy))) 
        #print(yy, fact_3, fact_pop * fact_2 * fact_3)
        return fact_pop * fact_2 * fact_3 

    # Scattering rate to mobility calculation
    def _mobility_calculator(self, inverse_scattering):  
        """
        This function calculates the sheet mobility from different scattering contributions.
        """
        # 1e4 is unit conversion from m^2 to cm^2
        # unit: cm^2 V^-1 S^-1
        return 1e4/(self.m0_by_e_ * inverse_scattering) if inverse_scattering else np.nan