import numpy as np
import pandas as pd
import scipy.integrate as integrate
from ._constants import *

## ==============================================================================
class _Mobility2DCarrier:
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
    
    Ref-3: Mondal et. al., Interplay of carrier density and mobility in Al-rich (Al,Ga)N-channel HEMTs: 
    Impact on high-power device performance potential. APL Electronic Devices 1, 026117 (2025)
    https://doi.org/10.1063/5.0277051
    '''
    
    def __init__(self):
        """
        Initiation function of the class _Mobility2DCarrier.
        
        Returns
        -------
        None.

        """
        self.eps_n_2d = self.eps_n

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def _calculate_sheet_mobility(self, n_2d=10, rms_roughness=0.1, corr_len=1, 
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
        return_sc_rates : float, optional 
            Return the scattering rates values.The default is False.

        Returns
        -------
        pandas dataframe of compositions and mobilities (unit: cm^2 V^-1 S^-1)
            Total (or individual contributions) sheet mobility. If return_sc_rates=True,
            then scattering rates (10^12 s^-1) and m_star_by_e (10^-12 V.m^-2.s^2) are also returned.

        """
        carrier_effective_mass = self.alloy_params_.get('carrier_effective_mass') 
        static_dielectric_constant = self.alloy_params_.get('static_dielectric_constant') 
        high_frequency_dielectric_constant = self.alloy_params_.get('high_frequency_dielectric_constant')
        lattice_a = self.alloy_params_.get('lattice_a0') 
        lattice_c = self.alloy_params_.get('lattice_c0') 
        sc_potential = self.alloy_params_.get('alloy_scattering_potential') 
        LA_velocity = self.alloy_params_.get('LA_phonon_velocity')
        mass_densitty = self.alloy_params_.get('mass_density')
        deformation_pot = self.alloy_params_.get('CB_deformation_potential')
        electromech_coupling_sqr = self.alloy_params_.get('electromechanical_coupling_const')
        POP_energy = self.alloy_params_.get('PO_phonon_energy')  
        isotropic_Poisson_ratio = self.alloy_params_.get('isotropic_Poisson_ratio')   
        n_ion_impurity = 0 # not implemented yet. 

        if isinstance(n_2d, int) or isinstance(n_2d, float):
            n_2d = [n_2d] * len(self.comps_)
        
        mobility = {}
        for ii in range(len(self.comps_)):
            #print(n_2d[ii])
            mobility[ii] = {'comp': f'{self.comps_[ii]:.3f}'}
            self._set_params(carrier_effective_mass[ii], static_dielectric_constant[ii], 
                             high_frequency_dielectric_constant[ii],
                             lattice_c[ii], lattice_a[ii], sc_potential[ii], self.comps_[ii],
                             n_2d[ii], rms_roughness, corr_len, n_dis, f_dis, n_ion_impurity, 
                             T, electromech_coupling_sqr[ii], deformation_pot[ii], 
                             mass_densitty[ii], LA_velocity[ii], POP_energy[ii],
                             isotropic_Poisson_ratio[ii])
            self._print_database_params()
            # mobility unit: cm^2 V^-1 S^-1
            if self.print_info is not None: print(f'- Composition: {self.comps_[ii]:.5f}')
            
            if return_sc_rates: mobility[ii]['m_star_by_e'] = self.m_star_by_e_
                
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
                    if return_sc_rates: 
                        mobility[ii]['DIS_sc'] = inv_sc
                    if self.mobility_model_ == 'v2':
                        inv_sc = self._inv_tau_dis_strain()
                        total_inv_sc += inv_sc
                        mobility[ii]['DIS_Strain'] = self._mobility_calculator(inv_sc)
                        if return_sc_rates: 
                            mobility[ii]['DIS_Strain_sc'] = inv_sc
                    
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
                    if self.print_info is not None: print('\t--- Calculating piezoelectric phonon effect mobility')
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
                    n_ion_impurity, T, K_square, E_D, mass_density, 
                    v_LA, E_pop, poisson_ratio):
        """
        This function sets the parameters for mobility calculations.
        """
        self._set_params_general(m_star, eps_s, eps_h, c_lattice, a_lattice, sc_potential, 
                                 n_dis, f_dis, n_ion_impurity, mass_density, v_LA, E_pop, 
                                 E_D, K_square, poisson_ratio, T)
        #==========================================
        self.comp_ = alloy_composition 
        self.n_2d_ = n_2d
        #==========================================
        self.corr_len_ = corr_len
        self.rms_roughness_ = rms_roughness
        #==========================================
        self._get_derived_params()

    def _get_derived_params(self):
        #-------------- derived parameters --------------
        self.m_star_by_eps_s = self.m_star_ / self.eps_s_ 
        self.m_star_by_eps_s_square = self.m_star_ / self.eps_s_ / self.eps_s_
        # n_2d = n_2d_ * 1e12 cm^-2
        # sqrt(2*pi_) = 2.5066282746310002
        self.k_F = 2.5066282746310002 * np.sqrt(self.n_2d_) # 1e6 cm^-1
        #e^2 * m0 / (2pi * eps_0 * h_bar^2) * 1e-2 = 377.94522518130917 * 1e6 cm^-1 
        #self.q_TF = 377.94522518130917 * tmp_ # 1e6 cm^-1 
        # 0.5*(e_mass*e_charge**2)/((2*pi_)**(3/2)*h_bar**2*eps_0)*1e-8 = 75.3891649487971
        self.q_TF_by_2k_F = 75.3891649487971 * self.m_star_by_eps_s / np.sqrt(self.n_2d_)
        #(33*e_charge**2*e_mass/(8*eps_0*h_bar**2) * 1e-2 * 1e12)**(1/3) = 2139.6573408935274 cm^-1
        self.b_ = 21.396573408935274*(self.n_2d_ * self.m_star_by_eps_s)**(1/3) # 1e6 cm^-1
        # np.sqrt(2*e_mass*e_charge/h_bar**2)*1e-2 = 51.23167219674931 1e6 cm^-1
        self.k_pop = 51.23167219674931*np.sqrt(self.m_star_*self.E_pop) # 1e6 cm^-1

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
            print(f'\t-- Mass density={self.mass_density_:.2f} | scattering potential={self.sc_potential_:.2f} eV | T={self.temp_:.1f} K')
            print(f'\t-- Interface rms roughness={self.rms_roughness_:.3f} nm | correlation length={self.corr_len_:.3f} nm')
            print(f'\t-- Dislocation density={self.n_dislocation_:.4f} nm^-2 | dislocation occupancy={self.f_dislocation_:.1f}')
            print(f'\t-- Electromechanical coupling coefficient={self.K_sqr:.5f} | deformation potential={self.E_d:.5f}')
            print(f'\t-- Longitudinal acoustic phonon velocity={self.v_LA:.2f} m/s | polar optical phonon energy={self.E_pop:.5f} eV')
            print(f'\t-- Fermi wave vector={self.k_F} | b={self.b_}')
            print('')

    def _form_factor(self, x, delta_2deg:bool=False, mode=None, numerator:bool=False):
        """
        This function calculates Fang-Howard form-factors.

        Parameters
        ----------
        x : float
            Scattering states. x=sing(theta/2), theta=scattering angle.
        delta_2deg : bool, optional
            If the 2DEG is perfect delta function 2DEG without any spread. For a 
            perfect 2DEG the FH form factors goes to 1. If true return 1. 
            The default is False.
        mode : string, optional ['IRF', 'DIS', 'DP', 'PE', 'POP']
            FW form factor which scattering mechanism. 
            The default is None. If None, None is returned.
        numerator : bool, optional
            Return Form factor for the numerator. The defalut is False.

        Returns
        -------
        float/None
            FW form factor.

        """
        # Fang-Howard form-factor
        if delta_2deg:
            return 1
        
        if mode in ['IRF','DIS','PE','DP']:
            # eta(u) = b/(b+2*k_f*u)
            eta = self.b_/(self.b_ + 2.0*self.k_F*x)
        else: #elif mode == 'POP':
            eta = self.b_/(self.b_ + self.k_pop)
        # F(eta) = eta^3    
        F_u_ = eta*eta*eta
        # G(eta) = (2*eta^3 + 3*eta^2 + 3*eta) / 8
        G_u_ = eta*(eta*(2*eta+3)+3)/8
        
        if mode == 'IRF':
            return G_u_
        elif mode == 'DIS' and self.mobility_model_ != 'v1':
            return F_u_**2 if numerator else G_u_
        elif mode in ['DP', 'PE']:
            return F_u_ if self.mobility_model_ == 'v1' else G_u_  
        elif mode == 'POP':
            return G_u_
        else:
            return 1
        return 1

    def _int_f_denomenator(self, x, mode=None):
        return ((x + self.q_TF_by_2k_F*self._form_factor(x, mode=mode))**2 * np.sqrt(1 - x**2))

    # ----- interface roughness ----------------
    def _inv_tau_ifr_f(self, x):
        return (x**4 * np.exp(-(self.corr_len_ * self.k_F * 0.1 * x)**2) / 
                self._int_f_denomenator(x, mode='IRF'))

    def _inv_tau_ifr_int(self):
        return integrate.quad(self._inv_tau_ifr_f, 0, 1)[0]

    def _inv_tau_ifr(self):
        # Since self.k_F propto sqrt(self.n_2d_); we must safe guard the scattering mechanism 
        # where division by either self.n_2d_ or self.k_F are done. 
        # In case of self.q_TF_by_2k_F, the division is done by k_F. This
        # term arrises in most scattering integration.
        # => this safe guard ensures we do not get any scattering in cases for
        # very low n_2d or no ally scarring for pure binary systems.
        if self.n_2d_ < self.eps_n_2d: return 0
        #*****************************************
        # (m0*e^4)/(8*h_bar^3*eps_0^2) * 1e-4 = 81.6046000430338 1e12 s^-1
        return 81.6046000430338 * self.m_star_by_eps_s_square \
                * (self.rms_roughness_ * self.corr_len_ * self.n_2d_)**2 \
                * self._inv_tau_ifr_int() # 1e12 s^-1
        
    # ----------- dislocation -------------------
    ## Scattering due to threading edge dislocation charge line 
    def _inv_tau_dis_f(self, x):
        return self._form_factor(x, mode='DIS', numerator=True) \
                /self._int_f_denomenator(x, mode='DIS')

    def _inv_tau_dis_int(self):
        return integrate.quad(self._inv_tau_dis_f, 0, 1)[0]

    def _inv_tau_dis(self):
        # Since self.k_F propto sqrt(self.n_2d_); we must safe guard the scattering mechanism 
        # where division by either self.n_2d_ or self.k_F are done. 
        # In case of self.q_TF_by_2k_F, the division is done by k_F. This
        # term arrises in most scattering integration.
        # => this safe guard ensures we do not get any scattering in cases for
        # very low n_2d or no ally scarring for pure binary systems.
        if self.n_2d_ < self.eps_n_2d: return 0
        #*****************************************
        # (m0*e^4)/(4*pi*h_bar^3*eps_0^2) * (1e8 / 1e6**4/ 1e-8**2) = 519511.0190323496 1e12 s^-1
        return 519511.0190323496 * self.m_star_by_eps_s_square \
                * self.n_dislocation_ * self.f_dislocation_**2  \
                * self._inv_tau_dis_int() / (self.k_F**4 * self.c_lp**2) # 1e12 s^-1
    
    ## Scattering due to strain field of threading edge dislocation line 
    def _inv_tau_dis_strain_f(self, x):
        return x**2*self._form_factor(x, mode='DIS', numerator=True)\
                /self._int_f_denomenator(x, mode='DIS')

    def _inv_tau_dis_strain_int(self):
        return integrate.quad(self._inv_tau_dis_strain_f, 0, 1)[0]
    
    def _inv_tau_dis_strain(self):
        # Since self.k_F propto sqrt(self.n_2d_); we must safe guard the scattering mechanism 
        # where division by either self.n_2d_ or self.k_F are done. 
        # In case of self.q_TF_by_2k_F, the division is done by k_F. This
        # term arrises in most scattering integration.
        # => this safe guard ensures we do not get any scattering in cases for
        # very low n_2d or no ally scarring for pure binary systems.
        if self.n_2d_ < self.eps_n_2d: return 0
        #*****************************************
        #(e_mass*e_charge**2)/(2*pi_*h_bar**3)*1e-4*1e-20 = 0.003173229123349822 1e12 s^-1
        # Burger's vector b_e = a_lp
        poisson_part = ((1-2*self.poisson_ratio)/(1-self.poisson_ratio))**2
        return 0.003173229123349822 * self.n_dislocation_ * self.m_star_ \
                * self.a_lp**2 * self.E_d**2 * poisson_part \
                * self._inv_tau_dis_strain_int() / self.k_F**2 # 1e12 s^-1

    # ----------- alloy disordered -------------------
    def _inv_tau_ado(self):
        # Since self.k_F propto sqrt(self.n_2d_); we must safe guard the scattering mechanism 
        # where division by either self.n_2d_ or self.k_F are done. 
        # In case of self.q_TF_by_2k_F, the division is done by k_F. This
        # term arrises in most scattering integration.
        # => this safe guard ensures we do not get any scattering in cases for
        # very low n_2d or no ally scarring for pure binary systems.
        if (self.comp_ < 1e-8) or (self.n_2d_ < self.eps_n_2d) or ((1-self.comp_)<1e-8): return 0
        #*****************************************
        #(3*e_mass*e_charge**2)/(16*h_bar**3)*1e6*1e-8**3*1e-4 = 0.37383724882773683 1e12 s^-1
        return 0.37383724882773683 * self.m_star_ * self.omega_0_ad * self.sc_potential_**2 \
                * self.comp_ * (1.0 - self.comp_) * self.b_ # 1e12 s^-1

    # ----------- deformation potential -------------------
    def _inv_tau_dp_f(self, x):
        return x**4/self._int_f_denomenator(x, mode='DP') 

    def _inv_tau_dp_int(self):
        return integrate.quad(self._inv_tau_dp_f, 0, 1)[0]
        
    def _inv_tau_dp(self):
        # Since self.k_F propto sqrt(self.n_2d_); we must safe guard the scattering mechanism 
        # where division by either self.n_2d_ or self.k_F are done. 
        # In case of self.q_TF_by_2k_F, the division is done by k_F. This
        # term arrises in most scattering integration.
        # => this safe guard ensures we do not get any scattering in cases for
        # very low n_2d or no ally scarring for pure binary systems.
        if self.n_2d_ < self.eps_n_2d: return 0
        #*****************************************
        # 3*(e_mass*e_charge**2*k_B)/(4*pi_*h_bar**3)*1e6*1e2 = 6571673.423885714 1e12 s^-1
        return 6571673.423885714 * (self.m_star_*self.E_d*self.E_d*self.temp_ \
                  *self.b_*self._inv_tau_dp_int()) / (self.mass_density_*self.v_LA*self.v_LA)  # 1e12 s^-1

    # ----------- piezoelectric -------------------
    def _inv_tau_pe_f(self, x):
        return x**3*self._form_factor(x, mode='PE')/self._int_f_denomenator(x, mode='PE')

    def _inv_tau_pe_int(self):
        return integrate.quad(self._inv_tau_pe_f, 0, 1)[0]
        
    def _inv_tau_pe(self):
        # Since self.k_F propto sqrt(self.n_2d_); we must safe guard the scattering mechanism 
        # where division by either self.n_2d_ or self.k_F are done. 
        # In case of self.q_TF_by_2k_F, the division is done by k_F. This
        # term arrises in most scattering integration.
        # => this safe guard ensures we do not get any scattering in cases for
        # very low n_2d or no ally scarring for pure binary systems.
        if self.n_2d_ < self.eps_n_2d: return 0
        #*****************************************
        # (e_mass*e_charge**2*k_B)/(pi_*eps_0*h_bar**3)*1e-2 * 1e-6 = 98.96143403667759 1e12 s^-1
        return 98.96143403667759 * (self.m_star_*self.K_sqr*self.temp_ \
                *self._inv_tau_pe_int())/(self.eps_s_*self.k_F) # 1e12 s^-1

    # ----------- polar optical phonon -------------------
    def _inv_tau_pop(self):
        # Since self.k_F propto sqrt(self.n_2d_); we must safe guard the scattering mechanism 
        # where division by either self.n_2d_ or self.k_F are done. 
        # In case of self.q_TF_by_2k_F, the division is done by k_F. This
        # term arrises in most scattering integration.
        # => this safe guard ensures we do not get any scattering in cases for
        # very low n_2d or no ally scarring for pure binary systems.
        if self.n_2d_ < self.eps_n_2d: return 0
        #*****************************************
        eps_star = 1/(1/self.eps_h_ - 1/self.eps_s_)
        fact_2 = np.sqrt(self.m_star_ * self.E_pop) * self._form_factor(None, mode='POP')/eps_star 
        # pi_*h_bar**2/(e_mass*k_B)*1e12 *1e4 = 27.77985128879875
        yy = 27.77985128879875*self.n_2d_/self.m_star_/self.temp_
        # e_charge/k_B = 11604.518121550082
        fact_3 = yy / ((np.exp(11604.518121550082*self.E_pop/self.temp_) - 1) * (1+yy-np.exp(-yy))) 
        # e_charge**2*np.sqrt(e_mass*e_charge)/(2*np.sqrt(2)*eps_0*h_bar**2) = 35210.68196468214 1e12 s^-1
        return 35210.68196468214 * fact_2 * fact_3 

    # Scattering rate to mobility calculation
    def _mobility_calculator(self, inverse_scattering):  
        """
        This function calculates the sheet mobility from different scattering contributions.
        """
        # 1e4 is unit conversion from m^2 to cm^2
        # self.m_star_by_e_ = 5.685630103565723 * self.m_star_ # 10^-12 V.m^-2.s^2
        # inverse_scattering is in 10^12 s^-1
        return 1e4/(self.m_star_by_e_ * inverse_scattering) if inverse_scattering else np.nan # cm^2 V^-1 S^-1
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def _calculate_figure_of_merit(self, n_2d, mobility,  
                                   temp:float=300,  mode:str='LFOM', 
                                   T_corect_bandgap:bool=False,
                                   direct_bandgap:bool=True, 
                                   indirect_bandgap:bool=False):
        """
        This function calculates the figure-of-merit (FOM). Available FOMs are
        LFOM: Lateral figure-of-merit
        
        Ref: J. L. Hudgins, G. S. Simin, E. Santi and M. A. Khan, 
        "An assessment of wide bandgap semiconductors for power devices," 
        in IEEE Transactions on Power Electronics, vol. 18, no. 3, pp. 907-914, 
        May 2003, doi: 10.1109/TPEL.2003.810840.
        
        direct_bandgap_critical_electric_field = 1.73e5*(bandgap_**2.5) # V/cm
        indirect_bandgap_critical_electric_field = 2.38e5*(bandgap_**2) # V/cm
        
        Units:
        bandgap_ => in eV.
        temp => in K
        n_2d => in 10^12 cm^-2
        E_cr => in V/cm
        e => 1.602176634e-19 C
        mobility (mu) => cm^2 V^-1 s^-1
        
        
        LFOM = e*n_2d*mu*E_cr^2 = 1.602e-19 C * 1e12 cm^-2 * cm^2 V^-1 s^-1 * V^2cm^-2
                                = 1.602e-7 CVs^-1cm^-2
                                = 1.602e-7 Wcm^-2    #1 watts = 1 coulombs*volt/second
                                = 1.602e-13 MW/cm^2
        
        Parameters
        ----------
        n_2d : 1D float array (unit: 10^12 cm^-2)
            Array containing carrier density data .
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
        assert mode in ['LFOM'], 'Requested mode is not implemented yet' 
        bandgap_ = self.alloy_params_.get('bandgap')
        bandgap_alpha_ = self.alloy_params_.get('bandgap_alpha')
        bandgap_beta_ = self.alloy_params_.get('bandgap_beta')
        if T_corect_bandgap:
            bandgap_ = self._apply_Varshni_T_correction_2_bandgap(bandgap_, temp=temp,
                                                                  bandgap_alpha=bandgap_alpha_,
                                                                  bandgap_beta=bandgap_beta_)
        if mode == 'LFOM':  
            if direct_bandgap:
                # 1.73**2*e_charge*1e10*1e12*1e-6 =  0.0047951544478985995
                e_times_e_cr_sqr_pre_fact = 0.0047951544478985995 
                e_cr_sqr_eg_pow = 5 # 2.5**2
            else: #elif indirect_bandgap:
                # 2.38**2*e_charge*1e10*1e12*1e-6 = 0.009075369325629598
                e_times_e_cr_sqr_pre_fact = 0.009075369325629598
                e_cr_sqr_eg_pow = 4 # 2**2
            return e_times_e_cr_sqr_pre_fact * n_2d * mobility \
                    * bandgap_**e_cr_sqr_eg_pow #MW/cm^2