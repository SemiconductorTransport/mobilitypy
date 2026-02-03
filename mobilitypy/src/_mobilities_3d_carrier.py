#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 16:15:33 2025

@author: badal.mondal
"""
import numpy as np
import pandas as pd
from ._Fermi_Dirac_integration import _FermiDiracInt
## ============================================================================
class _Mobility3DCarrier:
    '''
    The functions in this class calculates the mobility of 3D carrier gas.  
    The mobility models are implemented based on the following references.
    
    Note: Some of the equations in the references has prining mistakes. The mistakes
    are corrected in our implementation. 
    
    Ref-1: Rajan et al., Appl. Phys. Lett. 88, 042103 (2006) => alloy disorder, polar optical phonon
    Ref-2: DJ. and UKM., PRB 66, 241307(R) (2002) and DJ. et al., PRB 67, 153306 (2003)  => Dislocation
    Ref-3: Debdeep Jena's thesis, Chapter-6 APPENDIX, Sec. Three-dimensional carriers => Acoustic phonon 
    
    '''
    
    def __init__(self):
        """
        Initiation function of the class _Mobility3DCarrier.
        
        Returns
        -------
        None.

        """
        self.eps_n_3d = self.eps_n
        
    def _calculate_3d_mobility(self, n_3d=1, n_dis:float=1, f_dis:float=0.5, T:float=300):
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

        Returns
        -------
        pandas dataframe of compositions and mobilities (unit: cm^2 V^-1 S^-1).
            Total (or individual contributions) sheet mobility. 

        """
        self._set_params_general(self.alloy_params_.get('e_effective_mass'), 
                                 self.alloy_params_.get('static_dielectric_constant'),
                                 self.alloy_params_.get('high_frequency_dielectric_constant'), 
                                 self.alloy_params_.get('lattice_c0'), 
                                 self.alloy_params_.get('lattice_a0'), 
                                 self.alloy_params_.get('alloy_scattering_potential'), 
                                 n_dis, f_dis, 
                                 self.alloy_params_.get('mass_density'), 
                                 self.alloy_params_.get('LA_phonon_velocity'), 
                                 self.alloy_params_.get('PO_phonon_energy'), 
                                 self.alloy_params_.get('CB_deformation_potential'),
                                 self.alloy_params_.get('electromechanical_coupling_const'),
                                 T
                                 )
        # Remove small values for the n_3d to avoid 0-division
        self.n_3d_ = np.where(n_3d < self.eps_n_3d, np.nan, n_3d)        

        mobility = {}
        total_inv_mu = 0
        if self.only_total_mobility:
            if self.print_info is not None: print('\t-- Calculating only total mobility')
            if self.alloy_disordered_effect_: total_inv_mu += 1/self._alloy_disorder_mu()   
            if self.polar_optical_phonon_effect_: total_inv_mu += 1/self._pop_mu()
            if self.acoustic_phonon_effect_: total_inv_mu += 1/self._ac_dp_mu() # due to deformation_potential_effect
            if self.piezoelectric_effect_: total_inv_mu += self._mu_pz()
            if self.dislocation_effect_: total_inv_mu += 1/self._dis_mu()
            mobility['TOT'] = 1/total_inv_mu
        else:
            if self.alloy_disordered_effect_:
                if self.print_info is not None: print('\t-- Calculating alloy-disordered mobility')
                _mu_contrib = self._alloy_disorder_mu()
                total_inv_mu += 1/_mu_contrib
                mobility['AD'] = _mu_contrib
                    
            if self.polar_optical_phonon_effect_:
                if self.print_info is not None: print('\t-- Calculating polar optical phonon effect mobility')
                _mu_contrib = self._pop_mu()
                total_inv_mu += 1/_mu_contrib
                mobility['POP'] = _mu_contrib
 
            if self.acoustic_phonon_effect_:
                if self.print_info is not None: print('\t-- Calculating acoustic phonon effect mobility')
                _mu_contrib = self._ac_dp_mu()
                total_inv_mu += 1/_mu_contrib
                mobility['AP'] = _mu_contrib
                
            if self.piezoelectric_effect_:
                if self.print_info is not None: print('\t--- Calculating piezoelectric phonon effect mobility')
                _mu_contrib = self._mu_pz()
                total_inv_mu += 1/_mu_contrib
                mobility['PE'] = _mu_contrib
                                
            if self.dislocation_effect_:
                if self.print_info is not None: print('\t-- Calculating dislocation effect mobility')
                _mu_contrib = self._dis_mu()
                total_inv_mu += 1/_mu_contrib
                mobility['DIS'] = _mu_contrib
 
            if self.total_mobility_:
                if self.print_info is not None: print('\t-- Calculating total mobility')
                mobility['TOT'] = 1/total_inv_mu
                
        if self.print_info is not None: print(f'{"="*72}')
        return pd.DataFrame.from_dict(mobility, orient='index')
    
    def _ln_1p_exp_xi(self):
        """
        Zeroth order FD integral.
        """
        return _FermiDiracInt._cal_Fermi_Dirac_integral(self.n_3d_, self.m_star_, 
                                                        T = self.temp_, 
                                                        inv_half_FD_method = 
                                                        self.inverse_half_FD_method_,
                                                        FD_order = 'zero')
    ## Alloy disordered limited mobility
    def _alloy_disorder_mu(self, eps_den = 1e-8):        
        demoninator_ = self.m_star_*self.sc_potential_*self.sc_potential_*self.omega \
                       *self.n_3d*self.comp_*(1.0-self.comp_)          
        # Remove small values for both the n_3d and comp_ or (1-comp_)
        demoninator_ = np.where(demoninator_ < eps_den, np.nan, demoninator_)
        #fact_ad_3d = 21.16990563011839 # (2*e*h_bar*k_B)/(3*pi*m0*e^2) * 1e10 => cm^-2.K^-1.V-1.s-1
        return 21.16990563011839 * self.temp_ * self._ln_1p_exp_xi() / demoninator_ # cm^2.V^-1.s^-1
    
    ## Polar optical phonon limited mobility
    def _pop_mu(self):
        #4.0*pi*eps_0*h_bar*h_bar/(sqrt(2)*(e*m0)**(3/2))*1e4 = 0.1569266969277379 # cm^2V^-1s^-1
        eps_star = 1/((1/self.eps_h_) - (1/self.eps_s_))
        # e/k_B = 11604.518121550082
        exp_fact = np.exp(self.E_pop/self.temp_*11604.518121550082) - 1.0
        return 0.1569266969277379 * eps_star * exp_fact / self.m_star_ / np.sqrt(self.m_star_*self.E_pop)
    
    ## Deformation potential acoustic phonon limited mobility
    def _ac_dp_mu(self):
        #2*h_bar*e*1e-2/(3*pi*m0*e*e*1e18) = 1.53333002306295e-06 # cm^2V^-1s^-1
        numerator = self.mass_density_ * self.v_LA * self.v_LA * self._ln_1p_exp_xi() * 1.53333002306295e-06
        return numerator / (self.n_3d*self.m_star_*self.E_d*self.E_d)
    
    ## Piezoelectric phonon scattering limited mobility
    def _mu_pz(self):
        #16*k_B*eps_0/(3*pi*h_bar*e*1e18*1e2) = 0.12282713258060055 # cm^2V^-1s^-1K^-1
        FD_1 = _FermiDiracInt._cal_Fermi_Dirac_integral(self.n_3d_, self.m_star_, 
                                                        T = self.temp_, 
                                                        inv_half_FD_method = 
                                                        self.inverse_half_FD_method_,
                                                        FD_order = 'one',
                                                        use_numerical_integration =
                                                        self.use_numerical_FD_integration_)
        return self.temp_*self.eps_s_*FD_1*0.12282713258060055/(self.n_3d*self.K_sqr)
    
    ## Dislocation limited mobility
    def _dis_mu(self):
        #(3*pi_*pi_)**(1/3) = 3.0936677262801355e6 cm^-1
        fermi_wave_vect = 3.0936677262801355*(self.n_3d**(1/3)) # 1e6 cm^-1
        fermi_wave_vect_sq = fermi_wave_vect*fermi_wave_vect # 1e12 cm^-2
        # h_bar*h_bar/2/e_mass*1e12 = 6.10426432246119e-27 J^2s^2cm^-2kg^-1
        fermi_energy = 6.10426432246119*fermi_wave_vect_sq/self.m_star_ # 1e-27 J^2s^2cm^-2kg^-1
        # 2*eps_0/3/e_charge/e_charge*1e-43 = 2.2995173111302578e-17 cm^2
        TF_screening_len_sq = 2.2995173111302578*self.eps_s_*fermi_energy/self.n_3d # 1e-17 cm^2
        # k_f^2*lamda_TF^2 => 1e-5
        fact_12 = (1+4*fermi_wave_vect_sq*TF_screening_len_sq*1e-5)**(3/2)/TF_screening_len_sq**2
        # h_bar*h_bar*h_bar*eps_0*eps_0/e_charge**3/e_mass**2*1e22/1e10 = 26941.18974076079
        fact_11 = (self.eps_s_*self.c_lp/self.f_dislocation_/self.m_star_)**2 / self.n_dislocation_
        return fact_11 * fact_12 * 26941.18974076079