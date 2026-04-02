#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 16:15:33 2025

@author: badal.mondal
"""
import numpy as np
import pandas as pd
from ._Fermi_Dirac_integration import _FermiDiracInt
import scipy.integrate as integrate

## ============================================================================
class _Mobility3DCarrier:
    '''
    The functions in this class calculates the mobility of 3D carrier gas.  
    The mobility models are implemented based on the following references.
    
    Note: Some of the equations in the references has prining mistakes. The mistakes
    are corrected in our implementation. 
    
    '''
    
    def __init__(self):
        """
        Initiation function of the class _Mobility3DCarrier.
        
        Returns
        -------
        None.

        """
        self.eps_n_3d = self.eps_n
        
    def _calculate_3d_mobility(self, n_3d=1, n_dis:float=1, f_dis:float=0.5, 
                               n_ion_impurity:float=1, T:float=300):
        """
        This function calculates the sheet mobility from different scattering contributions.
        The mobility models are implemented based on the following references.
            
        The considered scattering mechanism are:
            Alloy disorder limited (AD)
            Threading dislocation mediated (DIS_TD)
            Piezoelectric phonon effect (PE)
            Acoustic deformation potential phonon (DP)
            Polar optical phonon (POP)
            Ionized impurity (ION_IMP)

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
        n_ion_impurity : float, optional (unit: 1e14 cm^-3)
            Ionized impurity density. The default is 1.
        T : float, optional (unit: K)
            Temperature at which mobility calculations will be done. 
            The default is 300K.

        Returns
        -------
        pandas dataframe of compositions and mobilities (unit: cm^2 V^-1 S^-1).
            Total (or individual contributions) sheet mobility. 

        """      
        #======================================================================
        self._set_params_general(self.alloy_params_.get('e_effective_mass'), 
                                 self.alloy_params_.get('static_dielectric_constant'),
                                 self.alloy_params_.get('high_frequency_dielectric_constant'), 
                                 self.alloy_params_.get('lattice_c0'), 
                                 self.alloy_params_.get('lattice_a0'), 
                                 self.alloy_params_.get('alloy_scattering_potential'), 
                                 n_dis, f_dis, n_ion_impurity,
                                 self.alloy_params_.get('mass_density'), 
                                 self.alloy_params_.get('LA_phonon_velocity'), 
                                 self.alloy_params_.get('PO_phonon_energy'), 
                                 self.alloy_params_.get('CB_deformation_potential'),
                                 self.alloy_params_.get('electromechanical_coupling_const'),
                                 self.alloy_params_.get('isotropic_Poisson_ratio'),
                                 T)
        #======================================================================
        # Remove small values for the n_3d to avoid 0-division
        self.n_3d_ = np.nan if (np.isscalar(n_3d) and n_3d < self.eps_n_3d) else\
            np.where(n_3d < self.eps_n_3d, np.nan, n_3d)    
        
        # Fermi eta = E_f/(k_B.T)
        self.eta_f_ = _FermiDiracInt._cal_eta_from_inv_FD(self.n_3d_, self.m_star_, T=self.temp_, 
                                                          method=self.inverse_half_FD_method_)
        #======================================================================    
        if self.td_dislocation_chg_effect_ or self.td_dislocation_strain_effect_:
            self._dis_facts()
        if self.alloy_disordered_effect_ or self.acoustic_phonon_effect_:
            self._ln_1p_exp_xi() # Calculates the FD oth order integral
        #======================================================================
        mobility = {}
        if self.alloy_disordered_effect_:
            if self.print_info is not None: print('\t-- Calculating alloy-disordered mobility')
            mobility['mu_AD'] = self._alloy_disorder_mu()
                
        if self.polar_optical_phonon_effect_:
            if self.print_info is not None: print('\t-- Calculating polar optical phonon effect mobility')
            # POP scattering does not depend on n_3d. For single comp and n_3d
            # array the return array shape would not match with other scattering 
            # mechanisms. This is to safe guard.
            if (not np.isscalar(self.n_3d_)) and (len(self.n_3d_) != len(self.comps_)):
                mobility['mu_POP'] = np.repeat(self._pop_mu(), len(self.n_3d_))
            else:
                mobility['mu_POP'] = self._pop_mu()

        if self.acoustic_phonon_effect_:
            if self.print_info is not None: print('\t-- Calculating acoustic phonon deformation potential effect mobility')
            mobility['mu_DP'] = self._ac_dp_mu()
            
        if self.piezoelectric_effect_:
            if self.print_info is not None: print('\t--- Calculating piezoelectric phonon effect mobility')
            mobility['mu_PE'] = self._mu_pz()
            
        if self.ionized_impurity_effect_: 
             if self.print_info is not None: print('\t-- Calculating ionized impurity limited mobility')
             mobility['mu_ION_IMP'] = self._ion_imp_mu()

        if self.td_dislocation_chg_effect_:
            if self.print_info is not None: print('\t-- Calculating charge line dislocation effect mobility')
            mobility['mu_DIS_TD_CHG'] = self._td_chg_dis_mu()
            
        if self.td_dislocation_strain_effect_: 
            if self.print_info is not None: print('\t-- Calculating dislocation strain field effect mobility')
            mobility['mu_DIS_TD_STR'] = self._td_str_dis_mu()
        #======================================================================
        MuDataframe = pd.DataFrame.from_dict(mobility)
        #======================================================================
        if self.total_mobility_:
            if self.print_info is not None: print('\t-- Calculating total mobility')
            #print(list(MuDataframe.keys()))
            MuDataframe['mu_TOT'] = 1/((1/MuDataframe).sum(axis=1, skipna=True))
        #======================================================================    
        if self.print_info is not None: print(f'{"="*72}')
        #======================================================================
        if self.only_total_mobility:
            return MuDataframe['mu_TOT']
        else:
            if self.td_dislocation_chg_effect_ and self.td_dislocation_strain_effect_:
                # Postprocessing: total DIS
                XX = ['mu_DIS_TD_CHG', 'mu_DIS_TD_STR']
                MuDataframe['mu_DIS_TD'] = 1/((1/MuDataframe[XX]).sum(axis=1, skipna=True))
                
            return MuDataframe
        #======================================================================
    
    def _ln_1p_exp_xi(self):
        """
        Zeroth order FD integral.
        """
        self.F0_eta = _FermiDiracInt._cal_Fermi_Dirac_integral(self.eta_f_, FD_order = 'zero')
        return 
    
    ## Alloy disordered limited mobility
    def _alloy_disorder_mu(self, eps_den = 1e-8):        
        demoninator_ = self.m_star_*self.sc_potential_*self.sc_potential_*self.omega_0_ad \
                       *self.n_3d_*self.comps_*(1.0-self.comps_)          
        # Remove small values for both the n_3d and comps_ or (1-comps_)
        demoninator_ = np.where(demoninator_ < eps_den, np.nan, demoninator_)
        #(2*e_charge*h_bar*k_B)/(3*pi_*e_mass*e_charge**2) * 1e10 = 21.16990563011839
        return 21.16990563011839 * self.temp_ * self.F0_eta / demoninator_ # cm^2.V^-1.s^-1
    
    ## Polar optical phonon limited mobility
    def _pop_mu(self):
        eps_star = (self.eps_h_*self.eps_s_)/(self.eps_s_-self.eps_h_)
        # e_charge/k_B = 11604.518121550082
        exp_fact = np.exp(self.E_pop/self.temp_*11604.518121550082) - 1.0
        #2.0*np.sqrt(2)*pi_*eps_0*h_bar**2/e_charge/np.sqrt(e_mass**3*e_charge)*1e4 = 0.1569266969277379 # cm^2V^-1s^-1
        return 0.1569266969277379 * eps_star * exp_fact / self.m_star_ / np.sqrt(self.m_star_*self.E_pop)
    
    ## Deformation potential acoustic phonon limited mobility
    def _ac_dp_mu(self):
        # 2*e_charge*h_bar/(3*pi_*e_mass*e_charge**2*1e20) = 1.53333002306295e-06 # cm^2V^-1s^-1
        numerator = self.mass_density_ * self.v_LA * self.v_LA * self.F0_eta * 1.53333002306295e-06
        return numerator / (self.n_3d_*self.m_star_*self.E_d*self.E_d)
    
    ## Piezoelectric phonon scattering limited mobility
    def _mu_pz(self):
        # 24*eps_0*e_mass*k_B**2/(h_bar**2*e_charge**2)*1e-24 = 0.00012925353328704564
        #xi_0 = 0.00012925353328704564*self.eps_s_*self.m_star_*self.temp_**2/self.n_3d_
        C_K0 = 1.0 #+(1/(1+xi_0))-2/xi_0*np.log(1.0+xi_0)
        
        if self.carrier_degenracy_limit_ == 'nondegenerate':
            # 16*np.sqrt(2*pi_)/3*h_bar**2*eps_0/(e_charge*e_mass**(3/2)*k_B**(1/2))*1e4 = 25.43338614569858
            return 25.43338614569858*self.eps_s_/(np.sqrt(self.m_star_**3*self.temp_)*self.K_sqr) #cm^2V^-1s^-1
        else:
            # 16*k_B*eps_0/(3*pi_*h_bar*e_charge*1e18*1e2) = 0.12282713258060055 # cm^2V^-1s^-1K^-1
            FD_1 = _FermiDiracInt._cal_Fermi_Dirac_integral(self.eta_f_, FD_order = 'one',
                                                            FD_int_approach=self.FD_int_approach_)
            return self.temp_*self.eps_s_*FD_1*0.12282713258060055/(self.n_3d_*self.K_sqr*C_K0) #cm^2V^-1s^-1
    
    ## Dislocation limited mobility
    def _dis_facts(self):
        # Only minimax_piecewise method is implemented for Fermi Diract 1/2, -1/2 integral.
        F_m_1h = _FermiDiracInt._cal_Fermi_Dirac_integral(self.eta_f_, FD_order = 'm_one_half',
                                                          FD_int_approach='minimax_piecewise')
        F_1h = _FermiDiracInt._cal_Fermi_Dirac_integral(self.eta_f_, FD_order = 'one_half',
                                                          FD_int_approach='minimax_piecewise')    
        self.F1hRatio = F_m_1h/F_1h
        # h_bar**2*e_charge**2/(8*e_mass*eps_0*k_B**2)*1e24 = 23210.19722793661
        self.dis_B_fact = 23210.19722793661*self.n_3d_*self.F1hRatio/(self.m_star_*self.eps_s_*self.temp_*self.temp_)
        return
    
    def _td_chg_dis_mu(self): 
        if self.carrier_degenracy_limit_ == 'degenerate':
            #4k_F^2 lambda^2 = 4*3**(1/3)*h_bar**2*pi_**(8/3)*eps_0/e_charge**2/e_mass*1e8 = 0.051430964880517044
            fact_12 = (1+0.051430964880517044*self.n_3d_**(1/3)*self.eps_s_/self.m_star_)**(3/2)
            # e_charge*3**(2/3)/(h_bar*pi_**(8/3))*1e-12 = 149.27327984905628 cm^2V^-1s^-1
            return 149.27327984905628 * self.c_lp * self.n_3d_**(2/3) * fact_12 \
                    / self.n_dislocation_ / self.f_dislocation_**2  
        elif self.carrier_degenracy_limit_ == 'nondegenerate':
            # 128*np.sqrt(2)/np.sqrt(pi_)*(eps_0**(3/2)*k_B/np.sqrt(e_mass)/e_charge**2*1e-16)=0.15163209085564022
            # 16*np.sqrt(2)*(eps_0**(3/2)*k_B/np.sqrt(e_mass)/e_charge**2*1e-16) = 0.03359511041974182
            return 0.15163209085564022 * self.c_lp*self.c_lp*self.temp_\
                    *np.sqrt(self.eps_s_*self.eps_s_*self.eps_s_*self.n_3d_/self.m_star_)\
                        / self.n_dislocation_ / self.f_dislocation_**2
        else:   
            I_eta = _FermiDiracInt._FD_dis_chg_Integration(self.eta_f_, self.dis_B_fact)            
            #8*np.sqrt(2)/pi_**2 * np.sqrt(e_mass*k_B**3)*eps_0/(e_charge*h_bar**2)*1e-28=0.027890781017309893
            return 0.027890781017309893*self.c_lp*self.c_lp*self.eps_s_*np.sqrt(self.m_star_*self.temp_**3)\
                        /(self.n_dislocation_*self.f_dislocation_**2)*self.F1hRatio*I_eta
                        
    def _td_str_dis_mu(self): 
        I_eta = _FermiDiracInt._FD_dis_str_Integration(self.eta_f_, self.dis_B_fact)     
        poisson_part = (1.0-self.poisson_ratio)/(1.0-2.0*self.poisson_ratio)
        #e_charge*np.sqrt(2*k_B/e_mass)/(3*pi_*pi_*eps_0)*1e12= 3364750.021017146
        return 3364750.021017146*np.sqrt(self.temp_/self.m_star_)*poisson_part*poisson_part\
                /(self.eps_s_*self.n_dislocation_*self.a_lp**2*self.E_d**2)*self.F1hRatio*I_eta 
                
    def _ion_imp_mu(self):
        # 24*eps_0*e_mass*k_B**2/(h_bar**2*e_charge**2)*1e-24 = 0.00012925353328704564
        xi_0 = 0.00012925353328704564*self.eps_s_*self.m_star_*self.temp_**2/self.n_3d_
        C_K0 = 1.0/(np.log(1.0+xi_0) - (xi_0/(1.0+xi_0)))
        #print(1+(1/(1+xi_0))-2/xi_0*np.log(1.0+xi_0))
        
        FD_2 = _FermiDiracInt._cal_Fermi_Dirac_integral(self.eta_f_, FD_order = 'two',
                                                        FD_int_approach=self.FD_int_approach_)
        
        # 128*e_mass*eps_0**2*k_B**3/(h_bar**3*e_charge**3)*1e-40 = 0.49875425045367205
        return 0.49875425045367205*self.m_star_*self.eps_s_*self.eps_s_*self.temp_**3*FD_2\
                /(self.n_ion_imp*self.n_3d_*C_K0)
        
    @staticmethod
    def _cal_elec_props_from_3DEC(n_3d, eps_s, m_star, pop_en, T, 
                                  inv_half_FD_method:str='minimax_piecewise',
                                  return_dis_ints:bool=False):
        """
    
        Parameters
        ----------
        n_3d : float or 1d array of float (unit: 1E18 cm^-3 )
            Volumetric carrier density.
        eps_s : float or 1d array of float (unit: epsilon_0)
            Static dielectic constants of the material. 
        m_star : float or 1d array of float (unit: m0)
            Carrier effective mass. 
        pop_en : float or 1d array of float (unit: eV)
            Polar optical phonon energy.
        T : float (unit: K)
            Temperature at which Fermi-Dirac integral calculations will be done. 
        inv_half_FD_method : str, optional [available: 'JD_approx', 'minimax_piecewise']
            The approximate method to calculate the scaled Fermi energy (E_f/k_BT) 
            using inverse Fermi-Dirac integral of order-1/2. The default is JD_approx.
            JD_approx : Joyce-Dixon approximation (APL 31, 354 (1977)).
            minimax_piecewise : minimax approximation (Applied Mathematics and 
                                                       Computation 259, 698 (2015))  
        return_dis_ints : bool, optional 
            Calculatd the n_3d  dependent integrals for dislocation related mobility.
            The default is False.
            
       Returns : tuple of lists/scalar 
       -------      
       If return_dis_ints=True returns the integrals only; (I_eta_chg, I_eta_str).
       Otherwise return (Scaled_Fermi_energy, Fermi_energy, Screening_wave_vector, 
       tau_c_by_tau_q_dis, Fermi_wave_vector, _pop_wave_vector)
       
       Scaled_Fermi_energy : list of float or 1d array of float list (unit: unitless)
           Fermi energy w.r.t conduction band w.r.t k_BT.
           general case: inverse Fermi-Dirac integral approach.
           degenerate case: assumes metallic ('degenerate') carriers.
           return : [general case, degenerate case]
       Fermi_energy : list of float or 1d array of float list (unit: eV)
           Fermi energy w.r.t conduction band. 
           general case: inverse Fermi-Dirac integral approach.
           degenerate case: assumes metallic ('degenerate') carriers.
           return : [general case, degenerate case]
       Screening_wave_vector : list of float or 1d array of float list (unit: 10^6 cm^-1)
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
        ## ========================= General ==================================
        scaled_Ef = _FermiDiracInt._cal_eta_from_inv_FD(n_3d, m_star, T=T, 
                                                        method=inv_half_FD_method)

        # Only minimax_piecewise method is implemented for Fermi Diract -1/2 integral.
        F_m_1h = _FermiDiracInt._cal_Fermi_Dirac_integral(scaled_Ef, FD_order = 'm_one_half',
                                                          FD_int_approach='minimax_piecewise')
        ## Integrals drom dislocation limited mobility
        if return_dis_ints:
            # Only minimax_piecewise method is implemented for Fermi Diract 1/2, -1/2 integral.
            F_1h = _FermiDiracInt._cal_Fermi_Dirac_integral(scaled_Ef, FD_order = 'one_half',
                                                            FD_int_approach='minimax_piecewise')    
            F1hRatio = F_m_1h/F_1h
            # h_bar**2*e_charge**2/(8*e_mass*eps_0*k_B**2)*1e24 = 23210.19722793661
            dis_B_fact = 23210.19722793661*n_3d*F1hRatio/(m_star*eps_s*T*T)
            I_eta_chg = _FermiDiracInt._FD_dis_chg_Integration(scaled_Ef,dis_B_fact)  
            I_eta_str = _FermiDiracInt._FD_dis_str_Integration(scaled_Ef,dis_B_fact)
            return (I_eta_chg, I_eta_str)
        
        Fermi_energy_Gen = 8.617333262145179e-05 * T * scaled_Ef.copy() # k_B J.K^-1 = k_B/e_charge eV.K^-1
        #np.sqrt(e_charge*e_charge/eps_0/k_B)*1e4 = 144.9086756309392
        #screening_wavevector_Gen = 144.9086756309392 * np.sqrt(self.n_3d_*F1hRatio/(self.eps_s_*self.temp_)) # 1e6 cm^-1
        #np.sqrt(e_charge**2*e_mass**(3/2)*k_B**(1/2)/(eps_0*(2*pi_*pi_*pi_)**
        # (1/2)*h_bar**3))*1e-2 = 10.070231429501227 # 1e6 cm^-1
        screening_wavevector_Gen = 10.070231429501227 * np.sqrt(np.sqrt(m_star**3*T)*F_m_1h/eps_s) # 1e6 cm^-1
        #
        
        ## =============== Non degenerate limit ===============================
        #np.sqrt(e_charge**2/eps_0/k_B)*1e10 = 144.9086756309392 # 1e6 cm^-1
        Debye_wave_vector = 144.9086756309392 * np.sqrt(n_3d/eps_s/T) # 1e6 cm^-1
        # 4*eps_0*e_mass*k_B*k_B/(h_bar**2*e_charge**2)*1e-24 = 2.1542255547840943e-05
        tau_c_by_tau_q_dis_ND = 1 + 2.1542255547840943e-05*eps_s*m_star*T*T/n_3d
         
        ## =================== Degenerate limit ===============================
        # (3*pi_*pi_)**(1/3) = 3.0936677262801355
        Fermi_wave_vector = 3.0936677262801355 * n_3d**(1/3) # 1e6 cm^-1
        # h_bar**2*1e4/2/e_mass/e_charge*1e12 = 0.0003809982110968585
        Fermi_energy_D = 0.0003809982110968585 * Fermi_wave_vector**2 / m_star # eV
        scaled_Fermi_Energy_D = Fermi_energy_D.copy() * 11604.518121550082 / T # eta_f = E_f * e / (k_B.T)
        # 10*np.sqrt(h_bar**2*eps_0*pi_**(4/3)/(3**(1/3)*e_charge**2*e_mass*1e6)) = 3.665292799274651e-08
        _Thomas_Fermi_screening_len = 0.03665292799274651 * np.sqrt(eps_s/m_star/(n_3d**(1/3))) # 1e-6 cm
        #2*eps_0*h_bar**2*pi_**(8/3)*3**(1/3)/e_charge/e_charge/e_mass*1e8 = 0.0257154824402585
        tau_c_by_tau_q_dis_D = 1 + 0.0257154824402585*eps_s*n_3d**(1/3)/m_star
        # np.sqrt(2*e_mass*e_charge/h_bar**2)*1e-2 = 51.23167223161843 # 1e6 cm^-1
        _pop_wave_vector = 51.23167223161843 * np.sqrt(m_star*pop_en)
        return ([scaled_Ef, scaled_Fermi_Energy_D], [Fermi_energy_Gen, Fermi_energy_D ], 
                [screening_wavevector_Gen, 1/_Thomas_Fermi_screening_len, Debye_wave_vector], 
                [tau_c_by_tau_q_dis_D, tau_c_by_tau_q_dis_ND], 
                Fermi_wave_vector, _pop_wave_vector)
    
    @staticmethod
    def _3dec_props(n_d, mu_d, position, eps_n_3d=1e-14, log_info=None):
        """
        This function calculates the effective/average properies of a 3D carrier distribution.

        Parameters
        ----------
        n_d : 1d numpy array of float (unit: 1E18 cm^-3 )
            The position dependent carrier density distribution.
        mu_d : 1d numpy array of float (unit: cm^2.V^-1.s^-1 )
            The position dependent carrier mobility distribution..
        position : 1d numpy array of float (unit: nm)
            The position array.
        eps_n_3d : float, optional (unit: 1e18 cm^-2)
            Carrier density below eps_n_3d will be considered as zero. 
            The default is 1e-14 1e18 cm^-2 == 1e4 cm^-2.
        log_info : string, optional [options: 'high','medium','low', None]
            Determines the level of log to be printed. The default is None.

        Returns
        -------
        IntegratedEdensity: float (unit: 1E13 cm^-2)
            Integrated carrier density.
        average_mu : tuple of two floats (unit: cm^2.V^-1.s^-1)
            First number is the effective/average mobility calculated using first 
            moment of carrier density w.r.t 1/mu. Second one is the average mobility 
            calculated using first moment of carrier density w.r.t mu.
        SheetResistance : float (unit: Ohm/sq)
            Effective sheet resistance.

        """
        n_d_ = 0 if (np.isscalar(n_d) and n_d < eps_n_3d) else np.where(n_d < eps_n_3d, 0, n_d)  
        IntegratedEdensity = integrate.trapezoid(n_d_, x=position) # 1e11 cm^-2
        
        ## Procedure-1: using inverse mobility
        density_mobility_ratio = n_d_/mu_d
        density_mobility_ratio[np.isnan(density_mobility_ratio)] = 0
        
        mobility_first_moment_nominator = integrate.trapezoid(density_mobility_ratio, x=position) 
        
        
        average_mu = IntegratedEdensity/mobility_first_moment_nominator # == 1.0/(mobility_first_moment_nominator/IntegratedEdensity)
        
        # Procedure-2: using mobility
        density_mobility_ = n_d_*mu_d
        density_mobility_[np.isnan(density_mobility_)] = 0
        density_weighted_mobility_ = integrate.trapezoid(density_mobility_, x=position)
        average_mu_ = density_weighted_mobility_/IntegratedEdensity
        
        # 1/(1e11*e_charge) = 62415090.744607635
        SheetResistance = 62415090.744607635/(average_mu*IntegratedEdensity) # Ohm/sq
       
        if log_info is not None:
           print(f'\to ave_es = {IntegratedEdensity*1e-2:0.2f} x 1E13 cm^-2')
           print(f'\to ave_mu = {average_mu:.2f} ({average_mu_:.2f}) cm^2.V^-1.s^-1')
           print(f'\to Rsh = {SheetResistance:.2f} Ohm/sq')
        return IntegratedEdensity*1e-2, (average_mu, average_mu_), SheetResistance