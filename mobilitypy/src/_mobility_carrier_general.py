#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 15:41:33 2025

@author: badal.mondal
"""

import numpy as np
from ._alloy_params import _AlloyParams

## ============================================================================
class _MobilityCarrier(_AlloyParams):
    '''
    The functions in this class sets general parameters for the mobility of nD carrier gas.  
    '''
    
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], 
                 pseudomorphic_strain:bool=False, substrate:str|float=None, 
                 alloy_crystal_structure:str='wz', use_mat_params:dict=None, 
                 alloy_type:str=None, print_log=None, eps_n=1e-10):
        """
        Initiation function of the class _MobilityCarrier.
        
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
        pseudomorphic_strain : bool, optional
            Whether to consider pseudomorphic strain.
            The default is False.
        substrate : string or float (unit: Angstrom), optional
            The substrate name (if string, warning: the name should be in the database) 
            or the substrate in-plane lattice parameter (if float, Angstrom unit).
            The default is None. Error will be raised if substrate=None and 
            pseudomorphic_strain=True.
        alloy_crystal_structure :  str, optional [options: 'WZ', 'ZB', 'DM']
            The crystal type of the materials. This will be considered when calculating
            parameters like Poisson ratio etc.
            Use following abbreviation name:
                for wurtzite use 'WZ' or 'wz'.
                for zincblende use 'ZB' or 'zb'.
                for diamond use 'DM' or 'dm'.
            The default is 'wz'. 
        use_mat_params : dict, optional
            To use different materials parameters from that given in the database.
            Simply join the binary names to construct the alloy name. e.g.,
            for binaries=['AlN', 'GaN'] the alloy name is 'AlNGaN' or 'GaNAlN'.
            Material parameters units should be same as in the database.
            e.g. use_mat_params = {'AlN': {'mass_density': 3000}}
        alloy_type : string (case sensitive), optional [options: 'CatAni']
            The alloy type name. Case sensitive. Only needed if alloy is of AxB1-xCxD1-y
            kind. Will be ignored for alloy of type AxB1-x, AxByC1-x-y, AxByCzD1-x-y-z etc.
            The default is None. 
        print_log : string, optional => ['high','medium','low', None]
            Determines the level of log to be printed. The default is None.
        eps_n : float, optional (unit: nm^-2 for 2DEG or 1e18 cm^-2 for 3DG)
            Carrier density below eps_n will be considered as zero. 
            For 2DEG: The default is 1e-10 nm^-2 == 1e4 cm^-2.
            For 3DEG: The default is 1e-14 1e18 cm^-2 == 1e4 cm^-2.

        Returns
        -------
        None.

        """
        if pseudomorphic_strain and (substrate is None):
            # This allows to return Error in the very begining without starting any calculations.
            raise ValueError('substrate tag can not be None when pseudomorphic_strain=True.')
        
        self.print_info = print_log
        if self.print_info is not None: self.print_info = self.print_info.lower()

        self.eps_n = eps_n
        _AlloyParams.__init__(self, compositions=compositions, binaries=binaries, 
                              alloy_crystal_structure=alloy_crystal_structure,
                              alloy_type=alloy_type)
        self._get_alloy_params(use_mat_params=use_mat_params)
        if pseudomorphic_strain: self._cal_pseudomorphic_strain(substrate)
        return
            
    def _set_params_general(self, m_star, eps_s, eps_h, c_lattice, a_lattice, sc_potential, 
                            n_dis, f_dis, n_ion_impurity, mass_density, v_LA, E_pop, 
                            E_D, K_square, poisson_ratio, T):
        """
        This function sets the parameters for mobility calculations.
        """
        self.m_star_ = m_star
        self.eps_s_ = eps_s
        self.eps_h_ = eps_h
        self.c_lp = c_lattice
        self.a_lp = a_lattice
        self.sc_potential_ = sc_potential
        self.n_dislocation_ = n_dis
        self.f_dislocation_ = f_dis
        self.n_ion_imp =n_ion_impurity
        self.mass_density_ = mass_density
        self.v_LA = v_LA
        self.E_pop = E_pop
        self.K_sqr = K_square
        self.E_d = E_D
        self.poisson_ratio = poisson_ratio
        self.temp_ = T if T > 1e-8 else 1e-5 # Make sure zero divison does not happen when T=0 is choosen
        self.omega_0_ad = self._cal_omega_0_ad(self.a_lp, self.a_lp, self.c_lp ) 
        # # m0 / e = 5.685630103565723*10^-12 V.m^-2.s^2
        self.m_star_by_e_ = 5.685630103565723 * self.m_star_ # 10^-12 V.m^-2.s^2
            
    @staticmethod
    def _calculate_sheet_resitance(carrier_density, mobility):
        """
        This function calculates the sheet resistance.
        
        1 coulomb/volt = 1 second/ohm
        1 ohm = 1 C^-1.V.s

        R = 1/(e * carrier_density * mu) ohm/square
          = 1/(1.602176634e-19 *carrier_density * mu C.cm^-2.cm^2.V^-1.S^-1) 
          = 6241509.0744607635/(carrier_density * mu) ohm/square

        Parameters
        ----------
        carrier_density : float/ndarray (unit: 10^12 cm^-2)
            Array containing carrier density data. 
        mobility : float/ndarray (unit: cm^2 V^-1 s^-1)
            Array containing mobility data.

        Returns
        -------
        float/ndarray (unit: ohm/square)
            Sheet resistance.

        """
        return 6241509.0744607635/(carrier_density * mobility)
    
    @staticmethod
    def _apply_Varshni_T_correction_2_bandgap(bandgap_0, temp:float=300, 
                                              bandgap_alpha:float=0, bandgap_beta:float=0):
        """
        This functions applies Varshni's formula for temperature correction to band gap.
        Eg(T) = Eg(T=0) - [aT^2/(T+b)]

        Parameters
        ----------
        bandgap_0 : 1D float array (unit: eV)
            Band gap values at 0K temperature.
        temp : float, optional (unit: K)
            Temperature in K. The default is 300K.
        bandgap_alpha : float, optional (unit: eV/K)
            Temperature correction coefficient alpha. The default is 0.
        bandgap_beta : float, optional (unit: K)
            Temperature correction coefficient beta. The default is 0.

        Returns
        -------
        1D float array (unit: eV)
            The temperature corrected band gap values.

        """
        return bandgap_0 - (bandgap_alpha*temp*temp/(temp+bandgap_beta))