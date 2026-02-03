#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 11:50:49 2026

@author: badal.mondal
"""

import scipy.integrate as integrate
import scipy.special as special
import numpy as np
## ============================================================================
        
class _FermiDiracInt:   
    def __init__(self):
        pass
    @classmethod
    def _Fukushima_iFD_half(cls, nu: float) -> float: 
        '''
        Taken from github.com/scott-maddox/fdint
        
        Double-precision minimax approximation to the inverse Fermi-Dirac
        integral of order k=1/2.
        
        T. Fukushima, "Precise and fast computation of inverse Fermi-Dirac
        integral of order 1/2 by minimax rational function approximation,"
        Applied Mathematics and Computation, vol. 259, pp. 698-707, May 2015.
        DOI: 10.1016/j.amc.2015.03.015
        '''
        t: float; z: float; y: float
        s: float; v: float; w: float        
        if(nu < 1.17683303804380831e0):
            t=nu*0.849738210666018375e0
            z=t*(156377.8333056294e0 +t*(48177.5705898287e0 +t*(5847.07218383812e0 \
            +t*(335.3978079672194e0 +t*7.84411868029912e0 ))))/(117762.02905535089e0 \
            +t*(-19007.26938370368e0 +t*(1376.2936928453140e0 +t*(-54.11372698481717e0 +t))))
            y=np.log(z)
        elif(nu < 3.82993088157949761e0):
            t=0.376917874490198033e0*nu -0.443569407329314587e0
            y=(489.140447310410217e0 +t*(5335.07269317261966e0 +t*(20169.0736140442509e0 \
            +t*(35247.8115595510907e0 +t*(30462.3668614714761e0 +t*(12567.9032426128967e0 \
            +t*(2131.86789357398657e0 +t*93.6520172085419439e0)))))))/(656.826207643060606e0 \
            +t*(4274.82831051941605e0 +t*(10555.7581310151498e0 +t*(12341.8742094611883e0 \
            +t*(6949.18854413197094e0 +t*(1692.19650634194002e0 +t*(129.221772991589751e0 +t)))))))
        elif(nu < 13.3854493161866553e0):
            t=0.104651569335924949e0*nu -0.400808277205416960e0
            y=(1019.84886406642351e0 +t*(9440.18255003922075e0 +t*(33947.6616363762463e0 \
            +t*(60256.7280980542786e0 +t*(55243.0045063055787e0 +t*(24769.8354802210838e0 \
            +t*(4511.77288617668292e0 +t*211.432806336150141e0)))))))/(350.502070353586442e0 \
            +t*(2531.06296201234050e0 +t*(6939.09850659439245e0 +t*(9005.40197972396592e0 \
            +t*(5606.73612994134056e0 +t*(1488.76634564005075e0 +t*(121.537028889412581e0 +t)))))))
        elif(nu < 53.2408277860982205e0):
            t=0.0250907164450825724e0*nu -0.335850513282463787e0
            y=(11885.8779398399498e0 +t*(113220.250825178799e0 +t*(408524.373881197840e0 \
            +t*(695674.357483475952e0 +t*(569389.917088505552e0 +t*(206433.082013681440e0 \
            +t*(27307.2535671974100e0 +t*824.430826794730740e0 )))))))/ (1634.40491220861182e0 \
            +t*(12218.1158551884025e0 +t*(32911.7869957793233e0 +t*(38934.6963039399331e0 \
            +t*(20038.8358438225823e0 +t*(3949.48380897796954e0 +t*(215.607404890995706e0 +t)))))))
        elif(nu < 188.411871723022843e0):
            t=0.00739803415638806339e0*nu -0.393877462475929313e0
            y=(11730.7011190435638e0 +t*(99421.7455796633651e0 +t*(327706.968910706902e0 \
            +t*(530425.668016563224e0 +t*(438631.900516555072e0 +t*(175322.855662315845e0 \
            +t*(28701.9605988813884e0 +t*1258.20914464286403e0)))))))/ (634.080470383026173e0 \
            +t*(4295.63159860265838e0 +t*(10868.5260668911946e0 +t*(12781.6871997977069e0 \
            +t*(7093.80732100760563e0 +t*(1675.06417056300026e0 +t*(125.750901817759662e0 +t)))))))
        else:
            v=nu**(-4.e0/3.e0)
            s=1080.13412050984017e0*v
            t=1.e0-s
            w=(1.12813495144821933e7 +t*(420368.911157160874e0 +t*(1689.69475714536117e0 +t))) / \
                (s*(6088.08350831295857e0 +t*(221.445236759466761e0 +t*0.718216708695397737e0 )))
            y=np.sqrt(w)
        return y
    
    @classmethod
    def _cal_eta_from_inv_FD(cls, n_d, m_star, T:float=300, method='JD_approx'):
        """
        Calculates scaled Fermi energy (E_f/kB.T) using inverse Fermi-Dirac integral
        of order-1/2.

        Parameters
        ----------
        n_d : 1D float array (unit: 1e18 cm^-3)
            Array containing carrier density data for compositions. Array size
            should be same as composition arrary. 
        m_star : ndarrary (unit: m0)
            Array containing effective mass for compositions.
        T : float, optional (unit: K)
            Temperature at which mobility calculations will be done. 
            The default is 300K.
        method : str, optional [available: 'JD_approx', 'minimax_piecewise']
            The approximate method to calculate the scalled Fermi energy (E_f/k_BT) 
            using inverse Fermi-Dirac integral of order-1/2. The default is JD_approx.
            JD_approx : Joyce-Dixon approximation (APL 31, 354 (1977)).
            minimax_piecewise : minimax approximation (Applied Mathematics and 
                                                       Computation 259, 698 (2015))
        Returns
        -------
        _eta : float/ndarray 
            The scaled Fermi energy eta_f = E_f/kB.T.

        """

        if method == 'JD_approx':
            # 2*((m0*k_B)**(3/2)) / (2*pi*h_bar**2)**(3/2) * 1e-6 => 1e15 K^(-3/2).cm^-3
            # Nc_3d = 4.829366089004772 * ((m_star*T)**(3/2)) # 1e15 cm^-3
            # n_by_Nc = n_d / Nc_3d * 1e3 # 1E18 cm^-3 / 1e15 cm^-3 = 1e3
            #         = (1/4.829366089004772 * 1e3)* n_d / ((m_star*T)**(3/2)) # cm^-3
            #         = 207.06651381777488 * n_d / ((m_star*T)**(3/2)) # cm^-3
            n_by_Nc = 207.06651381777488 * n_d / ((m_star*T)**(3/2))
            #----------------------------------------------------------------------
            # Joyce-Dixon approximation (APL 31, 354 (1977)) 
            # Upto cube term is considered
            return (np.log(n_by_Nc) + 0.3535534*n_by_Nc - 4.95009e-3*n_by_Nc*n_by_Nc
                    + 1.48386e-4*n_by_Nc*n_by_Nc*n_by_Nc)
        elif method == 'minimax_piecewise':
            # 2*((m0*k_B)**(3/2)) / (2*pi*h_bar**2)**(3/2) * 1e-6 => 1e15 K^(-3/2).cm^-3
            # Nc_3d = 4.829366089004772 * ((m_star*T)**(3/2)) # 1e15 cm^-3
            # N_3d = 2/sqrt(pi) * Nc_3d = 5.449356085110519*((m_star*T)**(3/2)) # 1e15 cm^-3
            # n_by_N = n_d / N_3d * 1e3 # 1E18 cm^-3 / 1e15 cm^-3 = 1e3
            #         = (1/5.449356085110519 * 1e3)* n_d / ((m_star*T)**(3/2)) # cm^-3
            #         = 183.5079199049476 * n_d / ((m_star*T)**(3/2)) # cm^-3
            n_by_N = 183.5079199049476 * n_d / ((m_star*T)**(3/2))
            if np.isscalar(n_by_N):
                return cls._Fukushima_iFD_half(n_by_N)
            else:
                # _Fukushima_iFD_half() is not vectorized. So vectorize the integrand function first.
                vec_Fukushima_iFD_half = np.vectorize(cls._Fukushima_iFD_half)
                return vec_Fukushima_iFD_half(n_by_N)
        else:
            raise ValueError(f'Requested {method} method is not implemeted yet. Contact developer.')
    
    @classmethod
    def _FD1_quad(cls, eta:float) -> float:
        """
        Integrand for F1(eta) = int_0_inf x/(1+exp(x-eta)) dx
        Use a numerically stable form: 1/(1+exp(x-eta)) = expit(eta-x)
        Numerial integration using scipy.integrate.quad over [0, inf)
        """
        _FD_func = lambda x, eta: x * special.expit(eta-x)
        return integrate.quad(_FD_func, 0, np.inf, args=(eta,))[0]
    
    @classmethod
    def _FD_integral_order_1(cls, eta_f, use_numerical_int:bool=False):
        """
        Compute the integration numerically, using scipy.quad or dilogarithm approach
        using spence function. The default is False. 
        """
        if use_numerical_int:
            if np.isscalar(eta_f):
                return cls._FD1_quad(eta_f)
            else:
                # quad is not vectorized. So vectorize the integrand function first.
                vec_FD1_quad = np.vectorize(cls._FD1_quad)
                return vec_FD1_quad(eta_f)
        return (-1) * special.spence(1+np.exp(eta_f))
       
    @classmethod
    def _cal_Fermi_Dirac_integral(cls, n_d, m_star, T:float=300, 
                                  inv_half_FD_method:str='JD_approx', 
                                  FD_order:str = 'zero', 
                                  use_numerical_integration:bool=False):
        """
        Calculates Fermi-Dirac integral.
        eta = E_f/(k_B.T)
        F_0(eta) = ln(1+e^eta) #=> analytical solution
        F_1(eta) = int_0_inf x/(1+exp(x-eta)) dx =>
            Use a numerically stable form: 1/(1+exp(x-eta)) = expit(eta-x)

        Parameters
        ----------
        n_d : 1D float array (unit: 1e18 cm^-3)
            Array containing carrier density data for compositions. Array size
            should be same as composition arrary. 
        m_star : ndarrary (unit: m0)
            Array containing effective mass for compositions.
        T : float, optional (unit: K)
            Temperature at which mobility calculations will be done. 
            The default is 300K.
        inv_half_FD_method : str, optional [available: 'JD_approx', 'minimax_piecewise']
            The approximate method to calculate the scalled Fermi energy (E_f/k_BT) 
            using inverse Fermi-Dirac integral of order-1/2. The default is JD_approx.
            JD_approx : Joyce-Dixon approximation (APL 31, 354 (1977)).
            minimax_piecewise : minimax approximation (Applied Mathematics and 
                                                       Computation 259, 698 (2015))
        FD_order : str, optional ['zero', 'one', 'two', ...]
            FD integral oder. The default is 'zero'.
        use_numerical_integration : bool, optional
            Compute the integration numerically, using scipy.quad. If False, 
            polylogaritm approach is used for FD_order > 1. For FD_order = 1 
            dilogarithm formulation is used.For FD_order=0, analytical solution 
            is used always.
            The default is False. 

        Returns
        -------
        float
            FD integral value.

        """
        # Fermi eta = E_f/(k_B.T)
        eta_f = cls._cal_eta_for_FD_integral(n_d, m_star, T=T, method=inv_half_FD_method)
        
        if FD_order == 'zero':
            #return np.log(1+np.exp(eta_f))
            # np.logaddexp(a,b) = np.log(exp(a)+exp(b)) => better stable in log operation
            return np.logaddexp(0,eta_f)
        if FD_order == 'one':
            cls._FD_integral_order_1(eta_f, use_numerical_int=use_numerical_integration)
        else:
            raise ValueError(f'{FD_order} FD integral is not implemented yet. Contact developer.')