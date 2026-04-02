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
        t: float; z: float; y: float; s: float; v: float; w: float     
        if np.isnan(nu): return nu
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
    def _Fukushima_FD_minus_one_half(cls, nu: float) -> float:
        '''
        Double precision rational minimax approximation of Fermi-Dirac integral 
        of order k=-1/2
        
        Taken from github.com/scott-maddox/fdint and
        https://www.researchgate.net/publication/267330765_xfdhtxt_Fortran_program_package_to_compute_Fermi-Dirac_integrals_of_half_integer_orders_k-92-72-52212_and_of_integer_orders_k01210
        
        T. Fukushima, "Precise and fast computation of Fermi–Dirac integral of 
        integer and half integer order by piecewise minimax rational approximation"
        Applied Mathematics and Computation, vol. 259, pp. 708-729, May 2015.
        DOI: https://doi.org/10.1016/j.amc.2015.03.009
        '''
        v: float; t: float; y: float; s: float; w: float  
        if np.isnan(nu): return nu
        if (nu < -2e0):
            v = np.exp(nu)
            t = v*7.38905609893065023e0
            y = v*(1.77245385090551603e0-v*(40641.4537510284430e0+t*(9395.7080940846442e0\
            +t*(649.96168315267301e0+t*(12.7972295804758967e0+t*0.00153864350767585460e0))))\
            /(32427.1884765292940e0+t*(11079.9205661274782e0+t*(1322.96627001478859e0\
            +t*(63.738361029333467e0+t)))))
        elif (nu < 0e0):
            s = -0.5e0*nu
            t = 1.e0 - s
            y = (272.770092131932696e0+t*(30.8845653844682850e0+t*(-6.43537632380366113e0\
            +t*(14.8747473098217879e0+t*(4.86928862842142635e0+t*(-1.53265834550673654e0\
            +t*(-1.02698898315597491e0+t*(-0.177686820928605932e0-t*0.00377141325509246441e0\
            ))))))))/(293.075378187667857e0+s*(305.818162686270816e0+s*(299.962395449297620e0\
            +s*(207.640834087494249e0+s*(92.0384803181851755e0+s*(37.0164914112791209e0\
            +s*(7.88500950271420583e0+s)))))))
        elif (nu < 2e0):
            t = 0.5e0*nu
            y = (3531.50360568243046e0+t*(6077.5339658420037e0+t*(6199.7700433981326e0\
            +t*(4412.78701919567594e0+t*(2252.27343092810898e0+t*(811.84098649224085e0\
            +t*(191.836401053637121e0+t*23.2881838959183802e0)))))))/(3293.83702584796268e0\
            +t*(1528.97474029789098e0+t*(2568.48562814986046e0+t*(925.64264653555825e0\
            +t*(574.23248354035988e0+t*(132.803859320667262e0+t*(29.8447166552102115e0\
            +t)))))))
        elif (nu < 5e0):
            t = 0.3333333333333333333e0 * (nu - 2.e0)
            y = (4060.70753404118265e0+t*(10812.7291333052766e0+t*(13897.5649482242583e0\
            +t*(10628.4749852740029e0+t*(5107.70670190679021e0+t*(1540.84330126003381e0\
            +t*(284.452720112970331e0+t*29.5214417358484151e0)))))))/(1564.58195612633534e0\
            +t*(2825.75172277850406e0+t*(3189.16066169981562e0+t*(1955.03979069032571e0\
            +t*(828.000333691814748e0+t*(181.498111089518376e0+t*(32.0352857794803750e0\
            +t)))))))
        elif(nu < 10e0):
            t = 0.2e0*nu-1.e0
            y=(1198.41719029557508e0+t*(3263.51454554908654e0+t*(3874.97588471376487e0\
            +t*(2623.13060317199813e0+t*(1100.41355637121217e0+t*(267.469532490503605e0\
            +t*(25.4207671812718340e0+t*0.389887754234555773e0\
            )))))))/(273.407957792556998e0+t*(595.918318952058643e0\
            +t*(605.202452261660849e0+t*(343.183302735619981e0+t*(122.187622015695729e0\
            +t*(20.9016359079855933e0+t))))))
        elif (nu < 20e0):
            t = 0.1e0*nu - 1.e0
            y = (9446.00169435237637e0+t*(36843.4448474028632e0+t*(63710.1115419926191e0\
            +t*(62985.2197361074768e0+t*(37634.5231395700921e0+t*(12810.9898627807754e0\
            +t*(1981.56896138920963e0+t*81.4930171897667580e0)))))))/(1500.04697810133666e0\
            +t*(5086.91381052794059e0+t*(7730.01593747621895e0+t*(6640.83376239360596e0\
            +t*(3338.99590300826393e0+t*(860.499043886802984e0+t*(78.8565824186926692e0\
            +t)))))))
        elif(nu < 40e0):
            t = 0.05e0*nu-1.e0
            y=(22977.9657855367223e0+t*(123416.616813887781e0+t*(261153.765172355107e0\
            +t*(274618.894514095795e0+t*(149710.718389924860e0+t*(40129.3371700184546e0\
            +t*(4470.46495881415076e0+t*132.684346831002976e0)))))))/(2571.68842525335676e0\
            +t*(12521.4982290775358e0+t*(23268.1574325055341e0+t*(20477.2320119758141e0\
            +t*(8726.52577962268114e0+t*(1647.42896896769909e0+t*(106.475275142076623e0\
            +t)))))))
        else:
            w = 1.e0/(nu*nu)
            t = 1600.e0 * w
            y = np.sqrt(nu)*2e0*(1.e0-w*(0.411233516712009968e0+t*(0.00110980410034088951e0\
            +t*(0.0000113689298990173683e0+t*(2.56931790679436797e-7+t*(9.97897786755446178e-9\
            +t*8.67667698791108582e-10))))))
        return 0.5641895835477563 * y # 1/(Gamma(1/2))*y
    
    @classmethod
    def _Fukushima_FD_one_half(cls, nu: float) -> float:
        '''
        Double precision rational minimax approximation of Fermi-Dirac integral 
        of order k=1/2
        
        Taken from github.com/scott-maddox/fdint and
        https://www.researchgate.net/publication/267330765_xfdhtxt_Fortran_program_package_to_compute_Fermi-Dirac_integrals_of_half_integer_orders_k-92-72-52212_and_of_integer_orders_k01210
        
        T. Fukushima, "Precise and fast computation of Fermi–Dirac integral of 
        integer and half integer order by piecewise minimax rational approximation"
        Applied Mathematics and Computation, vol. 259, pp. 708-729, May 2015.
        DOI: https://doi.org/10.1016/j.amc.2015.03.009
        '''
        t: float; y: float; s: float; w: float  
        if np.isnan(nu): return nu
        if(nu < -2.e0):
            s=np.exp(nu)
            t=s*7.38905609893065023e0
            y=s*(0.886226925452758014e0-s*(19894.4553386951666e0+t*(4509.64329955948557e0 \
            +t*(303.461789035142376e0+t*(5.7574879114754736e0 +t*0.00275088986849762610e0 \
            ))))/(63493.915041308052e0+t*(19070.1178243603945e0+t*(1962.19362141235102e0 \
            +t*(79.250704958640158e0+t)))))
        elif(nu < 0.e0):
            s=-0.5e0*nu
            t=1.e0-s
            y=(149.462587768865243e0+t*(22.8125889885050154e0+t*(-0.629256395534285422e0 \
            +t*(9.08120441515995244e0+t*(3.35357478401835299e0+t*(-0.473677696915555805e0 \
            +t*(-0.467190913556185953e0+t*(-0.0880610317272330793e0-t*0.00262208080491572673e0 \
            ))))))))/(269.94660938022644e0+s*(343.6419926336247e0+s*(323.9049470901941e0 \
            +s*(218.89170769294024e0+s*(102.31331350098315e0+s*(36.319337289702664e0 \
            +s*(8.3317401231389461e0+s)))))))
        elif(nu<2.e0):
            t=0.5e0*nu
            y=(71652.717119215557e0+t*(134954.734070223743e0+t*(153693.833350315645e0 \
            +t*(123247.280745703400e0+t*(72886.293647930726e0+t*(32081.2499422362952e0 \
            +t*(10210.9967337762918e0+t*(2152.71110381320778e0+t*232.906588165205042e0 \
            ))))))))/(105667.839854298798e0+t*(31946.0752989314444e0+t*(71158.788776422211e0 \
            +t*(15650.8990138187414e0+t*(13521.8033657783433e0+t*(1646.98258283527892e0 \
            +t*(618.90691969249409e0+t*(-3.36319591755394735e0+t))))))))
        elif(nu<5.e0):
            t=0.3333333333333333333e0*(nu-2.e0)
            y=(23744.8706993314289e0+t*(68257.8589855623002e0+t*(89327.4467683334597e0 \
            +t*(62766.3415600442563e0+t*(20093.6622609901994e0+t*(-2213.89084119777949e0 \
            +t*(-3901.66057267577389e0-t*948.642895944858861e0)))))))/(9488.61972919565851e0 \
            +t*(12514.8125526953073e0+t*(9903.44088207450946e0+t*(2138.15420910334305e0 \
            +t*(-528.394863730838233e0+t*(-661.033633995449691e0+t*(-51.4481470250962337e0 \
            +t)))))))
        elif(nu<10.e0):
            t=0.2e0*nu-1.e0
            y=(311337.452661582536e0+t*(1.11267074416648198e6+t*(1.75638628895671735e6 \
            +t*(1.59630855803772449e6+t*(910818.935456183774e0+t*(326492.733550701245e0 \
            +t*(65507.2624972852908e0+t*4809.45649527286889e0)))))))/(39721.6641625089685e0 \
            +t*(86424.7529107662431e0+t*(88163.7255252151780e0+t*(50615.7363511157353e0 \
            +t*(17334.9774805008209e0+t*(2712.13170809042550e0+t*(82.2205828354629102e0 \
            -t)))))))*0.999999999999999877e0
        elif(nu<20.e0):
            t=0.1e0*nu-1.e0
            y=(7.26870063003059784e6+t*(2.79049734854776025e7+t*(4.42791767759742390e7 \
            +t*(3.63735017512363365e7+t*(1.55766342463679795e7+t*(2.97469357085299505e6 \
            +t*154516.447031598403e0))))))/(340542.544360209743e0+t*(805021.468647620047e0 \
            +t*(759088.235455002605e0+t*(304686.671371640343e0+t*(39289.4061400542309e0 \
            +t*(582.426138126398363e0+t*(11.2728194581586028e0-t)))))))
        elif(nu<40.e0):
            t=0.05e0*nu-1.e0
            y=(4.81449797541963104e6+t*(1.85162850713127602e7+t*(2.77630967522574435e7 \
            +t*(2.03275937688070624e7+t*(7.41578871589369361e6+t*(1.21193113596189034e6 \
            +t*63211.9545144644852e0))))))/(80492.7765975237449e0+t*(189328.678152654840e0 \
            +t*(151155.890651482570e0+t*(48146.3242253837259e0+t*(5407.08878394180588e0 \
            +t*(112.195044410775577e0-t))))))
        else:
            w=1.e0/(nu*nu)
            s=1.e0-1600.e0*w
            y=nu*np.sqrt(nu)*0.666666666666666667e0*(1.e0+w \
            *(8109.79390744477921e0+s*(342.069867454704106e0+s*1.07141702293504595e0)) \
            /(6569.98472532829094e0+s*(280.706465851683809e0+s)))
        return 1.1283791670955126 * y # 1/(Gamma(3/2))*y
    
    @classmethod
    def _Fukushima_FD_one(cls, nu: float) -> float:
        '''
        Double precision rational minimax approximation of Fermi-Dirac integral 
        of order k=1
        
        Taken from github.com/scott-maddox/fdint and
        https://www.researchgate.net/publication/267330765_xfdhtxt_Fortran_program_package_to_compute_Fermi-Dirac_integrals_of_half_integer_orders_k-92-72-52212_and_of_integer_orders_k01210
        
        T. Fukushima, "Precise and fast computation of Fermi–Dirac integral of 
        integer and half integer order by piecewise minimax rational approximation"
        Applied Mathematics and Computation, vol. 259, pp. 708-729, May 2015.
        DOI: https://doi.org/10.1016/j.amc.2015.03.009
        '''
        t: float; y: float; absnu: float; s: float
        if np.isnan(nu): return nu
        absnu = -abs(nu)
        if(absnu < -2.e0):
            s=np.exp(absnu)
            t=s*7.38905609893065023e0
            y=s*(1.e0 -s*(22189.1070807945062e0 +t*(4915.92700908746777e0 \
            +t*(322.901386168881348e0+t*(5.9897442965804548e0+t*0.00397641173774375092e0 \
            ))))/(88756.428323178025e0+t*(25002.3197546553836e0+t*(2389.06277237306633e0 \
            +t*(88.376214553692756e0+t)))))
        elif(absnu <= 0.e0):
            s=-0.5e0*absnu
            t=1.e0-s
            y=(145.488167182330098e0+t*(251.392824471576922e0+t*(56.6537141912783024e0 \
            +t*(17.9918985363509694e0+t*(20.1369115558099802e0+t*(7.09659390228556164e0 \
            +t*(0.199701180197912643e0+t*(-0.403173132925886253e0-t*0.0792966701498222697e0 \
            ))))))))/(606.0757707716040e0+s*(374.1806357435014e0+s*(252.1367051536344e0 \
            +s*(27.2746245830016e0+s*(-61.57766112137513e0+s*(-53.72117554363975e0 \
            +s*(-25.678454878692950e0+s*(-7.1995819520154718e0-s))))))))
        if(nu > 0.e0):
            y=-y+1.64493406684822644e0+0.5e0*nu*nu
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
            n_by_Nc = 207.06651381777488 * n_d / (np.sqrt(m_star*T))**3
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
            n_by_N = 183.5079199049476 * n_d / (np.sqrt(m_star*T))**3
            if np.isscalar(n_by_N):
                return cls._Fukushima_iFD_half(n_by_N)
            else:
                # _Fukushima_iFD_half() is not vectorized. So vectorize the integrand function first.
                vec_Fukushima_iFD_half = np.vectorize(cls._Fukushima_iFD_half)
                return vec_Fukushima_iFD_half(n_by_N)
        else:
            raise ValueError(f'Requested {method} method is not implemeted yet. Contact developer.')
    
    @classmethod
    def _FD_dis_chg_Integration(cls, eta, B):
        """
        Use a numerically stable form: 1/(1+exp(x-eta)) = expit(eta-x)
        Numerial integration using scipy.integrate.quad over [0, inf)
        """
        _FD_func = lambda x, eta, B: (B+2*x)*np.sqrt(x+(x*x/B))*special.expit(eta-x)
        return integrate.quad_vec(_FD_func, 0, np.inf, args=(eta,B))[0]
    
    @classmethod
    def _FD_dis_str_Integral(cls, x, eta, B):
        """
        Use a numerically stable form: 1/(1+exp(x-eta)) = expit(eta-x)
        """
        y = np.sqrt(B/(B+x))
        return (5.0-5.0*y-(x/(B+x))*y)*x*np.sqrt(x)*special.expit(eta-x)/B/(1-y)**2
    
    @classmethod
    def _FD_dis_str_Integration(cls, eta_f, B):
        return integrate.quad_vec(cls._FD_dis_str_Integral, 0, np.inf, args=(eta_f,B))[0]
        
    @classmethod
    def _FD_integral_order_1(cls, eta_f, FD_integration_approach:str='minimax_piecewise'):
        """
        Compute the integration numerically, using scipy.quad or dilogarithm approach
        using spence function, or Fukushima's 'minimax_piecewise'. 
        The default is 'minimax_piecewise' 
        """
        if FD_integration_approach == 'num':
            """
            Integrand for F1(eta) = int_0_inf x/(1+exp(x-eta)) dx
            Use a numerically stable form: 1/(1+exp(x-eta)) = expit(eta-x)
            Numerial integration using scipy.integrate.quad over [0, inf)
            """
            _FD_func = lambda x, eta: x * special.expit(eta_f-x)
            return integrate.quad_vec(_FD_func, 0, np.inf, args=(eta_f,))[0]
        elif FD_integration_approach == 'polylog':
            return (-1) * special.spence(1.0+np.exp(eta_f))
        else:
            if np.isscalar(eta_f):
                return cls._Fukushima_FD_one(eta_f)
            else:
                vec_FD1_Fukushima = np.vectorize(cls._Fukushima_FD_one)
                return vec_FD1_Fukushima(eta_f)
            
    @classmethod
    def _FD_integral_order_2(cls, eta_f, FD_integration_approach:str='minimax_piecewise'):
        """
        Compute the integration numerically, using scipy.quad or dilogarithm approach
        using spence function, or Fukushima's 'minimax_piecewise'. 
        The default is 'minimax_piecewise' 
        """
        #if FD_integration_approach == 'num':
        """
        Use a numerically stable form: 1/(1+exp(x-eta)) = expit(eta-x)
        Numerial integration using scipy.integrate.quad over [0, inf)
        """
        _FD_func = lambda x, eta: x*x * special.expit(eta_f-x)
        return 0.5 * integrate.quad_vec(_FD_func, 0, np.inf, args=(eta_f,))[0]
            
    @classmethod
    def _FD_integral_order_m1h(cls, eta_f, FD_integration_approach:str='minimax_piecewise'):
        """
        Compute the integration numerically, using scipy.quad or dilogarithm approach
        using spence function, or Fukushima's 'minimax_piecewise'. 
        The default is 'minimax_piecewise' 
        """
        if FD_integration_approach == 'minimax_piecewise':
            if np.isscalar(eta_f):
                return cls._Fukushima_FD_minus_one_half(eta_f)
            else:
                vec_FDm1h_Fukushima = np.vectorize(cls._Fukushima_FD_minus_one_half)
                return vec_FDm1h_Fukushima(eta_f)
        else:
            raise ValueError(f'Only {FD_integration_approach} method is implemented for Fermi Diract -1/2 integral.')
    
    @classmethod
    def _FD_integral_order_1h(cls, eta_f, FD_integration_approach:str='minimax_piecewise'):
        """
        Compute the integration numerically, using scipy.quad or dilogarithm approach
        using spence function, or Fukushima's 'minimax_piecewise'. 
        The default is 'minimax_piecewise' 
        """
        if FD_integration_approach == 'minimax_piecewise':
            if np.isscalar(eta_f):
                return cls._Fukushima_FD_one_half(eta_f)
            else:
                vec_FD1h_Fukushima = np.vectorize(cls._Fukushima_FD_one_half)
                return vec_FD1h_Fukushima(eta_f)
        else:
            raise ValueError(f'Only {FD_integration_approach} method is implemented for Fermi Diract 1/2 integral.')
            
    @classmethod
    def _cal_Fermi_Dirac_integral(cls, eta_f, FD_order:str = 'zero', 
                                  FD_int_approach:str='minimax_piecewise'):
        """
        Calculates Fermi-Dirac integral.
        eta = E_f/(k_B.T)
        F_0(eta) = ln(1+e^eta) #=> analytical solution
        F_1(eta) = int_0_inf x/(1+exp(x-eta)) dx =>
            Use a numerically stable form: 1/(1+exp(x-eta)) = expit(eta-x)

        Parameters
        ----------
        eta_f : 1D float array (unit: uniless)
            Array containing carrier density data for compositions. Array size
            should be same as composition arrary. 
        FD_order : str, optional [options: 'zero', 'one', 'm_one_half', 'one_half']
            FD integral oder. The default is 'zero'.
        FD_int_approach : str, optional [options: 'num', 'minimax_piecewise', 'polylog']
            Compute the Fermi-Dirac integral. The default is minimax_piecewise.
            If num: calculated numerically, using scipy.quad. 
            If polylog: polylogaritm approach is used for FD_order > 1. For FD_order = 1 
            dilogarithm formulation is used.For FD_order=0, analytical solution 
            is used always.
            If minimax_piecewise: use Fukishima's minimax_piecewise approximation.

        Returns
        -------
        float
            FD integral value.

        """       
        if FD_order == 'zero':
            # np.log1p(x) = log(1 + x)
            return np.log1p(np.exp(eta_f))
            # np.logaddexp(a,b) = np.log(exp(a)+exp(b)) => better stable in log operation
            #return np.logaddexp(0,eta_f) # Produce warning when nan encounter 
        elif FD_order == 'm_one_half':
            return cls._FD_integral_order_m1h(eta_f, FD_integration_approach=FD_int_approach)
        elif FD_order == 'one_half':
            return cls._FD_integral_order_1h(eta_f, FD_integration_approach=FD_int_approach)
        elif FD_order == 'one':
            return cls._FD_integral_order_1(eta_f, FD_integration_approach=FD_int_approach)
        elif FD_order == 'two':
            return cls._FD_integral_order_2(eta_f, FD_integration_approach=FD_int_approach)
        else:
            raise ValueError(f'{FD_order} FD integral is not implemented yet. Contact developer.')