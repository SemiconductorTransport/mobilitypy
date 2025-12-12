import scipy as sc

## ============================================================================
# Universal constants
e_charge = sc.constants.e # C
pi_ = sc.constants.pi 
eps_0 = sc.constants.epsilon_0 # CV^-1m^-1
h_bar = sc.constants.hbar # J.s
e_mass = sc.constants.m_e # kg
k_B = sc.constants.Boltzmann # J K^-1
sqrt_3_by_2 = 0.8660254037844386 # sqrt(3)/2
## ============================================================================
## Constants for different 2DEG mobility equations
fact_q_TF = 37.79452249229504 # m0 * e^2 / (2pi * eps_0 * h_bar^2) * 1e-9 => nm^-1
fact_b = 9.931409618986013 # (33pi/4 * fact_q_TF)^(1/3) => nm^-1/3
fact_irf_dis = 6.528368003403906e+18 # (m0*e^4)/(h_bar^3*eps_0^2) => s^-1
fact_alloy = 3.7383724882773685e+15 # 3/16*m0*e^2/1e18 => nm s^-1
fact_phonon = 8.762231231847618e+10 # m0*e^2*k_b/(pi*h_bar^3) => kg.K^-1J^2s^-3
fact_pop_k0 = 5.123167219674931 # sqrt(2*m_0*e/h_bar^2)*1e9 => nm^-1
fact_pop_y = 2.777985128879875e+3 # #(pi*h_bar^2)*1e18/m_0/k_B
fact_pop = 1.8039021162385878e+17 # 1e-9*e^2*e*m_0/(2*eps_0*h_bar^3)
## ============================================================================
## Constants for different 3DEG mobility equations
#omega_0 * n_3d = X * 1e-24 cm^3 * Y * 1e18 cm^-3 = XY * 1e-6 
fact_ad_3d = 21.1699056589355 # (2*e*h_bar*k_B)/(3*pi*m0*e^2) * 1e10 => cm^-2.V-1.s-1
fact_Nc_3d = 4.829366079143929 # 2*((m0*k_B)**(3/2)) / (2*pi*h_bar**2)**(3/2) * 1e-6 => 1e15 cm^-3
fact_pop_3d = 0.156926697141818 # 4.0*pi*eps_0*h_bar*h_bar/(sqrt(2)*(e*m0)**(3/2))*1e4 => cm^2V^-1s^-1
fact_dp_3d = 1.5333300251501703e-6 # 2*h_bar*e*1e-2/(3*pi*m0*e*e*1e18) => cm^2V^-1s^-1
