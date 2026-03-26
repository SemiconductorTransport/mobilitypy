'''
The database contains the following information for each material. 
For ternary or higher-order materials, the parameter values represent 
the bowing coefficients.

Ref-1: Bassaler et. al., Adv. Electron. Mater., 11, 2400069 (2025)
Ref-1': Mondal et. al., APL Electronic Devices 1, 026117 (2025)
Ref-2: Pant et al, APL 117, 242105 (2020)
Ref-3: http://www.ioffe.ru/SVA/NSM/Semicond/index.html
Ref-4: Vurgaftman et. al., J. Appl. Phys. 94, 3675–3696 (2003)
Ref-5: Shimada, Jpn. J. Appl. Phys., Vol. 45, No. 12 (2006)
Ref-6: Dreyer et al. Phys. Rev. X 6, 021038 (2016)
Ref-7: Debdeep Jena. "Quantum physics of semiconductor materials and devices"", 
Chapter - scattering, mobility, and velocity saturation, page 555, USA, 2002. 

Units:
mass_density => kg/m3 # Ref-3
lattice_a0 => angstrom # Ref-4
lattice_c0 => angstrom # Ref-4
bandgap => eV # Ref-4
bandgap_alpha => eV/K # Ref-4
bandgap_beta => K # Ref-4
e_effective_mass => m0 # Ref-2
alloy_scattering_potential => eV # Ref-2
static_dielectric_constant => epsilon_0
high_frequency_dielectric_constant => epsilon_0
LA_phonon_velocity => m/s # Ref-7
TA_phonon_velocity => m/s # Ref-7
CB_deformation_potential => eV # Ref-7
PO_phonon_energy => eV
electromechanical_coupling_const => unitless
C_ij => GPa # Ref-4
e_ij => Piezoelectric constants # C/m^2 # e_15 from Ref-5; e33, e31 from Ref-6
'''
database = {
            # ========================== Binaries =============================
            # GaN, AlN
            'GaN': 
            {'mass_density': 6150,
             'lattice_a0': 3.189,
             'lattice_c0': 5.185,
             'bandgap': 3.510,
             'bandgap_alpha': 0.909e-3,
             'bandgap_beta': 830,
             'e_effective_mass': 0.20,
             'alloy_scattering_potential': 1.0,
             'static_dielectric_constant': 8.90,
             'high_frequency_dielectric_constant': 5.35,
             'LA_phonon_velocity': 7960,
             'TA_phonon_velocity': 5020,
             'CB_deformation_potential': 8.3,
             'PO_phonon_energy': 91.2e-3,
             #'electromechanical_coupling_const': 0.045,
             'C_11': 390, 'C_12': 145, 'C_13': 106, 'C_33': 398, 'C_44': 105,
             'e_31': -1.863, 'e_33': 1.020, 'e_15': -0.38 
             },
           'AlN': 
            {'mass_density': 3230,
             'lattice_a0': 3.112,
             'lattice_c0': 4.982,
             'bandgap': 6.25,
             'bandgap_alpha': 1.799e-3,
             'bandgap_beta': 1462,
             'e_effective_mass': 0.31,
             'alloy_scattering_potential': 1.8,
             'static_dielectric_constant': 8.50,
             'high_frequency_dielectric_constant': 4.60,
             'LA_phonon_velocity': 11270,
             'TA_phonon_velocity': 6220,
             'CB_deformation_potential': 9.5,
             'PO_phonon_energy': 99.0e-3,
             #'electromechanical_coupling_const': 0.106,
             'C_11': 396, 'C_12': 137, 'C_13': 108, 'C_33': 373, 'C_44': 116,
             'e_31': -2.027, 'e_33': 1.569, 'e_15': -0.41
             },
            # ========================== Ternaries ============================
            # AlGaN
           'AlGaN': # AlGaN == AlxGa(1-x)N, GaxAl(1-x)N  => Binaries=[AlN, GaN]
            {'mass_density': 0,
             'lattice_a0': 0,
             'lattice_c0': 0,
             'bandgap': 0.7,
             'bandgap_alpha': 0,
             'bandgap_beta': 0,
             'e_effective_mass': 0,
             'alloy_scattering_potential': -1.6,
             'static_dielectric_constant': 0,
             'high_frequency_dielectric_constant': 0,
             'LA_phonon_velocity': 0,
             'TA_phonon_velocity': 0,
             'CB_deformation_potential': 0,
             'PO_phonon_energy': 0,
             #'electromechanical_coupling_const': 0,
             'C_11': 0, 'C_12': 0, 'C_13': 0, 'C_33': 0, 'C_44': 0,
             'e_31': 0, 'e_33': 0, 'e_15': 0
             }
            }
