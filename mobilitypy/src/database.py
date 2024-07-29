'''
Units:
mass_density => kg/m3
lattice_a0 => angstrom
lattice_c0 => angstrom
bandgap => eV
e_effective_mass => m0
alloy_scattering_potential => eV
static_dielectric_constant => epsilon_0
high_frequency_dielectric_constant => epsilon_0
LA_phonon_velocity => m/s
TA_phonon_velocity => m/s
deformation_potential => eV
PO_phonon_energy => eV
electromechanical_coupling_const => unitless
'''
database = {'GaN': 
            {'mass_density': 6150,
             'lattice_a0': 3.189,
             'lattice_c0': 5.185,
             'bandgap': 3.43,
             'e_effective_mass': 0.20,
             'alloy_scattering_potential': 1.8,
             'static_dielectric_constant': 8.90,
             'high_frequency_dielectric_constant': 5.35,
             'LA_phonon_velocity': 6560,
             'TA_phonon_velocity': 2680,
             'deformation_potential': 8.3,
             'PO_phonon_energy': 91.2e-3,
             'electromechanical_coupling_const': 0.045},
           'AlN': 
            {'mass_density': 3230,
             'lattice_a0': 3.112,
             'lattice_c0': 4.982,
             'bandgap': 6.20,
             'e_effective_mass': 0.40,
             'alloy_scattering_potential': 1.8,
             'static_dielectric_constant': 8.50,
             'high_frequency_dielectric_constant': 4.60,
             'LA_phonon_velocity': 9060,
             'TA_phonon_velocity': 3700,
             'deformation_potential': 9.5,
             'PO_phonon_energy': 99.0e-3,
             'electromechanical_coupling_const': 0.106},
           'AlGaN': 
            {'mass_density': 0,
             'lattice_a0': 0,
             'lattice_c0': 0,
             'bandgap': 0.7,
             'e_effective_mass': 0,
             'alloy_scattering_potential': 0,
             'static_dielectric_constant': 0,
             'high_frequency_dielectric_constant': 0,
             'LA_phonon_velocity': 0,
             'TA_phonon_velocity': 0,
             'deformation_potential': 0,
             'PO_phonon_energy': 0,
             'electromechanical_coupling_const': 0}}