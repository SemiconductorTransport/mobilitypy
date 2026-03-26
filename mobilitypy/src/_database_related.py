#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 14:03:02 2026

@author: badal.mondal
"""
from .database import database

#%%# ==============================================================================
class _DataBase:
    database_units = {'mass_density': 'kg/m3',
                      'lattice_a0': 'angstrom',
                      'lattice_c0': 'angstrom',
                      'bandgap': 'eV',
                      'bandgap_alpha': 'eV/K',
                      'bandgap_beta': 'K',
                      'e_effective_mass': 'm0',
                      'alloy_scattering_potential': 'eV',
                      'static_dielectric_constant': 'epsilon_0',
                      'high_frequency_dielectric_constant': 'epsilon_0',
                      'LA_phonon_velocity': 'm/s',
                      'TA_phonon_velocity': 'm/s',
                      'CB_deformation_potential': 'eV',
                      'PO_phonon_energy': 'eV',
                      'electromechanical_coupling_const': 'unitless',
                      'C_ij': 'GPa',
                      'e_ij': 'C/m2'
                      }
    parameter_name_mapping = {'mass_density': 'Mass density',
                              'lattice_a0': 'a0 lattice parameter',
                              'lattice_c0': 'c0 lattice parameter',
                              'bandgap': 'Bandgap',
                              'bandgap_alpha': "Bandgap Varshni parameter alpha",
                              'bandgap_beta': 'Bandgap Varshni parameter beta',
                              'e_effective_mass': 'Electron effective mass',
                              'alloy_scattering_potential': 'Alloy scattering potential',
                              'static_dielectric_constant': 'Static dielectric constant',
                              'high_frequency_dielectric_constant': 'High frequency dielectric constant',
                              'LA_phonon_velocity': 'Longitudinal acoustic phonon velocity',
                              'TA_phonon_velocity': 'Transversal acoustic phonon velocity',
                              'CB_deformation_potential': 'Conduction band deformation potential',
                              'PO_phonon_energy': 'Polar optical phonon energy',
                              'electromechanical_coupling_const': 'Electromechanical coupling coefficient',
                              'C_ij': 'Elastic constants',
                              'e_ij': 'Piezoelectic constant'}
    def __init__(self):
        pass
    
    @classmethod
    def _print_database(cls, for_material=None):
        """
        This function prints the information of material parameters from the database.

        Parameters
        ----------
        for_material : string (case sensitive), optional
            The material name for which the parameters will be printed. 
            The name should match the name in database. If None, prints general
            information about the database.
            The default is None.

        Returns
        -------
        None.

        """
        
        database_info = '''
Ref-1: Bassaler et. al., Adv. Electron. Mater., 11, 2400069 (2025)
Ref-1': Mondal et. al., APL Electronic Devices 1, 026117 (2025)
Ref-2: Pant et al, APL 117, 242105 (2020)
Ref-3: http://www.ioffe.ru/SVA/NSM/Semicond/index.html
Ref-4: Vurgaftman et. al., J. Appl. Phys. 94, 3675–3696 (2003)
Ref-5: Shimada, Jpn. J. Appl. Phys., Vol. 45, No. 12 (2006)
Ref-6: Dreyer et al. Phys. Rev. X 6, 021038 (2016)
Ref-7: Debdeep Jena. "Quantum physics of semiconductor materials and devices"", 
Chapter - scattering, mobility, and velocity saturation, page 555, USA, 2002. 
        '''
        database_header_info = """
The database contains the following information for each material. 
For ternary or higher-order materials, the parameter values represent 
the bowing coefficients.
        """
        if for_material is None:
            print("""NOTE: Use 'print(DataBase().print_database(for_material=<material name> e.g. 'AlN'))'
to print parameter values for the specific material""")
            print(f'{database_header_info}{database_info}\nParameters (Name description => parameter name : unit):')  
            for key, value in cls.database_units.items():
                print(f'{cls.parameter_name_mapping.get(key)} => {key}: {value}')
        else:
            params_db = database.get(for_material).copy()
            if params_db:
                database_header_info = f"""
NOTE: For ternary or higher-order materials, the parameter values represent 
the bowing coefficients.
The database contains the following information for {for_material}: 
                """
                print(database_header_info)
                for key, value in params_db.items():
                    if key.startswith('C_'):
                        print(f'{cls.parameter_name_mapping.get("C_ij")} => {key}: {value} {cls.database_units.get("C_ij")}')
                    else:
                        print(f'{cls.parameter_name_mapping.get(key)} => {key}: {value} {cls.database_units.get(key)}')
            else:
                print(f'Error: "{for_material}" material does not exists in database yet. Contact developers.')
                
    def _update_database(self, for_material=None, with_new_database=None):
        """
        To update the material parameters in the database. 
        Not implemented yet. Contact developer.

        Parameters
        ----------
        for_material : string (case sensitive), optional
            The material name for which the parameters will be printed. 
            The name should match the name in database. If None, prints general
            information about the database.
            The default is None.
        with_new_database : dictionary, optional
            The material parameters. The parameter names should match in the database.
            Use 'print_database()' function from DataBase class for material parameter
            names in the database.
            The default is None.

        Returns
        -------
        None.

        """
        print_infot ="""
NOTE: Currently, you would not be able to update the database on the fly. You have to update
the database in the source code. If new material is needed contact the developers.
        """
        print(print_infot)
            
        
        