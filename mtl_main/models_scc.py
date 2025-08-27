""" This module contains the definition of the MOM-So systems used in the """

BARE_SINGLE_WIRE = {
    'name': 'BARED_SINGLE_WIRE',
    'type': 'scc',
    'note': 'An DE CONTI buried bare-wire. Ref.: Example 6a, p.157, 2023',
    'idx_ref_conductor': 0,
    0: {
        'line_id': 0,
        'conductor_name': 'soil',
        'line_type': 'return',
        'line_return': None,
        'conductivity': 0.01,
        'relative_permeability': 1.0,
        'relative_permittivity': 10.0,
        'relative_permittivity_out': 1.0,
    },
    1: {
        'line_id': 1,
        'conductor_name': 'core',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.0, - 2.0), # Negative numbers denotes cable depth above ground
        'radius': [0, 0.040],
        'conductivity': 1/1.934753E-8,
        'subconductors': None,
        'insulation': None,
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
}

PRYSMIAN_SINGLE_PHASE = {
    'name': 'PRYSMIAN_SINGLE_PHASE',
    'type': 'scc',
    'note': 'An Prysmian 138 kV SCC Cable',
    'idx_ref_conductor': 0,
    0: {
        'line_id': 0,
        'conductor_name': 'soil',
        'line_type': 'return',
        'line_return': None,
        'conductivity': 0.01,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
    },
    1: {
        'line_id': 1,
        'conductor_name': 'core',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.0, -1.7), # Negative numbers denotes cable depth above ground
        'radius': [0, 0.012975],
        'conductivity': 1/1.934753E-8,
        'subconductors': None,
        'insulation': {'name': 'primary_insulation',
                        'center_point': (0.0, -1.7),
                        'thickness': 0.02634,
                        'fourier_order': 0,
                        'relative_permittivity': 2.963538},
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
}

