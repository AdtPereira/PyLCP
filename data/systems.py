""" This module contains the definition of the MOM-So systems used in the """

SINGLE_OVERHEAD_DECONTI = [
    {
        'line_id': 0,
        'conductor_name': 'p',
        'line_type': 'active',
        'line_return': None,
        'center_point': (0, 10),
        'radius': [0, 0.01],
        'conductivity': 6.496E7,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    }
]

SINGLE_OVERHEAD_XUE = [
    {
        'line_id': 0,
        'conductor_name': 'p',
        'line_type': 'active',
        'line_return': None,
        'center_point': (0, 10),
        'radius': [0, 0.01],
        'conductivity': 1/1.68E-8,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    }
]

TL_138kV_OVERHEAD = [
    {
        'line_id': 0,
        'conductor_name': 'p',
        'line_type': 'active',
        'line_return': 3,
        'center_point': (0, 10),
        'radius': [0, 0.01],
        'conductivity': 6.496E7,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    },
    {
        'line_id': 1,
        'conductor_name': 'q',
        'line_type': 'active',
        'line_return': 3,
        'center_point': (0, 10),
        'radius': [0, 0.01],
        'conductivity': 6.496E7,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    },
    {
        'line_id': 2,
        'conductor_name': 'r',
        'line_type': 'active',
        'line_return': 3,
        'center_point': (0, 10),
        'radius': [0, 0.01],
        'conductivity': 6.496E7,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    },
    {
        'line_id': 3,
        'conductor_name': 'guard-wire',
        'line_type': 'return',
        'line_return': None,
        'center_point': (0, 0),
        'radius': [0, 0.01],
        'conductivity': 6.496E7,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'soil_relative_permeability': 1,
        'soil_relative_permittivity': 5,
        'soil_conductivity': 1/200,
        'fourier_order': 0,
    },
    {
        'line_id': 3,
        'conductor_name': 'soil',
        'line_type': 'return',
        'line_return': None,
        'center_point': None,
        'radius': None,
        'conductivity': 6.496E7,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    }
]

BIFILAR_S100_NP4 = [
    {
        'line_id': 0,
        'conductor_name': 'p',
        'line_type': 'active',
        'line_return': 1,
        'center_point': (-0.050, 10),
        'radius': [0, 0.01],
        'conductivity': 1/1.68E-8,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    },
    {
        'line_id': 1,
        'conductor_name': 'q',
        'line_type': 'return',
        'line_return': None,
        'center_point': (0.050, 10),
        'radius': [0, 0.01],
        'conductivity': 1/1.68E-8,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    }
]

BIFILAR_S25_NP0 = [
    {
        'line_id': 0,
        'conductor_name': 'p',
        'line_type': 'active',
        'line_return': 1,
        'center_point': (-0.0125, 10),
        'radius': [0, 0.01],
        'conductivity': 1/1.68E-8,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    },
    {
        'line_id': 1,
        'conductor_name': 'q',
        'line_type': 'return',
        'line_return': None,
        'center_point': (0.0125, 10),
        'radius': [0, 0.01],
        'conductivity': 1/1.68E-8,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 0,
    }
]

BIFILAR_S25_NP4 = [
    {
        'line_id': 0,
        'conductor_name': 'p',
        'line_type': 'active',
        'line_return': 1,
        'center_point': (-0.0125, 10),
        'radius': [0, 0.01],
        'conductivity': 1/1.68E-8,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 4,
    },
    {
        'line_id': 1,
        'conductor_name': 'q',
        'line_type': 'return',
        'line_return': None,
        'center_point': (0.0125, 10),
        'radius': [0, 0.01],
        'conductivity': 1/1.68E-8,
        'subconductors': None,
        'insulation': None,
        'relative_permeability': 1,
        'relative_permittivity': 1,
        'relative_permittivity_out': 1,
        'fourier_order': 4,
    }
]

MTL_MODELS = {
    'overhead': {
        'deConti': {'type': 'overhead_conductor', 'data': SINGLE_OVERHEAD_DECONTI},
        'xue': {'type': 'overhead_conductor', 'data': SINGLE_OVERHEAD_XUE},
    },
    'bifilar': {
        25: {
            0: {'type': 'bifilar_wires', 'data': BIFILAR_S25_NP0},
            4: {'type': 'bifilar_wires', 'data': BIFILAR_S25_NP4},
        },
        100: {
            4: {'type': 'bifilar_wires', 'data': BIFILAR_S100_NP4},
        },
    }
}
