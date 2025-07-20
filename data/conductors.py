""" This module contains the definition of the MOM-So systems used in the """

LINNET = { 'type': 'ACSR',
    'radius': {'value': [6.74E-3/2, 18.29E-3/2], 'unit': 'm'},
    'resistance_20': {'value': 0.1695, 'unit': 'ohm/km'},
    'resistance_65': {'value': 0.2002, 'unit': 'ohm/km'},
    'conductivity': None,
}

STEEL_3_8 = { 'type': 'Steel',
    'radius': {'value': [0, 9.52E-3/2], 'unit': 'm'},
    'resistance_20': {'value': 3.81, 'unit': 'ohm/km'},
    'resistance_65': {'value': 4.5815, 'unit': 'ohm/km'},
    'conductivity': None,
}
