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

'''
PRYSMIAN_138kV_CORE_SHEATH Cable Data (in meters):
1. CONDUTOR: Corda de cobre tipo circular compacta, de acordo com os requisitos da norma NBR NM 280 (classe 2). Seção nominal: 500 mm2
    Diâmetro nominal: 25,95 mm
2. ENFAIXAMENTO DO CONDUTOR:  Fita semicondutora contendo pó inchante e fita de nylon, ambas aplicadas helicoidalmente sobre o condutor.
    Diâmetro nominal: 26,73 mm
3.BLINDAGEM DO CONDUTOR: Camada extrudada de composto semicondutor à base de XLPE. Espessura nominal: 1,5 mm
    Diâmetro nominal: 29,73 mm
4. ISOLAÇÃO: Camada extrudada de polietileno reticulado (XLPE) Espessura nominal: 13,31 mm
    Diâmetro nominal: 60,35 mm
    Permitividade relativa nominal: 2.3
5. BLINDAGEM DA ISOLAÇÃO: Camada extrudada de composto semicondutor à base de XLPE. Espessura nominal: 1,5 mm
    Diâmetro nominal: 63,35 mm
6. ENFAIXAMENTO DA ISOLAÇÃO: Fita semicondutora contendo pó inchante, aplicada helicoidalmente sobre a blindagem da isolação.
    Diâmetro nominal: 64,63 mm
7. CAPA METÁLICA: Capa extrudada de liga de chumbo. Espessura nominal: 3,00 mm
    Diâmetro nominal: 70,63 mm
    Seção nominal: 637,4 mm2
8. COBERTURA: Camada extrudada de polietileno de alta densidade (HDPE) contendo aditivo de proteção contra térmitas e grafite em pó.
    Espessura nominal: 4,0 mm
    Diâmetro nominal: 78,63 mm

PROPRIEDADES ELÉTRICAS
1. TENSÃO EFICAZ ENTRE FASE E TERRA (kV): 79,69
2. TENSÃO EFICAZ ENTRE FASES (kV): 138
3. NÍVEL BÁSICO DE IMPULSO (NBI) (kV): 650
4. RESISTÊNCIA CC MÁXIMA DO CONDUTOR A 20º C (ohn/km): 0,0366
5. CAPACITÂNCIA (mF/km): 0,1805
'''

PRYSMIAN_138kV_CORE_CABLE = {
    'name': 'PRYSMIAN_SINGLE_PHASE',
    'type': 'scc',
    'note': 'An Prysmian 138 kV SCC Cable',
    'idx_ref_conductor': 0,
    0: {
        'line_id': 0,
        'conductor_name': 'soil',
        'line_type': 'return',
        'line_return': None,
        'conductivity': 0.001,
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
                        'type': 'XLPE',
                        'center_point': (0.0, -1.7),
                        'thickness': 0.02634,
                        'fourier_order': 0,
                        'relative_permittivity': 2.963538,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
}

PRYSMIAN_138kV_CORE_SHEATH = {
    'name': 'PRYSMIAN_SINGLE_PHASE',
    'type': 'scc',
    'note': 'An Prysmian 138 kV SCC Cable with sheath',
    'idx_ref_conductor': 0,
    0: {
        'line_id': 0,
        'conductor_name': 'soil',
        'line_type': 'return',
        'line_return': None,
        'conductivity': 0.001,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
    },
    1: {
        'line_id': 1,
        'conductor_name': 'core',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.0, -1.7),    # Negative numbers denotes cable depth above ground
        'radius': [0, 0.012975],
        'conductivity': 1/1.934753E-8,
        'subconductors': None,
        'insulation': {'name': 'primary_insulation',
                        'type': 'XLPE',
                        'center_point': (0.0, -1.7),
                        'thickness': 0.019340,
                        'fourier_order': 0,
                        'relative_permittivity': 2.963538,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
    2: {
        'line_id': 2,
        'conductor_name': 'sheath',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.0, -1.7),    # Negative numbers denotes cable depth above ground
        'radius': [0.032315, 0.035315],
        'conductivity': 1/2.2E-7,
        'subconductors': None,
        'insulation': {'name': 'secondary_insulation',
                        'type': 'HDPE',
                        'center_point': (0.0, -1.7),
                        'thickness': 0.004,
                        'fourier_order': 0,
                        'relative_permittivity': 2.3,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
}

XUE_FLAT_ARRANGEMENT = {
    'name': 'XUE_FLAT_ARRANGEMENT',
    'type': 'scc',
    'note': 'An XUE flat arrangement cable Fig. 4.18 (a) [Xue 2018]',
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
        'center_point': (0.0, -1.0),    # Negative numbers denotes cable depth above ground
        'radius': [0.0103, 0.019],
        'conductivity': 1/1.7E-8,
        'subconductors': None,
        'insulation': {'name': 'primary_insulation',
                        'type': 'XLPE',
                        'center_point': (0.0, -1.0),
                        'thickness': 0.0345 - 0.019,
                        'fourier_order': 0,
                        'relative_permittivity': 3.5,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
    2: {
        'line_id': 2,
        'conductor_name': 'sheath',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.0, -1.0),    # Negative numbers denotes cable depth above ground
        'radius': [0.0345, 0.0385],
        'conductivity': 1/2.1E-7,
        'subconductors': None,
        'insulation': {'name': 'secondary_insulation',
                        'type': 'HDPE',
                        'center_point': (0.0, -1.0),
                        'thickness': 0.0425 - 0.0385,
                        'fourier_order': 0,
                        'relative_permittivity': 4.0,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
    3: {
        'line_id': 3,
        'conductor_name': 'core',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.350, -1.0),    # Negative numbers denotes cable depth above ground
        'radius': [0.0103, 0.019],
        'conductivity': 1/1.7E-8,
        'subconductors': None,
        'insulation': {'name': 'primary_insulation',
                        'type': 'XLPE',
                        'center_point': (0.350, -1.0),
                        'thickness': 0.0345 - 0.019,
                        'fourier_order': 0,
                        'relative_permittivity': 3.5,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
    4: {
        'line_id': 4,
        'conductor_name': 'sheath',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.350, -1.0),    # Negative numbers denotes cable depth above ground
        'radius': [0.0345, 0.0385],
        'conductivity': 1/2.1E-7,
        'subconductors': None,
        'insulation': {'name': 'secondary_insulation',
                        'type': 'HDPE',
                        'center_point': (0.350, -1.0),
                        'thickness': 0.0425 - 0.0385,
                        'fourier_order': 0,
                        'relative_permittivity': 4.0,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
    5: {
        'line_id': 5,
        'conductor_name': 'core',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.700, -1.0),    # Negative numbers denotes cable depth above ground
        'radius': [0.0103, 0.019],
        'conductivity': 1/1.7E-8,
        'subconductors': None,
        'insulation': {'name': 'primary_insulation',
                        'type': 'XLPE',
                        'center_point': (0.700, -1.0),
                        'thickness': 0.0345 - 0.019,
                        'fourier_order': 0,
                        'relative_permittivity': 3.5,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
    6: {
        'line_id': 6,
        'conductor_name': 'sheath',
        'line_type': 'active',
        'line_return': 0,
        'center_point': (0.700, -1.0),    # Negative numbers denotes cable depth above ground
        'radius': [0.0345, 0.0385],
        'conductivity': 1/2.1E-7,
        'subconductors': None,
        'insulation': {'name': 'secondary_insulation',
                        'type': 'HDPE',
                        'center_point': (0.700, -1.0),
                        'thickness': 0.0425 - 0.0385,
                        'fourier_order': 0,
                        'relative_permittivity': 4.0,
                        'relative_permeability': 1.0,
                    },
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    },
}
