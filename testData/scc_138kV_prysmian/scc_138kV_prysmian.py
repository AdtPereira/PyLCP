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

import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt

# --- Import custom modules ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    from utils.case_utils import *
    from utils.comsol_data import ComsolPostProcessor
    from plotter.scc_models import PrysmianCableModels
    from .plot_config import PLOT_CONFIG
    from models.single_core_cable import SingleCoreCableModelGenerator
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters, PerUnitParameters
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    model = SingleCoreCableModelGenerator(__file__).underground_model()
    mtl_model = MulticonductorTransmissionLine(model)

    print("Loading COMSOL internal impedance results...")
    cmsl_processor = ComsolPostProcessor(__file__)
    scc_internal_cmsl = cmsl_processor.get_scc_internal_impedance_matrix_combined()

    pul_data = {
        'comsol': scc_internal_cmsl,
        'frequencies': np.logspace(0, 7, num=121),
        'logger_data': {
            'frequencies': [1e2, 1e4, 1e5],
            'scenario': 'ametani'
        },
        'scenarios': {
            'magalhaes_xue': {
                'mtl': mtl_model,
                'zg_form': 'magalhaes_xue'
            },
            'sunde': {
                'mtl': mtl_model,
                'zg_form': 'sunde'
            },
            'pollaczek': {
                'mtl': mtl_model,
                'zg_form': 'pollaczek'
            },
            'ametani': {
                'mtl': mtl_model,
                'zg_form': 'ametani'
            },
            'deconti': {
                'mtl': mtl_model,
                'zg_form': 'deconti'
            },
            'saad': {
                'mtl': mtl_model,
                'zg_form': 'saad'
            },
        }
    }

    print("Calculating internal parameters for all frequencies...")
    pul = InternalPerUnitParameters(mtl_model, pul_data['frequencies'])
    internal_matrices = pul.matrices()
    pul_data['internal_matrices'] = internal_matrices
    pul_data['internal_parameters'] = pul.parameters_hybrid()

    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        
        earth_return = pul.earth_return_parameters(value['zg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)

        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # Alias consumed by BasePlotter's generic engine (source='analytical'),
    # used by ground_return_impedance().
    pul_data['analytical'] = {'frequencies': pul_data['frequencies'], 'scenarios': pul_data['scenarios']}

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = PrysmianCableModels(__file__, pul_data, PLOT_CONFIG)
    plotter.internal_impedance_parameters(graph_key='core')
    plotter.internal_impedance_parameters(graph_key='sheath')
    plotter.internal_impedance_parameters(graph_key='core_sheath')
    plotter.internal_impedance_parameters(graph_key='internal_parameters')
    plotter.ground_return_impedance()
    GroundReturnMTLRepresentation(__file__, mtl_model, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()