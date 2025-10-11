# scc_flat_xue.py

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
import copy
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# --- Configure project root for module imports ---
try:
    os.system('cls' if os.name == 'nt' else 'clear')
    project_root = Path(__file__).resolve().parents[2]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root)) 
    print(f"Project root configured at: {project_root}")
    case_name = os.path.splitext(os.path.basename(__file__))[0]
    print(f"Case name identified as: '{case_name}'")
except IndexError:
    raise RuntimeError("Could not find project root. Ensure the directory structure is correct.")

# --- Import custom modules ---
try:
    from utils.case_utils import *
    from plotter.prysmian_models import PrysmianModels
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
    input_json = load_json_parameters(__file__, show_content=True)
    model_generator = SingleCoreCableModelGenerator(input_json)
    model = model_generator.underground_model(show_model=True)
    
    # --- Model setup ---
    mtl_model = MulticonductorTransmissionLine(model)

    # Define the calculation scenarios
    scenarios = {
        'magalhaes_xue': {'mtl': mtl_model, 'zg_form': 'magalhaes_xue'},
        'sunde': {'mtl': mtl_model, 'zg_form': 'sunde'},
        'pollaczek': {'mtl': mtl_model, 'zg_form': 'pollaczek'},
        'ametani': {'mtl': mtl_model, 'zg_form': 'ametani'},
        'deconti': {'mtl': mtl_model, 'zg_form': 'deconti'},
        'saad': {'mtl': mtl_model, 'zg_form': 'saad'},
    }

    # --- VECTORIZED CALCULATION ---
    pul_data = {'frequencies': np.logspace(1, 7, num=120)}

    # 1. Calculate internal parameters ONCE, as the cable geometry is the same for all scenarios.
    print("Calculating internal parameters for all frequencies...")
    #    This returns a dictionary of 3D matrices (e.g., shape (40, 6, 6)).
    internal = InternalPerUnitParameters(mtl_model, pul_data['frequencies'])
    pul_data['internal'] = internal.parameters_by_bessel()
    print("Internal parameters calculated.")

    # 2. Loop through scenarios to calculate ground-return effects.
    for key, value in scenarios.items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        pul_data[key] = pul.earth_return_parameters(zg_form=value['zg_form'])    
    print("Calculations completed for all scenarios.")

    # --- DISCRETE CALCULATION FOR LOGGING ---
    discrete_frequencies = [1e2, 1e4, 1e5]

    print("\nCalculating discrete points for logging...")
    internal = InternalPerUnitParameters(mtl_model, discrete_frequencies)
    pul = PerUnitParameters(mtl_model, discrete_frequencies)
    pul_data_discrete = pul.quasi_tem_approximation(
        internal.matrices(), zg_form='ametani', yg_form='ametani'
    )
    pul_data_discrete['frequencies'] = discrete_frequencies
    print("Discrete calculation for logging finished.")

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = PrysmianModels(pul_data, pul_data_discrete)
    plotter.core_parameters()
    plotter.sheath_parameters()   
    plotter.core_sheath_internal_impedance_matrix() 
    plotter.core_sheath_internal_parameters()
    plotter.ground_return_impedance()
    plotter.log_matricial_pul_parameters()
    GroundReturnMTLRepresentation(mtl_model, case_name, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()