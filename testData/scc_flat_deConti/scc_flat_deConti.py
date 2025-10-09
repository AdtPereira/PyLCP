# scc_flat_xue.py

import sys
import os
import time
import copy
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# --- Configure project root for module imports (sem alteração) ---
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
    from plotter.deConti_models import DeContiModels
    from models.scc import SingleCoreCableModelGenerator
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import PerUnitParameters
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    model_generator = SingleCoreCableModelGenerator(input_json)
    model = model_generator.generate_underground_model(show_model=True)
    
    # Define models for different physical scenarios
    flat_model = copy.deepcopy(model)
    mtl_model_a = MulticonductorTransmissionLine(flat_model)
    flat_model[0]['conductivity'] = 0.001
    mtl_model_b = MulticonductorTransmissionLine(flat_model)
    flat_model[0]['conductivity'] = 0.0001
    mtl_model_c = MulticonductorTransmissionLine(flat_model)

    # Define the calculation scenarios
    scenarios = {
        'p100':  {'mtl': mtl_model_a, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p1000': {'mtl': mtl_model_b, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p10000':{'mtl': mtl_model_c, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p100_deConti':  {'mtl': mtl_model_a, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        'p1000_deConti': {'mtl': mtl_model_b, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        'p10000_deConti':{'mtl': mtl_model_c, 'zg_form': 'deconti', 'yg_form': 'deconti'},
    }
    
    # --- VECTORIZED CALCULATION ---
    pul_data = {'frequencies': np.logspace(4, 7, num=80)}

    # 2. Loop through scenarios to calculate ground-return effects.
    for key, value in scenarios.items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        pul_data[key] = pul.earth_return_parameters(
            zg_form=value['zg_form'], yg_form=value['yg_form']
        )
    
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = DeContiModels(mtl_model_a, pul_data)
    plotter.fig_3(graph_form='norm_and_angle')
    plotter.fig_3(graph_form='resistance_and_inductance')
    # plotter.fig_6()
    # GroundReturnMTLRepresentation(mtl_model_a, case_name, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()