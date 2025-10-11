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
    from models.single_core_cable import SingleCoreCableModelGenerator
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
    model = model_generator.underground_model(show_model=True)
    
    # # --- Model setup (sem alteração) ---
    mtl_model = MulticonductorTransmissionLine(model)

    # Define the calculation scenarios
    scenarios = {
        'magalhaes_xue': {'mtl': mtl_model, 'zg_form': 'magalhaes_xue'},
        'sunde': {'mtl': mtl_model, 'zg_form': 'sunde'},
        'pollaczek': {'mtl': mtl_model, 'zg_form': 'pollaczek'},
        'ametani': {'mtl': mtl_model, 'zg_form': 'ametani'},
        'deconti': {'mtl': mtl_model, 'zg_form': 'deconti'},
        'saad': {'mtl': mtl_model, 'zg_form': 'saad'},
        'wedepohl': {'mtl': mtl_model, 'zg_form': 'wedepohl'},
    }
    
    # --- VECTORIZED CALCULATION ---
    pul_data = {'frequencies': np.logspace(-1, 8, num=400)}

    # 2. Loop through scenarios to calculate ground-return effects.
    for key, value in scenarios.items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        pul_data[key] = pul.earth_return_parameters(zg_form=value['zg_form'])
    
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = DeContiModels(pul_data)
    plotter.ground_return_impedance()
    GroundReturnMTLRepresentation(mtl_model, case_name, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()