# scc_flat_xue.py

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
    from plotter.xue_models import XueModels
    from models.scc import SingleCoreCableModelGenerator
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
    model = model_generator.generate_underground_model(show_model=True)
    
    # --- Model setup ---
    mtl_model_a = MulticonductorTransmissionLine(model)
    model_b_data = copy.deepcopy(model)
    model_b_data[0]['relative_permittivity'] = 20
    mtl_model_b = MulticonductorTransmissionLine(model_b_data)
    model_c_data = copy.deepcopy(model)
    model_c_data[0]['conductivity'] = 0.002 # rho = 500 Ohm.m
    mtl_model_c = MulticonductorTransmissionLine(model_c_data)

    scenarios = {
        'p100':      {'mtl': mtl_model_a, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p100_er20': {'mtl': mtl_model_b, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p500':      {'mtl': mtl_model_c, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p100_deconti':      {'mtl': mtl_model_a, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        'p100_er20_deconti': {'mtl': mtl_model_b, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        'p500_deconti':      {'mtl': mtl_model_c, 'zg_form': 'deconti', 'yg_form': 'deconti'},
    }
    
    # --- VECTORIZED CALCULATION ---
    pul_data = {'frequencies': np.logspace(3, 7, num=40)}

    # 1. Calculate internal parameters ONCE, as the cable geometry is the same for all scenarios.
    print("Calculating internal parameters for all frequencies...")
    internal = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
    print("Internal parameters calculated.")

    # 2. Loop through scenarios to calculate ground-return effects.
    for key, value in scenarios.items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        pul_data[key] = pul.quasi_tem_approximation(
            internal.internal_matrices(), zg_form=value['zg_form'], yg_form=value['yg_form']
        )
    
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = XueModels(pul_data)
    plotter.plot_fig419()
    plotter.plot_fig421()
    plotter.plot_fig423()
    GroundReturnMTLRepresentation(mtl_model_a, case_name, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()