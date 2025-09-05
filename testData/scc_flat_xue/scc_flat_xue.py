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
except IndexError:
    raise RuntimeError("Could not find project root. Ensure the directory structure is correct.")

# --- Import custom modules ---
try:
    from utils.case_utils import load_json_parameters
    from analyzer.plotter_scc_models import XueModels
    from models import single_core_cables as scc 
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.scc_vector import InternalPerUnitParameters, PerUnitParameters
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    model = scc.three_phase_flat_model(input_json)
    
    # --- Model setup (sem alteração) ---
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
    #    This returns a dictionary of 3D matrices (e.g., shape (40, 6, 6)).
    internal = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
    print("Internal parameters calculated.")

    # 2. Loop through scenarios to calculate ground-return effects.
    for key, value in scenarios.items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        pul_data[key] = pul.quasi_tem_pul(
            internal.internal_matrices(), zg_form=value['zg_form'], yg_form=value['yg_form']
        )
    
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = XueModels(pul_data)
    plotter.plot_fig419()
    plotter.plot_fig421()
    plotter.plot_fig423()
    MTLRepresentation(mtl_model_a, units='centimeter').ground_return_systems()
    plt.show()    

if __name__ == "__main__":
    main()