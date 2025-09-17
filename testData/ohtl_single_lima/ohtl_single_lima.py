# ohtl_single_xue.py

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
    from models import overhead_lines
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.overhead_lines import InternalPerUnitParameters, PerUnitParameters
    from plotter.lima_models import LimaModels
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    
    # Modelo Sólido (carregado diretamente do JSON original)
    print("Loading model from JSON...")
    model = overhead_lines.single_phase_model(input_json)

    # Define models for different physical scenarios
    model_a = copy.deepcopy(model)
    model_b = copy.deepcopy(model)
    model_c = copy.deepcopy(model)
    
    model_b[0]['conductivity'] = 1e+6 # rho = 1 uOhm.m
    model_c[0]['conductivity'] = 2e-4 # rho = 5000 Ohm.m
    
    mtl_model_a = MulticonductorTransmissionLine(model_a)
    mtl_model_b = MulticonductorTransmissionLine(model_b)
    mtl_model_c = MulticonductorTransmissionLine(model_c)

    # Define the calculation scenarios
    scenarios = {
        'carson':           {'mtl': mtl_model_a, 'zg_form': 'carson'},
        'sunde':            {'mtl': mtl_model_a, 'zg_form': 'sunde'},
        'nakagawa':         {'mtl': mtl_model_a, 'zg_form': 'nakagawa'},
        'quasi_tem':        {'mtl': mtl_model_a, 'zg_form': 'quasi_tem'},
        'quasi_tem_log':    {'mtl': mtl_model_a, 'zg_form': 'quasi_tem_log'},
        'nakagawa_1u':      {'mtl': mtl_model_b, 'zg_form': 'nakagawa'},
        'nakagawa_5k' :     {'mtl': mtl_model_c, 'zg_form': 'nakagawa'},
        'quasi_tem_1u':     {'mtl': mtl_model_b, 'zg_form': 'quasi_tem'},
        'quasi_tem_5k':     {'mtl': mtl_model_c, 'zg_form': 'quasi_tem'},
    }

    pul_data = {'frequencies': np.logspace(0, 10, num=200)}

    # Calculate internal parameters ONCE, as the cable geometry is the same for all scenarios.
    print("Calculating internal parameters for all frequencies...")
    internal = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
    zi = internal.approximations()['Zi_approx']
    print("Internal parameters calculated.")

    # Loop through scenarios to calculate ground-return effects.
    for key, value in scenarios.items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        pul_data[key] = pul.pul_matrices(zi, zg_form=value['zg_form'])

    print(f"End of simulations! Time spent: {(time.time() - st):.1f} seconds.\n")
    plotter = LimaModels(pul_data)
    plotter.propagation_constant_forms()
    plotter.propagation_constant_nakagawa()
    plotter.propagation_constant_quasi_tem()
    plotter.characteristic_impedance_matrix()
    plotter.characteristic_impedance_nakagawa()
    plotter.characteristic_impedance_quasi_tem()
    GroundReturnMTLRepresentation(mtl_model_a, case_name, units='millimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()