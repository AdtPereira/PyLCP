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
    from utils.case_utils import *
    from models import overhead_lines
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.overhead_lines import PerUnitParameters, InternalPerUnitParameters
    from analyzer.plotter_deConti_models import DeContiModels
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    print("Loading model from JSON...")
    model = overhead_lines.single_phase_model(input_json)
    
    # Cria os objetos de linha de transmissão para cada modelo
    mtl_model = MulticonductorTransmissionLine(model)
    discrete_frequencies = np.array([60, 100e3, 1e6])
    
    print("Calculating internal parameters for all frequencies...")
    internal = InternalPerUnitParameters(mtl_model, discrete_frequencies)
    zi = internal.approximations()['Zi_approx']
    print("Internal parameters calculated.")

    print("Calculating PUL parameters for all frequencies...")
    pul = PerUnitParameters(mtl_model, discrete_frequencies)
    pul_data = pul.pul_matrices(zi, zg_form='deri')
    pul_data['frequencies'] = discrete_frequencies
    print("PUL parameters calculated.")

    print(f"End of simulations! Time spent: {(time.time() - st):.1f} seconds.\n")
    plotter = DeContiModels(pul_data)
    plotter.log_matricial_pul_parameters(pul_data)
    MTLRepresentation(mtl_model, units='millimeter').ground_return_systems()
    plt.show()    

if __name__ == "__main__":
    main()