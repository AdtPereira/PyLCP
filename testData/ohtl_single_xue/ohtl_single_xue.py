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
    from utils.comsol_data import ComsolDataReader
    from models import overhead_lines
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.overhead_lines import InternalPerUnitParameters, PerUnitParameters
    from plotter.xue_models import XueModels
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

# --- Load COMSOL Data ---
COMSOL_DATA = {}
try:
    print(f"--- Instanciando ComsolDataReader para o caso '{case_name}' ---")
    reader = ComsolDataReader(project_root, case_name)
    COMSOL_DATA = reader.load_all_results()
    if COMSOL_DATA:
        reader.show_summary()
except FileNotFoundError as e:
    print(f"Aviso: Diretório de dados do COMSOL não encontrado. Detalhes: {e}")
except Exception as e:
    print(f"Ocorreu um erro ao carregar os dados do COMSOL: {e}")

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    
    # Modelo Sólido (carregado diretamente do JSON original)
    print("Loading model from JSON...")
    model = overhead_lines.single_phase_model(input_json)

    # Define models for different physical scenarios
    mtl_model = MulticonductorTransmissionLine(model)
    model_b_data = copy.deepcopy(model)
    model_c_data = copy.deepcopy(model)
    model_b_data[0]['relative_permittivity'] = 20
    model_c_data[0]['conductivity'] = 0.0005 # rho = 2000 Ohm.m
    
    # Define the calculation scenarios
    scenarios = {
        'p100':         {'zg_form': 'nakagawa', 'mtl': mtl_model},
        'p100_carson':  {'zg_form': 'carson',   'mtl': mtl_model},
        'p100_er20':    {'zg_form': 'nakagawa', 'mtl': MulticonductorTransmissionLine(model_b_data)},
        'p2000':        {'zg_form': 'nakagawa', 'mtl': MulticonductorTransmissionLine(model_c_data)},
        'p2000_carson': {'zg_form': 'carson',   'mtl': MulticonductorTransmissionLine(model_c_data)}
    }

    # --- VECTORIZED CALCULATION ---
    numerical_freqs = np.logspace(0, 6, num=51)

    pul_data = {
        'comsol': COMSOL_DATA,
        'frequencies': numerical_freqs}

    # Calculate internal parameters ONCE, as the cable geometry is the same for all scenarios.
    print("Calculating internal parameters for all frequencies...")
    internal = InternalPerUnitParameters(mtl_model, pul_data['frequencies'])
    zi = internal.approximations()['Zi_approx']
    print("Internal parameters calculated.")

    # Loop through scenarios to calculate ground-return effects.
    for key, value in scenarios.items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        pul_data[key] = pul.pul_matrices(zi, zg_form=value['zg_form'])

    print(f"End of simulations! Time spent: {(time.time() - st):.1f} seconds.\n")
    plotter = XueModels(pul_data)
    plotter.overhead_series_impedance_matrix()
    # plotter.plot_fig43()
    # plotter.plot_fig45()
    # plotter.plot_fig46()
    # plotter.plot_fig47()
    # plotter.plot_fig48()
    # GroundReturnMTLRepresentation(mtl_model, case_name, units='millimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()