# scc_flat_xue.py

import sys
import os
import time
import copy
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# --- Configure project root for module imports (unchanged) ---
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
    from utils.comsol_data import ComsolPostProcessor
    from models import overhead_lines
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.overhead_lines import InternalPerUnitParameters
    from plotter.deConti_models import InternalLinesModels
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    
    # --- 1. CREATION OF THE TWO MODELS (CORRECTED LOGIC) ---

    # Solid model (loaded directly from the original JSON)
    print("Loading SOLID model from JSON...")
    solid_model_dict = overhead_lines.single_phase_model(input_json)

    # Equivalent tubular model (created from the solid one)
    print("Creating equivalent TUBULAR model programmatically...")
    tubular_model_dict = copy.deepcopy(solid_model_dict)

    # Set the inner radius to 50% of the outer radius to create the tube
    # (this ratio can be adjusted as desired)
    outer_radius = tubular_model_dict[1]['radius'][1]
    tubular_model_dict[1]['radius'][0] = 0.5 * outer_radius

    # Create the transmission line objects for each model
    mtl_solid = MulticonductorTransmissionLine(solid_model_dict)
    mtl_tubular = MulticonductorTransmissionLine(tubular_model_dict)
    frequencies = np.logspace(0, 8, num=400)

    # Simulation for the solid model
    print("Running simulation for SOLID model...")
    pul_data = {'frequencies': frequencies}
    pul_data['internal'] = InternalPerUnitParameters(mtl_solid, frequencies).all_terms()

    # Simulation for the tubular model
    print("\nRunning simulation for TUBULAR model...")
    pul_data_tubular = {'frequencies': frequencies}
    pul_data_tubular['internal'] = InternalPerUnitParameters(mtl_tubular, frequencies).all_terms()

    print(f"End of simulations! Time spent: {(time.time() - st):.1f} seconds.\n")

    print("Loading COMSOL internal impedance results...")
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_internal = cmsl_processor.get_bare_and_hollow_wire_internal_impedance()

    plotter = InternalLinesModels(__file__, pul_data, pul_data_tubular, mtl_solid, model_tubular=mtl_tubular,
                                  comsol=cmsl_internal, autoSave=True)
    plotter.internal_solid_conductors()
    plotter.internal_impedance()
    plotter.internal_hollow_conductors()
    plotter.hollow_conductor_impedance()
    plotter.nahman_holt_comparison()
    plotter.internal_tubular_characteristics()
    GroundReturnMTLRepresentation(__file__, mtl_solid, units='millimeter').system_schematic(base_filename='system_schematic_solid')
    GroundReturnMTLRepresentation(__file__, mtl_tubular, units='millimeter').system_schematic(base_filename='system_schematic_tubular')
    plt.show()

if __name__ == "__main__":
    main()