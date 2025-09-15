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
    from models import single_core_cables as scc 
    from plotter.patel_models import PatelModels
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters
    from mom_so.quasi_static_green import QuasiStatic
    from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

# --- Load COMSOL Data ---
COMSOL_DATA = {}
try:
    COMSOL_DATA['core_sheath_return'] = load_comsol_results(__file__, comsol_tag='')
    COMSOL_DATA['core_exc'] = load_comsol_results(__file__, comsol_tag='_core_exc')
    COMSOL_DATA['sheath_exc'] = load_comsol_results(__file__, comsol_tag='_sheath_exc')
    print("COMSOL data loaded successfully.")

    # # Display the first few rows of the loaded data to verify
    # print("--- Data Head ---")
    # print(COMSOL_DATA['core_sheath_return'].head())

    # Display a concise summary of the DataFrame
    print("\n--- DataFrame Info core_sheath_return---")
    COMSOL_DATA['core_sheath_return'].info()
    
    print("\n--- DataFrame Info core_exc---")    
    COMSOL_DATA['core_exc'].info()

    print("\n--- DataFrame Info sheath_exc---")
    COMSOL_DATA['sheath_exc'].info()
except FileNotFoundError as e:
    print(f"Warning: COMSOL data file not found. Skipping comparison. Details: {e}")

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    model = scc.isolated_coaxial_cable(input_json, show_model=True)
    
    # --- Model setup (sem alteração) ---
    mtl_model = MulticonductorTransmissionLine(model)

    # --- VECTORIZED CALCULATION ---
    analytical_freqs = np.logspace(0, 6, num=200)
    numerical_freqs = np.logspace(0, 6, num=31)

    # Analytical Formulation (Ametani et al., 2015)
    print("Calculating internal parameters for all frequencies...")
    internal = InternalPerUnitParameters(mtl_model, analytical_freqs)
    
    # MoM-SO formulation (Patel, 2014)
    print("Vectorized numeric routine (MoM-SO)...")
    green_matrix = QuasiStatic(mtl_model).green_matrix()
    mom_so = HomogeneousLosslessMedium(mtl_model, numerical_freqs)
    post_processor = LosslessPostProcessing(mtl_model)
    z_partial_stack = mom_so.z_partial(green_matrix)    # Partial impedance matrix
    zs_stack = post_processor.z_total(z_partial_stack)  # Total series impedance matrix
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")

    # Populate the pul_data dictionary 
    pul_data = {
        'analytical': {
            "frequencies": analytical_freqs,
            "internal_parameters": {
                "bessel": internal.parameters_by_bessel(),
                "approximation": internal.parameters_approximation(),
                "hybrid": internal.parameters_hybrid(transition_frequency=1e5)
            },
            "internal_impedance_matrix": {
                "bessel": internal.internal_matrices(internal_form='bessel')['impedance_matrix'],
                "approximation": internal.internal_matrices(internal_form='approximation')['impedance_matrix']
            }
        },
        'numerical': {
            "frequencies": numerical_freqs,
            "partial_impedance_matrix": z_partial_stack,
            "series_impedance_matrix": zs_stack,
            "series_resistance_matrix": post_processor.rs_matrix(zs_stack),
            "series_inductance_matrix": post_processor.ls_matrix(zs_stack, numerical_freqs)
        },
        'comsol': {
            'core_sheath_return': COMSOL_DATA['core_sheath_return'],
            'core': COMSOL_DATA['core_exc'],
            'sheath': COMSOL_DATA['sheath_exc']
        }
    }

    plotter = PatelModels(pul_data)
    plotter.internal_impedance_elements()
    plotter.internal_impedance_matrix()
    MTLRepresentation(mtl_model, units='millimeter').isolated_coaxial_cables()
    plt.show()    

if __name__ == "__main__":
    main()