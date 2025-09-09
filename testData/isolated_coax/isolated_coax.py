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
    from plotter.patel_models import PatelModels
    from models import single_core_cables as scc 
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters
    from mom_so.quasi_static_green import QuasiStatic
    from mom_so.lossless_medium_vector import HomogeneousLosslessMedium, LosslessPostProcessing
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    model = scc.isolated_coaxial_cable(input_json, show_model=True)
    
    # --- Model setup (sem alteração) ---
    mtl_model = MulticonductorTransmissionLine(model)

    # --- VECTORIZED CALCULATION ---
    # Define the frequency arrays
    analytical_freqs = np.logspace(0, 6, num=200)
    numerical_freqs = np.logspace(0, 6, num=40)
    
    pul_data = {
        'frequencies': {
            'analytical': analytical_freqs,
            'numerical': numerical_freqs
        }
    }

    # 1. Calculate internal parameters ONCE, as the cable geometry is the same for all scenarios.
    print("Calculating internal parameters for all frequencies...")
    internal = InternalPerUnitParameters(mtl_model, analytical_freqs)
    pul_data['internal_parameters'] = internal.parameters_by_bessel()
    pul_data['internal_parameters'] = internal.parameters_approximation()
    print("Internal parameters calculated.")

    # 2. Perform numerical calculation for all frequencies at once.
    print("Iniciando rotina numérica vetorizada (MoM-SO)...")
    
    # The Green's matrix is frequency-independent
    green_matrix = QuasiStatic(mtl_model).green_matrix()
    
    # Instantiate the model ONCE with the full array of numerical frequencies
    mom_so = HomogeneousLosslessMedium(mtl_model, numerical_freqs)
    post_processor = LosslessPostProcessing(mtl_model)

    # Calculate partial impedance for all frequencies
    z_partial_stack = mom_so.z_partial(green_matrix)
    
    # Calculate total series impedance for all frequencies
    zs_stack = post_processor.z_total(z_partial_stack)
    
    # Calculate series resistance and inductance for all frequencies
    rs_stack = post_processor.rs_matrix(zs_stack)
    ls_stack = post_processor.ls_matrix(zs_stack, numerical_freqs)

    # 3. Populate the pul_data dictionary in the expected format for the plotter
    pul_data['numerical'] = {}
    for i, freq in enumerate(numerical_freqs):
        pul_data['numerical'][freq] = {
            'zs': zs_stack[i],
            'rs': rs_stack[i],
            'ls': ls_stack[i]
        }

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = PatelModels(pul_data)
    plotter.series_impedance_matrix()
    MTLRepresentation(mtl_model, units='millimeter').isolated_coaxial_cables()
    plt.show()    

if __name__ == "__main__":
    main()