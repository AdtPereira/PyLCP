import sys
import os
import time
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
    from utils.comsol_data import MergedComsolDataReader
    from models.single_core_cable import SingleCoreCableModelGenerator
    from mtl_main.source import MulticonductorTransmissionLine
    from mtl_main.graphics import IsolatedMTLRepresentation
    from .source import BifilarPULParameters as BifilarPul
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()
    model = SingleCoreCableModelGenerator(__file__).isolated_wires()
    mtl_model = MulticonductorTransmissionLine(model)

    # --- Model setup ---
    pul = BifilarPul(project_root, case_name, model)
    pul.run_analytical()
    pul.run_mom_so()

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    pul.show_header()
    pul.plot_impedance_results()
    pul.plot_partial_impedance_matrix()
    pul.print_impedance_matrix()
    IsolatedMTLRepresentation(__file__, mtl_model, units='millimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
