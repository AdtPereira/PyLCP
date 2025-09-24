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
    from models import isolated_wires
    from .source import BifilarCoatedWirePULParameters as BifilarPul
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

# --- Load COMSOL Data ---
COMSOL_DATA = {}
try:
    COMSOL_DATA['shunt_parameters'] = load_comsol_results(__file__, comsol_tag='_shunt_parameters')
    print("COMSOL data loaded successfully.")

    # Display the first few rows of the loaded data to verify
    print("--- Data Head ---")
    print(COMSOL_DATA['shunt_parameters'].head())

    # Display a concise summary of the DataFrame
    print("\n--- DataFrame Info shunt_parameters---")
    COMSOL_DATA['shunt_parameters'].info()    
    
except FileNotFoundError as e:
    print(f"Warning: COMSOL data file not found. Skipping comparison. Details: {e}")

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=False)
    model = isolated_wires.circular_conductor_wires(input_json, show_model=False)
    
    # --- Model setup ---
    pul = BifilarPul(project_root, case_name, model, SUM_MAX=18) 
    # pul.run_single_fortran()
    # pul.run_mom_methods()
    pul.run_srw_rates()
    # pul.run_convergence()

    # print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    # pul.show_header()    
    pul.plot_srw_rates(COMSOL_DATA)
    # pul.plot_capacitance_convergence()    
    plt.show()  

if __name__ == "__main__":
    main()