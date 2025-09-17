import sys
import os
import time
from pathlib import Path
import matplotlib.pyplot as plt

# --- Configure project root for module imports  ---
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
    from .source import BifilarBareWirePULParameters as BifilarPul
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    model = isolated_wires.circular_conductor_wires(input_json, show_model=True)
    
    # --- Model setup  ---
    pul = BifilarPul(project_root, case_name, model, SUM_MAX=18) 
    pul.run_single_fortran()
    pul.run_analytical()
    pul.run_fortran()
    pul.run_mom_methods(autoPlots=True)

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    pul.show_header()    
    plt.show()  

if __name__ == "__main__":
    main()