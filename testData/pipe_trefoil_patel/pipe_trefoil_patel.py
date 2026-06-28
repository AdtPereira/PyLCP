# root_folder/testData/patel_pipe_trefoil/patel_pipe_trefoil.py

import sys
import os
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt

# --- Configure project root for module imports ---
try:
    os.system('cls' if os.name == 'nt' else 'clear')
    # Assuming this script is in root_folder/testData/case_name/
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
    from models import pipe_type 
    from mtl_main.graphics import IsolatedMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """
    Main function to run the simulation and plotting.
    """
    st = time.time()    
    parameters = load_json_parameters(__file__, show_content=True)
    
    # 2. Generate the MODEL dictionary by calling the parametric function
    # The ** operator unpacks the dictionary into keyword arguments
    # MODEL = pipe_type.trefoil_symmetric_model(**parameters, show_model=True)
    MODEL = pipe_type.trefoil_asymmetric_model(**parameters, show_model=True)

    # The rest of the simulation runs as before
    frequency = {'Analytically': np.logspace(4, 7, num=80)}
    mtl_model = MulticonductorTransmissionLine(MODEL)

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    IsolatedMTLRepresentation(__file__, mtl_model, units='millimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()