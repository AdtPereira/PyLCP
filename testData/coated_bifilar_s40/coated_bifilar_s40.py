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
    from utils.comsol_data import ComsolDataReader
    from models import isolated_wires
    from .source import BifilarCoatedWirePULParameters as BifilarPul
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
    input_json = load_json_parameters(__file__, show_content=False)
    model = isolated_wires.circular_conductor_wires(input_json, show_model=False)
    
    # --- Model setup ---
    pul = BifilarPul(project_root, case_name, model, SUM_MAX=18) 
    pul.run_single_fortran()
    pul.run_mom_methods()
    pul.run_srw_rates()
    pul.run_convergence()

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    pul.show_header()    
    pul.plot_srw_rates(COMSOL_DATA)
    pul.plot_capacitance_convergence()    
    plt.show()  

if __name__ == "__main__":
    main()