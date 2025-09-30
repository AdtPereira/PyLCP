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
    from models import isolated_wires
    from .source import BifilarBareWirePULParameters as BifilarPul
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

# --- Read and process COMSOL data ---
COMSOL_DATA = {}
try:
    # 1. Definir o sufixo do arquivo
    txt_name = 'cmsl_surface_charge_density'
    file_to_read = Path(__file__).parent.parent / case_name / 'Results' / f"{txt_name}.txt"

    # 3. Criar uma instância da classe com o caminho construído
    print(f"--- Testing MergedComsolDataReader ---")
    print(f"Attempting to read file: {file_to_read}")
    reader = MergedComsolDataReader(str(file_to_read))
    df1, df2 = reader.df1, reader.df2
    COMSOL_DATA = {'curve_1': df1, 'curve_2': df2}
    
    print("\n--- Curve 1 Data ---")
    print("Columns identified:", reader.column_names)
    print(df1.info())
    print(df1.head())
    
    print("--- Curve 2 Data ---")
    print(df2.info())
    print(df2.head())
except (FileNotFoundError, ValueError) as e:
    print(f"Warning: COMSOL data file not found. Skipping comparison. Details: {e}")

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=False)
    model = isolated_wires.circular_conductor_wires(input_json, show_model=False)
    
    # --- Model setup ---
    pul = BifilarPul(project_root, case_name, model, SUM_MAX=18, comsol_data=COMSOL_DATA) 
    pul.run_single_fortran()
    pul.run_analytical()
    pul.run_fortran()
    pul.run_mom_methods(autoPlots=True)
    pul.run_mom_so()
    pul.run_srw_rates()
    pul.run_convergence()

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    pul.show_header()    
    pul.plot_resistance_results()
    pul.plot_inductance_results()
    pul.plot_capacitance_results()
    pul.plot_srw_rates()
    pul.plot_generalized_capacitance_convergence()    
    pul.plot_free_space_capacitance_convergence()
    plt.show()  

if __name__ == "__main__":
    main()