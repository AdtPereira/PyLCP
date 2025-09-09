# scc_flat_xue.py

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
    from models import overhead_lines
    from mtl_main.graphics import MTLRepresentation
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
    
    # --- 1. CRIAÇÃO DOS DOIS MODELOS (LÓGICA CORRIGIDA) ---
    
    # Modelo Sólido (carregado diretamente do JSON original)
    print("Loading SOLID model from JSON...")
    solid_model_dict = overhead_lines.single_phase_model(input_json)
    
    # Modelo Tubular Equivalente (criado a partir do sólido)
    print("Creating equivalent TUBULAR model programmatically...")
    tubular_model_dict = copy.deepcopy(solid_model_dict)
    
    # Define o raio interno como 50% do raio externo para criar o tubo
    # (pode ajustar essa proporção como desejar)
    outer_radius = tubular_model_dict[1]['radius'][1]
    tubular_model_dict[1]['radius'][0] = 0.5 * outer_radius
    
    # Cria os objetos de linha de transmissão para cada modelo
    mtl_solid = MulticonductorTransmissionLine(solid_model_dict)
    mtl_tubular = MulticonductorTransmissionLine(tubular_model_dict)
    frequencies = np.logspace(0, 8, num=400)
    
    # Simulação para o modelo sólido
    print("Running simulation for SOLID model...")
    pul_data = {'frequencies': frequencies}
    pul_data['internal'] = InternalPerUnitParameters(mtl_solid, frequencies).all_terms()

    # Simulação para o modelo tubular
    print("\nRunning simulation for TUBULAR model...")
    pul_data_tubular = {'frequencies': frequencies}
    pul_data_tubular['internal'] = InternalPerUnitParameters(mtl_tubular, frequencies).all_terms()

    print(f"End of simulations! Time spent: {(time.time() - st):.1f} seconds.\n")
    plotter = InternalLinesModels(pul_data, pul_data_tubular, mtl_solid)
    plotter.internal_solid_conductors()
    plotter.internal_impedance()
    plotter.nahman_holt_comparison()
    plotter.internal_tubular_characteristics()
    MTLRepresentation(mtl_solid, units='millimeter').ground_return_systems()
    MTLRepresentation(mtl_tubular, units='millimeter').ground_return_systems()
    plt.show()    

if __name__ == "__main__":
    main()