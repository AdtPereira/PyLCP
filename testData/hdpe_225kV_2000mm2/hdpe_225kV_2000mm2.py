# root_folder/testData/patel_pipe_trefoil/patel_pipe_trefoil.py

import sys
import os
import copy
import time
import numpy as np
from pathlib import Path
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
    from utils.comsol_data import ComsolDataReader
    from plotter.lafaia_models import LafaiaModels
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from models.hdpe import SingleCoreCableInHDPEModelGenerator
    from analytical_forms.single_core_cable import InternalPerUnitParameters
    from analytical_forms.single_core_cable import EquivalentRadiiSystems
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
    """
    Main function to run the simulation and plotting.
    """
    st = time.time()    
    
    # 1. Load parameters from the sibling .in.json file
    input_json = load_json_parameters(__file__, show_content=True)
    model_generator = SingleCoreCableInHDPEModelGenerator(input_json)
    
    model_0 = model_generator.eccentric_hdpe_enclosed_model(show_model=True)
    model_1 = model_generator.hdpe_ignored_model(show_model=True)
    mtl_0 = MulticonductorTransmissionLine(model_0)
    mtl_1 = MulticonductorTransmissionLine(model_1)

    # Equivalent Radii Systems (ERS) for shunt parameters
    print("Calculating equivalent radii systems for shunt parameters...")
    model_2 = copy.deepcopy(model_0)
    ers = EquivalentRadiiSystems(mtl_0)
    epsr_area = ers.equiv_rel_permittivity_epsr_area_weighted()
    r4, r7 = epsr_area['sheath_outer_radius'], epsr_area['sheath_enclosure_outer_radius']
    
    # Remove enclosure for Model Case 2
    model_2[2]['enclosure'] = None  
    model_2[2]['insulation']['thickness'] = r7 - r4
    model_2[2]['insulation']['relative_permittivity'] = epsr_area['equivalent_relative_permittivity']
    model_generator.show_model(model_2, title='Modified HDPE Model (Case 2)')
    mtl_2 = MulticonductorTransmissionLine(model_2)

    # Case 3.1: The outer radius r0 of the equivalent insulator is taken as the cable 
    # outer radius r5 in Fig. 1 (b) 
    print("Calculating equivalent radii systems for shunt parameters...")
    model_3 = copy.deepcopy(model_0)
    ers = EquivalentRadiiSystems(mtl_0)
    gmp = ers.equivalent_parameters_from_gmd()
    eps_a = gmp['equivalent_relative_permittivity']['case 3.1']
    model_3[2]['enclosure'] = None
    model_3[2]['insulation']['relative_permittivity'] = eps_a
    model_generator.show_model(model_3, title='Modified HDPE Model (Case 3)')
    mtl_3 = MulticonductorTransmissionLine(model_3)

    # --- VECTORIZED CALCULATION ---
    analytical_freqs = np.logspace(0, 6, num=31)

    # Analytical Formulation (Ametani et al., 2015)
    print("Calculating internal parameters for Model Case 1...")
    pul_1 = InternalPerUnitParameters(mtl_1, analytical_freqs) 
    print("Calculating internal parameters for Model Case 2...")
    pul_2 = InternalPerUnitParameters(mtl_2, analytical_freqs)
    print("Calculating internal parameters for Model Case 3.1...")
    pul_31 = InternalPerUnitParameters(mtl_3, analytical_freqs)

    # Populate the pul_data dictionary 
    pul_data = {
        0: {
            'analytical': None,
            'numerical': None,
            'comsol': COMSOL_DATA,
        },
        1: {
            'analytical': {
                "frequencies": analytical_freqs,
                "internal_parameters": pul_1.parameters_hybrid(transition_frequency=1e5),
                "internal_matrices": pul_1.matrices(internal_form='hybrid'),
            },
            'numerical': None,
            'comsol': None,
        },
        2: {
            'analytical': {
                "frequencies": analytical_freqs,
                "internal_parameters": pul_2.parameters_hybrid(transition_frequency=1e5),
                "internal_matrices": pul_2.matrices(internal_form='hybrid'),
            },
            'numerical': None,
            'comsol': None,
        },
        3: {
            'analytical': {
                "frequencies": analytical_freqs,
                "internal_parameters": pul_31.parameters_hybrid(transition_frequency=1e5),
                "internal_matrices": pul_31.matrices(internal_form='hybrid'),
            },
            'numerical': None,
            'comsol': None,
        },
    }

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = LafaiaModels(pul_data, case_name, autoSave=True)
    plotter.internal_impedance_matrix_js_method()
    plotter.internal_impedance_matrix_energy_method()
    plotter.internal_impedance_elements()
    plotter.internal_admittance_elements()

    # 1. Crie uma lista de configurações para cada esquemático
    schematic_configs = [
        {'mtl': mtl_0, 'filename': 'schematic_original_hdpe'},
        {'mtl': mtl_1, 'filename': 'schematic_ignored_hdpe'},
    ]

    # 2. Itere sobre a lista para gerar cada esquemático
    for config in schematic_configs:
        schematic = GroundReturnMTLRepresentation(
            config['mtl'], 
            case_name=case_name, 
            autoSave=True, 
            units='millimeter'
        )
        schematic.system_schematic(base_filename=config['filename'])
    plt.show()
    
if __name__ == "__main__":
    main()