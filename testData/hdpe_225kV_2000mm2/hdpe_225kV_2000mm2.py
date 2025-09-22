# root_folder/testData/patel_pipe_trefoil/patel_pipe_trefoil.py

import sys
import os
import copy
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
COMSOL_DATA_0 = {}
try:
    COMSOL_DATA_0['core_exc'] = load_comsol_results(__file__, comsol_tag='_core_exc')
    COMSOL_DATA_0['sheath_exc'] = load_comsol_results(__file__, comsol_tag='_sheath_exc')
    COMSOL_DATA_0['core_sheath'] = load_comsol_results(__file__, comsol_tag='_core_sheath')
    COMSOL_DATA_0['shunt_params'] = load_comsol_results(__file__, comsol_tag='_shunt_params')
    print("COMSOL data loaded successfully.")

    # Display the first few rows of the loaded data to verify
    print("--- Data Head ---")
    print(COMSOL_DATA_0['core_exc'].head())

    # Display a concise summary of the DataFrame
    print("\n--- DataFrame Info core_exc---")    
    COMSOL_DATA_0['core_exc'].info()
    print("\n--- DataFrame Info core_sheath---")
    COMSOL_DATA_0['core_sheath'].info() 
    print("\n--- DataFrame Info sheath_exc---")    
    COMSOL_DATA_0['sheath_exc'].info()
    print("\n--- DataFrame Info shunt_params---")    
    COMSOL_DATA_0['shunt_params'].info()
except FileNotFoundError as e:
    print(f"Warning: COMSOL data file not found. Skipping comparison. Details: {e}")

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
            'comsol': COMSOL_DATA_0,
        },
        1: {
            'analytical': {
                "frequencies": analytical_freqs,
                "internal_parameters": pul_1.parameters_hybrid(transition_frequency=1e5),
                "internal_matrices": pul_1.internal_matrices(internal_form='hybrid'),
            },
            'numerical': None,
            'comsol': None,
        },
        2: {
            'analytical': {
                "frequencies": analytical_freqs,
                "internal_parameters": pul_2.parameters_hybrid(transition_frequency=1e5),
                "internal_matrices": pul_2.internal_matrices(internal_form='hybrid'),
            },
            'numerical': None,
            'comsol': None,
        },
        3: {
            'analytical': {
                "frequencies": analytical_freqs,
                "internal_parameters": pul_31.parameters_hybrid(transition_frequency=1e5),
                "internal_matrices": pul_31.internal_matrices(internal_form='hybrid'),
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