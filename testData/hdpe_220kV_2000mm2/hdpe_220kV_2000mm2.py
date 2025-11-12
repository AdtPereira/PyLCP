import sys
import os
import copy
import time
import numpy as np
import matplotlib.pyplot as plt

# --- Import custom modules ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    from utils.case_utils import *
    from utils.comsol_data import ComsolPostProcessor
    from plotter.scc_plotter import HDPEPlotter
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from models.single_core_cable import SingleCoreCableModelGenerator
    from .plot_config import PLOT_CONFIG
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """
    Main function to run the simulation and plotting.
    """
    st = time.time()    
    model_generator = SingleCoreCableModelGenerator(__file__)    
    model_0 = model_generator.hdpe_shared_enclosed_model()
    mtl_0 = MulticonductorTransmissionLine(model_0)

    pul_data = {
        'frequencies': np.logspace(0, 6, num=31),
        'comsol': {
            'scenarios': {
                '1': {},
            },
        },
        'scenarios': {},
    }

    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_params = cmsl_processor.get_general_parameters('cmsl_coaxial_cable_impedance')
    pul_data['comsol'].update(cmsl_params)
    for key, value in pul_data['comsol']['scenarios'].items():
        print(f"  -> Processando COMSOL para: {key}")
        value['coaxial_cable_impedance'] = cmsl_processor.get_coaxial_cable_parameters()
        value['internal_impedance_matrix'] = cmsl_processor.get_internal_impedance_matrix()
        value['internal_impedance_elements'] = cmsl_processor.get_internal_impedance_elements()
    
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = HDPEPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=False)
    plotter.hdpe_internal_impedance_matrix()
    plotter.hdpe_internal_impedance_elements()
    GroundReturnMTLRepresentation(__file__, mtl_0, units='millimeter').system_schematic()
    plt.show()
    
if __name__ == "__main__":
    main()