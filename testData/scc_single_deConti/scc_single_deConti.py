# scc_flat_xue.py

import sys
import os
import time
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# --- Import custom modules ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    from utils.case_utils import *
    from plotter.deConti_models import DeContiModels
    from models.single_core_cable import SingleCoreCableModelGenerator
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import PerUnitParameters
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    model = SingleCoreCableModelGenerator(__file__).underground_model()
    mtl_model = MulticonductorTransmissionLine(model)

    # --- VECTORIZED CALCULATION ---
    pul_data = {
        'frequencies': np.logspace(-1, 8, num=400),
        'scenarios': {
            'magalhaes_xue': {
                'mtl': mtl_model,
                'zg_form': 'magalhaes_xue'
            },
            'sunde': {
                'mtl': mtl_model,
                'zg_form': 'sunde'
            },
            'pollaczek': {
                'mtl': mtl_model,
                'zg_form': 'pollaczek'
            },
            'ametani': {
                'mtl': mtl_model,
                'zg_form': 'ametani'
            },
            'deconti': {
                'mtl': mtl_model,
                'zg_form': 'deconti'
            },
            'saad': {
                'mtl': mtl_model,
                'zg_form': 'saad'
            },
            'wedepohl': {
                'mtl': mtl_model,
                'zg_form': 'wedepohl'
            },
        }
    }

    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        earth_return = pul.earth_return_parameters(value['zg_form'])
        value['earth_return_parameters'] = earth_return
    
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = DeContiModels(__file__, pul_data)
    plotter.ground_return_impedance()
    GroundReturnMTLRepresentation(__file__, mtl_model, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()