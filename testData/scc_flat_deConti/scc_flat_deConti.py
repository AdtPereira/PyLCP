# scc_flat_xue.py

import sys
import os
import time
import copy
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# --- Import custom modules ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    from utils.case_utils import *
    from utils.comsol_data import ComsolDataReader
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
    cmsl_reader = ComsolDataReader(__file__)
    model = SingleCoreCableModelGenerator(__file__).underground_model()
    
    # Define models for different physical scenarios
    flat_model = copy.deepcopy(model)
    mtl_model_a = MulticonductorTransmissionLine(flat_model)
    flat_model[0]['conductivity'] = 0.001
    mtl_model_b = MulticonductorTransmissionLine(flat_model)
    flat_model[0]['conductivity'] = 0.0001
    mtl_model_c = MulticonductorTransmissionLine(flat_model)

    # --- VECTORIZED CALCULATION ---
    pul_data = {
        'comsol': None, # cmsl_reader.data,
        'frequencies': np.logspace(4, 7, num=80),
        'scenarios': {
            'p100': {
                'mtl': mtl_model_a,
                'zg_form': 'magalhaes_xue',
                'yg_form': 'magalhaes_xue'
            },
            'p1000': {
                'mtl': mtl_model_b,
                'zg_form': 'magalhaes_xue',
                'yg_form': 'magalhaes_xue'
            },
            'p10000': {
                'mtl': mtl_model_c,
                'zg_form': 'magalhaes_xue',
                'yg_form': 'magalhaes_xue'
            },
            'p100_deConti': {
                'mtl': mtl_model_a,
                'zg_form': 'deconti',
                'yg_form': 'deconti'
            },
            'p1000_deConti': {
                'mtl': mtl_model_b,
                'zg_form': 'deconti',
                'yg_form': 'deconti'
            },
            'p10000_deConti': {
                'mtl': mtl_model_c,
                'zg_form': 'deconti',
                'yg_form': 'deconti'
            },
        }
    }

    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])        
        earth_return = pul.earth_return_parameters(value['zg_form'], value['yg_form'])
        value['earth_return_parameters'] = earth_return

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = DeContiModels(__file__, pul_data)
    plotter.fig_3(graph_form='norm_and_angle')
    plotter.fig_3(graph_form='resistance_and_inductance')
    plotter.fig_6()
    GroundReturnMTLRepresentation(__file__, mtl_model_a, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()