import sys
import os
import time
import copy
import numpy as np
import matplotlib.pyplot as plt

# --- Import custom modules ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    from utils.case_utils import *
    from utils.comsol_data import ComsolDataReader
    from plotter.xue_models import XueModels
    from models.single_core_cable import SingleCoreCableModelGenerator
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters, PerUnitParameters
    print("Modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    cmsl_reader = ComsolDataReader(__file__)
    model = SingleCoreCableModelGenerator(__file__).underground_model()
    
    # --- MTL setup ---
    mtl_model_a = MulticonductorTransmissionLine(model)
    model_b_data = copy.deepcopy(model)
    model_b_data[0]['relative_permittivity'] = 20
    mtl_model_b = MulticonductorTransmissionLine(model_b_data)
    model_c_data = copy.deepcopy(model)
    model_c_data[0]['conductivity'] = 0.002 # rho = 500 ohm.m
    mtl_model_c = MulticonductorTransmissionLine(model_c_data)

    scenarios = {
        'p100':      {'mtl': mtl_model_a, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p100_er20': {'mtl': mtl_model_b, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p500':      {'mtl': mtl_model_c, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
    }
    
    pul_data = {
        'comsol': cmsl_reader.data,
        'frequencies': np.logspace(3, 7, num=121)}
    
    print("Calculating internal parameters for all frequencies...")
    internal_pul = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
    pul_data['internal_matrices'] = internal_pul.matrices()
    print("Internal parameters calculated.")

    for key, value in scenarios.items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        pul_data[key] = pul.quasi_tem_approximation(
            pul_data['internal_matrices'], zg_form=value['zg_form'], yg_form=value['yg_form']
        )
    
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = XueModels(__file__, pul_data)
    plotter.plot_fig419()
    plotter.plot_fig421()
    plotter.plot_fig421a()
    plotter.plot_fig423()
    plotter.plot_fig425()
    plotter.plot_fig425a()
    GroundReturnMTLRepresentation(__file__, mtl_model_a, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()