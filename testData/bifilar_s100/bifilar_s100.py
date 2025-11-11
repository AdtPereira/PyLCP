import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt

# --- Import custom modules ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    from utils.case_utils import *
    from utils.comsol_data import ComsolPostProcessor
    from models.single_core_cable import SingleCoreCableModelGenerator
    from plotter.models_base import BasePlotter
    from mtl_main.graphics import IsolatedMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from mtl_main.graphics import IsolatedMTLRepresentation
    from analytical_forms.isolated_wires import WiresHomogeneousMedia
    from mom_so.quasi_static_green import QuasiStatic
    from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing
    from .plot_config import PLOT_CONFIG
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    model = SingleCoreCableModelGenerator(__file__).isolated_wires()
    mtl = MulticonductorTransmissionLine(model)

    pul_data = {
        'comsol': {
            'scenarios': {
                '1': {},
            },
        },
        'mom_so': {
            'frequencies': np.logspace(0, 5, num=26),
            'scenarios': {
                '1': {'mtl': mtl},
            },
        },
        'analytical': {
            'frequencies': np.logspace(0, 5, num=121),
            'scenarios': {
                '1': {'mtl': mtl},
            },
        },
    }

    print("Importing COMSOL data for cable model...")
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_params = cmsl_processor.get_general_parameters('cmsl_series_impedance_matrix')
    pul_data['comsol'].update(cmsl_params)
    for key, value in pul_data['comsol']['scenarios'].items():
        print(f"  -> Processando COMSOL para: {key}")
        
        series_impedance_terms = cmsl_processor.get_bifilar_data()
        value['partial_impedance_matrix'] = series_impedance_terms['partial_impedance_matrix']
        value['series_impedance_matrix'] = series_impedance_terms['series_impedance_matrix']

    print("\n===  Calculating per-unit-length parameters by Analytical Formulation (Ametani, 2015) ===")
    for key, value in pul_data['analytical']['scenarios'].items():
        print(f"  -> Calculating internal parameters for Model Case {key}...")        
        pul = WiresHomogeneousMedia(value['mtl'], pul_data['analytical']['frequencies'])
        series_parameters = pul.bifilar_series_impedance()
        value['series_impedance_matrix'] = series_parameters['series_impedance_matrix']
        value['high_frequency_limit'] = series_parameters['high_frequency_limit']

    print("\nCalculating per-unit-length parameters by MoM-SO (Patel, 2014)")
    for key, value in pul_data['mom_so']['scenarios'].items():
        print(f"  -> Calculating internal impedance matrix for Model Case {key}...")
        
        green_matrix = QuasiStatic(value['mtl']).green_matrix()
        mom_so = HomogeneousLosslessMedium(value['mtl'], pul_data['mom_so']['frequencies'])
        post_processor = LosslessPostProcessing(value['mtl'])
        z_partial = mom_so.z_partial(green_matrix)

        value['partial_impedance_matrix'] = z_partial
        value['series_impedance_matrix'] = post_processor.z_total(z_partial)

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = BasePlotter(__file__, pul_data, PLOT_CONFIG, autoSave=False)
    plotter.plot_graph(['partial_impedance_matrix'])
    plotter.plot_graph(['series_impedance_matrix'])
    IsolatedMTLRepresentation(__file__, mtl, units='millimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()