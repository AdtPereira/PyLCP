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
    from models.model_generator import IsolatedModels
    from plotter.models_base import BasePlotter
    from mtl_main.graphics import IsolatedMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters
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
    model = IsolatedModels(__file__).isolated_coaxial_cable()
    mtl = MulticonductorTransmissionLine(model)
    pul_data = mtl.get_pul_data_structure()

    print("Importing COMSOL data for coaxial cable model...")
    cmsl_processor = ComsolPostProcessor(__file__)
    for value in pul_data['comsol']['scenarios'].values():
        value['series_impedance_matrix'] = cmsl_processor.get_coaxial_series_impedance_matrix()
        value['coaxial_cable_parameters'] = cmsl_processor.get_coaxial_cable_parameters()
        value['internal_impedance_matrix'] = cmsl_processor.get_scc_internal_impedance_elements(excitation_type='core')

    print("\nCalculating per-unit-length parameters by Analytical Formulation (Ametani, 2015)")
    for value in pul_data['analytical']['scenarios'].values():        
        freq = pul_data['analytical']['frequencies']
        pul = InternalPerUnitParameters(value['mtl'], freq)
        parameters = pul.parameters_hybrid()['zcs']

        value['series_impedance_matrix'] = parameters['Zcs'][:, np.newaxis, np.newaxis]
        value['coaxial_cable_parameters'] = parameters
        value['internal_impedance_matrix'] = pul.matrices()['impedance_matrix']

    print("\nCalculating per-unit-length parameters by MoM-SO (Patel, 2014)")
    for value in pul_data['mom_so']['scenarios'].values():
        green_matrix = QuasiStatic(value['mtl']).green_matrix()
        mom_so = HomogeneousLosslessMedium(value['mtl'], pul_data['mom_so']['frequencies'])
        post_processor = LosslessPostProcessing(value['mtl'])
        z_partial = mom_so.z_partial(green_matrix)
        z_total = post_processor.z_total(z_partial)

        value['series_impedance_matrix'] = z_total
        value['coaxial_cable_parameters'] = {'Zcs': z_total[:, 0, 0]}
        value['internal_impedance_matrix'] = z_partial

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = BasePlotter(__file__, pul_data, PLOT_CONFIG, autoSave=False)
    graph_list = ['internal_impedance_matrix_js_method']
    # graph_list = ['series_impedance_matrix',
    #               'coaxial_cable_parameters',
    #               'internal_impedance_matrix_js_method',
    #               'internal_impedance_matrix_energy_method',
    #               'internal_impedance_matrix_comparison']
    plotter.plot_graph(graph_list)
    # IsolatedMTLRepresentation(__file__, mtl, units='millimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()