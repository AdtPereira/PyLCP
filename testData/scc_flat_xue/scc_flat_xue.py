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
    from utils.comsol_data import ComsolPostProcessor
    from plotter.scc_plotter import SCCPlotter
    from models.single_core_cable import SingleCoreCableModelGenerator
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters, PerUnitParameters
    from .plot_config import PLOT_CONFIG
    print("Modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    model = SingleCoreCableModelGenerator(__file__).underground_flat_model()
    mtl_model_a = MulticonductorTransmissionLine(model)
    
    # --- MTL setup ---
    model_b = copy.deepcopy(model)
    model_b[0]['relative_permittivity'] = 20
    mtl_model_b = MulticonductorTransmissionLine(model_b)
    model_c = copy.deepcopy(model)
    model_c[0]['conductivity'] = 0.002 # rho = 500 ohm.m
    mtl_model_c = MulticonductorTransmissionLine(model_c)

    pul_data = {
        'frequencies': np.logspace(4, 7, num=121),
        'comsol': {
            'scenarios': {
                'rho_g_100_epsr1_1_mf': {},
                'rho_g_100_epsr1_20_mf': {},
                'rho_g_500_epsr1_1_mf': {},
            },
        },
        'scenarios': {
            'p100_er1': {
                'mtl': mtl_model_a,
                'zg_form': 'magalhaes_xue',
                'yg_form': 'magalhaes_xue',
            },
            'p100_er20': {
                'mtl': mtl_model_b,
                'zg_form': 'magalhaes_xue',
                'yg_form': 'magalhaes_xue',
            },
            'p500_er1': {
                'mtl': mtl_model_c,
                'zg_form': 'magalhaes_xue',
                'yg_form': 'magalhaes_xue',
            },
            'p100_er1_vance': {
                'mtl': mtl_model_a,
                'zg_form': 'deconti',
                'yg_form': 'vance',
            },
            'p100_er20_vance': {
                'mtl': mtl_model_b,
                'zg_form': 'deconti',
                'yg_form': 'vance'
            },
            'p500_er1_vance': {
                'mtl': mtl_model_c,
                'zg_form': 'deconti',
                'yg_form': 'vance'
            },
            'p100_er1_deconti': {
                'mtl': mtl_model_a,
                'zg_form': 'deconti',
                'yg_form': 'deconti',
            },
            'p100_er20_deconti': {
                'mtl': mtl_model_b,
                'zg_form': 'deconti',
                'yg_form': 'deconti'
            },
            'p500_er1_deconti': {
                'mtl': mtl_model_c,
                'zg_form': 'deconti',
                'yg_form': 'deconti'
            },
        }
    }

    print("Construindo matrizes COMSOL...")
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_params = cmsl_processor.get_general_parameters('cmsl_ground_return_impedance')
    pul = InternalPerUnitParameters(mtl_model_a, cmsl_params['frequencies'])
    internal_matrices = pul.matrices(internal_form='hybrid')
    
    pul_data['comsol'].update(cmsl_params)
    pul_data['comsol']['internal_matrices'] = internal_matrices

    for key, value in pul_data['comsol']['scenarios'].items():
        print(f"  -> Processando COMSOL para: {key}")
        earth_return = cmsl_processor.get_earth_return_parameters(key)
        quasi_tem = cmsl_processor.get_quasi_tem_approx_matrices(internal_matrices, earth_return)
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem
    
    print("\nCalculating internal parameters for all frequencies...")
    pul = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
    internal_matrices = pul.matrices(internal_form='hybrid')
    pul_data['internal_matrices'] = internal_matrices

    print("\nCalculating per-unit-length parameters and quasi-TEM matrices for all scenarios...")
    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        earth_return = pul.earth_return_parameters(value['zg_form'], value['yg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem
    
    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = SCCPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=False)
    plotter.scc_earth_propagation_constant()
    plotter.scc_earth_return_impedance_matrix()
    plotter.scc_earth_return_admittance_matrix()
    plotter.scc_series_impedance_matrix(graph_key_list=['fig419', 'fig421a', 'fig421b'])
    plotter.scc_shunt_admittance_matrix(graph_key_list=['fig423', 'fig425a', 'fig425b'])
    GroundReturnMTLRepresentation(__file__, mtl_model_a, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
    