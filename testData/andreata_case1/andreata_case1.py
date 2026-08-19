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
    from utils.matlab_data import MatlabDataReader
    from plotter.scc_plotter import SCCPlotter
    from models.single_core_cable import SingleCoreCableModelGenerator
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import (
        InternalPerUnitParameters, PerUnitParameters, apply_semiconducting_layer_correction,
    )
    from .plot_config import PLOT_CONFIG
    print("Modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()
    cable_generator = SingleCoreCableModelGenerator(__file__)
    model = cable_generator.underground_flat_model()
    model = apply_semiconducting_layer_correction(model, cable_generator.core, cable_generator.sheath)
    mtl_model_a = MulticonductorTransmissionLine(model)
    
    # --- MTL setup ---
    model_b = copy.deepcopy(model)
    model_b[0]['relative_permittivity'] = 20
    mtl_model_b = MulticonductorTransmissionLine(model_b)
    
    model_c = copy.deepcopy(model)
    model_c[0]['conductivity'] = 0.002 # rho = 500 ohm.m
    mtl_model_c = MulticonductorTransmissionLine(model_c)

    pul_data = {
        'frequencies': np.logspace(-2, 7, num=90),
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

    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_processor.load_scc_earth_return_and_internal_scenarios(
        pul_data, internal_mtl_model=mtl_model_a, internal_form='approximation')

    print("Carregando dados de referência do MATLAB...")
    matlab_reader = MatlabDataReader(__file__, autoShow=False)
    matlab_data = matlab_reader.get_scc_scenario_data(
        prefix='andreata_case1',
        conductor_order=[0, 3, 1, 4, 2, 5],
    )
    matlab_data['frequencies'] = (
        matlab_data['frequencies'] if matlab_data['frequencies'] is not None else pul_data['frequencies'])
    pul_data['matlab'] = matlab_data

    print("\nCalculating internal parameters for all frequencies...")
    pul = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
    # internal_matrices = pul.matrices(internal_form='approximation')
    internal_matrices = pul.matrices()
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
    plotter.compare_complete_matrices(
        key_list=['self_impedance_phase_a_sheath',
                  'self_admittance_phase_a_sheath',
                  'earth_return_impedance_phase_a',
                  'earth_return_admittance_phase_a',
                  'earth_return_potential_coeff_phase_a'])
    plotter.compare_internal_matrices(
        key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])
    GroundReturnMTLRepresentation(__file__, mtl_model_a, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
    