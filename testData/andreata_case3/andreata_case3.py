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
    from analytical_forms.modal_analysis import ModalDecomposition, default_scc_roles
    from plotter.modal_plotter import ModalPropagationPlotter, print_modal_comparison_report
    from utils.passivity_check import check_pul_passivity, print_passivity_report
    from .plot_config import PLOT_CONFIG
    print("Modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()
    cable_generator = SingleCoreCableModelGenerator(__file__)
    model = cable_generator.flat_scc_with_ecc_cable_model()
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

    print("Loading MATLAB reference data...")
    matlab_reader = MatlabDataReader(__file__, autoShow=False)
    matlab_data = matlab_reader.get_scc_scenario_data(
        prefix='andreata_case3',
        conductor_order=[0, 3, 1, 4, 2, 5, 6],
    )
    matlab_data['frequencies'] = (
        matlab_data['frequencies'] if matlab_data['frequencies'] is not None else pul_data['frequencies'])
    pul_data['matlab'] = matlab_data

    # MATLAB reference for the modal-domain parameters (Andreata Ch. 5),
    # Configuration 3 (3 SCC + ECC -> 7 conductors -> 7 modes), 90-point
    # grid. MATLAB mode1..7 follow the verified Config. 3 order in
    # utils/matlab_data.py::MODAL_MODE_ORDER -- joined by label with the
    # pyLCP modes below (MODAL_CH5_DEVELOPMENT.md sec. 7.4).
    matlab_modal = matlab_reader.get_modal_scenario_data(prefix='andreata_case3', config_index=3)

    print("\nCalculating internal parameters for all frequencies...")
    pul = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
    internal_matrices = pul.matrices(internal_form='approximation')
    pul_data['internal_matrices'] = internal_matrices

    print("\nCalculating per-unit-length parameters and quasi-TEM matrices for all scenarios...")
    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        earth_return = pul.earth_return_parameters(value['zg_form'], value['yg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)
        quasi_tem = pul.propagation_matrices(quasi_tem)  # phase-domain gamma_v, gamma_i, Zc, Yc
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # ------------------------------------------------------------------ #
    # Modal-domain propagation characteristics -- Chapter 5 of Andreata  #
    # (Config. 3: 3 buried SCC + 1 ECC -> 7 conductors -> 7 modes).      #
    # Same earth-return formulation as Config. 1 (De Conti closed-form,  #
    # secs. 5.4/6.1). The 7 modes are labelled by                         #
    # ModalDecomposition._classify_7c (6 thesis modes + 'ecc').           #
    # ------------------------------------------------------------------ #
    base = pul_data['scenarios']['p100_er1_deconti']

    passivity = check_pul_passivity(pul_data['frequencies'], base['quasi_tem_matrices'])
    print_passivity_report(passivity, title="Config. 3 -- De Conti closed-form (rho=100)")

    print("\nModal decomposition (Config. 3, 3 SCC + ECC)...")
    modal = ModalDecomposition(
        pul_data['frequencies'],
        base['quasi_tem_matrices']['series_impedance_matrix'],
        base['quasi_tem_matrices']['shunt_admittance_matrix'],
        conductor_roles=default_scc_roles(num_cables=3, conductors_per_cable=2, num_ecc=1),
    )
    modal_data = modal.modal_parameters()
    pul_data['modal'] = modal_data
    _diag = modal_data['diagnostics']
    print(f"  mode labels: {modal_data['mode_labels']}")
    print(f"  max off-diag ratio Zm      : {_diag['max_offdiag_ratio_Zm']:.2e}")
    print(f"  max scalar-relation error  : {_diag['max_scalar_relation_error']:.2e}")
    print(f"  passive (Z', Y')           : {_diag['passivity']['passive']}")
    print_modal_comparison_report(modal_data, matlab_modal)

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")

    # plotter = SCCPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=True)
    # plotter.compare_complete_matrices(
    #     key_list=['mutual_impedance_phase_a_sheath_ecc',
    #               'self_impedance_ecc',
    #               'self_admittance_ecc',
    #               'earth_return_impedance_ecc',
    #               'earth_return_admittance_ecc',
    #               'earth_return_potential_coeff_ecc'])
    
    # plotter.compare_internal_matrices(
    #     key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])

    # Modal attenuation, phase velocity, |Z_cm| (7 modes) with the MATLAB
    # reference overlaid.
    ModalPropagationPlotter(__file__, pul_data['modal'],
                            config_name='Configuration 3', autoSave=True,
                            matlab_modal=matlab_modal).plot_all()

    GroundReturnMTLRepresentation(__file__, mtl_model_a, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
    