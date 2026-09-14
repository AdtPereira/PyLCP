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
    from utils.matlab_data import export_comsol_internal_matrices_to_mat
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

    print("Loading MATLAB reference data...")
    matlab_reader = MatlabDataReader(__file__, autoShow=False)
    matlab_data = matlab_reader.get_scc_scenario_data(
        prefix='andreata_case1',
        conductor_order=[0, 3, 1, 4, 2, 5],
    )
    matlab_data['frequencies'] = (
        matlab_data['frequencies'] if matlab_data['frequencies'] is not None else pul_data['frequencies'])
    pul_data['matlab'] = matlab_data

    # MATLAB reference for the modal-domain parameters (Andreata Ch. 5), sent
    # separately from the Z'/Y' dump above -- see MODAL_CH5_DEVELOPMENT.md.
    # `None` if the andreata_case1_{alpham,betam,velocm,Zcm,Ycm}.mat files are
    # not present; the modal plots/report below just skip the overlay then.
    matlab_modal = matlab_reader.get_modal_scenario_data(prefix='andreata_case1', config_index=1)

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
        quasi_tem = pul.propagation_matrices(quasi_tem)  # phase-domain gamma_v, gamma_i, Zc, Yc
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # ------------------------------------------------------------------ #
    # Modal-domain propagation characteristics -- Chapter 5 of Andreata  #
    # (Config. 1: 3 buried SCC -> 6 conductors -> 6 modes).              #
    #                                                                    #
    # Andreata secs. 5.4 / 6.1: the earth return uses the CLOSED-FORM     #
    # EXPRESSIONS of De Conti/Duarte/Alipio (2023) -- eqs. 4.59 (Z'_gjk)  #
    # and 4.63 (P'_gjk) -- NOT the Sommerfeld integrals (magalhaes_xue).  #
    # Scenario: rho = 100 Ohm.m, epsr = 1, zg_form = yg_form = 'deconti'. #
    # ------------------------------------------------------------------ #
    base = pul_data['scenarios']['p100_er1_deconti']

    # Passivity sanity check (assessment only -- Gustavsen 2008, eq. 3).
    passivity = check_pul_passivity(pul_data['frequencies'], base['quasi_tem_matrices'])
    print_passivity_report(passivity, title="Config. 1 -- De Conti closed-form (rho=100)")

    print("\nModal decomposition (Config. 1)...")
    modal = ModalDecomposition(
        pul_data['frequencies'],
        base['quasi_tem_matrices']['series_impedance_matrix'],
        base['quasi_tem_matrices']['shunt_admittance_matrix'],
        conductor_roles=default_scc_roles(
            mtl_model_a.num_sc_cables, mtl_model_a.num_conductors_per_scc),
    )
    modal_data = modal.modal_parameters()
    pul_data['modal'] = modal_data
    _diag = modal_data['diagnostics']
    print(f"  mode labels: {modal_data['mode_labels']}")
    print(f"  max off-diag ratio Zm      : {_diag['max_offdiag_ratio_Zm']:.2e}")
    print(f"  max scalar-relation error  : {_diag['max_scalar_relation_error']:.2e}")
    print(f"  classification similarity  : "
          f"{min(_diag['classification_similarity'].values()):.3f} (min over modes)")
    print(f"  passive (Z', Y')           : {_diag['passivity']['passive']}")
    print_modal_comparison_report(modal_data, matlab_modal)

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")

    # plotter = SCCPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=False)
    # plotter.compare_complete_matrices(
    #     key_list=['self_impedance_phase_a_sheath',
    #               'self_admittance_phase_a_sheath',
    #               'earth_return_impedance_phase_a',
    #               'earth_return_admittance_phase_a',
    #               'earth_return_potential_coeff_phase_a'])
    
    # plotter.compare_internal_matrices(
    #     key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])

    # Export the COMSOL/FEM internal matrices back in the MATLAB reference
    # format (see testData/andreata_common/COMSOL_TO_MATLAB_EXPORT.md).
    # strict_reference=False: the reference .mat is still on the off-standard
    # 90-point grid; the COMSOL sweep is the standard 91 points.
    export_comsol_internal_matrices_to_mat(
        __file__, prefix='andreata_case1', strict_reference=False)

    # Figs. 5.5 / 5.6 / 5.7 -- modal attenuation, phase velocity, |Z_cm| --
    # overlaid with the MATLAB reference (open circles) when available.
    ModalPropagationPlotter(__file__, pul_data['modal'],
                            config_name='Configuration 1', autoSave=True,
                            matlab_modal=matlab_modal).plot_all()

    # GroundReturnMTLRepresentation(__file__, mtl_model_a, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
    