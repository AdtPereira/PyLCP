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
        InternalPerUnitParameters, PerUnitParameters, EquivalentRadiiSystems,
        apply_semiconducting_layer_correction,
    )
    from .plot_config import PLOT_CONFIG
    print("Modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def _override_sheath_insulation(model: dict, thickness: float, relative_permittivity: float) -> dict:
    """
    Replaces each phase's physical duct ('enclosure') with an equivalent
    homogeneous outer insulation layer on the 'sheath' conductor -- the
    ERS/GMD technique used to fold the HDPE duct's dielectric effect into
    the standard core+sheath internal-parameter formulas (which have no
    notion of a duct layer on their own). Applied identically to the 3 SCC
    phases; the ECC conductor (conductor_name='ecc') is untouched, since
    the duct-modeling scenarios never represent the ECC as being inside a
    duct (see README.md, Limitations).
    """
    for conductor in model.values():
        if isinstance(conductor, dict) and conductor.get('conductor_name') == 'sheath':
            conductor['enclosure'] = None
            conductor['insulation']['thickness'] = thickness
            conductor['insulation']['relative_permittivity'] = relative_permittivity
    return model

def validate_against_case3_reference(pul_data: dict) -> None:
    """
    Checks that scenario '1' ("Underground", duct ignored) agrees with the
    andreata_case3 MATLAB reference -- same physical geometry of 3 SCC cables
    + 1 ECC in a flat arrangement, without a duct. It does not validate the
    duct effect itself (which has no counterpart in case3); it only serves as a
    code cross-check for the part of the pipeline the two cases share
    (generation of the 'flat_scc_with_ecc_cable_model' model,
    InternalPerUnitParameters, heterogeneous PerUnitParameters).
    """
    case3 = pul_data['matlab']['scenarios'].get('case3_no_duct')
    if case3 is None or case3.get('series_impedance_matrix') is None:
        print("  Warning: andreata_case3 MATLAB reference unavailable; cross-validation skipped.")
        return

    def _max_rel_error(analytical, matlab):
        analytical, matlab = np.asarray(analytical), np.asarray(matlab)
        scale = max(np.max(np.abs(matlab)), 1e-30)
        return float(np.max(np.abs(analytical - matlab)) / scale)

    quasi_tem = pul_data['scenarios']['1']['quasi_tem_matrices']
    zs_err = _max_rel_error(quasi_tem['series_impedance_matrix'], case3['series_impedance_matrix'])
    ysh_err = _max_rel_error(quasi_tem['shunt_admittance_matrix'], case3['shunt_admittance_matrix'])

    print("\n" + "=" * 70)
    print("  Cross-validation: scenario '1' (duct ignored) vs. andreata_case3 MATLAB")
    print("  (same physical geometry -- 3 SCC cables + 1 ECC in a flat arrangement, no duct)")
    print("=" * 70)
    print(f"  Zs  (full series impedance):  max relative error (inf norm) = {zs_err:.2%}")
    print(f"  Ysh (full shunt admittance):  max relative error (inf norm) = {ysh_err:.2%}")
    print("=" * 70 + "\n")

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()
    cable_generator = SingleCoreCableModelGenerator(__file__)

    # --- Base model for EquivalentRadiiSystems: 3 SCC in individual HDPE
    # ducts, homogeneous, without ECC (same model_0 as andreata_case2) --
    # sufficient because the duct's dielectric effect does not depend on the
    # presence of the ECC inside it. EquivalentRadiiSystems requires a flat
    # 'scc' dict (a single cross-section), incompatible with the heterogeneous
    # model below -- see README.md, Limitations section.
    model_ers_base = cable_generator.flat_hdpe_enclosed_model()
    # 'flat_hdpe_enclosed_model' reads 'type' from the JSON (see Item 0 of the
    # andreata_common standardization); since the single andreata_case4.json
    # describes the case as a whole ('scc-flat-hdpe-ecc', for the heterogeneous
    # model below), here the type must be forced back to the homogeneous one --
    # this model never has an ECC, it only serves EquivalentRadiiSystems.
    model_ers_base['type'] = 'hdpe'
    model_ers_base = apply_semiconducting_layer_correction(model_ers_base, cable_generator.core, cable_generator.sheath)
    mtl_ers_base = MulticonductorTransmissionLine(model_ers_base)

    # --- Full physical model (real duct on the 3 phases + ECC sharing the
    # phase C duct) -- used only for the schematic. It never goes through
    # InternalPerUnitParameters/PerUnitParameters: the SCC+ECC geometry sharing
    # a duct is non-canonical, with no analytical formulation available in
    # pyLCP (see README.md, Limitations section).
    model_schematic = cable_generator.flat_hdpe_enclosed_with_shared_ecc_model()
    model_schematic = apply_semiconducting_layer_correction(model_schematic, cable_generator.core, cable_generator.sheath)
    mtl_schematic = MulticonductorTransmissionLine(model_schematic)

    # --- Model 1: "Underground" -- ignores the duct entirely (identical to
    # the physics of andreata_case3: 3 SCC + ECC, no duct). ---
    model_1 = cable_generator.flat_scc_with_ecc_cable_model()
    model_1 = apply_semiconducting_layer_correction(model_1, cable_generator.core, cable_generator.sheath)
    mtl_1 = MulticonductorTransmissionLine(model_1)

    # --- Equivalent Radii Systems (ERS) for the shunt parameters ---
    print("Calculating equivalent radii systems for shunt parameters...")
    ers = EquivalentRadiiSystems(mtl_ers_base)
    epsr_area = ers.equiv_rel_permittivity_epsr_area_weighted()
    r4, r7 = epsr_area['sheath_outer_radius'], epsr_area['sheath_enclosure_outer_radius']

    # --- Model 2: replaces the duct with an area-weighted equivalent
    # insulation, applied to the 3 SCC sheaths of the heterogeneous model (the
    # ECC is never affected -- see _override_sheath_insulation). ---
    model_2 = _override_sheath_insulation(
        copy.deepcopy(model_1), thickness=r7 - r4, relative_permittivity=epsr_area['equivalent_relative_permittivity'])
    mtl_2 = MulticonductorTransmissionLine(model_2)

    # --- Model 3: GMD case 3.1 -- r0 of the equivalent insulator = r5 (cable
    # outer radius). ---
    print("Calculating GMD-based equivalent parameters for shunt parameters...")
    gmp = ers.equivalent_parameters_from_gmd()
    eps_a = gmp['equivalent_relative_permittivity']['case 3.1']
    model_3 = _override_sheath_insulation(
        copy.deepcopy(model_1), thickness=model_1[2]['insulation']['thickness'], relative_permittivity=eps_a)
    mtl_3 = MulticonductorTransmissionLine(model_3)

    pul_data = {
        'frequencies': np.logspace(-2, 7, num=90),
        'comsol': {
            'scenarios': {
                'rho_g_100_epsr1_1_mf': {},
            },
        },
        'scenarios': {
            '1': {'mtl': mtl_1, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
            '2': {'mtl': mtl_2, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
            '3': {'mtl': mtl_3, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        }
    }

    print("Loading MATLAB reference data...")
    matlab_reader = MatlabDataReader(__file__, autoShow=False)
    matlab_data = matlab_reader.get_scc_scenario_data(
        prefix='andreata_case4',
        conductor_order=[0, 3, 1, 4, 2, 5, 6],
    )

    # Auxiliary validation reference: the andreata_case3 MATLAB data describes
    # the same system of 3 SCC cables + ECC in a flat arrangement, without a
    # duct -- it must agree with scenario '1' ("Underground", duct ignored) of
    # this case. Point MatlabDataReader at case3's 'Results' directory (same
    # prefix / conductor order it uses itself).
    print("Loading andreata_case3 MATLAB reference data (validation, duct ignored)...")
    case3_script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'andreata_case3', 'andreata_case3.py')
    case3_matlab_reader = MatlabDataReader(case3_script_path, autoShow=False)
    case3_matlab_data = case3_matlab_reader.get_scc_scenario_data(
        prefix='andreata_case3',
        conductor_order=[0, 3, 1, 4, 2, 5, 6],
    )
    matlab_data['scenarios']['case3_no_duct'] = case3_matlab_data['scenarios']['measured']

    matlab_data['frequencies'] = (
        matlab_data['frequencies'] if matlab_data['frequencies'] is not None else
        case3_matlab_data['frequencies'] if case3_matlab_data['frequencies'] is not None else
        pul_data['frequencies'])
    pul_data['matlab'] = matlab_data

    # COMSOL internal impedance data already available (Results/cmsl_internal_
    # impedance_matrix.txt, 3-conductor simulation -- core/sheath/ECC --
    # specific to this case; see ComsolPostProcessor.get_scc_internal_
    # impedance_matrix_combined). Earth return and internal admittance do not
    # have a file yet -- the block stays "guarded" for those two parts (it only
    # emits warnings), same pattern as andreata_case3.
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_processor.load_scc_earth_return_and_internal_scenarios(
        pul_data, internal_mtl_model=mtl_3, internal_form='approximation')

    print("\nCalculating per-unit-length parameters and quasi-TEM matrices for all scenarios...")
    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")

        # Internal parameters (core + sheath + ECC): they depend on the duct
        # model (Underground/ERS/GMD), but not on the soil.
        pul_internal = InternalPerUnitParameters(value['mtl'], pul_data['frequencies'])
        internal_matrices = pul_internal.matrices(internal_form='hybrid')
        value['internal_matrices'] = internal_matrices

        # Earth return + full quasi-TEM matrices. NOTE: the earth-return
        # formulation (Zg/Yg) does not model the HDPE duct -- it sees only the
        # position / outer radius of each cable, never the presence of the duct
        # itself. For that reason there is no "analytical solution with a duct"
        # here: the three scenarios below are all approximations, and none of
        # them models the non-canonical ECC-inside-the-duct coupling (see
        # README.md, Limitations).
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        earth_return = pul.earth_return_parameters(value['zg_form'], value['yg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # Canonical internal matrix for 'internal_impedance_matrix'/'internal_admittance_matrix'
    # (see plot_config.py): uses scenario '3' (GMD case 3.1), the most complete
    # of the three duct models -- same choice already made in andreata_case2.
    pul_data['internal_matrices'] = pul_data['scenarios']['3']['internal_matrices']

    validate_against_case3_reference(pul_data)

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = SCCPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=True)
    # plotter.compare_complete_matrices(
    #     key_list=['core_self_impedance',
    #               'mutual_impedance_core_sheath',
    #               'sheath_self_impedance',
    #               'earth_return_impedance_phase_a',
    #               'self_impedance_ecc',
    #               'self_admittance_ecc',
    #               'earth_return_impedance_ecc',
    #               'earth_return_potential_coeff_ecc'])
    plotter.compare_internal_matrices(
        key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])
    GroundReturnMTLRepresentation(__file__, mtl_schematic, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
