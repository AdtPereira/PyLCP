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
        InternalParametersFromFEM, build_pul_matrices,
    )
    from analytical_forms.modal_analysis import ModalDecomposition, default_scc_roles
    from plotter.modal_plotter import ModalPropagationPlotter, print_modal_comparison_report
    from utils.passivity_check import check_pul_passivity, print_passivity_report
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

def _validate_fem_hybrid_vs_matlab(pul_data: dict) -> None:
    """Checks the 'fem' scenario (FEM-hybrid: Zi/Yi from COMSOL for all 3
    ducts + analytical ground return with the pipe radius) against the
    MATLAB reference of andreata_case4 itself -- which is Andreata's FEM
    solution for Configuration 4. Same pattern as andreata_case2's function
    of the same name."""
    fem = pul_data['scenarios'].get('fem')
    if fem is None:
        print("  Warning: 'fem' scenario unavailable (missing COMSOL data); "
              "FEM-hybrid vs. MATLAB validation skipped.")
        return
    mat = pul_data['matlab']['scenarios'].get('measured', {})
    qt = fem['quasi_tem_matrices']
    er = fem['earth_return_parameters']

    def _err(a, b):
        if a is None or b is None:
            return None
        a, b = np.asarray(a), np.asarray(b)
        # pyLCP's raw Zg/Pg are (Nf, 4, 4) per ground-return OBJECT (A, B, C,
        # ECC); MATLAB's are (Nf, 7, 7) per CONDUCTOR. Reduce both to their
        # [0, 0] scalar (phase A's self ground-return term) when the shapes
        # don't match -- same reconciliation as andreata_case2's function.
        if a.shape[1:] != b.shape[1:]:
            a = a[:, :1, :1]
            b = b[:, :1, :1]
        return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-30))

    rows = [
        ("Z'  (full series)", qt.get('series_impedance_matrix'), mat.get('series_impedance_matrix')),
        ("Y'  (full shunt)", qt.get('shunt_admittance_matrix'), mat.get('shunt_admittance_matrix')),
        ("Zg  (earth return)", er.get('impedance_matrix'), mat.get('earth_return_impedance_matrix')),
        ("Pg  (potential coeff.)", er.get('potential_coefficient'), mat.get('earth_return_potential_coefficient_matrix')),
    ]
    print("\n" + "=" * 70)
    print("  Validation: scenario 'fem' (FEM-hybrid) vs. MATLAB (Andreata's FEM, Config. 4)")
    print("=" * 70)
    for label, a, b in rows:
        err = _err(a, b)
        print(f"  {label:24s}  (reference unavailable)" if err is None else
              f"  {label:24s}  max rel. error (inf norm) = {err:.2%}")
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
    # phase C duct). Used for the schematic, AND as the geometry of the
    # 'fem' scenario below (ground-return position/outer-radius only --
    # `SingleCoreCableWithECCStrategy` already resolves the pipe outer
    # radius via the shared `_cable_external_geometry` helper, Part A of
    # HYBRID_PIPELINE_PLAN.md). It never goes through
    # InternalPerUnitParameters (no *analytical* internal formulation exists
    # for SCC+ECC sharing a duct, see README.md, Limitations) -- the 'fem'
    # scenario instead supplies Zi/Pi externally, from COMSOL.
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
            # De Conti/Duarte/Alipio (2023) closed-form ground return (Andreata
            # secs. 5.4/6.1) on the GMD case 3.1 duct model (the most complete
            # of the three -- same choice as the canonical internal matrix
            # below) -- used as the base scenario for the modal decomposition,
            # matching the earth-return formulation of andreata_case1/2/3.
            '3_deconti': {'mtl': mtl_3, 'zg_form': 'deconti', 'yg_form': 'deconti'},
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

    # MATLAB reference for the modal-domain parameters (Andreata Ch. 5),
    # Configuration 4 (3 SCC + ECC, individual HDPE ducts -> 7 conductors ->
    # 7 modes). Not sent yet by the external developer as of this writing --
    # returns None until the andreata_case4_{alpham,betam,velocm,Zcm,Ycm}.mat
    # files show up (see utils/matlab_data.py::get_modal_scenario_data), at
    # which point the report/overlay below start comparing automatically.
    # Mode *labels* stay unset either way (None): classification
    # (ModalDecomposition._classify) only covers the 6-conductor configs so
    # far -- see MODAL_CH5_DEVELOPMENT.md, known limitation 4.
    matlab_modal = matlab_reader.get_modal_scenario_data(prefix='andreata_case4', config_index=4)

    # COMSOL internal impedance data already available (Results/cmsl_internal_
    # impedance_matrix.txt, 3-conductor simulation -- core/sheath/ECC --
    # specific to this case; see ComsolPostProcessor.get_scc_internal_
    # impedance_matrix_combined). Earth return and internal admittance do not
    # have a file yet -- the block stays "guarded" for those two parts (it only
    # emits warnings), same pattern as andreata_case3.
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_processor.load_scc_earth_return_and_internal_scenarios(
        pul_data, internal_mtl_model=mtl_3, internal_form='approximation')

    # --- 'fem' scenario: FEM-hybrid internal (Zi/Yi from COMSOL, all 3
    # ducts) + analytical ground return (De Conti closed-form, pipe outer
    # radius) -- see HYBRID_PIPELINE_PLAN.md sec. 9: "Config. 4: internal =
    # FEM (3x3 block of the shared pipe + 2x2 for the others)". Phases A/B
    # sit alone in their ducts, so their Zi/Yi are the SAME FEM measurement
    # already used by andreata_case2 (bare SCC inside an HDPE duct, no ECC).
    # Phase C's duct (which also hosts the ECC) has its own dedicated
    # 3-conductor COMSOL model (core/sheath/ECC, all mutually coupled) --
    # andreata_case4's own cmsl_internal_impedance_matrix.txt /
    # cmsl_internal_admittance_charge_method.txt, already loaded above via
    # `cmsl_processor`.
    print("Loading andreata_case2 FEM internal data (bare duct, phases A/B)...")
    case2_script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'andreata_case2', 'andreata_case2.py')
    cmsl_case2 = ComsolPostProcessor(case2_script_path, autoShow=False)
    _zi_ab = cmsl_case2.get_scc_internal_impedance_matrix_combined()
    _c_ab = cmsl_case2.get_scc_internal_capacitance_matrix_combined()
    _zi_c_ecc = cmsl_processor.get_scc_internal_impedance_matrix_combined()
    _c_c_ecc = cmsl_processor.get_scc_internal_capacitance_matrix_combined()

    if _zi_ab is None or _c_ab is None or _zi_c_ecc is None or _c_c_ecc is None:
        print("  Warning: COMSOL internal data missing for the FEM-hybrid 'fem' scenario "
              "(needs andreata_case2's + andreata_case4's own internal impedance/admittance "
              "files) -- 'fem' scenario skipped.")
    else:
        zi_ab = _zi_ab['scenarios']['measured']['impedance_matrix']    # (Nf_ab, 2, 2)
        freq_ab = _zi_ab['frequencies']
        cap_ab = _c_ab['scenarios']['measured']['capacitance_matrix']  # (Nf_ab', 2, 2)

        # Slice the populated [core_C, sheath_C, ECC] sub-block out of
        # case4's own global 7x7 (indices 0/1/6 -- see
        # get_scc_internal_impedance_matrix_combined's docstring); it is
        # dense (real core-ECC/sheath-ECC coupling), unlike the
        # zero-elsewhere block-diagonal assumption used for phases A/B.
        _sel = [0, 1, 6]
        zi_full = _zi_c_ecc['scenarios']['measured']['impedance_matrix']
        zi_c_ecc = zi_full[:, _sel, :][:, :, _sel]                     # (Nf_c, 3, 3)
        freq_c_ecc = _zi_c_ecc['frequencies']
        c_full = _c_c_ecc['scenarios']['measured']['capacitance_matrix']
        cap_c_ecc = c_full[:, _sel, :][:, :, _sel]                     # (Nf_c', 3, 3)

        fem_internal = InternalParametersFromFEM.from_component_blocks(
            pul_data['frequencies'],
            component_blocks=[
                {'zi': zi_ab, 'zi_frequencies': freq_ab, 'capacitance': cap_ab},          # phase A
                {'zi': zi_ab, 'zi_frequencies': freq_ab, 'capacitance': cap_ab},          # phase B (identical duct, no ECC)
                {'zi': zi_c_ecc, 'zi_frequencies': freq_c_ecc, 'capacitance': cap_c_ecc}, # phase C + ECC
            ],
            # Ground-return sees 4 objects (A, B, C, ECC), not 3 -- see
            # InternalParametersFromFEM.from_component_blocks's docstring.
            block_sizes=[2, 2, 2, 1],
        )
        pul_data['scenarios']['fem'] = {
            'mtl': mtl_schematic, 'internal_source': 'fem', 'fem_internal': fem_internal,
            'zg_form': 'deconti', 'yg_form': 'deconti',
        }

    print("\nCalculating per-unit-length parameters and quasi-TEM matrices for all scenarios...")
    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")

        if value.get('internal_source') == 'fem':
            built = build_pul_matrices(
                value['mtl'], pul_data['frequencies'],
                internal_source='fem', fem_internal=value['fem_internal'],
                zg_form=value['zg_form'], yg_form=value['yg_form'])
            value['internal_matrices'] = built['internal_matrices']
            value['earth_return_parameters'] = built['earth_return_parameters']
            value['quasi_tem_matrices'] = built['quasi_tem_matrices']
            continue

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
        quasi_tem = pul.propagation_matrices(quasi_tem)  # phase-domain gamma_v, gamma_i, Zc, Yc
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # Canonical internal matrix for 'internal_impedance_matrix'/'internal_admittance_matrix'
    # (see plot_config.py): uses scenario '3' (GMD case 3.1), the most complete
    # of the three duct models -- same choice already made in andreata_case2.
    pul_data['internal_matrices'] = pul_data['scenarios']['3']['internal_matrices']

    validate_against_case3_reference(pul_data)
    _validate_fem_hybrid_vs_matlab(pul_data)

    # ------------------------------------------------------------------ #
    # Modal-domain propagation characteristics -- Chapter 5 of Andreata  #
    # (Config. 4: 3 SCC + 1 ECC in individual HDPE ducts -> 7 conductors  #
    # -> 7 modes). Runs on the 'fem' scenario (FEM-hybrid, no Lafaia GMD  #
    # -- same choice as andreata_case2), falling back to '3_deconti' (GMD #
    # case 3.1 + De Conti ground return) if the COMSOL data needed for    #
    # 'fem' is unavailable. Mode classification/labelling for this        #
    # 7-conductor (ECC) layout is not implemented yet (see                #
    # MODAL_CH5_DEVELOPMENT.md, known limitation 4) -- ModalDecomposition #
    # still runs and produces alpha_m/v_m/|Z_cm| per mode, just without   #
    # semantic names; the MATLAB overlay/report activate automatically    #
    # once both the reference data and a 7-mode classifier are in place.  #
    # ------------------------------------------------------------------ #
    base_key = 'fem' if 'fem' in pul_data['scenarios'] else '3_deconti'
    base = pul_data['scenarios'][base_key]

    passivity = check_pul_passivity(pul_data['frequencies'], base['quasi_tem_matrices'])
    print_passivity_report(passivity, title=f"Config. 4 -- scenario '{base_key}'")

    print("\nModal decomposition (Config. 4, 3 SCC + ECC, HDPE ducts)...")
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
    #     key_list=['core_self_impedance',
    #               'mutual_impedance_core_sheath',
    #               'sheath_self_impedance',
    #               'earth_return_impedance_phase_a',
    #               'self_impedance_ecc',
    #               'self_admittance_ecc',
    #               'earth_return_impedance_ecc',
    #               'earth_return_potential_coeff_ecc'])
    
    # plotter.compare_internal_matrices(
    #     key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])

    # Modal attenuation, phase velocity, |Z_cm| (7 modes, unlabelled -- see
    # note above) -- overlaid with the MATLAB reference once it is available.
    ModalPropagationPlotter(__file__, pul_data['modal'],
                            config_name='Configuration 4', autoSave=True,
                            matlab_modal=matlab_modal).plot_all()

    GroundReturnMTLRepresentation(__file__, mtl_schematic, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
