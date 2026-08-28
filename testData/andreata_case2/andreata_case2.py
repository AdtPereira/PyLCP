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
    from plotter.modal_plotter import ModalPropagationPlotter
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
    notion of a duct layer on their own). Applied identically to all three
    phases, since the duct is identical and concentric on each of them.
    """
    for conductor in model.values():
        if isinstance(conductor, dict) and conductor.get('conductor_name') == 'sheath':
            conductor['enclosure'] = None
            conductor['insulation']['thickness'] = thickness
            conductor['insulation']['relative_permittivity'] = relative_permittivity
    return model

def compare_capacitance(pul_data: dict, c_comsol: dict = None) -> None:
    """
    Print a comparison table of C_11 and C_22 (capacitance matrix diagonal,
    phase A) for the three duct-modeling scenarios plus the COMSOL reference
    (energy method). Both values are frequency-independent, expressed in
    uF/km. Phase A occupies indices [0, 1] of the capacitance matrix
    regardless of the number of phases, since conductor ordering is
    core/sheath per cable (np.kron(identity(N), Zij)).
    """
    rows = {}
    labels = {'1': 'Underground (ignoring HDPE)', '2': 'Area-weighted ERS', '3': 'GMD case 3.1'}
    for key, label in labels.items():
        cap = pul_data['scenarios'][key]['internal_matrices']['capacitance_matrix']
        rows[label] = {'c11': cap[0, 0], 'c22': cap[1, 1]}

    col_w = 26
    print("\n" + "=" * 62)
    print("  Per-Unit-Length Capacitance, phase A [uF/km]")
    print("=" * 62)
    print(f"  {'Model':<{col_w}} {'C_11':>12}  {'C_22':>12}")
    print("-" * 62)
    for label, vals in rows.items():
        print(f"  {label:<{col_w}} {vals['c11'] * 1e9:>12.4f}  {vals['c22'] * 1e9:>12.4f}")
    if c_comsol is not None:
        print("-" * 62)
        print(f"  {'COMSOL (energy method)':<{col_w}} {c_comsol['c11'] * 1e9:>12.4f}  {c_comsol['c22'] * 1e9:>12.4f}")
    print("=" * 62 + "\n")

def validate_against_case1_reference(pul_data: dict) -> None:
    """
    Confere que o cenário '1' ("Underground", duto ignorado) concorda com a
    referência MATLAB do andreata_case1 -- mesma geometria física de 3 cabos
    SCC em arranjo plano, sem duto. Não valida o efeito do duto em si (que
    não tem contraparte no case1); serve só como checagem cruzada de código
    para a parte do pipeline que os dois casos compartilham (geração do
    modelo 'underground_flat_model', InternalPerUnitParameters,
    PerUnitParameters).
    """
    case1 = pul_data['matlab']['scenarios'].get('case1_no_duct')
    if case1 is None or case1.get('series_impedance_matrix') is None:
        print("  Aviso: referência MATLAB do andreata_case1 indisponível; validação cruzada ignorada.")
        return

    def _max_rel_error(analytical, matlab):
        analytical, matlab = np.asarray(analytical), np.asarray(matlab)
        scale = max(np.max(np.abs(matlab)), 1e-30)
        return float(np.max(np.abs(analytical - matlab)) / scale)

    quasi_tem = pul_data['scenarios']['1']['quasi_tem_matrices']
    zs_err = _max_rel_error(quasi_tem['series_impedance_matrix'], case1['series_impedance_matrix'])
    ysh_err = _max_rel_error(quasi_tem['shunt_admittance_matrix'], case1['shunt_admittance_matrix'])

    print("\n" + "=" * 70)
    print("  Validação cruzada: cenário '1' (duto ignorado) vs. MATLAB do andreata_case1")
    print("  (mesma geometria física -- 3 cabos SCC em arranjo plano, sem duto)")
    print("=" * 70)
    print(f"  Zs  (impedância série completa):  erro relativo máx. (norma inf) = {zs_err:.2%}")
    print(f"  Ysh (admitância shunt completa):  erro relativo máx. (norma inf) = {ysh_err:.2%}")
    print("=" * 70 + "\n")

def _validate_fem_hybrid_vs_matlab(pul_data: dict) -> None:
    """Confere o cenário 'fem' (FEM-híbrido: Zi/Yi do COMSOL + retorno pela
    terra analítico com raio do tubo) contra a referência MATLAB do próprio
    andreata_case2 -- que é a solução FEM de Andreata para a Configuração 2.
    Esta é a validação de fato do caminho FEM-híbrido (substitui o GMD de
    Lafaia)."""
    mat = pul_data['matlab']['scenarios'].get('measured', {})
    fem = pul_data['scenarios']['fem']['quasi_tem_matrices']
    er = pul_data['scenarios']['fem']['earth_return_parameters']

    def _err(a, b):
        a, b = np.asarray(a), np.asarray(b)
        return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-30))

    rows = [
        ("Z'  (série completa)", fem.get('series_impedance_matrix'), mat.get('series_impedance_matrix')),
        ("Y'  (shunt completa)", fem.get('shunt_admittance_matrix'), mat.get('shunt_admittance_matrix')),
        ("Zg  (retorno terra)", er.get('impedance_matrix'), mat.get('earth_return_impedance_matrix')),
        ("Pg  (coef. potencial)", er.get('potential_coefficient'), mat.get('earth_return_potential_coefficient_matrix')),
    ]
    print("\n" + "=" * 70)
    print("  Validação: cenário 'fem' (FEM-híbrido) vs. MATLAB (FEM de Andreata, Config. 2)")
    print("=" * 70)
    for label, a, b in rows:
        if a is None or b is None:
            print(f"  {label:24s}  (referência indisponível)")
            continue
        # Zg/Pg do pyLCP são (Nf,3,3) por cabo; MATLAB é (Nf,6,6) por condutor.
        if a.shape[1:] != b.shape[1:]:
            a = a[:, :1, :1]
            b = b[:, :1, :1]
        print(f"  {label:24s}  erro rel. máx. (norma inf) = {_err(a, b):.2%}")
    print("=" * 70 + "\n")


def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()
    cable_generator = SingleCoreCableModelGenerator(__file__)

    # --- Modelo 0: geometria física real (Configuração 2 / Figura 5.2) ---
    # Usado para o esquemático e como base geométrica do ERS/GMD (fornece
    # r4..r7 via EquivalentRadiiSystems). Seu retorno à terra e seus
    # parâmetros internos NÃO são usados diretamente nos cenários -- ver
    # nota de retorno à terra abaixo.
    model_0 = cable_generator.flat_hdpe_enclosed_model()
    model_0 = apply_semiconducting_layer_correction(model_0, cable_generator.core, cable_generator.sheath)
    mtl_0 = MulticonductorTransmissionLine(model_0)

    # --- Modelo 1: "Underground" -- ignora o duto por completo ---
    model_1 = cable_generator.underground_flat_model()
    model_1 = apply_semiconducting_layer_correction(model_1, cable_generator.core, cable_generator.sheath)
    mtl_1 = MulticonductorTransmissionLine(model_1)

    # --- Equivalent Radii Systems (ERS) para os parâmetros shunt ---
    print("Calculating equivalent radii systems for shunt parameters...")
    ers = EquivalentRadiiSystems(mtl_0)
    epsr_area = ers.equiv_rel_permittivity_epsr_area_weighted()
    r4, r7 = epsr_area['sheath_outer_radius'], epsr_area['sheath_enclosure_outer_radius']

    # --- Modelo 2: substitui o duto por uma isolação equivalente ponderada por área ---
    model_2 = _override_sheath_insulation(
        copy.deepcopy(model_0), thickness=r7 - r4, relative_permittivity=epsr_area['equivalent_relative_permittivity'])
    mtl_2 = MulticonductorTransmissionLine(model_2)

    # --- Modelo 3: GMD case 3.1 -- r0 do isolante equivalente = r5 (raio externo do cabo) ---
    print("Calculating GMD-based equivalent parameters for shunt parameters...")
    gmp = ers.equivalent_parameters_from_gmd()
    eps_a = gmp['equivalent_relative_permittivity']['case 3.1']
    model_3 = _override_sheath_insulation(
        copy.deepcopy(model_0), thickness=model_0[2]['insulation']['thickness'], relative_permittivity=eps_a)
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
            # FEM-híbrido: Zi/Yi do COMSOL (geometria excêntrica exata, ar +
            # tubo HDPE) + retorno pela terra analítico com o raio externo do
            # tubo (mtl_0 -> _cable_external_geometry). Sem o GMD de Lafaia.
            'fem': {'mtl': mtl_0, 'internal_source': 'fem',
                    'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        }
    }

    print("Carregando dados de referência do MATLAB...")
    matlab_reader = MatlabDataReader(__file__, autoShow=False)
    matlab_data = matlab_reader.get_scc_scenario_data(
        prefix='andreata_case2',
        conductor_order=[0, 3, 1, 4, 2, 5],
    )

    # Referência auxiliar de validação: os dados MATLAB do andreata_case1
    # descrevem o mesmo sistema de 3 cabos SCC em arranjo plano, sem duto --
    # devem concordar com o cenário '1' ("Underground", duto ignorado) deste
    # caso. Aponta o MatlabDataReader para o diretório 'Results' do case1
    # (mesmo prefixo/ordem de condutores que ele próprio usa).
    print("Carregando dados de referência do MATLAB do andreata_case1 (validação, duto ignorado)...")
    case1_script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'andreata_case1', 'andreata_case1.py')
    case1_matlab_reader = MatlabDataReader(case1_script_path, autoShow=False)
    case1_matlab_data = case1_matlab_reader.get_scc_scenario_data(
        prefix='andreata_case1',
        conductor_order=[0, 3, 1, 4, 2, 5],
    )
    matlab_data['scenarios']['case1_no_duct'] = case1_matlab_data['scenarios']['measured']

    matlab_data['frequencies'] = (
        matlab_data['frequencies'] if matlab_data['frequencies'] is not None else
        case1_matlab_data['frequencies'] if case1_matlab_data['frequencies'] is not None else
        pul_data['frequencies'])
    pul_data['matlab'] = matlab_data

    # Os parâmetros internos usados para compor a matriz quasi-TEM do COMSOL
    # vêm do modelo GMD case 3.1 (a variante mais completa das três) -- só o
    # retorno à terra é, de fato, medido no COMSOL. A impedância interna
    # combinada (segundo bloco de load_scc_earth_return_and_internal_scenarios)
    # é a mesma medição legada de cabo único reaproveitada do estudo
    # monofásico anterior hdpe_300mm2 -- válida aqui porque a seção
    # transversal núcleo+blindagem+duto é idêntica em cada fase; só o retorno
    # à terra muda com o número de fases. Fornece a referência 'measured' de
    # 'internal_impedance_matrix' (ver plot_config.py).
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_processor.load_scc_earth_return_and_internal_scenarios(
        pul_data, internal_mtl_model=mtl_3, internal_form='approximation')

    # Bloco interno FEM (por cabo, 2x2) para o cenário 'fem': Zi já com a
    # indutância externa do ar + tubo (magnetodinâmica), C nodal já com o
    # dielétrico ar + HDPE (método de carga). Ambos da geometria excêntrica
    # exata da Config. 2.
    _zi_fem = cmsl_processor.get_scc_internal_impedance_matrix_combined()
    _c_fem = cmsl_processor.get_scc_internal_capacitance_matrix_combined()
    fem_internal = InternalParametersFromFEM(
        pul_data['frequencies'],
        zi=_zi_fem['scenarios']['measured']['impedance_matrix'],
        zi_frequencies=_zi_fem['frequencies'],
        capacitance=_c_fem['scenarios']['measured']['capacitance_matrix'],
        num_cables=mtl_0.num_sc_cables,
    )

    print("\nCalculating per-unit-length parameters and quasi-TEM matrices for all scenarios...")
    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")

        if value.get('internal_source') == 'fem':
            # Zi/Yi do COMSOL + retorno pela terra analítico (raio do tubo).
            built = build_pul_matrices(
                value['mtl'], pul_data['frequencies'],
                internal_source='fem', fem_internal=fem_internal,
                zg_form=value['zg_form'], yg_form=value['yg_form'])
            value['internal_matrices'] = built['internal_matrices']
            value['earth_return_parameters'] = built['earth_return_parameters']
            value['quasi_tem_matrices'] = built['quasi_tem_matrices']
            continue

        # Cenários 1/2/3: interno analítico (com o modelo de duto bare/ERS/GMD).
        pul_internal = InternalPerUnitParameters(value['mtl'], pul_data['frequencies'])
        internal_matrices = pul_internal.matrices(internal_form='hybrid')
        value['internal_matrices'] = internal_matrices

        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        earth_return = pul.earth_return_parameters(value['zg_form'], value['yg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)
        quasi_tem = pul.propagation_matrices(quasi_tem)  # phase-domain gamma_v, gamma_i, Zc, Yc
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # Validação cruzada Z'/Y' FEM-híbrido vs. MATLAB (referência FEM de Andreata).
    _validate_fem_hybrid_vs_matlab(pul_data)

    # ------------------------------------------------------------------ #
    # Modal-domain propagation characteristics -- Chapter 5 of Andreata  #
    # (Config. 2: 3 SCC em dutos HDPE individuais -> 6 condutores -> 6   #
    # modos). Roda sobre o cenario 'fem' (FEM-híbrido, sem o GMD de       #
    # Lafaia). Figs. 5.8 / 5.9 / 5.10.                                    #
    # ------------------------------------------------------------------ #
    base = pul_data['scenarios']['fem']
    passivity = check_pul_passivity(pul_data['frequencies'], base['quasi_tem_matrices'])
    print_passivity_report(passivity, title="Config. 2 -- cenario FEM-hibrido")

    print("Modal decomposition (Config. 2, FEM-hibrido)...")
    modal = ModalDecomposition(
        pul_data['frequencies'],
        base['quasi_tem_matrices']['series_impedance_matrix'],
        base['quasi_tem_matrices']['shunt_admittance_matrix'],
        conductor_roles=default_scc_roles(mtl_0.num_sc_cables, mtl_0.num_conductors_per_scc),
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

    # Matriz interna canônica para 'internal_impedance_matrix' (ver
    # plot_config.py): agora vem do FEM-híbrido.
    pul_data['internal_matrices'] = pul_data['scenarios']['fem']['internal_matrices']

    c_comsol = None
    if 'cmsl_shunt_params' in cmsl_processor.cmsl_reader.data:
        c_comsol = cmsl_processor.get_shunt_capacitance_elements()
    compare_capacitance(pul_data, c_comsol)
    validate_against_case1_reference(pul_data)

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = SCCPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=True)
    plotter.compare_internal_matrices(
            key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])

    # Comparação FEM-híbrido (pyLCP) vs. MATLAB (FEM de Andreata) -- equivalente
    # aos gráficos do andreata_case1, agora com o duto HDPE modelado de fato.
    plotter.compare_complete_matrices(
        key_list=['self_impedance_phase_a_sheath',
                  'self_admittance_phase_a_sheath',
                  'earth_return_impedance_phase_a',
                  'earth_return_potential_coeff_phase_a'])

    # Figs. 5.8 / 5.9 / 5.10 -- modal attenuation, phase velocity, |Z_cm|
    ModalPropagationPlotter(__file__, pul_data['modal'],
                            config_name='Configuracao 2', autoSave=True).plot_all()

    GroundReturnMTLRepresentation(__file__, mtl_0, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
