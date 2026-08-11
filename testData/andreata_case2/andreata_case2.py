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

    print("\nCalculating per-unit-length parameters and quasi-TEM matrices for all scenarios...")
    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")

        # Parâmetros internos (núcleo + blindagem): dependem do modelo de
        # duto (bare/ERS/GMD), mas não do solo.
        pul_internal = InternalPerUnitParameters(value['mtl'], pul_data['frequencies'])
        internal_matrices = pul_internal.matrices(internal_form='hybrid')
        value['internal_matrices'] = internal_matrices

        # Retorno à terra + matrizes quasi-TEM completas. NOTA: a formulação
        # de retorno à terra (Zg/Yg) não modela o duto de HDPE -- ela enxerga
        # apenas a posição/raio externo de cada cabo (via d_matrix/D_matrix),
        # nunca a presença do duto em si. Por isso não existe aqui uma
        # "solução analítica com duto": os três cenários abaixo são todos
        # aproximações (a diferença entre eles está apenas no raio externo
        # efetivo usado no termo próprio de imagem -- ver
        # `_cable_distance_matrices` -- que ERS/GMD estendem até r7). A
        # validação de fato do efeito do duto no retorno à terra fica a
        # cargo da comparação com COMSOL/MATLAB.
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        earth_return = pul.earth_return_parameters(value['zg_form'], value['yg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # Matriz interna canônica para 'internal_impedance_matrix' (ver
    # plot_config.py): usa o cenário '3' (GMD case 3.1), o mais completo dos
    # três modelos de duto -- mesma escolha já feita para os parâmetros
    # internos do COMSOL de retorno à terra acima.
    pul_data['internal_matrices'] = pul_data['scenarios']['3']['internal_matrices']

    c_comsol = None
    if 'cmsl_shunt_params' in cmsl_processor.cmsl_reader.data:
        c_comsol = cmsl_processor.get_shunt_capacitance_elements()
    compare_capacitance(pul_data, c_comsol)
    validate_against_case1_reference(pul_data)

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = SCCPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=True)
    # plotter.compare_complete_matrices(
    #     key_list=['core_self_impedance',
    #               'mutual_impedance_core_sheath',
    #               'sheath_self_impedance',
    #               'earth_return_impedance_phase_a'])
    plotter.compare_internal_matrices(
            key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])
    GroundReturnMTLRepresentation(__file__, mtl_0, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
