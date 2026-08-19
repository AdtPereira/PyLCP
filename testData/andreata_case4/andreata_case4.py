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
    duct (see README.md, Limitações).
    """
    for conductor in model.values():
        if isinstance(conductor, dict) and conductor.get('conductor_name') == 'sheath':
            conductor['enclosure'] = None
            conductor['insulation']['thickness'] = thickness
            conductor['insulation']['relative_permittivity'] = relative_permittivity
    return model

def validate_against_case3_reference(pul_data: dict) -> None:
    """
    Confere que o cenário '1' ("Underground", duto ignorado) concorda com a
    referência MATLAB do andreata_case3 -- mesma geometria física de 3 cabos
    SCC + 1 ECC em arranjo plano, sem duto. Não valida o efeito do duto em
    si (que não tem contraparte no case3); serve só como checagem cruzada de
    código para a parte do pipeline que os dois casos compartilham
    (geração do modelo 'flat_scc_with_ecc_cable_model', InternalPerUnitParameters,
    PerUnitParameters heterogêneos).
    """
    case3 = pul_data['matlab']['scenarios'].get('case3_no_duct')
    if case3 is None or case3.get('series_impedance_matrix') is None:
        print("  Aviso: referência MATLAB do andreata_case3 indisponível; validação cruzada ignorada.")
        return

    def _max_rel_error(analytical, matlab):
        analytical, matlab = np.asarray(analytical), np.asarray(matlab)
        scale = max(np.max(np.abs(matlab)), 1e-30)
        return float(np.max(np.abs(analytical - matlab)) / scale)

    quasi_tem = pul_data['scenarios']['1']['quasi_tem_matrices']
    zs_err = _max_rel_error(quasi_tem['series_impedance_matrix'], case3['series_impedance_matrix'])
    ysh_err = _max_rel_error(quasi_tem['shunt_admittance_matrix'], case3['shunt_admittance_matrix'])

    print("\n" + "=" * 70)
    print("  Validação cruzada: cenário '1' (duto ignorado) vs. MATLAB do andreata_case3")
    print("  (mesma geometria física -- 3 cabos SCC + 1 ECC em arranjo plano, sem duto)")
    print("=" * 70)
    print(f"  Zs  (impedância série completa):  erro relativo máx. (norma inf) = {zs_err:.2%}")
    print(f"  Ysh (admitância shunt completa):  erro relativo máx. (norma inf) = {ysh_err:.2%}")
    print("=" * 70 + "\n")

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()
    cable_generator = SingleCoreCableModelGenerator(__file__)

    # --- Modelo base para EquivalentRadiiSystems: 3 SCC em dutos HDPE
    # individuais, homogêneo, sem ECC (mesmo model_0 de andreata_case2) --
    # suficiente porque o efeito dielétrico do duto não depende da presença
    # do ECC dentro dele. EquivalentRadiiSystems exige um dict 'scc' achatado
    # (uma seção transversal só), incompatível com o modelo heterogêneo
    # abaixo -- ver README.md, seção Limitações.
    model_ers_base = cable_generator.flat_hdpe_enclosed_model()
    # 'flat_hdpe_enclosed_model' lê 'type' do JSON (ver Item 0 da padronização
    # de andreata_common); como o único andreata_case4.json descreve o caso
    # como um todo ('scc-flat-hdpe-ecc', para o modelo heterogêneo abaixo),
    # aqui o tipo precisa ser forçado de volta para o homogêneo -- este
    # modelo nunca tem ECC, só serve para EquivalentRadiiSystems.
    model_ers_base['type'] = 'hdpe'
    model_ers_base = apply_semiconducting_layer_correction(model_ers_base, cable_generator.core, cable_generator.sheath)
    mtl_ers_base = MulticonductorTransmissionLine(model_ers_base)

    # --- Modelo físico completo (duto real nas 3 fases + ECC compartilhando
    # o duto da fase C) -- usado só para o esquemático. Nunca passa por
    # InternalPerUnitParameters/PerUnitParameters: a geometria SCC+ECC
    # compartilhando um duto é não-canônica, sem formulação analítica
    # disponível no pyLCP (ver README.md, seção Limitações).
    model_schematic = cable_generator.flat_hdpe_enclosed_with_shared_ecc_model()
    model_schematic = apply_semiconducting_layer_correction(model_schematic, cable_generator.core, cable_generator.sheath)
    mtl_schematic = MulticonductorTransmissionLine(model_schematic)

    # --- Modelo 1: "Underground" -- ignora o duto por completo (idêntico à
    # física de andreata_case3: 3 SCC + ECC, sem duto). ---
    model_1 = cable_generator.flat_scc_with_ecc_cable_model()
    model_1 = apply_semiconducting_layer_correction(model_1, cable_generator.core, cable_generator.sheath)
    mtl_1 = MulticonductorTransmissionLine(model_1)

    # --- Equivalent Radii Systems (ERS) para os parâmetros shunt ---
    print("Calculating equivalent radii systems for shunt parameters...")
    ers = EquivalentRadiiSystems(mtl_ers_base)
    epsr_area = ers.equiv_rel_permittivity_epsr_area_weighted()
    r4, r7 = epsr_area['sheath_outer_radius'], epsr_area['sheath_enclosure_outer_radius']

    # --- Modelo 2: substitui o duto por uma isolação equivalente ponderada
    # por área, aplicada às 3 bainhas SCC do modelo heterogêneo (o ECC nunca
    # é afetado -- ver _override_sheath_insulation). ---
    model_2 = _override_sheath_insulation(
        copy.deepcopy(model_1), thickness=r7 - r4, relative_permittivity=epsr_area['equivalent_relative_permittivity'])
    mtl_2 = MulticonductorTransmissionLine(model_2)

    # --- Modelo 3: GMD case 3.1 -- r0 do isolante equivalente = r5 (raio
    # externo do cabo). ---
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

    print("Carregando dados de referência do MATLAB...")
    matlab_reader = MatlabDataReader(__file__, autoShow=False)
    matlab_data = matlab_reader.get_scc_scenario_data(
        prefix='andreata_case4',
        conductor_order=[0, 3, 1, 4, 2, 5, 6],
    )

    # Referência auxiliar de validação: os dados MATLAB do andreata_case3
    # descrevem o mesmo sistema de 3 cabos SCC + ECC em arranjo plano, sem
    # duto -- devem concordar com o cenário '1' ("Underground", duto
    # ignorado) deste caso. Aponta o MatlabDataReader para o diretório
    # 'Results' do case3 (mesmo prefixo/ordem de condutores que ele próprio
    # usa).
    print("Carregando dados de referência do MATLAB do andreata_case3 (validação, duto ignorado)...")
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

    # Dado COMSOL de impedância interna já disponível (Results/cmsl_internal_
    # impedance_matrix.txt, simulação de 3 condutores -- núcleo/blindagem/ECC
    # -- específica deste caso; ver ComsolPostProcessor.get_scc_internal_
    # impedance_matrix_combined). Retorno à terra e admitância interna ainda
    # não têm arquivo -- o bloco permanece "guardado" para essas duas partes
    # (só emite avisos), mesmo padrão de andreata_case3.
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_processor.load_scc_earth_return_and_internal_scenarios(
        pul_data, internal_mtl_model=mtl_3, internal_form='approximation')

    print("\nCalculating per-unit-length parameters and quasi-TEM matrices for all scenarios...")
    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")

        # Parâmetros internos (núcleo + blindagem + ECC): dependem do modelo
        # de duto (Underground/ERS/GMD), mas não do solo.
        pul_internal = InternalPerUnitParameters(value['mtl'], pul_data['frequencies'])
        internal_matrices = pul_internal.matrices(internal_form='hybrid')
        value['internal_matrices'] = internal_matrices

        # Retorno à terra + matrizes quasi-TEM completas. NOTA: a formulação
        # de retorno à terra (Zg/Yg) não modela o duto de HDPE -- ela enxerga
        # apenas a posição/raio externo de cada cabo, nunca a presença do
        # duto em si. Por isso não existe aqui uma "solução analítica com
        # duto": os três cenários abaixo são todos aproximações, e nenhum
        # deles modela o acoplamento não-canônico ECC-dentro-do-duto (ver
        # README.md, Limitações).
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        earth_return = pul.earth_return_parameters(value['zg_form'], value['yg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # Matriz interna canônica para 'internal_impedance_matrix'/'internal_admittance_matrix'
    # (ver plot_config.py): usa o cenário '3' (GMD case 3.1), o mais completo
    # dos três modelos de duto -- mesma escolha já feita em andreata_case2.
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
