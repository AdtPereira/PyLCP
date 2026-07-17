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
    model = cable_generator.underground_model()
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
            # 'p100_er20': {
            #     'mtl': mtl_model_b,
            #     'zg_form': 'magalhaes_xue',
            #     'yg_form': 'magalhaes_xue',
            # },
            # 'p500_er1': {
            #     'mtl': mtl_model_c,
            #     'zg_form': 'magalhaes_xue',
            #     'yg_form': 'magalhaes_xue',
            # },
            # 'p100_er1_vance': {
            #     'mtl': mtl_model_a,
            #     'zg_form': 'deconti',
            #     'yg_form': 'vance',
            # },
            # 'p100_er20_vance': {
            #     'mtl': mtl_model_b,
            #     'zg_form': 'deconti',
            #     'yg_form': 'vance'
            # },
            # 'p500_er1_vance': {
            #     'mtl': mtl_model_c,
            #     'zg_form': 'deconti',
            #     'yg_form': 'vance'
            # },
            # 'p100_er1_deconti': {
            #     'mtl': mtl_model_a,
            #     'zg_form': 'deconti',
            #     'yg_form': 'deconti',
            # },
            # 'p100_er20_deconti': {
            #     'mtl': mtl_model_b,
            #     'zg_form': 'deconti',
            #     'yg_form': 'deconti'
            # },
            # 'p500_er1_deconti': {
            #     'mtl': mtl_model_c,
            #     'zg_form': 'deconti',
            #     'yg_form': 'deconti'
            # },
        }
    }

    print("Construindo matrizes COMSOL...")
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_params = cmsl_processor.get_general_parameters('cmsl_ground_return_impedance')
    if cmsl_params is not None:
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
    else:
        print("  Aviso: Processamento COMSOL ignorado (dados não disponíveis).")
        pul_data['comsol'] = {}

    print("Carregando dados de referência do MATLAB...")
    matlab_reader = MatlabDataReader(__file__, autoShow=False)
    matlab_freq = matlab_reader.data.get('andreata_frequency_range')

    # O MATLAB exporta os condutores agrupados por tipo: [core_A, core_B, core_C,
    # sheath_A, sheath_B, sheath_C]. O pyLCP monta suas matrizes (interna e
    # quasi-TEM) agrupadas por cabo: [core_A, sheath_A, core_B, sheath_B,
    # core_C, sheath_C] (ver np.kron(np.identity(N), Zij) em
    # InternalPerUnitParameters.matrices()). Sem essa reordenação, M[p,q]
    # (pyLCP) e Z[p,q] (MATLAB) apontam para pares de condutores fisicamente
    # diferentes para os mesmos índices (p, q).
    matlab_to_pylcp_order = [0, 3, 1, 4, 2, 5]

    def _reorder_matlab_matrix(matrix):
        if matrix is None:
            return None
        return matrix[:, matlab_to_pylcp_order, :][:, :, matlab_to_pylcp_order]

    matlab_internal_z = _reorder_matlab_matrix(matlab_reader.data.get('andreata_internal_impedance_matrix'))
    matlab_series_z = _reorder_matlab_matrix(matlab_reader.data.get('andreata_series_impedance_matrix'))
    matlab_shunt_y = _reorder_matlab_matrix(matlab_reader.data.get('andreata_shunt_admittance_matrix'))
    matlab_internal_y = _reorder_matlab_matrix(matlab_reader.data.get('andreata_internal_admittance_matrix'))

    # A matriz de coeficiente de potencial de retorno à terra (Pg) só existe a nível de
    # cabo/fase (3x3: A, B, C) -- ela não distingue núcleo de bainha, já que apenas o
    # condutor mais externo de cada cabo "enxerga" o retorno pela terra. O MATLAB a
    # exporta como um 6x6 redundante (bloco 2x2 de Pg repetido), por isso basta extrair
    # o bloco 3x3 superior-esquerdo -- sem necessidade da permutação núcleo/bainha.
    matlab_earth_return_pg_full = matlab_reader.data.get('andreata_earth_return_potential_coefficient_matrix')
    matlab_earth_return_pg = matlab_earth_return_pg_full[:, :3, :3] if matlab_earth_return_pg_full is not None else None

    # Mesmo raciocínio para a impedância de retorno à terra (Zg): 3x3 por fase,
    # exportada pelo MATLAB como 6x6 redundante -- basta o bloco 3x3 superior-esquerdo.
    matlab_earth_return_zg_full = matlab_reader.data.get('andreata_earth_return_impedance_matrix')
    matlab_earth_return_zg = matlab_earth_return_zg_full[:, :3, :3] if matlab_earth_return_zg_full is not None else None

    pul_data['matlab'] = {
        'frequencies': matlab_freq.flatten() if matlab_freq is not None else pul_data['frequencies'],
        'scenarios': {
            'measured': {
                'internal_impedance_matrix': matlab_internal_z,
                'series_impedance_matrix': matlab_series_z,
                'shunt_admittance_matrix': matlab_shunt_y,
                'internal_admittance_matrix': matlab_internal_y,
                'earth_return_potential_coefficient_matrix': matlab_earth_return_pg,
                'earth_return_impedance_matrix': matlab_earth_return_zg,
            },
        },
    }

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
        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = SCCPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=True)
    plotter.scc_series_impedance_matrix(graph_key_list=['series_impedance_all_scenarios'])
    plotter.scc_shunt_admittance_matrix(graph_key_list=['shunt_admittance_all_scenarios'])
    plotter.scc_series_impedance_internal_vs_matlab('internal_impedance_core_sheath')
    plotter.scc_shunt_admittance_internal_vs_matlab('internal_admittance_core_sheath')
    plotter.scc_earth_return_potential_coefficient_matrix()
    plotter.scc_earth_return_impedance_matrix()
    GroundReturnMTLRepresentation(__file__, mtl_model_a, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()
    