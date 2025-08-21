""" 
Este script executa simulações de impedância de linhas de transmissão coaxiais
usando tanto uma abordagem analítica (formulação de Patel) quanto uma abordagem numérica
(Método dos Momentos - MoM-SO). Ele gera gráficos comparativos dos resultados
obtidos por ambas as metodologias.

Este arquivo é parte do projeto PyLCP, que é um pacote Python para análise de linhas de transmissão.

REFERENCES:
[1] PATEL, Utkarsh R. A Surface Admittance Approach For Fast Calculation of the 
    Series Impedance of Cables Including Skin, Proximity, and Ground Return Effects.
    2014. University of Toronto, Graduate Department of The Edward S. Rogers Sr. 
    Department of Electrical & Computer Engineering. 

[2] U. R. Patel, B. Gustavsen and P. Triverio, "An Equivalent Surface Current Approach
    for the Computation of the Series Impedance of Power Cables with Inclusion of Skin
    and Proximity Effects," in IEEE Transactions on Power Delivery, vol. 28, no. 4, pp.
    2474-2482, Oct. 2013, doi: 10.1109/TPWRD.2013.2267098.

[3] U. R. Patel, B. Gustavsen and P. Triverio, "Application of the MoM-SO Method for 
    Accurate Impedance Calculation of Single-Core Cables Enclosed by a Conducting Pipe," 
    Proc. International Conference on Power Systems Transients (IPST 2013), Vancouver, 
    Canada July 18-20, 2013. https://www.ipstconf.org/Proc_IPST2013.php

[4] A. Ametani, "A General Formulation of Impedance and Admittance of Cables," in IEEE
    Transactions on Power Apparatus and Systems, vol. PAS-99, no. 3, pp. 902-910, May
    1980, doi: 10.1109/TPAS.1980.319718.

[5] A. Ametani, "Wave Propagation Characteristics of Cables," in IEEE Transactions on
    Power Apparatus and Systems, vol. PAS-99, no. 2, pp. 499-505, March 1980, 
    doi: 10.1109/TPAS.1980.319685.
"""
import os
import sys
import time
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# RAIZ DO PROJETO E DIRETÓRIOS
os.system('cls' if os.name == 'nt' else 'clear')
try:
    script_dir = Path(__file__).resolve().parent
    print(f"Script directory: {script_dir}")
    project_root = script_dir.parents[1]
    print(f"Project root: {project_root}")
    sys.path.append(str(project_root))
    print("Caminhos do projeto configurados com sucesso.")
except IndexError:
    raise FileNotFoundError(
        "Não foi possível encontrar a raiz do projeto. "
        "Certifique-se de que o script está em 'examples/coated_wires'."
    )

# IMPORTAÇÕES DOS MÓDULOS E MODELO DE DADOS
try:
    from mtl_main.models_wires import COAXIAL_CABLE as MODEL
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_formulation.isolated_wires import CoaxialCable, Ametani
    from mom_so.quasi_static_green import QuasiStatic
    from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

def run_analytical_simulation(mtl, frequencies):
    """
    Executa a simulação analítica da impedância da linha de transmissão.

    Args:
        mtl_config (dict): Dicionário de configuração da linha de transmissão.
        frequencies (np.ndarray): Array de frequências para a análise.

    Returns:
        tuple: Uma tupla contendo três listas: impedâncias série,
               resistências de alta frequência e indutâncias externas.
    """
    print("Iniciando rotina analítica...")
    analytical_data = {}    
    coaxial = CoaxialCable(mtl)
    ametani = Ametani(mtl)

    for freq in frequencies:
        zs = coaxial.pul_parameters(freq)
        l_ext = coaxial.external_inductance()
        z11, z12, z22 = ametani.impedance_two_layered_conductor(freq)

        analytical_data[freq] = {
            'zs': zs,
            'le': l_ext,
            'rs': np.real(zs),
            'ls': np.imag(zs) / (2 * np.pi * freq),
            'r11': np.real(z11),
            'r12': np.real(z12),
            'r22': np.real(z22),
            'l11': np.imag(z11) / (2 * np.pi * freq),
            'l12': np.imag(z12) / (2 * np.pi * freq),
            'l22': np.imag(z22) / (2 * np.pi * freq),
        }

    return analytical_data

def run_momso_simulation(mtl, frequencies):
    """
    Executa a simulação da impedância usando o Método dos Momentos (MoM-SO).

    Args:
        mtl_config (dict): Dicionário de configuração da linha de transmissão.
        frequencies (np.ndarray): Array de frequências para a análise.
        green_mode (GreenFunctionMode): O modo de cálculo para a função de Green.

    Returns:
        list: Uma lista contendo as impedâncias série totais calculadas via MoM.
    """
    print("Iniciando rotina numérica (MoM-SO)...")
    green_matrix = QuasiStatic(mtl).green_matrix()
    post_processor = LosslessPostProcessing(mtl)
    momso_data = {}

    for freq in frequencies:
        mom_so = HomogeneousLosslessMedium(mtl, freq)
        z_partial = mom_so.z_partial(green_matrix)
        zs = post_processor.z_total(z_partial)
        momso_data[freq] = {
            'zs': zs,
            'rs': post_processor.rs_matrix(zs),
            'ls': post_processor.ls_matrix(zs, freq)
        }

    return momso_data

def plot_results(freqs, analytical_data, mom_so_data):
    """
    Gera e exibe os gráficos dos resultados da simulação de forma flexível,
    organizados em subplots. Esta versão é refatorada para maior clareza e
    manutenibilidade.

    Args:
        freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
        analytical_data (dict): Dicionário com os resultados da simulação analítica.
        mom_so_data (dict): Dicionário com os resultados da simulação MoM-SO.
    """
    print("Gerando gráficos...")

    def _configure_subplot(ax, ylabel, data_to_plot, y_lim, yscale='log'):
        for plot_params in data_to_plot:
            frequencies, values = plot_params['data']
            plot_type = plot_params.get('type', 'plot')
            label = plot_params.get('label')

            if plot_type == 'scatter':
                ax.scatter(frequencies, values, label=label, facecolors='none', edgecolors='k', marker='o')
            else:
                ax.plot(
                    frequencies, values, label=label,
                    color=plot_params.get('color', 'k'),
                    linestyle=plot_params.get('linestyle', '-')
                )

        ax.set_xscale('log')
        ax.set_yscale(yscale)
        ax.set_xlim(1E0, 1E6)
        ax.set_ylim(y_lim)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(False)

    def _build_plot_data(param_prefix):
        """Função auxiliar que usa a configuração explícita para construir os dados."""
        data_list = []
        param_prefix_upper = param_prefix.upper()

        for p_def in plot_definitions:
            data_key = p_def[f'{param_prefix}_key'] # Acessa a chave correta: 'r_key' ou 'l_key'
            if data_key is None:
                continue
            plot_frequencies = freqs.get(p_def['freq_key'])            
            label = p_def.get('label') or p_def.get('label_template', '').format(P=param_prefix_upper)
            plot_dict = {
                'type': p_def['type'],
                'data': (plot_frequencies, processed_data[data_key]),
                'label': label
            }
            plot_dict.update(p_def.get('style', {}))
            data_list.append(plot_dict)
        return data_list

    # --- 1. Processamento de Dados com List Comprehensions ---
    R_FACTOR = 1E3
    L_FACTOR = 1e6
    analytical_values = analytical_data.values()
    mom_so_values = mom_so_data.values()

    processed_data = {
        'rs':     [d['rs'] * R_FACTOR for d in analytical_values],
        'ls':     [d['ls'] * L_FACTOR for d in analytical_values],
        'le':     [d['le'] * L_FACTOR for d in analytical_values],
        'r11':    [d['r11'] * R_FACTOR for d in analytical_values],
        'l11':    [d['l11'] * L_FACTOR for d in analytical_values],
        'r12':    [d['r12'] * R_FACTOR for d in analytical_values],
        'l12':    [d['l12'] * L_FACTOR for d in analytical_values],
        'r22':    [d['r22'] * R_FACTOR for d in analytical_values],
        'l22':    [d['l22'] * L_FACTOR for d in analytical_values],
        'rs_mom': [d['rs'][0, 0] * R_FACTOR for d in mom_so_values],
        'ls_mom': [d['ls'][0, 0] * L_FACTOR for d in mom_so_values],
        'req':    [(d['r11']-2*d['r12']+d['r22']) * R_FACTOR for d in analytical_values],
        'leq':    [(d['l11']-2*d['l12']+d['l22']) * L_FACTOR for d in analytical_values],
    }

    plot_definitions = [
        {'r_key': 'rs_mom', 'l_key': 'ls_mom', 'freq_key': 'mom', 'type': 'scatter', 'label': 'MoM-SO'},
        {'r_key': 'rs',     'l_key': 'ls',     'freq_key': 'ana', 'type': 'plot',    'label': 'Analytical'},
        {'r_key': None,     'l_key': 'le',     'freq_key': 'ana', 'type': 'plot',    'label': 'Asymptotic',            'style': {'color': 'k', 'linestyle': '-.'}},
        {'r_key': 'r11',    'l_key': 'l11',    'freq_key': 'ana', 'type': 'plot',    'label_template': '${P}_{{11}}$', 'style': {'color': 'k', 'linestyle': ':'}},
        {'r_key': 'r12',    'l_key': 'l12',    'freq_key': 'ana', 'type': 'plot',    'label_template': '${P}_{{12}}$', 'style': {'color': 'g', 'linestyle': ':'}},
        {'r_key': 'r22',    'l_key': 'l22',    'freq_key': 'ana', 'type': 'plot',    'label_template': '${P}_{{22}}$', 'style': {'color': 'b', 'linestyle': ':'}},
        {'r_key': 'req',    'l_key': 'leq',    'freq_key': 'ana', 'type': 'plot',    'label_template': '${P}_{{11}}-2{P}_{{12}}+{P}_{{22}}$', 'style': {'color': 'r', 'linestyle': ':'}},
    ]
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    resistance_data = _build_plot_data('r')
    inductance_data = _build_plot_data('l')    
    _configure_subplot(axes[0], r'Series Resistance p.u.l. ($\Omega$/km)', resistance_data, y_lim=(1E-2, 1E1))
    _configure_subplot(axes[1], 'Series Inductance p.u.l. (mH/km)', inductance_data, y_lim=(0.0, 0.2), yscale='linear')
    plt.tight_layout()

if __name__ == "__main__":
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    st = time.time()
    os.system('cls' if os.name == 'nt' else 'clear')
    print("Iniciando cálculos da impedância p.u.l. do cabo coaxial...")

    # Rotinas Analítica e MoM-SO
    mtl_model = MulticonductorTransmissionLine(MODEL)
    FREQUENCY = {'ana': np.logspace(0, 5.9, num=200), 'mom': np.logspace(0, 5.9, num=30)}    
    analytical_data = run_analytical_simulation(MODEL, FREQUENCY['ana'])
    momso_data = run_momso_simulation(MODEL, FREQUENCY['mom'])

    print(f"\nRotinas de cálculo finalizadas! Tempo de simulação: {(time.time() - st):.2f} segundos.")
    plot_results(FREQUENCY, analytical_data, momso_data)
    MTLRepresentation(mtl_model, units='millimeter').isolated_coaxial_cables()
    plt.show()