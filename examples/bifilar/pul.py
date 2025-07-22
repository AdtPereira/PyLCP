"""
Reprodução da Figura 4.9 do livro "Introduction to Electromagnetic Compatibility".

Este script calcula e plota a capacitância por unidade de comprimento de uma
linha de transmissão bifilar usando as fórmulas exata e aproximada. O objetivo
é recriar a Figura 4.9, que compara essas duas formulações em função da razão
entre a separação dos condutores e seu raio (s/r_w).

Estrutura:
- Define-se uma faixa de valores para a razão s/r_w.
- Para cada valor, uma geometria de linha bifilar é criada.
- A classe Bifilar é usada para calcular as capacitâncias.
- Os resultados são plotados usando Matplotlib.
"""

import os
import sys
import time
import numpy as np
from enum import Enum
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH para importação de módulos.
# ATENÇÃO: Esta é uma solução frágil. O ideal é instalar o projeto como um pacote.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from lossless_systems.bifilar_line import Bifilar

from data.systems import MTL_MODELS
from data.mtl import MulticonductorTransmissionLine
from data.graph import GraphicRepresentation as graph

from mom_so.green import QuasiStatic
from mom_so.patel import HomogeneousLosslessMedium, LosslessPostProcessing
from mom_so.utils import clear_screen

# --- Configurações da Simulação ---
# MTL_MODELS['bifilar'][separation(mm)][fourier_order]
MTL = MTL_MODELS['bifilar'][25][4]
FREQUENCY_RANGE = {
    'ana': np.logspace(0, 6, num=200),
    'mom': np.logspace(0, 6, num=30)
}

class GreenFunctionMode(Enum):
    """Define o modo de cálculo para a função de Green."""
    ANALYTICAL = 'Analytically'
    NUMERICAL = 'Numerically'


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
    bifilar_analytical = Bifilar(mtl)
    analytical_data = {}
    
    for freq in frequencies:
        z_s, r_hf, l_ext = bifilar_analytical.series_impedance(freq)
        analytical_data[freq] = {
            'zs': z_s,
            'rs': np.real(z_s),
            'ls': np.imag(z_s) / (2 * np.pi * freq),
            'rhf': r_hf,
            'le': l_ext
        }

    return analytical_data


def run_momso_simulation(mtl, frequencies, green_mode=GreenFunctionMode.ANALYTICAL):
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
    green_matrix = QuasiStatic(mtl).g_tanaka(mode=green_mode.value)
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


def plot_results(freqs, analytical, mom_so):
    """
    Gera e exibe os gráficos dos resultados da simulação de forma flexível,
    organizados em subplots.

    Args:
        mtl (dict): Dicionário de configuração da linha de transmissão.
        freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
        analytical (dict): Dicionário com os resultados da simulação analítica.
        mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
    """
    print("Gerando gráficos...")

    def _configure_subplot(ax, ylabel, data_to_plot, ref_data=None, yscale='log'):
        """
        Função auxiliar para configurar um único subplot.

        Args:
            ax (matplotlib.axes.Axes): O eixo do subplot a ser configurado.
            title (str): Título do subplot.
            ylabel (str): Rótulo do eixo Y.
            data_to_plot (dict): Dados principais para plotagem.
            ref_data (tuple, optional): Dados de referência para plotagem.
        """
        # Itera sobre os dados para plotagem
        for label, (frequencies, values) in data_to_plot.items():
            if frequencies is not None and values is not None:
                if label == 'MoM-SO':
                    ax.scatter(frequencies, values, label=label, facecolors='none', edgecolors='k', marker='o')
                else:
                    ax.plot(frequencies, values, label=label, color='k', linestyle='-')

        # Plota os dados de referência, se existirem
        if ref_data:
            frequencies, values, label = ref_data
            if frequencies is not None and values is not None:
                ax.plot(frequencies, values, 'k--', label=label)

        # Configurações do subplot
        ax.set_xscale('log')
        ax.set_yscale(yscale)
        ax.set_xlim(1E0, 1E6)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(False)

    # Fatores de conversão de unidade
    R_FACTOR = 1000  # de Ohm/m para Ohm/km
    L_FACTOR = 1e6   # de H/m para mH/km

    # Extrai as frequências
    ana_freqs = freqs.get('ana')
    mom_freqs = freqs.get('mom')

    # 1. Inicialize listas vazias para armazenar os resultados
    ana_r, r_hf, ana_l, l_ext = [], [], [], []
    mom_r, mom_l = [], []

    # 2. Processe os dados analíticos em um único laço
    for data in analytical.values():
        ana_r.append(data['rs'][0, 0] * R_FACTOR)
        r_hf.append(data['rhf'][0, 0] * R_FACTOR)
        ana_l.append(data['ls'][0, 0] * L_FACTOR)
        l_ext.append(data['le'][0, 0] * L_FACTOR)

    # 3. Processe os dados do MoM-SO em um laço separado
    for data in mom_so.values():
        mom_r.append(data['rs'][0, 0] * R_FACTOR)
        mom_l.append(data['ls'][0, 0] * L_FACTOR)

    resistance_data = {'Analytical': (ana_freqs, ana_r), 'MoM-SO': (mom_freqs, mom_r)}
    inductance_data = {'Analytical': (ana_freqs, ana_l), 'MoM-SO': (mom_freqs, mom_l)}
    hf_resistance_ref = (ana_freqs, r_hf, 'Resistência HF (Analítica)')
    external_inductance_ref = (ana_freqs, l_ext, 'Indutância Externa (Analítica)')

    # Cria a figura com subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Configura cada subplot
    _configure_subplot(axes[0], r'Series Resistance p.u.l. ($\Omega$/km)',
                       resistance_data, hf_resistance_ref)
    
    _configure_subplot(axes[1], 'Series Inductance p.u.l. (mH/km)',
                       inductance_data, external_inductance_ref, yscale='linear')
    
    plt.tight_layout()
    plt.show()


def run_capacitance_simulation(s_rw_ratios, wire_radius):
    """
    Executa a simulação de capacitância para diferentes geometrias.

    Args:
        s_rw_ratios (np.ndarray): Array com as razões s/r_w a serem analisadas.
        wire_radius (float): O raio fixo do condutor a ser usado na simulação.

    Returns:
        tuple: Uma tupla contendo duas listas com os valores de capacitância
               exata e aproximada em F/m.
    """
    print("Iniciando simulação de capacitância...")
    capacitances_exact = []
    capacitances_approx = []

    for ratio in s_rw_ratios:
        separation = ratio * wire_radius

        # Cria a configuração da linha bifilar dinamicamente para cada razão.
        # Os condutores são posicionados em (0, 0) e (separation, 0).
        mtl_config = {
            'type': 'bifilar_line',
            'data': [
                {
                    'center_point': (0.0, 0.0),
                    'radius': [0, wire_radius],
                    'conductivity': 5.8E7, 'relative_permeability': 1,
                    'relative_permittivity': 1, 'relative_permittivity_out': 1,
                    'fourier_order': 4
                },
                {
                    'center_point': (separation, 0.0),
                    'radius': [0, wire_radius],
                    'conductivity': 5.8E7, 'relative_permeability': 1,
                    'relative_permittivity': 1, 'relative_permittivity_out': 1,
                    'fourier_order': 4
                }
            ]
        }

        # Instancia a classe e calcula as capacitâncias
        capacitance_data = Bifilar(mtl_config).capacitance()

        capacitances_exact.append(capacitance_data['exact'])
        capacitances_approx.append(capacitance_data['approximate'])

    print("Simulação finalizada.")
    return (capacitances_exact, capacitances_approx)


def plot_capacitance_comparison(s_rw_ratios, sim_results):
    """
    Gera o gráfico comparativo das capacitâncias, similar à Figura 4.9.

    Args:
        s_rw_ratios (np.ndarray): Array com as razões s/r_w (eixo X).
        sim_results (tuple): Tupla com as listas de resultados de capacitância.
    """
    print("Gerando gráfico da capacitância...")
    exact_caps, approx_caps = sim_results

    # Converte de F/m para pF/m para corresponder ao eixo Y da figura
    exact_caps_pF = np.array(exact_caps) * 1e12
    approx_caps_pF = np.array(approx_caps) * 1e12

    # Configuração da plotagem para imitar a Figura 4.9
    plt.figure(figsize=(8, 6))
    plt.plot(s_rw_ratios, exact_caps_pF, 'k-', label='Exact')
    plt.plot(s_rw_ratios, approx_caps_pF, 'k--', label='Approximate')

    plt.xlabel('Ratio of separation to wire radius, s/r$_w$')
    plt.ylabel('Per-unit-length capacitance (pF/m)')
    plt.title('Figure 4.9: A comparison of the exact and approximate formulas \n for the per-unit-length capacitance of two wires [1].')
    plt.xlim(2, 8)
    plt.ylim(10, 90)
    plt.legend()
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.show()


def main():
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    clear_screen()
    print("Iniciando cálculos da impedância p.u.l. ...")
    start_time = time.time()

    try:
        # 1. Exibe a geometria da linha (em uma figura separada)
        graph(MTL).wires_and_cables(line_type='bifilar')

        # 1. Rotina Analítica
        analytical_data = run_analytical_simulation(MTL, FREQUENCY_RANGE['ana'])

        # 2. Rotina MoM-SO
        momso_data = run_momso_simulation(MTL, FREQUENCY_RANGE['mom'], GreenFunctionMode.ANALYTICAL)

        # 3. Medição de tempo
        elapsed_time = time.time() - start_time
        print(f"\nRotinas de cálculo finalizadas! Tempo de simulação: {elapsed_time:.2f} segundos.")

        # 5. Geração e exibição dos resultados com a função revisada
        plot_results(FREQUENCY_RANGE, analytical_data, momso_data)

        # 6. Simulação de Capacitância
        WIRE_RADIUS = 0.010
        S_RW_RATIOS = np.linspace(2, 8, num=100)
        capacitance_results = run_capacitance_simulation(S_RW_RATIOS, WIRE_RADIUS)
        plot_capacitance_comparison(S_RW_RATIOS, capacitance_results)

    except Exception as e:
        print(f"\nOcorreu um erro durante a execução do script: {e}")
        print("Verifique as configurações de entrada e as dependências do projeto.")


if __name__ == "__main__":
    main()