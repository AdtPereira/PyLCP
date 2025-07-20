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

from data.systems import MULTICONDUCTOR_TRANSMISSION_LINE
from data.mtl import MulticonductorTransmissionLine
from data.graph import GraphicRepresentation as graph

from mom_so.green import QuasiStatic
from mom_so.patel import HomogeneousLosslessMedium, LosslessPostProcessing
from mom_so.utils import clear_screen

# --- Configurações da Simulação ---
BIFILAR_LINE_INDEX = 2
FREQUENCY_RANGE_ANA = np.logspace(0, 7, num=200)
FREQUENCY_RANGE_MOM = np.logspace(0, 7, num=30)
MTL_CONFIG = MULTICONDUCTOR_TRANSMISSION_LINE[BIFILAR_LINE_INDEX]

class GreenFunctionMode(Enum):
    """Define o modo de cálculo para a função de Green."""
    ANALYTICAL = 'Analytically'
    NUMERICAL = 'Numerically'


def run_analytical_simulation(mtl_config, frequencies):
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
    bifilar_analytical = Bifilar(mtl_config)
    series_impedances = []
    hf_resistances = []
    external_inductances = []
    
    for freq in frequencies:
        z_s, r_hf, l_ext = bifilar_analytical.series_impedance(freq)
        series_impedances.append(z_s)
        hf_resistances.append(r_hf)
        external_inductances.append(l_ext)
        
    return series_impedances, hf_resistances, external_inductances


def run_momso_simulation(mtl_config, frequencies, green_mode=GreenFunctionMode.ANALYTICAL):
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
    green_matrix = QuasiStatic(mtl_config).g_tanaka(mode=green_mode.value)
    post_processor = LosslessPostProcessing(mtl_config)
    series_impedances_mom = []

    for freq in frequencies:
        mom_so = HomogeneousLosslessMedium(mtl_config, freq)
        z_partial = mom_so.z_partial(green_matrix)
        series_impedances_mom.append(post_processor.z_total(z_partial))
        
    return series_impedances_mom


def plot_results(mtl_config, frequencies, analytical_results, mom_results):
    """
    Gera e exibe os gráficos dos resultados da simulação.

    Args:
        mtl_config (dict): Dicionário de configuração da linha de transmissão.
        frequencies (tuple): Tupla com os arrays de frequência para cada método.
        analytical_results (tuple): Resultados da simulação analítica.
        mom_results (list): Resultados da simulação MoM.
    """
    print("Gerando gráficos...")
    freq_analytical, freq_mom = frequencies
    impedance_analytical, hf_resistance, external_inductance = analytical_results

    # Exibe a geometria da linha
    graph(mtl_config).wires_and_cables(line_type='bifilar')

    # Parâmetros para os gráficos
    plot_data = [0, mtl_config['data'][0]['fourier_order'], MulticonductorTransmissionLine(mtl_config).distance_matrices[0]]
    
    simulation_results = [impedance_analytical, mom_results]
    
    # Plota a resistência série
    Bifilar(mtl_config).plot_series_resistance((freq_analytical, freq_mom), simulation_results, hf_resistance, plot_data)
    
    # Plota a indutância série
    Bifilar(mtl_config).plot_series_inductance((freq_analytical, freq_mom), simulation_results, external_inductance, plot_data)


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
        # # 1. Rotina Analítica
        # zi_ana, res_hf, ind_externa = run_analytical_simulation(
        #     MTL_CONFIG, FREQUENCY_RANGE_ANA
        # )

        # # 2. Rotina MoM-SO
        # zi_momso = run_momso_simulation(
        #     MTL_CONFIG, FREQUENCY_RANGE_MOM, GreenFunctionMode.ANALYTICAL
        # )

        # # 3. Medição de tempo
        # elapsed_time = time.time() - start_time
        # print(f"\nRotinas de cálculo finalizadas! Tempo de simulação: {elapsed_time:.2f} segundos.")

        # # 4. Geração e exibição dos resultados
        # plot_results(
        #     MTL_CONFIG,
        #     (FREQUENCY_RANGE_ANA, FREQUENCY_RANGE_MOM),
        #     (zi_ana, res_hf, ind_externa),
        #     zi_momso
        # )

        # 5. Simulação de Capacitância
        S_RW_RATIOS = np.linspace(2, 8, num=100)
        WIRE_RADIUS = 0.01  # Raio do fio em metros
        capacitance_results = run_capacitance_simulation(S_RW_RATIOS, WIRE_RADIUS)
        plot_capacitance_comparison(S_RW_RATIOS, capacitance_results)

    except Exception as e:
        print(f"\nOcorreu um erro durante a execução do script: {e}")
        print("Verifique as configurações de entrada e as dependências do projeto.")


if __name__ == "__main__":
    main()