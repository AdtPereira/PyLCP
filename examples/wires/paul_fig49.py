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
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH para importação de módulos.
# ATENÇÃO: Esta é uma solução frágil. O ideal é instalar o projeto como um pacote.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from lossless_systems.wires_homogeneous_media import WiresHomogeneousMedia
from mom_so.utils import *
from mom_so.green import QuasiStatic
from mom_so.patel import HomogeneousLosslessMedium


def run_simulation(s_rw_ratios, wire_radius):
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
    results = {'exact': [], 'approximate': [], 'mom-so': []}

    for ratio in s_rw_ratios:
        separation = ratio * wire_radius
        frequency = 0

        # Cria a configuração da linha bifilar dinamicamente para cada razão.
        # Os condutores são posicionados em (0, 0) e (separation, 0).
        mtl = {
            'type': 'bifilar_line',
            'data': [
                {
                    'line_id': 0,
                    'conductor_name': 'p',
                    'line_type': 'active',
                    'line_return': 1,
                    'center_point': (0.0, 0.0),
                    'radius': [0, wire_radius],
                    'conductivity': 1/1.68E-8,
                    'subconductors': None,
                    'insulation': None,
                    'relative_permeability': 1,
                    'relative_permittivity': 1,
                    'relative_permittivity_out': 1,
                    'fourier_order': 4,
                },
                {
                    'line_id': 1,
                    'conductor_name': 'q',
                    'line_type': 'return',
                    'line_return': None,
                    'center_point': (separation, 0.0),
                    'radius': [0, wire_radius],
                    'conductivity': 1/1.68E-8,
                    'subconductors': None,
                    'insulation': None,
                    'relative_permeability': 1,
                    'relative_permittivity': 1,
                    'relative_permittivity_out': 1,
                    'fourier_order': 4,
                }
            ]
        }

        # Instancia analítica
        wires = WiresHomogeneousMedia(mtl)
        pul_bifilar = wires.bifilar_pul_inductance_and_capacitance()
        results['exact'].append(pul_bifilar['capacitante']['exact'])
        results['approximate'].append(pul_bifilar['capacitante']['approximate'])

        # Instancia MoM-SO
        green_matrix = QuasiStatic(mtl).g_tanaka()
        mom_so = HomogeneousLosslessMedium(mtl, frequency)
        gc = mom_so.generalized_capacitance_matrix(green_matrix)
        c = mom_so.maxwellian_capacitance_matrix(gc)
        # post_processor = LosslessPostProcessing(mtl)
        # z_partial = mom_so.z_partial(green_matrix)
        # zs = post_processor.z_total(z_partial)
        # ls = post_processor.ls_matrix(zs, frequency)
        # from scipy.constants import mu_0, epsilon_0
        # results['mom-so'].append(mu_0 * epsilon_0 / ls[0, 0])
        results['mom-so'].append(np.real(c[0,0]))

    print("\nSimulação finalizada.")
    return results


def plot_comparison(s_rw_ratios, sim_results):
    """
    Gera o gráfico comparativo das capacitâncias, similar à Figura 4.9.

    Args:
        s_rw_ratios (np.ndarray): Array com as razões s/r_w (eixo X).
        sim_results (tuple): Tupla com as listas de resultados de capacitância.
    """
    print("Gerando gráfico da capacitância...")
       
    C_FACTOR = 1e12 # Converte de F/m para pF/m
    exact_pF = np.array(sim_results['exact']) * C_FACTOR
    approx_pF = np.array(sim_results['approximate']) * C_FACTOR
    mom_so_pF = np.array(sim_results['mom-so']) * C_FACTOR

    plt.figure(figsize=(8, 6))
    plt.plot(s_rw_ratios, exact_pF, 'k-', label='Exact')
    plt.plot(s_rw_ratios, approx_pF, 'k--', label='Approximate')
    plt.scatter(s_rw_ratios, mom_so_pF, color='r', marker='o', s=10, label=r'MoM-SO (Generalized Cap. Concept)')
    plt.xlabel('Ratio of separation to wire radius, s/r$_w$')
    plt.ylabel('Per-unit-length capacitance (pF/m)')
    # plt.title('Figure 4.9: A comparison of the exact and approximate formulas \n for the per-unit-length capacitance of two wires [1].')
    plt.xlim(2, 8)
    plt.ylim(10, 90)
    plt.legend()
    plt.grid(True, linestyle='--', linewidth=0.5)


def main():
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    clear_screen()
    print("Iniciando cálculos da impedância p.u.l. ...")
    start_time = time.time()

    WIRE_RADIUS = 0.010
    S_RW_RATIOS = np.linspace(2.1, 8, num=120)
    capacitance_results = run_simulation(S_RW_RATIOS, WIRE_RADIUS)
    
    # 3. Medição de tempo
    elapsed_time = time.time() - start_time
    print(f"\nRotinas de cálculo finalizadas! Tempo de simulação: {elapsed_time:.2f} segundos.")
    
    # 7. Geração e exibição do gráfico comparativo
    plot_comparison(S_RW_RATIOS, capacitance_results)


if __name__ == "__main__":
    main()
    plt.show()
