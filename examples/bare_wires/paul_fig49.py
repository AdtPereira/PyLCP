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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from mom_so.utils import *
from mom_so.quasi_static_green import QuasiStatic
from mom_so.lossless_medium import HomogeneousLosslessMedium

from mtl_data.models import *
from analytical_forms.bare_wires import WiresHomogeneousMedia
from mom.bare_wire_systems import MulticonductorBareWireSystems

# --- Configurações da Linha Bifilar ---
MTL = BIFILAR_BARE_WIRE_S21 

def run_analytical(mtl, s_rw_ratios):
    """
    Executa a simulação de capacitância para diferentes geometrias.

    Args:
        s_rw_ratios (np.ndarray): Array com as razões s/r_w a serem analisadas.
        wire_radius (float): O raio fixo do condutor a ser usado na simulação.

    Returns:
        tuple: Uma tupla contendo duas listas com os valores de capacitância
               exata e aproximada em F/m.
    """
    print("Iniciando instância analítica de cálculos...")
    results = {'exact': [], 'approximate': []}

    for ratio in s_rw_ratios:
        # Cria a configuração da linha bifilar dinamicamente para cada razão.
        separation = ratio * mtl[0]['radius'][1]  
        mtl[1]['center_point'] = (separation, 0.0)        

        wires = WiresHomogeneousMedia(mtl)
        pul_bifilar = wires.bifilar_pul_inductance_and_capacitance()
        results['exact'].append(pul_bifilar['capacitance']['exact'])
        results['approximate'].append(pul_bifilar['capacitance']['approximate'])

    print("\nSimulação analítica finalizada.")
    return results


def run_numerical(mtl, s_rw_ratios):
    """
    Executa a simulação de capacitância para diferentes geometrias.

    Args:
        s_rw_ratios (np.ndarray): Array com as razões s/r_w a serem analisadas.
        wire_radius (float): O raio fixo do condutor a ser usado na simulação.

    Returns:
        tuple: Uma tupla contendo duas listas com os valores de capacitância
               exata e aproximada em F/m.
    """
    print("Iniciando instância numérica de cálculos...")
    results = {'mom-so': [], 'mom': []}

    for ratio in s_rw_ratios:
        # Cria a configuração da linha bifilar dinamicamente para cada razão.
        separation = ratio * mtl[0]['radius'][1]  
        mtl[1]['center_point'] = (separation, 0.0)       

        # Instância MoM-SO
        green_matrix = QuasiStatic(mtl).green_matrix()
        mom_so = HomogeneousLosslessMedium(mtl, frequency=0)
        gc = mom_so.generalized_capacitance_matrix(green_matrix)
        c = mom_so.maxwellian_capacitance_matrix(gc).item()
        results['mom-so'].append(np.real(c))

        # Instância MoM
        mom_wires = MulticonductorBareWireSystems(mtl)
        mom_wires.run_simulation()
        results['mom'].append(mom_wires.C_maxwellian)

    print("\nSimulação finalizada.")
    return results


def plot_comparison(s_rw_ratios, analytical_data, numerical_data):
    """
    Gera o gráfico comparativo das capacitâncias, similar à Figura 4.9.

    Args:
        s_rw_ratios (np.ndarray): Array com as razões s/r_w (eixo X).
        sim_results (tuple): Tupla com as listas de resultados de capacitância.
    """
    print("Gerando gráfico da capacitância...")
       
    C_FACTOR = 1e12 # Converte de F/m para pF/m
    exact_pF = np.array(analytical_data['exact']) * C_FACTOR
    approx_pF = np.array(analytical_data['approximate']) * C_FACTOR
    mom_so_pF = np.array(numerical_data['mom-so']) * C_FACTOR
    mom_pF = np.array(numerical_data['mom']) * C_FACTOR

    plt.figure(figsize=(8, 5))
    plt.plot(s_rw_ratios['ana'], exact_pF, 'k-', label='Exact')
    plt.plot(s_rw_ratios['ana'], approx_pF, 'k--', label='Approximate')
    
    plt.scatter(s_rw_ratios['mom'], mom_so_pF, color='r', marker='o', s=8, label='MoM-SO')
    plt.plot(s_rw_ratios['mom'], mom_pF, label='MoM', linestyle='none', marker='d', markersize=6, fillstyle='none', markeredgecolor='black')
    
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
    start_time = time.time()

    S_RW_RATIOS = {'ana': np.linspace(2.1, 8, num=300), 'mom': np.linspace(2.1, 8, num=40)}
    data_ana = run_analytical(MTL, S_RW_RATIOS['ana'])
    data_mom = run_numerical(MTL, S_RW_RATIOS['mom'])

    # 7. Geração e exibição do gráfico comparativo
    elapsed_time = time.time() - start_time
    print(f"\nRotinas de cálculo finalizadas! Tempo de simulação: {elapsed_time:.2f} segundos.")
    plot_comparison(S_RW_RATIOS, data_ana, data_mom)


if __name__ == "__main__":
    main()
    plt.show()
