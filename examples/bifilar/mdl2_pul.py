"""
Análise Comparativa da Impedância por Unidade de Comprimento (p.u.l.) de uma Linha de Transmissão Bifilar.

Este script calcula e compara a impedância série por unidade de comprimento de uma
linha de transmissão de dois fios (bifilar). A análise é realizada através de duas abordagens:
1.  Um método analítico clássico.
2.  Um método numérico baseado no Método dos Momentos com Ortogonalidade de Superfície (MoM-SO).

Os resultados de resistência e indutância em função da frequência são plotados para
comparação visual das duas metodologias.

Dependências do Projeto:
- numpy
- As classes e dados dos módulos:
  - data.systems
  - data.multiconductor
  - data.graph
  - lossless_medium.pul_parameters
  - mom_so.*

Execução:
Para executar o script, execute o seguinte comando no terminal, a partir do diretório raiz do projeto:
$ python caminho/para/mdl2_pul_refatorado.py

NOTA:
Este script utiliza uma manipulação de `sys.path` para localizar os módulos do projeto.
Para uma solução mais robusta, recomenda-se transformar o projeto em um pacote Python
instalável (usando `setup.py` ou `pyproject.toml`).

REFERÊNCIAS:
[1] Preencher com as referências relevantes.
"""

import os
import sys
import time
import numpy as np
from enum import Enum

# Adiciona a raiz do projeto ao PYTHONPATH para importação de módulos.
# ATENÇÃO: Esta é uma solução frágil. O ideal é instalar o projeto como um pacote.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from data.systems import MULTICONDUCTOR_TRANSMISSION_LINE
from data.mtl import MulticonductorTransmissionLine
from data.graph import GraphicRepresentation as graph
from lossless_systems.bifilar_line import Bifilar
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


def plot_results(mtl, frequencies, analytical_results, mom_results):
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
    graph(mtl).wires_and_cables(line_type='bifilar')

    # Parâmetros para os gráficos
    d = MulticonductorTransmissionLine(mtl).D_pq
    plot_data = [0, mtl['data'][0]['fourier_order'], d]
    
    simulation_results = [impedance_analytical, mom_results]
    
    # Plota a resistência série
    Bifilar(mtl).plot_series_resistance((freq_analytical, freq_mom), simulation_results, hf_resistance, plot_data)
    
    # Plota a indutância série
    Bifilar(mtl).plot_series_inductance((freq_analytical, freq_mom), simulation_results, external_inductance, plot_data)


def main():
    """
    Função principal para orquestrar a análise, cálculo e visualização dos resultados.
    """
    clear_screen()
    print("Iniciando cálculos da impedância p.u.l. ...")
    start_time = time.time()

    try:
        # 1. Rotina Analítica
        zi_ana, res_hf, ind_externa = run_analytical_simulation(
            MTL_CONFIG, FREQUENCY_RANGE_ANA
        )

        # 2. Rotina MoM-SO
        zi_momso = run_momso_simulation(
            MTL_CONFIG, FREQUENCY_RANGE_MOM, GreenFunctionMode.ANALYTICAL
        )

        # 3. Medição de tempo
        elapsed_time = time.time() - start_time
        print(f"\nRotinas de cálculo finalizadas! Tempo de simulação: {elapsed_time:.2f} segundos.")

        # 4. Geração e exibição dos resultados
        plot_results(
            MTL_CONFIG,
            (FREQUENCY_RANGE_ANA, FREQUENCY_RANGE_MOM),
            (zi_ana, res_hf, ind_externa),
            zi_momso
        )

    except Exception as e:
        print(f"\nOcorreu um erro durante a execução do script: {e}")
        print("Verifique as configurações de entrada e as dependências do projeto.")


if __name__ == "__main__":
    main()