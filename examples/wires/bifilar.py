import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH para importação de módulos.
# ATENÇÃO: Esta é uma solução frágil. O ideal é instalar o projeto como um pacote.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from data.models import MTL_MODELS
from lossless_systems.wires_homogeneous_media import WiresHomogeneousMedia

from mom_so.utils import *
from mom_so.green import QuasiStatic
from mom_so.mtl_graphics import MTLRepresentation
from mom_so.patel import HomogeneousLosslessMedium, LosslessPostProcessing


# --- Configurações da Simulação ---
#[num_conductor][separation(mm)][fourier_order]
MTL = MTL_MODELS['wires'][2][100][4]
MTL = MTL_MODELS['wires'][2][25][0]
MTL = MTL_MODELS['wires'][2][25][4]
FREQUENCY_RANGE = {'ana': np.logspace(0, 6, num=200), 'mom': np.logspace(0, 6, num=30)}


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
    wires = WiresHomogeneousMedia(mtl)
    pul_bifilar = wires.bifilar_pul_inductance_and_capacitance()
    
    for freq in frequencies:
        z_s, r_hf = wires.bifilar_pul_series_impedance(freq)
        analytical_data[freq] = {
            'rs': np.real(z_s),
            'ls': np.imag(z_s) / (2 * np.pi * freq),
            'rhf': r_hf,
            'le_exact': pul_bifilar['indutância']['exact'],
            'le_approx': pul_bifilar['indutância']['approximate'],
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
    green_matrix = QuasiStatic(mtl).g_tanaka()
    post_processor = LosslessPostProcessing(mtl)
    momso_data = {}

    for freq in frequencies:
        mom_so = HomogeneousLosslessMedium(mtl, freq)
        z_partial = mom_so.z_partial(green_matrix)
        zs = post_processor.z_total(z_partial)
        momso_data[freq] = {
            'rs': post_processor.rs_matrix(zs),
            'ls': post_processor.ls_matrix(zs, freq)
        }

    return momso_data


def plot_results(freqs, analytical_data, mom_so_data):
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

    def _configure_subplot(ax, ylabel, data_to_plot, yscale='log'):
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
        for key, data in data_to_plot.items():
            frequencies, values = data['data']
            label = data['label']

            # Verifica se as frequências e valores estão disponíveis
            if frequencies is not None and values is not None:
                if key == 'MoM-SO':
                    ax.scatter(frequencies, values, label=label, facecolors='none', edgecolors='k', marker='o')
                elif key == 'Analytical':
                    ax.plot(frequencies, values, label=label, color='k', linestyle='--')
                elif key == 'Approximate':
                    ax.plot(frequencies, values, label=label, color='k', linestyle=':')
                elif key == 'Exactly':
                    ax.plot(frequencies, values, label=label, color='k', linestyle='-')

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

    # 1. Inicialize listas vazias para armazenar os resultados
    rs, r_hf, ls, le_exact, le_approx = [], [], [], [], []
    rs_mom, ls_mom = [], []

    # 2. Processe os dados analíticos em um único laço
    for data in analytical_data.values():
        rs.append(data['rs'][0, 0] * R_FACTOR)
        ls.append(data['ls'][0, 0] * L_FACTOR)
        r_hf.append(data['rhf'][0, 0] * R_FACTOR)
        le_exact.append(data['le_exact'] * L_FACTOR)
        le_approx.append(data['le_approx'] * L_FACTOR)

    # 3. Processe os dados do MoM-SO em um laço separado
    for data in mom_so_data.values():
        rs_mom.append(data['rs'][0, 0] * R_FACTOR)
        ls_mom.append(data['ls'][0, 0] * L_FACTOR)

    resistance_data = {
        'MoM-SO': {'data': (freqs.get('mom'), rs_mom), 'label': 'MoM-SO'},
        'Analytical': {'data': (freqs.get('ana'), rs), 'label':'$R_i$'},
        'Approximate': {'data': (freqs.get('ana'), r_hf), 'label': '$R_{HF}$'},
    }

    inductance_data = {
        'MoM-SO': {'data': (freqs.get('mom'), ls_mom), 'label': 'MoM-SO'},
        'Analytical': {'data': (freqs.get('ana'), ls), 'label': '$\ell_s$'},
        'Approximate': {'data': (freqs.get('ana'), le_approx), 'label': '$\ell_{e}$ (Approx.)'},
        'Exactly': {'data': (freqs.get('ana'), le_exact), 'label': '$\ell_{e,bifilar}$ (Exactly)'},
    }

    # Cria a figura com subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Configura cada subplot
    _configure_subplot(axes[0], r'Series Resistance p.u.l. ($\Omega$/km)', resistance_data)    
    _configure_subplot(axes[1], 'Series Inductance p.u.l. (mH/km)', inductance_data, yscale='linear')    
    plt.tight_layout()
    plt.show()


def main():
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    clear_screen()
    print("Iniciando cálculos da impedância p.u.l. ...")
    start_time = time.time()

    try:
        # 1. Exibe a geometria da linha (em uma figura separada)
        MTLRepresentation(MTL).wires_and_cables()

        # 1. Rotina Analítica
        analytical_data = run_analytical_simulation(MTL, FREQUENCY_RANGE['ana'])

        # 2. Rotina MoM-SO
        momso_data = run_momso_simulation(MTL, FREQUENCY_RANGE['mom'])

        # 3. Medição de tempo
        elapsed_time = time.time() - start_time
        print(f"\nRotinas de cálculo finalizadas! Tempo de simulação: {elapsed_time:.2f} segundos.")

        # 5. Geração e exibição dos resultados com a função revisada
        plot_results(FREQUENCY_RANGE, analytical_data, momso_data)

    except Exception as e:
        print(f"\nOcorreu um erro durante a execução do script: {e}")
        print("Verifique as configurações de entrada e as dependências do projeto.")


if __name__ == "__main__":
    main()