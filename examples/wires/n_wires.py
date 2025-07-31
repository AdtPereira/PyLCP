import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0

# Adiciona a raiz do projeto ao PYTHONPATH para importação de módulos.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from data.models import MTL_MODELS
from lossless_systems.wires_homogeneous_media import WiresHomogeneousMedia

from mom_so.utils import *
from mom_so.mom import TwoBareWireSystem
from mom_so.mom import MulticonductorBareWireSystems
from mom_so.green import QuasiStatic
from mom_so.mtl_graphics import MTLRepresentation
from mom_so.patel import HomogeneousLosslessMedium, LosslessPostProcessing

# --- Configurações da Simulação ---
MTL = MTL_MODELS['wires'][2][21][8]
FREQUENCY_RANGE = {'ana': np.logspace(0, 6, num=200), 'mom': np.logspace(0, 6, num=30)}

def run_n_wires_mom():
    """
    Executa a simulação clássica do Método dos Momentos (MoM) para a linha de transmissão bifilar.

    Returns:
        BifilarMoM: Instância do objeto BifilarMoM configurado.
    """
    print("\n=== Rotina numérica MoM MulticonductorBareWireSystems (Collocation Method) ===")
    mom_wires = MulticonductorBareWireSystems(MTL)
    mom_wires.run_simulation()
    mom_wires.print_results()
    mom_wires.plot_charge_density()
    mom_wires.plot_collocation_points()
    mom_wires.plot_harmonic_coefficients()
    MulticonductorBareWireSystems.plot_convergence_rates(MTL, nf_max=20)
    return mom_wires


def run_analytical(mtl, frequencies):
    """
    Executa a simulação analítica da impedância da linha de transmissão.

    Args:
        mtl_config (dict): Dicionário de configuração da linha de transmissão.
        frequencies (np.ndarray): Array de frequências para a análise.

    Returns:
        tuple: Uma tupla contendo três listas: impedâncias série,
               resistências de alta frequência e indutâncias externas.
    """
    print("\n=== Rotina Analítica ===")
    analytical_data = {}
    wires = WiresHomogeneousMedia(mtl)
    le_wires = wires.n_wires_inductance_matrix()
    pul_bifilar = wires.bifilar_pul_inductance_and_capacitance()
    
    for freq in frequencies:
        z_s, r_hf = wires.bifilar_pul_series_impedance(freq)
        analytical_data[freq] = {
            'ls': np.imag(z_s) / (2 * np.pi * freq),
            'le_wires': le_wires,
            'cap_wires': wires.n_wires_capacitance_matrix(le_wires),
            'le_exact': pul_bifilar['indutância']['exact'],
            'le_bifilar': pul_bifilar['indutância']['approximate'],
            'cap_exact': pul_bifilar['capacitante']['exact'],
            'cap_approx': pul_bifilar['capacitante']['approximate'],
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
    print("\n=== Rotina numérica MoM-SO ===")
    green_matrix = QuasiStatic(mtl).g_tanaka()
    post_processor = LosslessPostProcessing(mtl)
    momso_data = {}

    for freq in frequencies:
        mom_so = HomogeneousLosslessMedium(mtl, freq)
        zs = post_processor.z_total(mom_so.z_partial(green_matrix))
        general_cap = mom_so.generalized_capacitance_matrix(green_matrix)
        maxwell_cap = mom_so.maxwellian_capacitance_matrix(general_cap)
        momso_data[freq] = {
            'zs': zs,
            'rs': post_processor.rs_matrix(zs),
            'ls': post_processor.ls_matrix(zs, freq),
            'c': maxwell_cap,
        }

    if green_matrix.shape[0] < 6:
        print(f"\nGreen's Matrix (Dim: {green_matrix.shape}):\n{green_matrix}")
        print(f"\n2*pi*e0*G:\n{- 2 * np.pi * epsilon_0 * np.real(green_matrix)}")
    
    print(f"\nMoM-SO Generalized Capacitance Matrix (Dim: {general_cap.shape}):\n{np.real(general_cap)}")
    print(f"\nMoM-SO Bifilar Capacitance: {np.real(maxwell_cap.item()) * 1E12:.4f} pF/m")
    
    return momso_data


def plot_results(freqs, analytical_data, mom_so_data, mom_data):
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
                    ax.scatter(frequencies, values, label=label, facecolors='k', edgecolors='k', marker='o', s=12)
                elif key == 'MoM':
                    ax.scatter(frequencies, values, label=label, facecolors='g', marker='x', s=12)
                elif key == 'Analytical':
                    ax.plot(frequencies, values, label=label, color='k', linestyle='--')
                elif key == 'Wires':
                    ax.plot(frequencies, values, label=label, color='k', linestyle='-.')
                elif key == 'Bifilar':
                    ax.plot(frequencies, values, label=label, color='r', linestyle=':')
                elif key == 'Exact':
                    ax.plot(frequencies, values, label=label, color='b', linestyle='-')

        # Configurações do subplot
        ax.set_xscale('log')
        ax.set_yscale(yscale)
        ax.set_xlim(1E0, 1E6)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(False)

    # Fatores de conversão de unidade
    C_FACTOR = 1e12  # de F/m para nF/km
    L_FACTOR = 1e6   # de H/m para mH/km

    # 1. Inicialize listas vazias para armazenar os resultados
    ls, le, le_approx, le_wires, ls_mom = [], [], [], [], []
    cap_approx, cap_exact, cap_mom_so, cap_mom, cap_wires = [], [], [], [], []

    # 2. Processe os dados analíticos em um único laço
    for data in analytical_data.values():
        ls.append(data['ls'][0, 0] * L_FACTOR)
        le.append(data['le_exact'] * L_FACTOR)    
        cap_approx.append(data['cap_approx'] * C_FACTOR)
        cap_exact.append(data['cap_exact'] * C_FACTOR)
        cap_wires.append(data['cap_wires'][0, 0] * C_FACTOR)
        le_wires.append(data['le_wires'][0, 0] * L_FACTOR)
        le_approx.append(data['le_bifilar'] * L_FACTOR)

    # 3. Processe os dados do MoM-SO em um laço separado
    for data in mom_so_data.values():
        ls_mom.append(data['ls'][0, 0] * L_FACTOR)
        cap_mom_so.append(np.real(data['c'][0, 0]) * C_FACTOR)
        cap_mom.append(mom_data.C_maxwellian * C_FACTOR)

    inductance_data = {
        'MoM-SO': {'data': (freqs.get('mom'), ls_mom), 'label': 'MoM-SO'},
        'Analytical': {'data': (freqs.get('ana'), ls), 'label': '$\ell_s$'},
        'Exact': {'data': (freqs.get('ana'), le), 'label': '$\ell_{e,bifilar}$ (Exact)'},
        'Wires': {'data': (freqs.get('ana'), le_wires), 'label': r'$\ell_{e,n+1 \; wires}$'},
        'Bifilar': {'data': (freqs.get('ana'), le_approx), 'label': '$\ell_{e,bifilar}$ (Approx.)'},
    }

    capacitante_data = {
        'MoM-SO': {'data': (freqs.get('mom'), cap_mom_so), 'label': 'MoM-SO'},
        'MoM': {'data': (freqs.get('mom'), cap_mom), 'label': 'MoM (Collocation Method)'},
        'Exact': {'data': (freqs.get('ana'), cap_exact), 'label': r'$c_{bifilar}$ (Exact)'},
        'Wires': {'data': (freqs.get('ana'), cap_wires), 'label': r'$c_{n+1 \; wires}$'},
        'Bifilar': {'data': (freqs.get('ana'), cap_approx), 'label': r'$c_{bifilar}$ (Approx.)'},
    }

    # Cria a figura com subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    _configure_subplot(axes[0], 'Series Inductance p.u.l. (mH/km)', inductance_data, yscale='linear')    
    _configure_subplot(axes[1], 'Capacitance p.u.l. (nF/km)', capacitante_data, yscale='linear')
    plt.tight_layout()


def main():
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    clear_screen()
    print("Iniciando cálculos da impedância p.u.l. ...")
    start_time = time.time()

    MTLRepresentation(MTL).wires_and_cables()
    wires_mom_data = run_n_wires_mom()
    analytical_data = run_analytical(MTL, FREQUENCY_RANGE['ana'])
    momso_data = run_momso_simulation(MTL, FREQUENCY_RANGE['mom'])
    
    print(f"\nRotinas de cálculo finalizadas em {(time.time()-start_time):.2f} segundos.")
    plot_results(FREQUENCY_RANGE, analytical_data, momso_data, wires_mom_data)


if __name__ == "__main__":
    main()
    plt.show()