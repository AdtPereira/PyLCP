import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0

# Adiciona a raiz do projeto ao PYTHONPATH para importação de módulos.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from mtl_data.models import *
from mtl_data.graphics import MTLRepresentation
from mom.bare_wire_systems import MulticonductorBareWireSystems
from analytical_forms.bare_wires import WiresHomogeneousMedia

from mom_so.utils import *
from mom_so.quasi_static_green import QuasiStatic
from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing

# --- Configurações da Linha Bifilar ---
MTL = BIFILAR_BARE_WIRE_S21
FREQUENCY_RANGE = {'ana': np.logspace(0, 6, num=200), 'mom': np.logspace(0, 6, num=30)}

def run_mom(autoPlots=False):
    """
    Executa a simulação clássica do Método dos Momentos (MoM) para a linha de transmissão bifilar.

    Returns:
        BifilarMoM: Instância do objeto BifilarMoM configurado.
    """
    print("\n==============================================================")
    print("=== MoM MulticonductorBareWireSystems (Collocation Method) ===")
    print("==============================================================")
    
    mom_wires = MulticonductorBareWireSystems(MTL)
    mom_wires.run_simulation()
    mom_wires.print_results()

    if autoPlots:
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
    print("\n=============================")
    print("=== Analytical Processing ===")
    print("=============================")
    
    analytical_data = {}
    bifilar_wires = WiresHomogeneousMedia(mtl)
    le_wires = bifilar_wires.n_wires_inductance_matrix()
    pul_bifilar = bifilar_wires.bifilar_pul_inductance_and_capacitance()
    
    for freq in frequencies:
        z_s, r_hf = bifilar_wires.bifilar_pul_series_impedance(freq)
        analytical_data[freq] = {
            'rhf': r_hf,
            'rs': np.real(z_s),
            'ls': np.imag(z_s) / (2 * np.pi * freq),
            'le_wires': le_wires,
            'cap_wires': bifilar_wires.n_wires_capacitance_matrix(le_wires),
            'le_exact': pul_bifilar['inductance']['exact'],
            'le_bifilar': pul_bifilar['inductance']['approximate'],
            'cap_exact': pul_bifilar['capacitance']['exact'],
            'cap_approx': pul_bifilar['capacitance']['approximate'],
        }
    return analytical_data


def run_mom_so(mtl, frequencies):
    """
    Executa a simulação da impedância usando o Método dos Momentos (MoM-SO).

    Args:
        mtl_config (dict): Dicionário de configuração da linha de transmissão.
        frequencies (np.ndarray): Array de frequências para a análise.
        green_mode (GreenFunctionMode): O modo de cálculo para a função de Green.

    Returns:
        list: Uma lista contendo as impedâncias série totais calculadas via MoM.
    """
    print("\n========================================")
    print("=== MoM-SO HomogeneousLosslessMedium ===")
    print("========================================")

    green_matrix = QuasiStatic(mtl).green_matrix()
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
    
    print(f"\nMoM-SO Generalized Capacitance Matrix (F/m):\n{np.real(general_cap)}")
    print(f"\nMoM-SO Bifilar Capacitance: {np.real(maxwell_cap.item()) * 1E12:.4f} pF/m")
    
    return momso_data


def _configure_plot(ax, ylabel, data_to_plot, yscale='log'):
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
                ax.scatter(frequencies, values, label=label, facecolors='k', marker='o', s=12)
            elif key == 'MoM':
                ax.plot(frequencies, values, label=label, linestyle='none', marker='d', markersize=8, fillstyle='none', markeredgecolor='black')
            elif key == 'Analytical':
                ax.plot(frequencies, values, label=label, color='k', linestyle='--')
            elif key == 'Wires':
                ax.plot(frequencies, values, label=label, color='k', linestyle='-.')
            elif key == 'Bifilar':
                ax.plot(frequencies, values, label=label, color='r', linestyle=':')
            elif key == 'Exactly':
                ax.plot(frequencies, values, label=label, color='k', linestyle='-', linewidth=1.0)
            elif key == 'Approximate':
                ax.plot(frequencies, values, label=label, color='k', linestyle=':')

    # Configurações do subplot
    ax.set_xscale('log')
    ax.set_yscale(yscale)
    ax.set_xlim(1E0, 1E6)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel(ylabel)
    ax.legend()
    ax.grid(False)


def plot_resistance_results(freq, analytical_data, mom_so_data):
    """
    Gera e exibe os gráficos dos resultados da simulação de forma flexível,
    organizados em subplots.

    Args:
        mtl (dict): Dicionário de configuração da linha de transmissão.
        freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
        analytical (dict): Dicionário com os resultados da simulação analítica.
        mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
    """
    R_FACTOR = 1e3   # de Ohm/m para Ohm/km
    rs, rhf, rs_mom = [], [], []

    # 2. Processe os dados analíticos em um único laço
    for data in analytical_data.values():
        rs.append(data['rs'] * R_FACTOR)
        rhf.append(data['rhf'] * R_FACTOR)

    # 3. Processe os dados do MoM-SO em um laço separado
    for data in mom_so_data.values():
        rs_mom.append(data['rs'][0, 0] * R_FACTOR)

    resistance_data = {
        'MoM-SO': {'data': (freq.get('mom'), rs_mom), 'label': 'MoM-SO'},
        'Analytical': {'data': (freq.get('ana'), rs), 'label':'$R_i$'},
        'Approximate': {'data': (freq.get('ana'), rhf), 'label': '$R_{HF}$'},
    }

    # Cria a figura com subplots
    fig, ax = plt.subplots(figsize=(8, 5))
    _configure_plot(ax, r'Series Resistance p.u.l. ($\Omega$/km)', resistance_data)    
    plt.tight_layout()


def plot_inductance_results(freq, analytical_data, mom_so_data):
    """
    Gera e exibe os gráficos dos resultados da simulação de forma flexível,
    organizados em subplots.

    Args:
        mtl (dict): Dicionário de configuração da linha de transmissão.
        freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
        analytical (dict): Dicionário com os resultados da simulação analítica.
        mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
    """
    L_FACTOR = 1e6   # de H/m para mH/km
    ls, le, le_approx, le_wires, ls_mom_so = [], [], [], [], []

    # 2. Processe os dados analíticos em um único laço
    for data in analytical_data.values():
        ls.append(data['ls'] * L_FACTOR)
        le.append(data['le_exact'] * L_FACTOR)    
        le_wires.append(data['le_wires'][0, 0] * L_FACTOR)
        le_approx.append(data['le_bifilar'] * L_FACTOR)

    # 3. Processe os dados do MoM-SO em um laço separado
    for data in mom_so_data.values():
        ls_mom_so.append(data['ls'][0, 0] * L_FACTOR)

    inductance_data = {
        'MoM-SO':       {'data': (freq.get('mom'), ls_mom_so), 'label': 'MoM-SO'},
        'Analytical':   {'data': (freq.get('ana'), ls), 'label': '$\ell_s$'},
        'Exactly':      {'data': (freq.get('ana'), le), 'label': '$\ell_{e,bifilar}$ (Exact)'},
        'Wires':        {'data': (freq.get('ana'), le_wires), 'label': r'$\ell_{e,n+1 \; wires}$'},
        'Bifilar':      {'data': (freq.get('ana'), le_approx), 'label': '$\ell_{e,bifilar}$ (Approx.)'},
    }

    # Cria a figura com subplots
    fig, ax = plt.subplots(figsize=(8, 5))
    _configure_plot(ax, 'Series Inductance p.u.l. (mH/km)', inductance_data, yscale='linear')
    plt.tight_layout()


def plot_capacitance_results(freq, analytical_data, mom_so_data, mom_data):
    """
    Gera e exibe os gráficos dos resultados da simulação de forma flexível,
    organizados em subplots.

    Args:
        mtl (dict): Dicionário de configuração da linha de transmissão.
        freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
        analytical (dict): Dicionário com os resultados da simulação analítica.
        mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
    """
    C_FACTOR = 1e12  # de F/m para nF/km
    cap_approx, cap_exact, cap_mom_so, cap_mom, cap_wires = [], [], [], [], []

    # 2. Processe os dados analíticos em um único laço
    for data in analytical_data.values():
        cap_approx.append(data['cap_approx'] * C_FACTOR)
        cap_exact.append(data['cap_exact'] * C_FACTOR)
        cap_wires.append(data['cap_wires'][0, 0] * C_FACTOR)

    # 3. Processe os dados do MoM-SO em um laço separado
    for data in mom_so_data.values():
        cap_mom_so.append(np.real(data['c'][0, 0]) * C_FACTOR)
        cap_mom.append(mom_data.C_maxwellian * C_FACTOR)

    capacitante_data = {
        'MoM-SO':   {'data': (freq.get('mom'), cap_mom_so), 'label': 'MoM-SO'},
        'MoM':      {'data': (freq.get('mom'), cap_mom), 'label': 'MoM'},
        'Exactly':  {'data': (freq.get('ana'), cap_exact), 'label': r'$c_{bifilar}$ (Exact)'},
        'Wires':    {'data': (freq.get('ana'), cap_wires), 'label': r'$c_{n+1 \; wires}$'},
        'Bifilar':  {'data': (freq.get('ana'), cap_approx), 'label': r'$c_{bifilar}$ (Approx.)'},
    }

    # Cria a figura com subplots
    fig, ax = plt.subplots(figsize=(8, 5))
    _configure_plot(ax, 'Capacitance p.u.l. (nF/km)', capacitante_data, yscale='linear')
    plt.tight_layout()


def main():
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    clear_screen()
    print("Iniciando cálculos da impedância p.u.l. ...")
    start_time = time.time()

    MTLRepresentation(MTL).bare_and_coated_wires()
    mom_data = run_mom(autoPlots=False)
    analytical_data = run_analytical(MTL, FREQUENCY_RANGE['ana'])
    mom_so_data = run_mom_so(MTL, FREQUENCY_RANGE['mom'])
    
    print(f"\nRotinas de cálculo finalizadas em {(time.time()-start_time):.2f} segundos.")
    plot_resistance_results(FREQUENCY_RANGE, analytical_data, mom_so_data)
    plot_inductance_results(FREQUENCY_RANGE, analytical_data, mom_so_data)
    plot_capacitance_results(FREQUENCY_RANGE, analytical_data, mom_so_data, mom_data)


if __name__ == "__main__":
    main()
    plt.show()