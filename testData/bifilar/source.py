import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any

from utils.case_utils import *
from mtl_main.graphics import IsolatedMTLRepresentation
from mtl_main.source import MulticonductorTransmissionLine
from mtl_paul.py_fortran import FortranRunner
from analytical_forms.isolated_wires import WiresHomogeneousMedia
from mom_so.quasi_static_green import QuasiStatic
from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing
from mom.bare_wire_systems import BareWireMoMSolver
from plotter.mom_models import MoMVisualizer

class BifilarBareWirePULParameters:
    """
    The ConvergenceAnalyzer class is a tool designed to perform and visualize a convergence analysis 
    for the electrical parameters of multiconductor transmission lines (MTLs). 
    
    It systematically runs simulations with an increasing number of Fourier coefficients to observe 
    how the calculated inductance and capacitance values stabilize. The class facilitates a comparison
    between a Fortran-based simulation (RIBBON.FOR) and a Python-based Method of Moments (MoM) implementation.
    """
    def __init__(self, project_root: Path, case_name: str, mtl: Dict[str, Any], comsol_data: Dict[str, pd.DataFrame] = {}):
        """
        Inicializa o analisador de convergência.

        Args:
            mtl_config (Dict[str, Any]): Dicionário com a configuração do modelo MTL.
            nf_max (int): Número máximo de coeficientes/ordem harmônica para testar.
        """

        assert len([key for key in mtl.keys() if isinstance(key, int)]) == 2, "A linha bifilar deve conter exatamente dois condutores."

        self.project_root = project_root
        self.case_name = case_name
        self.comsol_data = comsol_data
        self.mtl_copy = copy.deepcopy(mtl)
        self.freq_range = {'ana': np.logspace(0, 6, num=200),   'mom': np.logspace(0, 6, num=31)}
        self.srw_ratios = {'ana': np.linspace(2.1, 8, num=300), 'mom': np.linspace(2.1, 8, num=20)}
        self.N = len([key for key in mtl.keys() if isinstance(key, int)])
        self.results_df = None

        # Extrai parâmetros e prepara o executor do Fortran
        self._bifilar_analytical_solution()

        # Parâmetros de dados
        self.srw_data = {}
        self.srw_mum_data = {}
        self.analytical_data = {}
        self.mom_collocation_data = {}
        self.mom_galerkin_data = {}
        self.mom_so_data = {}
        self.ribbon_data = {}

        # Parâmetros adicionais
        self.c_factor = 1e12  # F/m to nF/km
        self.l_factor = 1e6   # H/m to mH/km
        self.r_factor = 1e3   # Ohm/m to Ohm/km
        self.figsize = (12, 5)
        self.pt1 = 63
        self.pt2 = 10

        # Parâmetros de plotagem
        self.plot_params = {
            'linestyles': [':', '-.', '--', '-', ':', '-.', '--'],
            'markers': ['o', 's', '^', 'd', 'v', '<', '>'],
            'colors': ['black', 'gray', 'lightgray', 'darkgray', 'dimgray', 'silver', 'gainsboro']
        }

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.case_name, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

    def _bifilar_analytical_solution(self):
        """Calcula a solução analítica para fios nus como referência."""
        R = self.mtl_copy[0]['radius'][1]
        D = np.linalg.norm(np.array(self.mtl_copy[0]['center_point']) - np.array(self.mtl_copy[1]['center_point']))
        self.DR_ratio = D/R
        from scipy.constants import epsilon_0
        self.analytical_bifilar_capacitance = (np.pi * epsilon_0) / np.arccosh(0.5*self.DR_ratio)

    def _configure_plot_appearance(self, ax, ylabel, data_to_plot, yscale='log'):
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
            freq, value = data['data']
            label = data['label']
            if freq is not None and value is not None:
                if key == 'mom-so':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=3, zorder=2)
                elif key == 'mom':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='s', markersize=10, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'ribbon':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=8, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'analytical':
                    ax.plot(freq, value, label=label, color='k', linestyle='--', linewidth=1.0)
                elif key == 'wires':
                    ax.plot(freq, value, label=label, color='k', linestyle='-.', linewidth=1.0)
                elif key == 'bifilar':
                    ax.plot(freq, value, label=label, color='r', linestyle=':', linewidth=1.0)
                elif key == 'exactly':
                    ax.plot(freq, value, label=label, color='k', linestyle='-', linewidth=1.0)
                elif key == 'approx':
                    ax.plot(freq, value, label=label, color='k', linestyle=':', linewidth=1.0)

        ax.set_xscale('log')
        ax.set_yscale(yscale)
        ax.set_xlim(1E0, 1E6)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(False)

    def _configure_plot_matrix_appearance(self, ax, ylabel, data_to_plot, yscale='log', ylim=None):
        """
        Função auxiliar para configurar um subplot para os elementos da matriz de resistência (R11, R12, R22).

        Args:
            ax (matplotlib.axes.Axes): O eixo do subplot a ser configurado.
            ylabel (str): Rótulo do eixo Y.
            data_to_plot (dict): Dados para plotagem (espera chaves 'r11', 'r12', 'r22').
            yscale (str, optional): Escala do eixo Y ('log' or 'linear').
        """
        
        # Usar os parâmetros de plotagem definidos na classe para consistência
        ls = self.plot_params['linestyles']
        mk = self.plot_params['markers']
        cl = self.plot_params['colors']

        styles = {
            '11': {'color': cl[0], 'linestyle': ls[3], 'marker': mk[0], 'markersize': 6, 'fillstyle': 'none'}, # Preto, Sólido, Círculo
            '22': {'color': cl[1], 'linestyle': ls[2], 'marker': mk[1], 'markersize': 8, 'fillstyle': 'none'}, # Cinza, Tracejado, Quadrado
            '12': {'color': cl[2], 'linestyle': ls[0], 'marker': mk[2], 'markersize': 7, 'fillstyle': 'none'}, # Cinza claro, Pontilhado, Triângulo
        }

        for key, data_dict in data_to_plot.items():
            freq, value = data_dict['data']
            label = data_dict['label']
            
            if freq is not None and value is not None:
                if key in styles:
                    style = styles[key]
                    ax.plot(freq, value, 
                            label=label, 
                            color=style['color'], 
                            linestyle=style['linestyle'], 
                            marker=style['marker'],
                            markersize=style['markersize'],
                            fillstyle=style['fillstyle'],
                            zorder=2)
                else:
                    ax.plot(freq, value, label=label)

        ax.set_xscale('log')
        ax.set_yscale(yscale)
        ax.set_xlim(1E0, 1E6)
        if ylim is not None:
            ax.set_ylim(ylim)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(False)
        
    def show_header(self):
        """Exibe o cabeçalho do script."""
        print("\n")
        print("="*self.pt2 + " BIFILAR BARE-WIRE RIBBON CABLE SIMULATION " + "="*self.pt2)
        print(f"Project: {self.project_root}")
        print(f"Model: {self.mtl_copy['name']}")
        print(f"D/R = {self.DR_ratio:.3f}. Fourier Order (k) = {self.mtl_copy[0]['fourier_order']}.")
        print(f"RIBBON Fourier Coef./cond. (NF) = {self.mtl_copy[0]['fourier_order']+1}.")
        print(f"PYTHON Fourier Coef./cond. (NF) = {2*self.mtl_copy[0]['fourier_order']+1}.")
        print(f"Exact Bifilar Bare Wire Capacitance: {self.analytical_bifilar_capacitance * 1E12:.4f} pF/m")
        print("="*self.pt1)

    def run_analytical(self):
        """
        Executa a simulação analítica da impedância da linha de transmissão.

        Args:
            mtl_config (dict): Dicionário de configuração da linha de transmissão.
            frequencies (np.ndarray): Array de frequências para a análise.

        Returns:
            tuple: Uma tupla contendo três listas: impedâncias série,
                resistências de alta frequência e indutâncias externas.
        """
        print("\n==============         Analytical Processing       =============")

        bifilar_wires = WiresHomogeneousMedia(self.mtl_copy)
        le_wires = bifilar_wires.n_wires_inductance_matrix()
        c_wires = bifilar_wires.n_wires_capacitance_matrix(le_wires)
        pul_bifilar = bifilar_wires.bifilar_pul_inductance_and_capacitance()

        for freq in self.freq_range['ana']:
            z_s, r_hf = bifilar_wires.bifilar_pul_series_impedance(freq)
            self.analytical_data[freq] = {
                'rhf':          self.r_factor * r_hf,
                'rs':           self.r_factor * np.real(z_s),
                'ls':           self.l_factor * np.imag(z_s) / (2 * np.pi * freq),
                'le_wires':     self.l_factor * le_wires.item(),
                'le_exact':     self.l_factor * pul_bifilar['inductance']['exact'],
                'le_bifilar':   self.l_factor * pul_bifilar['inductance']['approximate'],
                'c_exact':      self.c_factor * pul_bifilar['capacitance']['exact'],
                'c_approx':     self.c_factor * pul_bifilar['capacitance']['approximate'],
                'c_wires':      self.c_factor * c_wires.item(),
            }

    def run_mom_so(self):
        """
        Executa a simulação da impedância usando o Método dos Momentos (MoM-SO).

        Args:
            mtl_config (dict): Dicionário de configuração da linha de transmissão.
            frequencies (np.ndarray): Array de frequências para a análise.
            green_mode (GreenFunctionMode): O modo de cálculo para a função de Green.

        Returns:
            list: Uma lista contendo as impedâncias série totais calculadas via MoM.
        """
        print("\n=============   MoM-SO HomogeneousLosslessMedium   =============")

        model = MulticonductorTransmissionLine(self.mtl_copy)
        freq = self.freq_range['mom']

        # The Green's matrix is frequency-independent
        green_matrix = QuasiStatic(model).green_matrix()

        # Instantiate the model ONCE with the full array of numerical frequencies
        mom_so = HomogeneousLosslessMedium(model, freq)
        post_processor = LosslessPostProcessing(model)

        # Calculate partial impedance for all frequencies
        z_partial = mom_so.z_partial(green_matrix)

        # Calculate total series impedance for all frequencies
        zs_stack = post_processor.z_total(z_partial)

        # Calculate series resistance and inductance for all frequencies
        rs_stack = post_processor.rs_matrix(zs_stack)
        ls_stack = post_processor.ls_matrix(zs_stack, freq)

        for i, freq in enumerate(freq):
            self.mom_so_data[freq] = {
                'zp': z_partial[i, :, :],
                'zs': zs_stack[i].item(),
                'rs': rs_stack[i].item() * self.r_factor,
                'ls': ls_stack[i].item() * self.l_factor,
            }

    def plot_impedance_results(self):
        """
        Gera e exibe os gráficos dos resultados da simulação de forma flexível,
        organizados em subplots.

        Args:
            mtl (dict): Dicionário de configuração da linha de transmissão.
            freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
            analytical (dict): Dicionário com os resultados da simulação analítica.
            mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
        """
        freq_mom = self.freq_range.get('mom')
        freq_ana = self.freq_range.get('ana')
        data_mom = self.mom_so_data.values()
        data_ana = self.analytical_data.values()

        resistance_data = {
            'mom-so':     {'data': (freq_mom, [data['rs']  for data in data_mom]), 'label': r'MoM-SO'},
            'analytical': {'data': (freq_ana, [data['rs']  for data in data_ana]), 'label': r'$R_i$'},
            'approx':     {'data': (freq_ana, [data['rhf'] for data in data_ana]), 'label': r'$R_{HF}$'},
        }

        inductance_data = {
            'mom-so':     {'data': (freq_mom, [data['ls']       for data in data_mom]), 'label': 'MoM-SO'},
            'exactly':    {'data': (freq_ana, [data['le_exact'] for data in data_ana]), 'label': r'$\ell_e$'},
            'analytical': {'data': (freq_ana, [data['ls']       for data in data_ana]), 'label': r'$\ell_s$'},
        }
        
        fig1, ax1 = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax1, r'Series Resistance p.u.l. ($\Omega$/km)', resistance_data)

        fig2, ax2 = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax2, 'Series Inductance p.u.l. (mH/km)', inductance_data, yscale='linear')
        
        save_figure(fig1, self.results_dir, base_filename='pul_series_resistance')
        save_figure(fig2, self.results_dir, base_filename='pul_series_inductance')
        plt.tight_layout()

    def plot_partial_impedance_matrix(self):
        """
        Gera e exibe os gráficos dos resultados da simulação da matriz de 
        resistência (R11, R12, R22).
        """
        freq = self.freq_range.get('mom')
        w = 2 * np.pi * freq
        zp_matrices = [data['zp'] for data in self.mom_so_data.values()]

        # Extrai os componentes da matriz (assumindo matriz 2x2)
        r11 = [np.real(zp[0, 0]) * self.r_factor for zp in zp_matrices]
        r12 = [np.real(zp[0, 1]) * self.r_factor for zp in zp_matrices]
        r22 = [np.real(zp[1, 1]) * self.r_factor for zp in zp_matrices]
        x11 = [np.imag(zp[0, 0]) * self.r_factor for zp in zp_matrices]
        x12 = [np.imag(zp[0, 1]) * self.r_factor for zp in zp_matrices]
        z22 = [np.imag(zp[1, 1]) * self.r_factor for zp in zp_matrices]
        l11 = np.array(x11) / w
        l12 = np.array(x12) / w 
        l22 = np.array(z22) / w

        resistance_matrix_data = {
            '11': {'data': (freq, r11), 'label': '$R_{11}$'},
            '22': {'data': (freq, r22), 'label': '$R_{22}$'},
            '12': {'data': (freq, r12), 'label': '$R_{12} = R_{21}$'},
        }

        reactance_matrix_data = {
            '11': {'data': (freq, x11), 'label': '$L_{11}$'},
            '22': {'data': (freq, z22), 'label': '$L_{22}$'},
            '12': {'data': (freq, x12), 'label': '$L_{12} = L_{21}$'},
        }

        inductance_matrix_data = {
            '11': {'data': (freq, l11), 'label': '$L_{11}$'},
            '22': {'data': (freq, l22), 'label': '$L_{22}$'},
            '12': {'data': (freq, l12), 'label': '$L_{12} = L_{21}$'},
        }

        fig1, ax1 = plt.subplots(figsize=self.figsize)
        fig2, ax2 = plt.subplots(figsize=self.figsize)
        fig3, ax3 = plt.subplots(figsize=self.figsize)

        self._configure_plot_matrix_appearance(
            ax1, 
            r'Partial Resistance Matrix p.u.l. ($\Omega$/km)', 
            resistance_matrix_data,
            yscale='log',
            ylim=(1e-2, 1e1)
        )

        self._configure_plot_matrix_appearance(
            ax2, 
            r'Partial Reactance Matrix p.u.l. ($\Omega$/km)', 
            reactance_matrix_data,
            yscale='log',
        )

        self._configure_plot_matrix_appearance(
            ax3, 
            r'Partial Inductance Matrix p.u.l. (mH/km)', 
            inductance_matrix_data,
            yscale='linear',
        )

        save_figure(fig1, self.results_dir, base_filename='pul_series_resistance_matrix')
        save_figure(fig2, self.results_dir, base_filename='pul_series_reactance_matrix')
        save_figure(fig3, self.results_dir, base_filename='pul_series_inductance_matrix')
        plt.tight_layout()

    def print_impedance_matrix(self):
        """
        Imprime os resultados da matriz de impedância (Zp) no terminal,
        frequência a frequência, formatado como no exemplo da figura.
        """
        print("\n--- Impedance Matrix (Zp) Results ---\n")
        header = f"{'freq (Hz)':<12} {'z11 (Ω/m)':<30} {'z12 (Ω/m)':<30} {'z22 (Ω/m)':<30}"
        print(header)
        print("-" * len(header))

        freq = self.freq_range.get('mom')
        zp_matrices = [data['zp'] for data in self.mom_so_data.values()]

        if (freq is None or freq.size == 0) or not zp_matrices:
            print("No impedance matrix data available to display.")
            return

        for i, freq in enumerate(freq):
            if i < len(zp_matrices):
                z11 = format_complex_number(zp_matrices[i][0, 0])
                z12 = format_complex_number(zp_matrices[i][0, 1])
                z22 = format_complex_number(zp_matrices[i][1, 1])
                print(f"{freq:<12.2E} {z11:<30} {z12:<30} {z22:<30}")
            else:
                print(f"Warning: No impedance matrix found for frequency {freq} Hz.")
                