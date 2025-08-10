import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any

from mtl_paul.py_fortran import FortranRunner
from mom.bare_wire_systems import MulticonductorBareWireSystems
from analytical_forms.bare_wires import WiresHomogeneousMedia
from mom_so.quasi_static_green import QuasiStatic
from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing


class BifilarBareWirePULParameters():
    """
    Encapsula a lógica para executar e analisar o estudo de convergência
    de capacitância, comparando MoM Python e Fortran.
    """
    def __init__(self, project_root: Path, mtl: Dict[str, Any]):
        """
        Inicializa o analisador de convergência.

        Args:
            mtl_config (Dict[str, Any]): Dicionário com a configuração do modelo MTL.
            nf_max (int): Número máximo de coeficientes/ordem harmônica para testar.
        """

        assert len([key for key in mtl.keys() if isinstance(key, int)]) == 2, "A linha bifilar deve conter exatamente dois condutores."

        self.project_root = project_root
        self.mtl_copy = copy.deepcopy(mtl)
        self.freq_range = {'ana': np.logspace(0, 6, num=200), 'mom': np.logspace(0, 6, num=30)}
        self.srw_ratios = {'ana': np.linspace(2.1, 8, num=300), 'mom': np.linspace(2.1, 8, num=40)}

        # Extrai parâmetros e prepara o executor do Fortran
        self._bifilar_analytical_solution()

        # Parâmetros de dados
        self.analytical_data = {}
        self.mom_data = {}
        self.mom_so_data = {}
        self.ribbon_data = {}

        # Parâmetros adicionais
        self.c_factor = 1e12  # F/m to nF/km
        self.l_factor = 1e6   # H/m to mH/km
        self.r_factor = 1e3   # Ohm/m to Ohm/km
        self.pt1 = 63
        self.pt2 = 10


    def _prepare_fortran_runner(self, mtl):
        """Prepara os parâmetros e o executor para a simulação Fortran."""

        # Obtenha o dicionário 'sheath' de forma segura.
        #    Se 'sheath' não existir ou for None, use um dicionário vazio {} como fallback.
        refIdx = mtl['idx_ref_conductor']
        sheath_dict = mtl[refIdx].get('sheath') or {}

        self.fortran_base_params = {
            'N':    2,
            'NF':   mtl[refIdx]['fourier_order'] + 1,
            'IREF': refIdx,
            'RW':   mtl[refIdx]['radius'][1],
            'TD':   sheath_dict.get('thickness', 0.0),
            'ER':   sheath_dict.get('relative_permittivity', 1.0),
            'S':    np.linalg.norm(np.array(mtl[0]['center_point']) - np.array(mtl[1]['center_point'])),
        }

        fortran_exe_path = self.project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.runner = FortranRunner(exe_path=str(fortran_exe_path), silent=True)


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


    def show_header(self):
        """Exibe o cabeçalho do script."""
        print("\n")
        print("="*self.pt2 + " BIFILAR BARE-WIRE RIBBON CABLE SIMULATION " + "="*self.pt2)
        print(f"Project: {self.project_root}")
        print(f"Model: {self.mtl_copy['name']}")
        print(f"D/R = {self.DR_ratio}. Fourier Order (k) = {self.mtl_copy[0]['fourier_order']}.")
        print(f"RIBBON Fourier Coef./cond. (NF) = {self.mtl_copy[0]['fourier_order']+1}.")
        print(f"PYTHON Fourier Coef./cond. (NF) = {2*self.mtl_copy[0]['fourier_order']+1}.")
        print(f"Exact Bifilar Bare Wire Capacitance: {self.analytical_bifilar_capacitance * 1E12:.4f} pF/m")
        print("="*self.pt1)


    def run_fortran(self, displayTerminal: bool = True):
        """Executa uma simulação única para um valor específico de k."""
        self._prepare_fortran_runner(self.mtl_copy)
        self.runner.run_fortran(self.fortran_base_params)

        self.ribbon_data = {
            freq: {'c': self.c_factor * self.runner.C0_matrix[0,0], 
                   'le': self.l_factor * self.runner.L_matrix[0,0]} for freq in self.freq_range['mom']}

        if displayTerminal:
            print("\n")
            print("="*self.pt2 + "                 RIBBON.FOR                " + "="*self.pt2)
            if self.runner.L_matrix is not None:
                print("\nMatriz de Indutância Externa (Le):")
                print(self.runner.L_matrix)

                print("\nMatriz de Capacitância (C):")
                print(self.runner.C_matrix)

                print("\nMatriz de Capacitância no Vácuo (C0):")
                print(self.runner.C0_matrix)

                print("\nMatriz de Capacitância Generalizada (CGEN):")
                print(self.runner.CGEN_matrix)
            else:
                print("\nNenhum resultado foi analisado. Verifique os logs de erro.")


    def run_py_mom(self, autoPlots=False):
        """
        Executa a simulação clássica do Método dos Momentos (MoM) para a linha de transmissão bifilar.

        Returns:
            BifilarMoM: Instância do objeto BifilarMoM configurado.
        """
        print("\n============== pyMoM MulticonductorBareWireSystems =============")

        mom_wires = MulticonductorBareWireSystems(self.mtl_copy)
        mom_wires.run_simulation()
        mom_wires.print_results()

        self.mom_data = {freq: {'c': self.c_factor * mom_wires.C_maxwellian.item()} for freq in self.freq_range['mom']}

        if autoPlots:
            mom_wires.plot_charge_density()
            mom_wires.plot_collocation_points()
            mom_wires.plot_harmonic_coefficients()
            MulticonductorBareWireSystems.plot_convergence_rates(self.mtl_copy, nf_max=20)


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

        green_matrix = QuasiStatic(self.mtl_copy).green_matrix()
        post_processor = LosslessPostProcessing(self.mtl_copy)

        for freq in self.freq_range['mom']:
            mom_so = HomogeneousLosslessMedium(self.mtl_copy, freq)
            zs = post_processor.z_total(mom_so.z_partial(green_matrix))
            general_cap = mom_so.generalized_capacitance_matrix(green_matrix)
            maxwell_cap = mom_so.maxwellian_capacitance_matrix(general_cap)
            
            self.mom_so_data[freq] = {
                'zs':   zs,
                'rs':   self.r_factor * post_processor.rs_matrix(zs).item(),
                'ls':   self.l_factor * post_processor.ls_matrix(zs, freq).item(),
                'c':    self.c_factor * np.real(maxwell_cap.item()),
            }

        if green_matrix.shape[0] < 6:
            from scipy.constants import epsilon_0
            print(f"\nGreen's Matrix (Dim: {green_matrix.shape}):\n{green_matrix}")
            print(f"\n2*pi*e0*G:\n{- 2 * np.pi * epsilon_0 * np.real(green_matrix)}")

        print(f"\nGreen's Matrix Shape: {green_matrix.shape}.")
        print(f"\nMoM-SO Generalized Capacitance Matrix (F/m): \n{np.real(general_cap)}")
        print(f"\nMoM-SO Bifilar Capacitance: \n{np.real(maxwell_cap.item()) * 1E12:.4f} pF/m")


    def srw_rates_analytical(self):
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

        for ratio in self.srw_ratios['ana']:
            # Cria a configuração da linha bifilar dinamicamente para cada razão.
            separation = ratio * self.mtl_copy[0]['radius'][1]
            self.mtl_copy[1]['center_point'] = (separation, 0.0)

            wires = WiresHomogeneousMedia(self.mtl_copy)
            le_wires = wires.n_wires_inductance_matrix()
            pul_bifilar = wires.bifilar_pul_inductance_and_capacitance()

            self.analytical_data[ratio] = {
                'le_wires':     self.l_factor * le_wires,
                'le_exact':     self.l_factor * pul_bifilar['inductance']['exact'],
                'le_bifilar':   self.l_factor * pul_bifilar['inductance']['approximate'],
                'c_exact':      self.c_factor * pul_bifilar['capacitance']['exact'],
                'c_approx':     self.c_factor * pul_bifilar['capacitance']['approximate'],
                'c_wires':      self.c_factor * wires.n_wires_capacitance_matrix(le_wires),
            }


    def srw_rates_py_mom(self):
        """
        Executa a simulação clássica do Método dos Momentos (MoM) para a linha de transmissão bifilar.

        Returns:
            BifilarMoM: Instância do objeto BifilarMoM configurado.
        """
        print("\n==============             RIBBON.FOR              =============")
        print("\n============== pyMoM MulticonductorBareWireSystems =============")
        print("\n==============  MoM-SO HomogeneousLosslessMedium   =============")

        # Cria a configuração da linha bifilar dinamicamente para cada razão.
        for ratio in self.srw_ratios['mom']:
            separation = ratio * self.mtl_copy[0]['radius'][1]
            self.mtl_copy[1]['center_point'] = (separation, 0.0)

            mom_wires = MulticonductorBareWireSystems(self.mtl_copy)
            mom_wires.run_simulation()
            
            self._prepare_fortran_runner(self.mtl_copy)
            self.runner.run_fortran(self.fortran_base_params)

            green_matrix = QuasiStatic(self.mtl_copy).green_matrix()
            mom_so = HomogeneousLosslessMedium(self.mtl_copy, frequency=0)
            gen_cap = mom_so.generalized_capacitance_matrix(green_matrix)
            capacitance = mom_so.maxwellian_capacitance_matrix(gen_cap)

            self.mom_so_data[ratio] = {'c': self.c_factor * np.real(capacitance.item())}
            self.ribbon_data[ratio] = {'c': self.c_factor * self.runner.C0_matrix[0, 0]} 
            self.mom_data[ratio]    = {'c': self.c_factor * mom_wires.C_maxwellian.item()}


    def plot_resistance_results(self):
        """
        Gera e exibe os gráficos dos resultados da simulação de forma flexível,
        organizados em subplots.

        Args:
            mtl (dict): Dicionário de configuração da linha de transmissão.
            freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
            analytical (dict): Dicionário com os resultados da simulação analítica.
            mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
        """
        resistance_data = {
            'mom-so':       {'data': (self.freq_range.get('mom'), [data['rs']  for data in self.mom_so_data.values()]),     'label': 'MoM-SO'},
            'analytical':   {'data': (self.freq_range.get('ana'), [data['rs']  for data in self.analytical_data.values()]), 'label': '$R_i$'},
            'approx':       {'data': (self.freq_range.get('ana'), [data['rhf'] for data in self.analytical_data.values()]), 'label': '$R_{HF}$'},
        }

        fig, ax = plt.subplots(figsize=(8, 5))
        self._configure_plot_appearance(ax, r'Series Resistance p.u.l. ($\Omega$/km)', resistance_data)
        plt.tight_layout()


    def plot_inductance_results(self):
        """
        Gera e exibe os gráficos dos resultados da simulação de forma flexível,
        organizados em subplots.

        Args:
            mtl (dict): Dicionário de configuração da linha de transmissão.
            freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
            analytical (dict): Dicionário com os resultados da simulação analítica.
            mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
        """
        inductance_data = {
            'ribbon':       {'data': (self.freq_range.get('mom'), [data['le']           for data in self.ribbon_data.values()]),     'label': r'$\ell_{e,RIBBON.FOR}$'},
            'mom-so':       {'data': (self.freq_range.get('mom'), [data['ls']           for data in self.mom_so_data.values()]),     'label': 'MoM-SO'},
            'wires':        {'data': (self.freq_range.get('ana'), [data['le_wires']     for data in self.analytical_data.values()]), 'label': r'$\ell_{e,n+1 \; wires}$'},
            'exactly':      {'data': (self.freq_range.get('ana'), [data['le_exact']     for data in self.analytical_data.values()]), 'label': r'$\ell_{e,Bifilar (Exactly)}$'},
            'bifilar':      {'data': (self.freq_range.get('ana'), [data['le_bifilar']   for data in self.analytical_data.values()]), 'label': r'$\ell_{e,Bifilar (Approx.)}$'},
            'analytical':   {'data': (self.freq_range.get('ana'), [data['ls']           for data in self.analytical_data.values()]), 'label': '$\ell_s$'},
        }

        fig, ax = plt.subplots(figsize=(8, 5))
        self._configure_plot_appearance(ax, 'Series Inductance p.u.l. (mH/km)', inductance_data, yscale='linear')
        plt.tight_layout()


    def plot_capacitance_results(self):
        """
        Gera e exibe os gráficos dos resultados da simulação de forma flexível,
        organizados em subplots.

        Args:
            mtl (dict): Dicionário de configuração da linha de transmissão.
            freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
            analytical (dict): Dicionário com os resultados da simulação analítica.
            mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
        """
        capacitante_data = {
            'mom':      {'data': (self.freq_range.get('mom'), [data['c']        for data in self.mom_data.values()]),        'label': 'MoM'},
            'ribbon':   {'data': (self.freq_range.get('mom'), [data['c']        for data in self.ribbon_data.values()]),     'label': 'RIBBON.FOR'},
            'mom-so':   {'data': (self.freq_range.get('mom'), [data['c']        for data in self.mom_so_data.values()]),     'label': 'MoM-SO'},
            'wires':    {'data': (self.freq_range.get('ana'), [data['c_wires']  for data in self.analytical_data.values()]), 'label': 'n+1 wires'},
            'exactly':  {'data': (self.freq_range.get('ana'), [data['c_exact']  for data in self.analytical_data.values()]), 'label': 'Bifilar (Exactly)'},
            'bifilar':  {'data': (self.freq_range.get('ana'), [data['c_approx'] for data in self.analytical_data.values()]), 'label': 'Bifilar (Approx.)'},
        }

        fig, ax = plt.subplots(figsize=(8, 5))
        self._configure_plot_appearance(ax, 'Capacitance p.u.l. (nF/km)', capacitante_data, yscale='linear')
        plt.tight_layout()


    def plot_paul_fig49(self):
        """
        Gera e exibe os gráficos dos resultados da simulação de forma flexível,
        organizados em subplots.

        Args:
            mtl (dict): Dicionário de configuração da linha de transmissão.
            freqs (dict): Dicionário contendo os arrays de frequência para cada simulação.
            analytical (dict): Dicionário com os resultados da simulação analítica.
            mom_so (dict): Dicionário com os resultados da simulação MoM-SO.
        """

        capacitante_data = {
            'mom':      {'data': (self.srw_ratios.get('mom'), [data['c']        for data in self.mom_data.values()]),        'label': 'MoM'},
            'ribbon':   {'data': (self.srw_ratios.get('mom'), [data['c']        for data in self.ribbon_data.values()]),     'label': 'RIBBON.FOR'},
            'mom-so':   {'data': (self.srw_ratios.get('mom'), [data['c']        for data in self.mom_so_data.values()]),     'label': 'MoM-SO'},
            'exactly':  {'data': (self.srw_ratios.get('ana'), [data['c_exact']  for data in self.analytical_data.values()]), 'label': 'Exactly'},
            'approx':   {'data': (self.srw_ratios.get('ana'), [data['c_approx'] for data in self.analytical_data.values()]), 'label': 'Approx.'},
        }

        fig, ax = plt.subplots(figsize=(8, 5))
        for key, data in capacitante_data.items():
            freq, value = data['data']
            label = data['label']
            if freq is not None and value is not None:
                if key == 'ribbon':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=8, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'mom':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='s', markersize=10, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'mom-so':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=3, zorder=2)
                elif key == 'approx':
                    ax.plot(freq, value, label=label, color='k', linestyle='--', linewidth=1.0, zorder=1)
                elif key == 'exactly':
                    ax.plot(freq, value, label=label, color='k', linestyle=':', linewidth=1.0, zorder=1)

        ax.set_xlabel('Ratio of separation to wire radius, s/r$_w$')
        ax.set_ylabel('Per-unit-length capacitance (pF/m)')
        ax.set_xlim(2, 8)
        ax.set_ylim(10, 90)
        ax.legend()
        ax.grid(True, linestyle='--', linewidth=0.5)
        plt.tight_layout()


class BifilarBareWireConvergence():
    """
    The ConvergenceAnalyzer class is a tool designed to perform and visualize a convergence analysis 
    for the electrical parameters of multiconductor transmission lines (MTLs). 
    
    It systematically runs simulations with an increasing number of Fourier coefficients to observe 
    how the calculated inductance and capacitance values stabilize. The class facilitates a comparison
    between a Fortran-based simulation (RIBBON.FOR) and a Python-based Method of Moments (MoM) implementation.
    """
    def __init__(self, project_root: Path, mtl: Dict[str, Any], SUM_MAX: int = 10):
        """
        Inicializa o analisador de convergência.

        Args:
            mtl_config (Dict[str, Any]): Dicionário com a configuração do modelo MTL.
            nf_max (int): Número máximo de coeficientes/ordem harmônica para testar.
        """
        self.project_root = project_root
        self.mtl_copy = copy.deepcopy(mtl)
        self.N = len([key for key in mtl.keys() if isinstance(key, int)])

        self.sum_max = SUM_MAX
        self.results_df = None
        
        # Parâmetros adicionais
        self.c_factor = 1e12  # F/m to nF/km
        self.l_factor = 1e6   # H/m to mH/km
        self.r_factor = 1e3   # Ohm/m to Ohm/km
        self.pt1 = 63
        self.pt2 = 10

        # Parâmetros de plotagem
        self.plot_params = {
            'linestyles': [':', '-.', '--', '-', ':', '-.', '--'],
            'markers': ['o', 's', '^', 'd', 'v', '<', '>'],
            'colors': ['black', 'gray', 'lightgray', 'darkgray', 'dimgray', 'silver', 'gainsboro']
        }


    def _prepare_fortran_runner(self, mtl):
        """Prepara os parâmetros e o executor para a simulação Fortran."""

        # Obtenha o dicionário 'sheath' de forma segura.
        #    Se 'sheath' não existir ou for None, use um dicionário vazio {} como fallback.
        refIdx = mtl['idx_ref_conductor']
        sheath_dict = mtl[refIdx].get('sheath') or {}

        self.fortran_base_params = {
            'N':    len([key for key in mtl.keys() if isinstance(key, int)]),
            'NF':   mtl[refIdx]['fourier_order'] + 1,
            'IREF': refIdx,
            'RW':   mtl[refIdx]['radius'][1],
            'TD':   sheath_dict.get('thickness', 0.0),
            'ER':   sheath_dict.get('relative_permittivity', 1.0),
            'S':    np.linalg.norm(np.array(mtl[0]['center_point']) - np.array(mtl[1]['center_point'])),
        }
        
        fortran_exe_path = self.project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.runner = FortranRunner(exe_path=str(fortran_exe_path), silent=False)
        self.runner_silent = FortranRunner(exe_path=str(fortran_exe_path), silent=True)


    def _extract_matrix_element(self, column_name: str, row: int, col: int) -> pd.Series:
        """
        Extrai e processa um elemento específico de uma coluna de matrizes no DataFrame de resultados.
        O fator de conversão (para indutância ou capacitância) é determinado automaticamente
        com base no nome da coluna.
        """
        sign = 1.0
        
        # Decide qual fator de conversão usar com base no nome da coluna
        if column_name.startswith('L'):
            factor = self.l_factor
        elif column_name.startswith('C'):
            factor = self.c_factor
            if row != col:
                sign = -1.0
        else:
            # Lança um erro se a coluna não for de Indutância ('L') ou Capacitância ('C')
            raise ValueError(f"Não foi possível determinar o fator de conversão para a coluna: '{column_name}'")
        
        extractor = lambda matrix: (
            sign * matrix[row, col] * factor
            if isinstance(matrix, np.ndarray)
            else np.nan
        )
        return self.results_df[column_name].apply(extractor)


    def run_convergence(self):
        """Executa o laço de convergência para ambas as simulações e armazena os resultados."""
        print(f"\nRunning Convergence Rate until k = {self.sum_max}!")
        
        results = []
        for k in range(0, self.sum_max):
            temp_mtl = copy.deepcopy(self.mtl_copy)
            for key in temp_mtl.keys():
                if isinstance(key, int):
                    temp_mtl[key]['fourier_order'] = k  
                    if temp_mtl['type'] == 'coated_wires':
                        temp_mtl[key]['sheath']['fourier_order'] = k

            # === Fortran Instance ===
            self._prepare_fortran_runner(temp_mtl)
            self.runner_silent.run_fortran(self.fortran_base_params)

            # === MoM Instance ===
            bare_wire_mtl = copy.deepcopy(self.mtl_copy)
            bare_wire_mtl['type'] = 'bare_wires'
            for key in bare_wire_mtl.keys():
                if isinstance(key, int):
                    bare_wire_mtl[key]['fourier_order'] = k
                    bare_wire_mtl[key]['sheath'] = None

            mom_bare_wires = MulticonductorBareWireSystems(bare_wire_mtl)
            mom_bare_wires.run_simulation()

            # === MoM-SO Instance ===
            green_matrix = QuasiStatic(bare_wire_mtl).green_matrix()
            mom_so = HomogeneousLosslessMedium(bare_wire_mtl, frequency=0)
            general_cap = mom_so.generalized_capacitance_matrix(green_matrix)
            maxwell_cap = mom_so.maxwellian_capacitance_matrix(general_cap)
            
            # Coleta de resultados
            results.append({
                'k': k,
                'L (RIBBON.FOR)':       self.runner_silent.L_matrix if self.fortran_base_params is not None else np.nan,
                'C (RIBBON.FOR)':       self.runner_silent.C_matrix if self.runner_silent.C_matrix is not None else np.nan,
                'C0 (RIBBON.FOR)':      self.runner_silent.C0_matrix if self.runner_silent.C0_matrix is not None else np.nan,
                'C0 (BARE-WIRE.PY)':    mom_bare_wires.C_maxwellian if mom_bare_wires.C_maxwellian is not None else np.nan,
                'C0 (MOM-SO.PY)':       np.real(maxwell_cap) if maxwell_cap is not None else np.nan,
                'CGEN (RIBBON.FOR)':    self.runner_silent.CGEN_matrix if self.runner_silent.CGEN_matrix is not None else np.nan,
                'CGEN (BARE-WIRE.PY)':  mom_bare_wires.C_generalized if mom_bare_wires.C_generalized is not None else np.nan,
                'CGEN (MOM-SO.PY)':     np.real(general_cap) if general_cap is not None else np.nan,
            })
            print(f"  Complete for k = {k}.")

        self.results_df = pd.DataFrame(results).set_index('k')


    def run_single_fortran_simulation(self):
        """Executa uma simulação única para um valor específico de k."""
        
        self._prepare_fortran_runner(self.mtl_copy)
        self.runner.run_fortran(self.fortran_base_params)

        dist1, dist2 = 59, 8
        print("\n" + "="*dist1)
        print("="*dist2 + " THREE-WIRE RIBBON CABLE SYSTEM SIMULATION " + "="*dist2)
        print("="*dist1)
        
        print("\n")
        print("="*dist2 + "                 RIBBON.FOR                " + "="*dist2)
        if self.runner.L_matrix is not None:
            print("\nMatriz de Indutância Externa (Le):")
            print(self.runner.L_matrix)
            
            print("\nMatriz de Capacitância (C):")
            print(self.runner.C_matrix)

            print("\nMatriz de Capacitância no Vácuo (C0):")
            print(self.runner.C0_matrix)

            print("\nMatriz de Capacitância Generalizada (CGEN):")
            print(self.runner.CGEN_matrix)
        else:
            print("\nNenhum resultado foi analisado. Verifique os logs de erro.")


    def plot_paul_fig514a(self):
        """
        Gera o gráfico de convergência da indutância a partir dos resultados armazenados,
        replicando a figura de referência.
        """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_convergence()'.")
            return

        # Garante que a simulação foi executada completamente
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Extrai os componentes da matriz de indutância 'L' e converte para µH/m.
        try:
            # Passa o fator de conversão correto (self.l_factor)
            l11 = self._extract_matrix_element('L (RIBBON.FOR)', row=0, col=0)
            l22 = self._extract_matrix_element('L (RIBBON.FOR)', row=1, col=1)
            l12 = self._extract_matrix_element('L (RIBBON.FOR)', row=0, col=1)
        except (TypeError, IndexError) as e:
            # Corrige as mensagens de erro para o contexto de indutância
            print(f"Erro ao extrair elementos da matriz de indutância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'L (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))
        
        fortran_nf_axis = self.results_df.index + 1
        ax.plot(fortran_nf_axis, l22, color='k', marker='o', linestyle='-',  label='$L_{22}$')
        ax.plot(fortran_nf_axis, l11, color='k', marker='o', linestyle=':',  label='$L_{11}$', markerfacecolor='white', markeredgecolor='k')
        ax.plot(fortran_nf_axis, l12, color='k', marker='s', linestyle='-.', label='$L_{12}$')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Inductance (µH/m)', fontsize=12)
        
        # Adiciona uma pequena margem (padding) aos limites do eixo x para melhor visualização
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))        
        ax.set_ylim(0.2, 1.1)
        ax.set_yticks(np.arange(0.2, 1.2, 0.1))        
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=3, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()


    def plot_paul_fig514b(self):
        """
        Gera o gráfico de convergência a partir dos resultados armazenados,
        replicando a figura de referência.
        """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_convergence()'.")
            return

        # Garante que a simulação foi executada completamente
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Extrai os componentes da matriz de indutância 'L' e converte para µH/m.
        try:
            # Passa o fator de conversão correto (self.l_factor)
            c11 = self._extract_matrix_element('C (RIBBON.FOR)', row=0, col=0)
            c22 = self._extract_matrix_element('C (RIBBON.FOR)', row=1, col=1)
            c12 = self._extract_matrix_element('C (RIBBON.FOR)', row=0, col=1)
        except (TypeError, IndexError) as e:
            # Corrige as mensagens de erro para o contexto de indutância
            print(f"Erro ao extrair elementos da matriz de indutância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'L (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        fortran_nf_axis = self.results_df.index + 1
        ax.plot(fortran_nf_axis, c11, color='k', marker='o', linestyle='-',  label='$C_{11}$')
        ax.plot(fortran_nf_axis, c22, color='k', marker='o', linestyle=':',  label='$C_{22}$', markerfacecolor='white', markeredgecolor='k')
        ax.plot(fortran_nf_axis, c12, color='k', marker='s', linestyle='-.', label='$C_{12}$')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_yticks(np.arange(15, 40, 5))
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_ylim(15, 40)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=3, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()


    def plot_paul_fig514c(self):
        """
        Gera o gráfico de convergência a partir dos resultados armazenados,
        replicando a figura de referência.
        """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_convergence()'.")
            return

        # Garante que a simulação foi executada completamente
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Extrai os componentes da matriz de indutância 'L' e converte para µH/m.
        try:
            # Passa o fator de conversão correto (self.l_factor)
            c0_11 = self._extract_matrix_element('C0 (RIBBON.FOR)', row=0, col=0)
            c0_22 = self._extract_matrix_element('C0 (RIBBON.FOR)', row=1, col=1)
            c0_12 = self._extract_matrix_element('C0 (RIBBON.FOR)', row=0, col=1)
        except (TypeError, IndexError) as e:
            # Corrige as mensagens de erro para o contexto de indutância
            print(f"Erro ao extrair elementos da matriz de indutância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'L (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados com o estilo da figura de referência
        fortran_nf_axis = self.results_df.index + 1
        ax.plot(fortran_nf_axis, c0_11, color='k', marker='o', linestyle='-',  label='$C0_{11}$')
        ax.plot(fortran_nf_axis, c0_22, color='k', marker='o', linestyle=':',  label='$C0_{22}$', markerfacecolor='white', markeredgecolor='k')
        ax.plot(fortran_nf_axis, c0_12, color='k', marker='s', linestyle='-.', label='$C0_{12}$')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Bare-Wire Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_yticks(np.arange(10, 24, 2))
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_ylim(10, 24)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=3, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()


    def plot_paul_fig514d(self):
        """ Gera o gráfico de convergência a partir dos resultados armazenados, replicando a figura de referência. """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_convergence()'.")
            return

        # Garante que a simulação foi executada completamente
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Extrai os componentes da matriz de indutância 'L' e converte para µH/m.
        try:
            # Passa o fator de conversão correto (self.l_factor)
            cgen_00 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=0, col=0)
            cgen_11 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=1, col=1)
            cgen_22 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=2, col=2)
            cgen_01 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=0, col=1)
            cgen_02 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=0, col=2)
            cgen_12 = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=1, col=2)

        except (TypeError, IndexError) as e:
            # Corrige as mensagens de erro para o contexto de indutância
            print(f"Erro ao extrair elementos da matriz de indutância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'L (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados com o estilo da figura de referência
        fortran_nf_axis = self.results_df.index + 1
        ax.plot(fortran_nf_axis, cgen_00, color='k', marker='o', label='$CGEN_{00}$', linewidth=1.0, zorder=0, markersize=3, linestyle=':')
        ax.plot(fortran_nf_axis, cgen_11, color='k', marker='o', label='$CGEN_{11}$', linewidth=1.0, zorder=1, markersize=3, linestyle='-')
        ax.plot(fortran_nf_axis, cgen_22, color='k', marker='o', label='$CGEN_{22}$', linewidth=1.0, zorder=1, markersize=8, linestyle=':', fillstyle='none')
        ax.plot(fortran_nf_axis, cgen_01, color='k', marker='s', label='$CGEN_{01}$', linewidth=1.0, zorder=0, markersize=4, linestyle='--')
        ax.plot(fortran_nf_axis, cgen_12, color='k', marker='s', label='$CGEN_{12}$', linewidth=1.0, zorder=1, markersize=8, linestyle='--', fillstyle='none')
        ax.plot(fortran_nf_axis, cgen_02, color='k', marker='d', label='$CGEN_{02}$', linewidth=1.0, zorder=0, markersize=6, linestyle='-.')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Generalized Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_yticks(np.arange(0, 40, 5))
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_ylim(0, 40)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=6, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()


    def plot_bifilar_generalized_capacitance_matrix(self):
        """ Gera o gráfico de convergência a partir dos resultados armazenados, replicando a figura de referência. """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Extrai os componentes da matriz de capacitância 'C0' e converte para pF/m.
        # Por convenção, os elementos Cgen_ij (i != j) são negativos. 
        try:
            c11_ribbon = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=0, col=0)
            c22_ribbon = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=1, col=1)
            c12_ribbon = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=0, col=1)
            c11_pymom  = self._extract_matrix_element('CGEN (BARE-WIRE.PY)', row=0, col=0)
            c22_pymom  = self._extract_matrix_element('CGEN (BARE-WIRE.PY)', row=1, col=1)
            c12_pymom  = self._extract_matrix_element('CGEN (BARE-WIRE.PY)', row=0, col=1)
        except (TypeError, IndexError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'CGEN' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados do FORTRAN
        fortran_nf_axis = self.results_df.index + 1
        ax.plot(fortran_nf_axis, c11_ribbon, color='k', marker='o', label='$C_{11} (.FOR)$', linewidth=1.0, markersize=3, linestyle=':')
        ax.plot(fortran_nf_axis, c22_ribbon, color='k', marker='o', label='$C_{22} (.FOR)$', linewidth=1.0, markersize=8, linestyle=':', fillstyle='none')
        ax.plot(fortran_nf_axis, c12_ribbon, color='k', marker='^', label='$C_{12} (.FOR)$', linewidth=1.0, markersize=5, linestyle=':')

        # Plotagem dos dados do Python com o eixo x corrigido
        python_nf_axis = 2 * self.results_df.index + 1
        ax.plot(python_nf_axis, c11_pymom, color='gray', marker='s', label='$C_{11} (.PY)$', linewidth=1.0, markersize=3, linestyle='-.', fillstyle='none')
        ax.plot(python_nf_axis, c22_pymom, color='gray', marker='s', label='$C_{22} (.PY)$', linewidth=1.0, markersize=9, linestyle='-.', fillstyle='none')
        ax.plot(python_nf_axis, c12_pymom, color='gray', marker='^', label='$C_{12} (.PY)$', linewidth=1.0, markersize=9, linestyle='-.', fillstyle='none')

        ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=12)
        ax.set_ylabel('Generalized Capacitance (pF/m)', fontsize=12)
        ax.set_xlim(1, self.sum_max)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_ylim(20, 100)
        ax.set_yticks(np.arange(20, 100, 10))
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=6, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()


    def plot_generalized_capacitance_matrix(self):
        """
        Gera o gráfico de convergência da capacitância generalizada, com subplots
        separados para os termos CGEN_00 e CGEN_01.
        """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_convergence()'.")
            return

        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Extração de dados de forma programática.
        # Esta parte permanece flexível para extrair todos os elementos,
        # mesmo que apenas alguns sejam plotados.
        ribbon, mom, mom_so = {}, {}, {}
        try:
            for i in range(self.N):
                for j in range(self.N):
                    ribbon[f'c_{i}{j}'] =   self._extract_matrix_element('CGEN (RIBBON.FOR)', row=i, col=j)
                    mom_so[f'c_{i}{j}'] =   self._extract_matrix_element('CGEN (MOM-SO.PY)', row=i, col=j)
                    mom[f'c_{i}{j}'] =      self._extract_matrix_element('CGEN (BARE-WIRE.PY)', row=i, col=j)
        except (TypeError, IndexError, ValueError, KeyError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se as matrizes 'CGEN' foram populadas com as dimensões corretas.")
            return

        plt.style.use('default')
        # Cria uma figura com dois subplots (1 linha, 2 colunas)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        # fig.suptitle('')

        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1
        max_nf_fortran = fortran_nf_axis.max()
        mask_py = python_nf_axis <= max_nf_fortran
        
        # --- Subplot 1: CGEN_00 ---
        ax1.plot(fortran_nf_axis, ribbon['c_00'], label='RIBBON.FOR',
                color=self.plot_params['colors'][0], marker=self.plot_params['markers'][0], linestyle=self.plot_params['linestyles'][0])
        ax1.plot(python_nf_axis[mask_py], mom['c_00'][mask_py], label='MoM.PY',
                color=self.plot_params['colors'][1], marker=self.plot_params['markers'][1], linestyle=self.plot_params['linestyles'][1], fillstyle='none')
        ax1.plot(python_nf_axis[mask_py], mom_so['c_00'][mask_py], label='MoM-SO.PY',
                color=self.plot_params['colors'][2], marker=self.plot_params['markers'][2], linestyle=self.plot_params['linestyles'][2], fillstyle='none')

        # --- Subplot 2: CGEN_12 ---
        ax2.plot(fortran_nf_axis, ribbon['c_01'], label='RIBBON.FOR',
                color=self.plot_params['colors'][0], marker=self.plot_params['markers'][0], linestyle=self.plot_params['linestyles'][0])
        ax2.plot(python_nf_axis[mask_py], mom['c_01'][mask_py], label='MoM.PY',
                color=self.plot_params['colors'][1], marker=self.plot_params['markers'][1], linestyle=self.plot_params['linestyles'][1], fillstyle='none')
        ax2.plot(python_nf_axis[mask_py], mom_so['c_01'][mask_py], label='MoM-SO.PY',
                color=self.plot_params['colors'][2], marker=self.plot_params['markers'][2], linestyle=self.plot_params['linestyles'][2], fillstyle='none')

        # Configuração dos eixos para ambos os subplots
        for ax in [ax1, ax2]:
            ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=11)
            ax.set_xlim(0.8, max_nf_fortran + 0.2)
            ax.set_xticks(np.arange(1, max_nf_fortran + 1, 1))
            ax.tick_params(top=True, right=True, direction='in', which='both')
            ax.legend(loc='lower right', fontsize=10)
            ax.grid(False)

        # Configurações específicas por subplot
        ax1.set_ylabel('Generalized Capacitance Matrix, $CGEN$ (pF/m)', fontsize=11)
        ax1.set_title('Auto-Capacitance Term $CGEN_{00}$')
        ax2.set_title('Mutual Capacitance Term $CGEN_{01}$')
        plt.tight_layout(rect=[0, 0, 1, 0.96])


    def plot_free_space_capacitance_matrix(self):
        """ Gera o gráfico de convergência da capacitância do espaço livro, adaptando-se ao número de condutores do sistema. """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Extração de dados de forma programática
        ribbon, mom, mom_so = {}, {}, {}
        try:
            for i in range(self.N-1):
                for j in range(self.N-1):
                    ribbon[f'c_{i}{j}'] =   self._extract_matrix_element('C0 (RIBBON.FOR)', row=i, col=j)
                    mom_so[f'c_{i}{j}'] =   self._extract_matrix_element('C0 (MOM-SO.PY)', row=i, col=j)
                    mom[f'c_{i}{j}'] =      self._extract_matrix_element('C0 (BARE-WIRE.PY)', row=i, col=j)

        except (TypeError, IndexError, ValueError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))
        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1
        max_nf_fortran = fortran_nf_axis.max()
        mask_py = python_nf_axis <= max_nf_fortran
        
        legend_items_count = 0
        for i in range(self.N-1):
            legend_items_count += 3
            idx = i % len(self.plot_params['markers'])
            
            # Plot da diagonal principal
            ax.plot(fortran_nf_axis, ribbon[f'c_{i}{i}'], label='RIBBON.FOR',
                    color=self.plot_params['colors'][0], marker=self.plot_params['markers'][idx], linestyle=self.plot_params['linestyles'][0])

            ax.plot(python_nf_axis[mask_py], mom[f'c_{i}{i}'][mask_py], label='MoM.PY',
                    color=self.plot_params['colors'][1], marker=self.plot_params['markers'][idx], linestyle=self.plot_params['linestyles'][1], fillstyle='none')

            ax.plot(python_nf_axis[mask_py], mom_so[f'c_{i}{i}'][mask_py], label='MoM-SO.PY',
                    color=self.plot_params['colors'][2], marker=self.plot_params['markers'][idx], linestyle=self.plot_params['linestyles'][2], fillstyle='none')

            # Plot dos elementos fora da diagonal
            for j in range(i + 1, self.N-1):
                legend_items_count += 3
                ijdx = (i+j) % len(self.plot_params['markers'])
                ax.plot(fortran_nf_axis, ribbon[f'c_{i}{j}'], label=f'$C_{{{i+1}{j+1}}} (RIBBON.FOR)$',
                        color=self.plot_params['colors'][3], marker=self.plot_params['markers'][ijdx], linestyle=self.plot_params['linestyles'][0])

                ax.plot(python_nf_axis[mask_py], mom[f'c_{i}{j}'][mask_py], label=f'$C_{{{i+1}{j+1}}} (MoM.PY)$',
                        color=self.plot_params['colors'][4], marker=self.plot_params['markers'][ijdx], linestyle=self.plot_params['linestyles'][1], fillstyle='none')

                ax.plot(python_nf_axis[mask_py], mom_so[f'c_{i}{j}'][mask_py], label=f'$C_{{{i+1}{j+1}}} (MoM-SO.PY)$',
                        color=self.plot_params['colors'][5], marker=self.plot_params['markers'][ijdx], linestyle=self.plot_params['linestyles'][2], fillstyle='none')

        ax.set_title('Auto-Capacitance Term $C0_{11}$')
        ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=12)
        ax.set_ylabel('Free-Space Capacitance Matrix, $C0$ (pF/m)', fontsize=12)
        ax.set_xlim(0.8, max_nf_fortran + 0.2)
        ax.set_xticks(np.arange(1, max_nf_fortran + 1, 1))
        ax.tick_params(top=True, right=True, direction='in', which='both')
        # ax.legend(loc='upper center', fontsize=9, ncol=legend_items_count, bbox_to_anchor=(0.5, 1.1), fancybox=True)
        ax.legend(loc='lower right', fontsize=10)
        ax.grid(False)
        plt.tight_layout()