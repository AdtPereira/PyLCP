import copy
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any

from mtl_paul.py_fortran import FortranRunner
from mom.bare_wire_systems import MulticonductorBareWireSystems
from analytical_forms.bare_wires import WiresHomogeneousMedia
from mom_so.quasi_static_green import QuasiStatic
from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing

class CompareBifilarPULParameters():
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
            'NF':   mtl[refIdx]['fourier_order'],
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
        print(f"D/R = {self.DR_ratio} --- NF = {self.mtl_copy[0]['fourier_order']} per conductor.")
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

        self.mom_data = {freq: {'c': self.c_factor * mom_wires.C_maxwellian} for freq in self.freq_range['mom']}

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
            self.ribbon_data[ratio] = {'c': self.c_factor * self.runner.C0_matrix[0,0]} 
            self.mom_data[ratio]    = {'c': self.c_factor * mom_wires.C_maxwellian}


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
