import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any

from mtl_main.utils import *
from mtl_main.graphics import IsolatedMTLRepresentation
from mtl_main.source import MulticonductorTransmissionLine
from mtl_paul.py_fortran import FortranRunner
from tulip.py_tulip import PyTulip
from analytical_forms.isolated_wires import WiresHomogeneousMedia
from mom.coated_wire_systems import MulticonductorCoatedWireSystems

class RibbonCoatedCablesPULParameters:
    """
    Encapsula a lógica para executar e analisar o estudo de convergência
    de capacitância, comparando MoM Python e Fortran.
    """
    def __init__(self, project_root: Path, case_name: str, mtl: Dict[str, Any], SUM_MAX: int = 10):
        """
        Inicializa o analisador de convergência.

        Args:
            mtl_config (Dict[str, Any]): Dicionário com a configuração do modelo MTL.
            nf_max (int): Número máximo de coeficientes/ordem harmônica para testar.
        """

        assert len([key for key in mtl.keys() if isinstance(key, int)]) > 1, "A linha deve conter mais de um condutor."

        self.project_root = project_root
        self.case_name = case_name
        self.mtl_copy = copy.deepcopy(mtl)
        self.freq_range = {'ana': np.logspace(0, 6, num=200), 'mom': np.logspace(0, 6, num=30)}
        self.srw_ratios = {'ana': np.linspace(4.0, 10.0, num=300), 'mom': np.linspace(4.0, 10.0, num=40)}

        # Parâmetros de dados
        self.srw_data = {}
        self.srw_mum_data = {}
        self.analytical_data = {}
        self.mom_data = {}
        self.mom_so_data = {}
        self.ribbon_data = {}
        self.tulip_runner = None

        # Parâmetros adicionais
        self.c_factor = 1e12  # F/m to nF/km
        self.l_factor = 1e6   # H/m to mH/km
        self.r_factor = 1e3   # Ohm/m to Ohm/km
        self.figsize = (12, 5)
        self.pt1 = 63
        self.pt2 = 10

        self.N = len([key for key in mtl.keys() if isinstance(key, int)])
        self.sum_max = SUM_MAX
        self.results_df = None
        
        # Parâmetros de plotagem
        self.plot_params = {
            'linestyles': [':', '-.', '--', '-', ':', '-.', '--'],
            'markers': ['o', 'd', 's', '^', 'v', '<', '>'],
            'colors': ['black', 'gray', 'lightgray', 'gainsboro']
        }

        self.R = self.mtl_copy[0]['radius'][1]
        self.D = np.linalg.norm(np.array(self.mtl_copy[0]['center_point']) - np.array(self.mtl_copy[1]['center_point']))
        self.swr_ratio = self.D / self.R

    def _prepare_fortran_runner(self, mtl):
        """Prepara os parâmetros e o executor para a simulação Fortran."""

        # Obtenha o dicionário 'insulation' de forma segura.
        #    Se 'insulation' não existir ou for None, use um dicionário vazio {} como fallback.
        refIdx = mtl['idx_ref_conductor']
        insulation_dict = mtl[refIdx].get('insulation') or {}

        self.fortran_base_params = {
            'N':    len([key for key in mtl.keys() if isinstance(key, int)]),
            'NF':   mtl[refIdx]['fourier_order'] + 1,
            'IREF': refIdx + 1,
            'RW':   mtl[refIdx]['radius'][1],
            'TD':   insulation_dict.get('thickness', 0.0),
            'ER':   insulation_dict.get('relative_permittivity', 1.0),
            'S':    np.linalg.norm(np.array(mtl[0]['center_point']) - np.array(mtl[1]['center_point'])),
        }

        fortran_exe_path = self.project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.runner = FortranRunner(exe_path=str(fortran_exe_path), silent=True)

    def _prepare_tulip_runner(self, case_name: str = "three_wires_ribbon"):
        """Prepara o executor para a simulação Tulip."""
        EXE_PATH = r"C:\git\tulip\pulmtln-build\rls\bin\Release\pulmtln.exe"
        print(f"--- Preparing to run simulation for case: '{case_name}' ---")
        
        try:
            self.tulip_runner = PyTulip(
                executable_path=EXE_PATH,
                case_name=case_name,
                silent=False,
                overwrite_mesh=False
            )

        except (FileNotFoundError, NotADirectoryError) as e:
            print(f"\nConfiguration Error: {e}")
        except Exception as e:
            print(f"\nAn unexpected error occurred: {e}")
        
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
        print("="*self.pt2 + " BIFILAR COATED-WIRE RIBBON CABLE SIMULATION " + "="*self.pt2)
        print(f"Project: {self.project_root}")
        print(f"Model: {self.mtl_copy['name']}")
        print(f"D/R = {self.swr_ratio:.3f}. Fourier Order (k) = {self.mtl_copy[0]['fourier_order']}.")
        print(f"RIBBON Fourier Coef./cond. (NF) = {self.mtl_copy[0]['fourier_order']+1}.")
        print(f"PYTHON Fourier Coef./cond. (NF) = {2*self.mtl_copy[0]['fourier_order']+1}.")
        print("="*self.pt1)

    def run_single_fortran(self, mtl: dict = None, displayTerminal: bool = True):
        """Executa uma simulação única para um valor específico de k."""

        if mtl is None:
            self._prepare_fortran_runner(self.mtl_copy)
        else:
            self._prepare_fortran_runner(mtl)

        self._prepare_fortran_runner(self.mtl_copy)
        self.runner.run_fortran(self.fortran_base_params)

        print("\n")
        print("="*self.pt2 + "                 RIBBON.FOR                " + "="*self.pt2)
        if self.runner.CAP_matrix is not None:
            if self.runner.A_matrix.shape[0] < 6:
                matrix_viewer(self.runner.A_matrix,     "Block Matrix A")
                matrix_viewer(self.runner.B_matrix,     "Block Matrix B")
                matrix_viewer(self.runner.C_matrix,     "Block Matrix C")
                matrix_viewer(self.runner.D_matrix,     "Block Matrix D")
            else:
                print(f"\nBlock matrices Shape: {self.runner.A_matrix.shape}.")

            matrix_viewer(self.runner.IND_matrix,   "External Inductance Matrix, Le (H/m)")
            matrix_viewer(self.runner.CAP_matrix,   "Capacitance Matrix, C (F/m)")
            matrix_viewer(self.runner.CAP0_matrix,  "Capacitance Matrix in Vacuum, C0 (F/m)")
            matrix_viewer(self.runner.CGEN0_matrix, "Free-Space Generalized Capacitance Matrix, CGEN0 (F/m)")
            matrix_viewer(self.runner.CGEN_matrix,  "Generalized Capacitance Matrix, CGEN (F/m)")
        else:
            print("\nNo results found.")

    def run_fortran(self, mtl: dict = None, displayTerminal: bool = True):
        """Executa uma simulação única para um valor específico de k."""

        if mtl is None:
            self._prepare_fortran_runner(self.mtl_copy)
        else:
            self._prepare_fortran_runner(mtl)

        self._prepare_fortran_runner(self.mtl_copy)
        self.runner.run_fortran(self.fortran_base_params)

        self.ribbon_data = {
            freq: {'c': self.c_factor * self.runner.CAP0_matrix[0, 0], 
                   'le': self.l_factor * self.runner.IND_matrix[0, 0]} for freq in self.freq_range['mom']}

    def run_single_tulip(self):
        """
        Executa a simulação usando o PyTulip e processa os resultados.
        """
        print("\n============  pyTulip Simulation  ===========")
        if self.tulip_runner is None:
            self._prepare_tulip_runner()

        self.tulip_runner.run()
        
        if not self.tulip_runner.output_filepath.is_file():
            print(f"Error: Output file '{self.tulip_runner.output_filepath.name}' was not generated.")
            return
        
        print("\n--- Accessing Results ---")
        if self.tulip_runner.C_matrix is not None:
            matrix_viewer(self.tulip_runner.C_matrix, "Capacitance Matrix (C)")

        if self.tulip_runner.L_matrix is not None:
            matrix_viewer(self.tulip_runner.L_matrix, "Inductance Matrix (L)")

    def run_mom_methods(self):
        """
        Executa a simulação clássica do Método dos Momentos (MoM) para a linha de transmissão bifilar.

        Returns:
            BifilarMoM: Instância do objeto BifilarMoM configurado.
        """
        print("\n============  pyMoM TwoCoatedWireSystem (Bare-Wire)  ===========")
        bare_wire_mtl = copy.deepcopy(self.mtl_copy)
        bare_wire_mtl['type'] = 'bare_wires'
        for key in bare_wire_mtl.keys():
            if isinstance(key, int):
                bare_wire_mtl[key]['insulation'] = None

        bare_wires_model = MulticonductorTransmissionLine(bare_wire_mtl)
        bare_wires = MulticonductorCoatedWireSystems(bare_wires_model)
        bare_wires.run_simulation()
        bare_wires.print_results()

        print("\n==============      pyMoM TwoCoatedWireSystem      =============")
        coated_wires_model = MulticonductorTransmissionLine(self.mtl_copy)
        IsolatedMTLRepresentation(coated_wires_model, self.case_name, units='millimeter').system_schematic()
        coated_wires = MulticonductorCoatedWireSystems(coated_wires_model)
        coated_wires.run_simulation()
        coated_wires.print_results()
        coated_wires.plot_collocation_points()

        self.mom_data = {
            freq: {'c_bare_wire': self.c_factor * bare_wires.C_maxwellian,
                    'c_coated_wire': self.c_factor * coated_wires.C_maxwellian} for freq in self.freq_range['mom']}

    def srw_rates(self):
        """
        Executa a simulação analítica da impedância da linha de transmissão.

        Args:
            mtl_config (dict): Dicionário de configuração da linha de transmissão.
            frequencies (np.ndarray): Array de frequências para a análise.

        Returns:
            tuple: Uma tupla contendo três listas: impedâncias série,
                resistências de alta frequência e indutâncias externas.
        """
        print("\n==============         SRW RATES EVALUATION        =============")

        # Cria a configuração da linha bifilar dinamicamente para cada razão.
        for ratio in self.srw_ratios['ana']:
            temp_mtl = copy.deepcopy(self.mtl_copy)
            separation = ratio * temp_mtl[0]['radius'][1]
            temp_mtl[1]['center_point'] = (separation, 0.0)

            wires = WiresHomogeneousMedia(temp_mtl)
            le_wires = wires.n_wires_inductance_matrix()
            pul_bifilar = wires.bifilar_pul_inductance_and_capacitance()

            self.srw_data[ratio] = {
                'le_wires':     self.l_factor * le_wires,
                'le_exact':     self.l_factor * pul_bifilar['inductance']['exact'],
                'le_bifilar':   self.l_factor * pul_bifilar['inductance']['approximate'],
                'c_exact':      self.c_factor * pul_bifilar['capacitance']['exact'],
                'c_approx':     self.c_factor * pul_bifilar['capacitance']['approximate'],
                'c_wires':      self.c_factor * wires.n_wires_capacitance_matrix(le_wires),
            }

        for ratio in self.srw_ratios['mom']:
            temp_mtl = copy.deepcopy(self.mtl_copy)
            separation = ratio * temp_mtl[0]['radius'][1]
            temp_mtl[1]['center_point'] = (separation, 0.0)

            # === Fortran RIBBON Instance ===            
            self._prepare_fortran_runner(temp_mtl)
            self.runner.run_fortran(self.fortran_base_params)

            # === MoM TwoCoatedWireSystem Instance ===
            mom_coated_model = MulticonductorTransmissionLine(temp_mtl)
            mom_coated = MulticonductorCoatedWireSystems(mom_coated_model)
            mom_coated.run_simulation()
            
            temp_mtl['type'] = 'bare_wires'
            for key in temp_mtl.keys():
                if isinstance(key, int):
                    temp_mtl[key]['insulation'] = None

            mom_bare_model = MulticonductorTransmissionLine(temp_mtl)
            mom_bare = MulticonductorCoatedWireSystems(mom_bare_model)
            mom_bare.run_simulation()

            self.srw_mum_data[ratio] = {
                'ribbon-c':     self.c_factor * self.runner.CAP_matrix.item(),
                'ribbon-c0':    self.c_factor * self.runner.CAP0_matrix.item(),
                'mom-c':        self.c_factor * mom_coated.C_maxwellian.item(),
                'mom-c0':       self.c_factor * mom_bare.C_maxwellian.item(),
            }

    def run_convergence(self):
        """Executa o laço de convergência para ambas as simulações e armazena os resultados."""
        print(f"\nRunning Convergence Rate until k = {self.sum_max}!")

        if self.tulip_runner is None:
            self._prepare_tulip_runner()
            self.tulip_runner.run()
        
        results = []
        for k in range(0, self.sum_max):
            temp_mtl = copy.deepcopy(self.mtl_copy)
            for key in temp_mtl.keys():
                if isinstance(key, int):
                    temp_mtl[key]['fourier_order'] = k  
                    temp_mtl[key]['insulation']['fourier_order'] = k

            # === Fortran RIBBON Instance ===
            self._prepare_fortran_runner(temp_mtl)
            self.runner.run_fortran(self.fortran_base_params)

            # === MoM TwoCoatedWireSystem Instance ===
            mom_coated_model = MulticonductorTransmissionLine(temp_mtl)
            mom_coated = MulticonductorCoatedWireSystems(mom_coated_model)
            mom_coated.run_simulation()

            results.append({
                'k': k,
                'L (RIBBON.FOR)':    self.runner.IND_matrix if self.fortran_base_params is not None else np.nan,
                'C (RIBBON.FOR)':    self.runner.CAP_matrix if self.fortran_base_params is not None else np.nan,
                'C0 (RIBBON.FOR)':   self.runner.CAP0_matrix if self.fortran_base_params is not None else np.nan,
                'CGEN (RIBBON.FOR)': self.runner.CGEN0_matrix if self.fortran_base_params is not None else np.nan,
                'C (MoM.PY)':        mom_coated.C_maxwellian if mom_coated.C_maxwellian is not None else np.nan,
                'C (SEMBA-TULIP)':   self.tulip_runner.C_matrix if self.tulip_runner.C_matrix is not None else np.nan,
                'L (SEMBA-TULIP)':   self.tulip_runner.L_matrix if self.tulip_runner.L_matrix is not None else np.nan,
            })
            print(f"  Complete for k = {k}.")

        self.results_df = pd.DataFrame(results).set_index('k')

    def plot_inductance_data(self):
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
            'wires':        {'data': (self.freq_range.get('ana'), [data['le_wires']     for data in self.srw_data.values()]), 'label': r'$\ell_{e,n+1 \; wires}$'},
            'exactly':      {'data': (self.freq_range.get('ana'), [data['le_exact']     for data in self.srw_data.values()]), 'label': r'$\ell_{e,Bifilar (Exactly)}$'},
            'bifilar':      {'data': (self.freq_range.get('ana'), [data['le_bifilar']   for data in self.srw_data.values()]), 'label': r'$\ell_{e,Bifilar (Approx.)}$'},
            'analytical':   {'data': (self.freq_range.get('ana'), [data['ls']           for data in self.srw_data.values()]), 'label': '$\ell_s$'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax, 'Series Inductance p.u.l. (mH/km)', inductance_data, yscale='linear')
        plt.tight_layout()

    def plot_capacitance_data(self):
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
            'wires':    {'data': (self.freq_range.get('ana'), [data['c_wires']  for data in self.srw_data.values()]), 'label': 'n+1 wires'},
            'exactly':  {'data': (self.freq_range.get('ana'), [data['c_exact']  for data in self.srw_data.values()]), 'label': 'Bifilar (Exactly)'},
            'bifilar':  {'data': (self.freq_range.get('ana'), [data['c_approx'] for data in self.srw_data.values()]), 'label': 'Bifilar (Approx.)'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        self._configure_plot_appearance(ax, 'Capacitance p.u.l. (nF/km)', capacitante_data, yscale='linear')
        plt.tight_layout()

    def plot_srw_rates(self):
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
            'ribbon-c':     {'data': (self.srw_ratios.get('mom'), [data['ribbon-c']     for data in self.srw_mum_data.values()]), 'label': '$C$ (RIBBON.FOR)'},
            'mom-c':        {'data': (self.srw_ratios.get('mom'), [data['mom-c']        for data in self.srw_mum_data.values()]), 'label': '$C$ (MoM.PY)'},
            'ribbon-c0':    {'data': (self.srw_ratios.get('mom'), [data['ribbon-c0']    for data in self.srw_mum_data.values()]), 'label': '$C_0$ (RIBBON.FOR)'},
            'mom-c0':       {'data': (self.srw_ratios.get('mom'), [data['mom-c0']       for data in self.srw_mum_data.values()]), 'label': '$C_0$ (MoM.PY)'},
            'exactly':      {'data': (self.srw_ratios.get('ana'), [data['c_exact']      for data in self.srw_data.values()]),     'label': 'Exactly'},
            'approx':       {'data': (self.srw_ratios.get('ana'), [data['c_approx']     for data in self.srw_data.values()]),     'label': 'Approx.'},
        }

        fig, ax = plt.subplots(figsize=self.figsize)
        for key, data in capacitante_data.items():
            freq, value = data['data']
            label = data['label']
            if freq is not None and value is not None:
                if key == 'ribbon-c':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='s', markersize=5)
                elif key == 'ribbon-c0':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=4)
                elif key == 'mom-c':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='s', markersize=9, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'mom-c0':
                    ax.plot(freq, value, label=label, color='k', linestyle='none', marker='o', markersize=9, fillstyle='none', markeredgecolor='k', zorder=2)
                elif key == 'approx':
                    ax.plot(freq, value, label=label, color='k', linestyle='--', linewidth=1.0, zorder=1)
                elif key == 'exactly':
                    ax.plot(freq, value, label=label, color='k', linestyle=':', linewidth=1.0, zorder=1)

        ax.set_xlabel('Ratio of separation to wire radius, s/r$_w$')
        ax.set_ylabel('Per-unit-length capacitance (pF/m)')
        ax.set_xlim(4, 8)
        ax.legend()
        ax.grid(True, linestyle='--', linewidth=0.5)
        plt.tight_layout()
        
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
            l11_ribbon = self._extract_matrix_element('L (RIBBON.FOR)', row=0, col=0)
            l22_ribbon = self._extract_matrix_element('L (RIBBON.FOR)', row=1, col=1)
            l12_ribbon = self._extract_matrix_element('L (RIBBON.FOR)', row=0, col=1)
            l11_tulip  = self._extract_matrix_element('L (SEMBA-TULIP)', row=0, col=0)
            l22_tulip  = self._extract_matrix_element('L (SEMBA-TULIP)', row=1, col=1)
            l12_tulip  = self._extract_matrix_element('L (SEMBA-TULIP)', row=0, col=1)
        except (TypeError, IndexError) as e:
            # Corrige as mensagens de erro para o contexto de indutância
            print(f"Erro ao extrair elementos da matriz de indutância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'L (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)

        for_nf_axis = self.results_df.index + 1
        ax.plot(for_nf_axis, l22_ribbon, color='k', marker='o', linestyle=':', label='$L_{22}$ (RIBBON.FOR)')
        ax.plot(for_nf_axis, l11_ribbon, color='k', marker='o', linestyle=':', label='$L_{11}$ (RIBBON.FOR)', markerfacecolor='white', markeredgecolor='k')
        ax.plot(for_nf_axis, l12_ribbon, color='k', marker='s', linestyle=':', label='$L_{12}$ (RIBBON.FOR)')
        ax.plot(for_nf_axis, l22_tulip, color='gray', linestyle='-.', label='(SEMBA-TULIP)')
        ax.plot(for_nf_axis, l11_tulip, color='gray', linestyle='-.')
        ax.plot(for_nf_axis, l12_tulip, color='gray', linestyle='-.')

        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Inductance (µH/m)', fontsize=12)
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))        
        ax.set_yticks(np.arange(0.2, 1.2, 0.1))        
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=4, bbox_to_anchor=(0.5, 1.12), fancybox=True)
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
            c11_ribbon = self._extract_matrix_element('C (RIBBON.FOR)', row=0, col=0)
            c22_ribbon = self._extract_matrix_element('C (RIBBON.FOR)', row=1, col=1)
            c12_ribbon = self._extract_matrix_element('C (RIBBON.FOR)', row=0, col=1)
            c11_tulip  = self._extract_matrix_element('C (SEMBA-TULIP)', row=0, col=0)
            c22_tulip  = self._extract_matrix_element('C (SEMBA-TULIP)', row=1, col=1)
            c12_tulip  = self._extract_matrix_element('C (SEMBA-TULIP)', row=0, col=1)
        except (TypeError, IndexError) as e:
            # Corrige as mensagens de erro para o contexto de indutância
            print(f"Erro ao extrair elementos da matriz de indutância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'L (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)

        for_nf_axis = self.results_df.index + 1
        ax.plot(for_nf_axis, c11_ribbon, color='k', marker='o', linestyle=':', label='$C_{11}$ (RIBBON.FOR)')
        ax.plot(for_nf_axis, c22_ribbon, color='k', marker='o', linestyle=':', label='$C_{22}$ (RIBBON.FOR)', markerfacecolor='white', markeredgecolor='k')
        ax.plot(for_nf_axis, c12_ribbon, color='k', marker='s', linestyle=':', label='$C_{12}$ (RIBBON.FOR)')
        ax.plot(for_nf_axis, c11_tulip, color='gray', linestyle='-.',  label='(SEMBA-TULIP)')
        ax.plot(for_nf_axis, c22_tulip, color='gray', linestyle='-.')
        ax.plot(for_nf_axis, c12_tulip, color='gray', linestyle='-.')

        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_yticks(np.arange(15, 40, 5))
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=4, bbox_to_anchor=(0.5, 1.12), fancybox=True)
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
        fig, ax = plt.subplots(figsize=self.figsize)

        # Plotagem dos dados com o estilo da figura de referência
        fortran_nf_axis = self.results_df.index + 1
        ax.plot(fortran_nf_axis, c0_11, color='k', marker='o', linestyle='-',  label='$C0_{11}$')
        ax.plot(fortran_nf_axis, c0_22, color='k', marker='o', linestyle=':',  label='$C0_{22}$', markerfacecolor='white', markeredgecolor='k')
        ax.plot(fortran_nf_axis, c0_12, color='k', marker='s', linestyle='-.', label='$C0_{12}$')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Bare-Wire Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.sum_max + 1, 1))
        ax.set_xlim(0.8, self.sum_max + 0.2)
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
        fig, ax = plt.subplots(figsize=self.figsize)

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
        ax.set_xlim(0.8, self.sum_max + 0.2)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper center', fontsize=9, markerscale=1.0, ncol=6, bbox_to_anchor=(0.5, 1.12), fancybox=True)
        plt.tight_layout()

    def plot_generalized_capacitance_convergence(self):
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
        except (TypeError, IndexError, ValueError, KeyError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se as matrizes 'CGEN' foram populadas com as dimensões corretas.")
            return

        plt.style.use('default')
        # Cria uma figura com dois subplots (1 linha, 2 colunas)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=True)
        # fig.suptitle('')

        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1
        max_nf_fortran = fortran_nf_axis.max()
        mask_py = python_nf_axis <= max_nf_fortran
        
        # --- Subplot 1: CGEN_00 ---
        ax1.plot(fortran_nf_axis, ribbon['c_00'], label='RIBBON.FOR',
                color=self.plot_params['colors'][0], marker=self.plot_params['markers'][0], linestyle=self.plot_params['linestyles'][0])

        # --- Subplot 2: CGEN_12 ---
        ax2.plot(fortran_nf_axis, ribbon['c_01'], label='RIBBON.FOR',
                color=self.plot_params['colors'][0], marker=self.plot_params['markers'][0], linestyle=self.plot_params['linestyles'][0])

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

    def plot_dielectric_coated_capacitance_convergence(self):
        """ Gera o gráfico de convergência da capacitância do espaço livro, adaptando-se ao número de condutores do sistema. """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Extração de dados de forma programática
        ribbon, mom =  {}, {}

        try:
            for i in range(self.N-1):
                for j in range(self.N-1):
                    ribbon[f'c_{i}{j}'] = self._extract_matrix_element('C (RIBBON.FOR)', row=i, col=j)
                    mom[f'c_{i}{j}'] = self._extract_matrix_element('C (MoM.PY)', row=i, col=j)

        except (TypeError, IndexError, ValueError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)
        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1
        max_nf_fortran = fortran_nf_axis.max()
        mask_py = python_nf_axis <= max_nf_fortran
        
        for i in range(self.N-1):
            for j in range(i, self.N-1):
                ax.plot(fortran_nf_axis, ribbon[f'c_{i}{j}'], markersize=4, zorder = 1,
                        label=f'$C_{{{i+1}{j+1}}}$ (RIBBON.FOR)', 
                        color=self.plot_params['colors'][i+j],
                        marker=self.plot_params['markers'][i+j],
                        linestyle=self.plot_params['linestyles'][i+j])

                ax.plot(python_nf_axis[mask_py], mom[f'c_{i}{j}'][mask_py], markersize=10, fillstyle='none', zorder = 2,
                        label=f'$C_{{{i+1}{j+1}}}$ (MoM.PY)', 
                        color=self.plot_params['colors'][i+j],
                        marker=self.plot_params['markers'][i+j],
                        linestyle=self.plot_params['linestyles'][i+j])
            
        ax.set_title('Three Ribbon Cable Capacitance Convergence')
        ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=12)
        ax.set_ylabel('Capacitance (pF/m)', fontsize=12)
        ax.set_xlim(0.8, max_nf_fortran + 0.2)
        ax.set_xticks(np.arange(1, max_nf_fortran + 1, 1))
        ax.tick_params(top=True, right=True, direction='in', which='both')
        ax.legend(loc='best', fontsize=10)
        ax.grid(False)
        plt.tight_layout()
    