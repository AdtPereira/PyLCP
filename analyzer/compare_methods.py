import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any

from mtl_data.mtl import MulticonductorTransmissionLine as MTL
from mtl_paul.py_fortran import FortranRunner
from mom.bare_wire_systems import MulticonductorBareWireSystems
from analytical_forms.bare_wires import WiresHomogeneousMedia
from mom_so.quasi_static_green import QuasiStatic
from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing

class PULparameters(MTL):
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
        super().__init__(mtl)
        self.project_root = project_root
        self.mtl_data = copy.deepcopy(mtl)
        self.freq_range = {'ana': np.logspace(0, 6, num=200), 'mom': np.logspace(0, 6, num=30), 'mom-so': np.logspace(0, 6, num=30)}

        # Extrai parâmetros e prepara o executor do Fortran
        self._prepare_fortran_runner()
        self._calculate_analytical_solution()

        # Parâmetros adicionais
        self.c_factor = 1e12  # F/m para nF/km
        self.dist1 = 63
        self.dist2 = 10

    def _prepare_fortran_runner(self):
        """Prepara os parâmetros e o executor para a simulação Fortran."""
        center_0 = np.array(self.mtl[0]['center_point'])
        center_1 = np.array(self.mtl[1]['center_point'])
        
        # Obtenha o dicionário 'sheath' de forma segura. 
        #    Se 'sheath' não existir ou for None, use um dicionário vazio {} como fallback.
        ref_conductor = self.mtl[self.idx_ref]
        sheath_dict = ref_conductor.get('sheath') or {}

        # 3. Construa seu dicionário de parâmetros usando as variáveis seguras
        self.fortran_base_params = {
            'N':    len(self.mtl),
            'NF':   self.NF,
            'IREF': self.idx_ref,
            'RW':   ref_conductor['radius'][1],
            'TD':   sheath_dict.get('thickness', 0.0),
            'ER':   sheath_dict.get('relative_permittivity', 1.0),
            'S':    np.linalg.norm(center_0 - center_1)
        }
        
        fortran_exe_path = self.project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.runner = FortranRunner(exe_path=str(fortran_exe_path), silent=True)
        self.D_over_R = self.fortran_base_params['S'] / self.fortran_base_params['RW']

    def _calculate_analytical_solution(self):
        """Calcula a solução analítica para fios nus como referência."""
        R = self.fortran_base_params['RW']
        D = self.fortran_base_params['S']
        from scipy.constants import epsilon_0
        self.analytical_bare_capacitance = (np.pi * epsilon_0) / np.arccosh(D / (2 * R))

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
            frequencies, values = data['data']
            label = data['label']

            # Verifica se as frequências e valores estão disponíveis
            if frequencies is not None and values is not None:
                if key == 'MoM-SO':
                    ax.scatter(frequencies, values, label=label, facecolors='k', marker='o', s=11)
                elif key == 'MoM':
                    ax.plot(frequencies, values, label=label, linestyle='none', marker='d', markersize=8, fillstyle='none', markeredgecolor='black')
                elif key == 'ribbon':
                    ax.plot(frequencies, values, label=label, linestyle='none', marker='o', markersize=13, fillstyle='none', markeredgecolor='blue')
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

    def show_header(self):
        """Exibe o cabeçalho do script."""
        pt1, pt2 = self.dist1, self.dist2
        print("\n")
        print("="*pt2 + " BIFILAR BARE-WIRE RIBBON CABLE SIMULATION " + "="*pt2)
        print(f"Projeto: {self.project_root.name}")
        print(f"Modelo: {self.mtl_data['type']}")
        print("="*pt1)

    def run_fortran(self, displayTerminal: bool = True):
        """Executa uma simulação única para um valor específico de k."""
        
        self.runner.run_fortran(self.fortran_base_params)
        
        if displayTerminal:
            pt1, pt2 = self.dist1, self.dist2
            print("\n")
            print("="*pt2 + "                 RIBBON.FOR                " + "="*pt2)
            if self.runner.L_matrix is not None:
                print("\nMatriz de Indutância (L):")
                print(self.runner.L_matrix)

                print("\nMatriz de Capacitância (C):")
                print(self.runner.C_matrix)

                print("\nMatriz de Capacitância no Vácuo (C0):")
                print(self.runner.C0_matrix)

                print("\nMatriz de Capacitância Generalizada (CGEN):")
                print(self.runner.CGEN_matrix)
            else:
                print("\nNenhum resultado foi analisado. Verifique os logs de erro.")

    def run_mom(self, autoPlots=False):
        """
        Executa a simulação clássica do Método dos Momentos (MoM) para a linha de transmissão bifilar.

        Returns:
            BifilarMoM: Instância do objeto BifilarMoM configurado.
        """
        print("\n============== pyMoM MulticonductorBareWireSystems =============")
        
        self.mom_data = MulticonductorBareWireSystems(self.mtl_data)
        self.mom_data.run_simulation()
        self.mom_data.print_results()

        if autoPlots:
            self.mom_data.plot_charge_density()
            self.mom_data.plot_collocation_points()
            self.mom_data.plot_harmonic_coefficients()
            MulticonductorBareWireSystems.plot_convergence_rates(self.mtl_data, nf_max=20)
    
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
        
        self.ana_data = {}
        bifilar_wires = WiresHomogeneousMedia(self.mtl_data)
        le_wires = bifilar_wires.n_wires_inductance_matrix()
        pul_bifilar = bifilar_wires.bifilar_pul_inductance_and_capacitance()
        
        for freq in self.freq_range['ana']:
            z_s, r_hf = bifilar_wires.bifilar_pul_series_impedance(freq)
            self.ana_data[freq] = {
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

        print(f"\nExact Bifilar Bare Wire Capacitance: {pul_bifilar['capacitance']['exact'] * 1E12:.4f} pF/m")

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

        green_matrix = QuasiStatic(self.mtl_data).green_matrix()
        post_processor = LosslessPostProcessing(self.mtl_data)
        self.mom_so_data = {}

        for freq in self.freq_range['mom-so']:
            mom_so = HomogeneousLosslessMedium(self.mtl_data, freq)
            zs = post_processor.z_total(mom_so.z_partial(green_matrix))
            general_cap = mom_so.generalized_capacitance_matrix(green_matrix)
            maxwell_cap = mom_so.maxwellian_capacitance_matrix(general_cap)
            self.mom_so_data[freq] = {
                'zs': zs,
                'rs': post_processor.rs_matrix(zs),
                'ls': post_processor.ls_matrix(zs, freq),
                'c': maxwell_cap,
            }

        if green_matrix.shape[0] < 6:
            from scipy.constants import epsilon_0
            print(f"\nGreen's Matrix (Dim: {green_matrix.shape}):\n{green_matrix}")
            print(f"\n2*pi*e0*G:\n{- 2 * np.pi * epsilon_0 * np.real(green_matrix)}")
        
        print(f"\nMoM-SO Generalized Capacitance Matrix (F/m):\n{np.real(general_cap)}")
        print(f"\nMoM-SO Bifilar Capacitance: {np.real(maxwell_cap.item()) * 1E12:.4f} pF/m")
        
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
        R_FACTOR = 1e3   # de Ohm/m para Ohm/km
        rs, rhf, rs_mom = [], [], []

        # 2. Processe os dados analíticos em um único laço
        for data in self.ana_data.values():
            rs.append(data['rs'] * R_FACTOR)
            rhf.append(data['rhf'] * R_FACTOR)

        # 3. Processe os dados do MoM-SO em um laço separado
        for data in self.mom_so_data.values():
            rs_mom.append(data['rs'][0, 0] * R_FACTOR)

        resistance_data = {
            'MoM-SO':       {'data': (self.freq_range.get('mom-so'), rs_mom),   'label': 'MoM-SO'},
            'Analytical':   {'data': (self.freq_range.get('ana'), rs),          'label':'$R_i$'},
            'Approximate':  {'data': (self.freq_range.get('ana'), rhf),         'label': '$R_{HF}$'},
        }

        # Cria a figura com subplots
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
        L_FACTOR = 1e6   # de H/m para mH/km
        ls, le_exact, le_approx, le_wires = [], [], [], []
        ls_mom_so, le_ribbon = [], []

        # 2. Processe os dados analíticos em um único laço
        for data in self.ana_data.values():
            ls.append(data['ls'] * L_FACTOR)
            le_exact.append(data['le_exact'] * L_FACTOR)    
            le_wires.append(data['le_wires'][0, 0] * L_FACTOR)
            le_approx.append(data['le_bifilar'] * L_FACTOR)

        # 3. Processe os dados do MoM-SO em um laço separado
        for data in self.mom_so_data.values():
            ls_mom_so.append(data['ls'][0, 0] * L_FACTOR)
            le_ribbon.append(self.runner.L_matrix[0, 0] * L_FACTOR)

        inductance_data = {
            'MoM-SO':       {'data': (self.freq_range.get('mom-so'), ls_mom_so),    'label': 'MoM-SO'},
            'Analytical':   {'data': (self.freq_range.get('ana'), ls),              'label': '$\ell_s$'},
            'ribbon':       {'data': (self.freq_range.get('mom'), le_ribbon),       'label': r'$\ell_{e,RIBBON.FOR}$'},
            'Wires':        {'data': (self.freq_range.get('ana'), le_wires),        'label': r'$\ell_{e,n+1 \; wires}$'},
            'Exactly':      {'data': (self.freq_range.get('ana'), le_exact),        'label': r'$\ell_{e,bifilar (Exactly)}$'},
            'Bifilar':      {'data': (self.freq_range.get('ana'), le_approx),       'label': r'$\ell_{e,bifilar (Approx.)}$'},
        }

        # Cria a figura com subplots
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
        cap_approx, cap_exact, cap_wires = [], [], []
        cap_mom_so, cap_mom, cap_ribbon = [], [], []

        # 2. Processe os dados analíticos em um único laço
        for data in self.ana_data.values():
            cap_approx.append(data['cap_approx'] * self.c_factor)
            cap_exact.append(data['cap_exact'] * self.c_factor)
            cap_wires.append(data['cap_wires'][0, 0] * self.c_factor)

        # 3. Processe os dados do MoM-SO em um laço separado
        for data in self.mom_so_data.values():
            cap_ribbon.append(self.runner.C0_matrix[0, 0] * self.c_factor)
            cap_mom_so.append(np.real(data['c'][0, 0]) * self.c_factor)
            cap_mom.append(self.mom_data.C_maxwellian * self.c_factor)

        capacitante_data = {
            'MoM-SO':   {'data': (self.freq_range.get('mom-so'), cap_mom_so),   'label': 'MoM-SO'},
            'ribbon':   {'data': (self.freq_range.get('mom'), cap_ribbon),      'label': 'RIBBON.FOR'},
            'MoM':      {'data': (self.freq_range.get('mom'), cap_mom),         'label': 'MoM'},
            'Wires':    {'data': (self.freq_range.get('ana'), cap_wires),       'label': r'$c_{n+1 \; wires}$'},
            'Exactly':  {'data': (self.freq_range.get('ana'), cap_exact),       'label': r'$c_{bifilar (Exactly)}$'},
            'Bifilar':  {'data': (self.freq_range.get('ana'), cap_approx),      'label': r'$c_{bifilar (Approx.)}$'},
        }

        # Cria a figura com subplots
        fig, ax = plt.subplots(figsize=(8, 5))
        self._configure_plot_appearance(ax, 'Capacitance p.u.l. (nF/km)', capacitante_data, yscale='linear')
        plt.tight_layout()

