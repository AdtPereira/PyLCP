import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any

from mtl_paul.py_fortran import FortranRunner
from mom.bare_wire_systems import MulticonductorBareWireSystems

class ConvergenceAnalyzer():
    """
    Encapsula a lógica para executar e analisar o estudo de convergência
    de capacitância, comparando MoM Python e Fortran.
    """
    def __init__(self, project_root: Path, mtl: Dict[str, Any], NF_MAX: int = 10):
        """
        Inicializa o analisador de convergência.

        Args:
            mtl_config (Dict[str, Any]): Dicionário com a configuração do modelo MTL.
            nf_max (int): Número máximo de coeficientes/ordem harmônica para testar.
        """
        self.project_root = project_root
        self.mtl_copy = copy.deepcopy(mtl)
        self.nf_max = NF_MAX
        self.results_df = None
        
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
            'N':    len([key for key in mtl.keys() if isinstance(key, int)]),
            'NF':   mtl[refIdx]['fourier_order'],
            'IREF': refIdx,
            'RW':   mtl[refIdx]['radius'][1],
            'TD':   sheath_dict.get('thickness', 0.0),
            'ER':   sheath_dict.get('relative_permittivity', 1.0),
            'S':    np.linalg.norm(np.array(mtl[0]['center_point']) - np.array(mtl[1]['center_point'])),
        }
        
        fortran_exe_path = self.project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.runner = FortranRunner(exe_path=str(fortran_exe_path), silent=False)
        self.runner_silent = FortranRunner(exe_path=str(fortran_exe_path), silent=True)


    def run_convergence(self):
        """Executa o laço de convergência para ambas as simulações e armazena os resultados."""
        print(f"\nRunning Convergence Rate until k = {self.nf_max}!")
        
        results = []
        for NF in range(1, self.nf_max + 1):
            temp_mtl = copy.deepcopy(self.mtl_copy)
            for key in temp_mtl.keys():
                if isinstance(key, int):
                    temp_mtl[key]['fourier_order'] = NF
                    if temp_mtl['type'] == 'coated_wires':
                        temp_mtl[key]['sheath']['fourier_order'] = NF
            
            # Fortran Instance
            self._prepare_fortran_runner(temp_mtl)
            self.runner_silent.run_fortran(self.fortran_base_params)

            # MoM Instance
            bare_wire_mtl = copy.deepcopy(self.mtl_copy)
            bare_wire_mtl['type'] = 'bare_wires'
            for key in bare_wire_mtl.keys():
                if isinstance(key, int):
                    bare_wire_mtl[key]['fourier_order'] = NF
                    bare_wire_mtl[key]['sheath'] = None

            mom_bare_wires = MulticonductorBareWireSystems(bare_wire_mtl)
            mom_bare_wires.run_simulation()

            # Coleta de resultados
            results.append({
                'NF': NF,
                'L (RIBBON.FOR)': self.runner_silent.L_matrix if self.fortran_base_params is not None else np.nan,
                'C (RIBBON.FOR)': self.runner_silent.C_matrix if self.runner_silent.C_matrix is not None else np.nan,
                'C0 (RIBBON.FOR)': self.runner_silent.C0_matrix if self.runner_silent.C0_matrix is not None else np.nan,
                'CGEN (RIBBON.FOR)': self.runner_silent.CGEN_matrix if self.runner_silent.CGEN_matrix is not None else np.nan,
                'C0 (BARE_WIRE.PY)': mom_bare_wires.C_maxwellian if mom_bare_wires.C_maxwellian is not None else np.nan
            })

            print(f"  Complete for NF = {NF}.")            

        self.results_df = pd.DataFrame(results).set_index('NF')


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
        Gera o gráfico de convergência a partir dos resultados armazenados,
        replicando a figura de referência.
        """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        self.nf_max = self.results_df.index.max()

        # Extrai os componentes da matriz de capacitância 'L' e converte para uH/m.
        try:
            l11 = self.results_df['L (RIBBON.FOR)'].apply(lambda x: x[0, 0] * self.l_factor if isinstance(x, np.ndarray) else np.nan)
            l22 = self.results_df['L (RIBBON.FOR)'].apply(lambda x: x[1, 1] * self.l_factor if isinstance(x, np.ndarray) else np.nan)
            l12 = self.results_df['L (RIBBON.FOR)'].apply(lambda x: x[0, 1] * self.l_factor if isinstance(x, np.ndarray) else np.nan)
        except (TypeError, IndexError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados com o estilo da figura de referência
        ax.plot(self.results_df.index, l22, color='k', marker='o', linestyle='-',  label='$L_{22}$')
        ax.plot(self.results_df.index, l11, color='k', marker='o', linestyle=':',  label='$L_{11}$', markerfacecolor='white', markeredgecolor='k')
        ax.plot(self.results_df.index, l12, color='k', marker='s', linestyle='-.', label='$L_{12}$')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Inductance (µH/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.nf_max + 1, 1))
        ax.set_yticks(np.arange(0.2, 1.2, 0.1))
        ax.set_xlim(0.8, self.nf_max + 0.2)
        ax.set_ylim(0.2, 1.1)
        ax.set_title('')
        ax.grid(False)
        ax.legend()        
        plt.tight_layout()


    def plot_paul_fig514b(self):
        """
        Gera o gráfico de convergência a partir dos resultados armazenados,
        replicando a figura de referência.
        """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        self.nf_max = self.results_df.index.max()

        # Extrai os componentes da matriz de capacitância 'C' e converte para pF/m.
        # Por convenção, os elementos Cij (i != j) são negativos. 
        try:
            c11 = self.results_df['C (RIBBON.FOR)'].apply(lambda x: +x[0, 0] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            c22 = self.results_df['C (RIBBON.FOR)'].apply(lambda x: +x[1, 1] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            c12 = self.results_df['C (RIBBON.FOR)'].apply(lambda x: -x[0, 1] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
        except (TypeError, IndexError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados com o estilo da figura de referência
        ax.plot(self.results_df.index, c11, color='k', marker='o', linestyle='-',  label='$C_{11}$')
        ax.plot(self.results_df.index, c22, color='k', marker='o', linestyle=':',  label='$C_{22}$', markerfacecolor='white', markeredgecolor='k')
        ax.plot(self.results_df.index, c12, color='k', marker='s', linestyle='-.', label='$C_{12}$')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.nf_max + 1, 1))
        ax.set_yticks(np.arange(15, 40, 5))
        ax.set_xlim(0.8, self.nf_max + 0.2)
        ax.set_ylim(15, 40)
        ax.set_title('')
        ax.grid(False)
        ax.legend()        
        plt.tight_layout()


    def plot_paul_fig514c(self):
        """
        Gera o gráfico de convergência a partir dos resultados armazenados,
        replicando a figura de referência.
        """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        self.nf_max = self.results_df.index.max()

        # Extrai os componentes da matriz de capacitância 'C0' e converte para pF/m.
        # Por convenção, os elementos C0ij (i != j) são negativos. 
        try:
            c0_11 = self.results_df['C0 (RIBBON.FOR)'].apply(lambda x: +x[0, 0] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            c0_22 = self.results_df['C0 (RIBBON.FOR)'].apply(lambda x: +x[1, 1] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            c0_12 = self.results_df['C0 (RIBBON.FOR)'].apply(lambda x: -x[0, 1] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
        except (TypeError, IndexError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C0 (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados com o estilo da figura de referência
        ax.plot(self.results_df.index, c0_11, color='k', marker='o', linestyle='-',  label='$C0_{11}$')
        ax.plot(self.results_df.index, c0_22, color='k', marker='o', linestyle=':',  label='$C0_{22}$', markerfacecolor='white', markeredgecolor='k')
        ax.plot(self.results_df.index, c0_12, color='k', marker='s', linestyle='-.', label='$C0_{12}$')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Bare-Wire Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.nf_max + 1, 1))
        ax.set_yticks(np.arange(10, 24, 2))
        ax.set_xlim(0.8, self.nf_max + 0.2)
        ax.set_ylim(10, 24)
        ax.set_title('')
        ax.grid(False)
        ax.legend()        
        plt.tight_layout()


    def plot_paul_fig514d(self):
        """ Gera o gráfico de convergência a partir dos resultados armazenados, replicando a figura de referência. """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        self.nf_max = self.results_df.index.max()

        # Extrai os componentes da matriz de capacitância 'C0' e converte para pF/m.
        # Por convenção, os elementos Cgen_ij (i != j) são negativos. 
        try:
            cgen_11 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: +x[0, 0] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            cgen_22 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: +x[1, 1] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            cgen_33 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: +x[2, 2] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            cgen_12 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: -x[0, 1] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            cgen_13 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: -x[0, 2] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            cgen_23 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: -x[1, 2] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
        except (TypeError, IndexError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C0 (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados com o estilo da figura de referência
        ax.plot(self.results_df.index, cgen_11, color='k', marker='o', label='$CGEN_{11}$', linewidth=1.0, zorder=0, markersize=3, linestyle=':')
        ax.plot(self.results_df.index, cgen_22, color='k', marker='o', label='$CGEN_{22}$', linewidth=1.0, zorder=1, markersize=3, linestyle='-')
        ax.plot(self.results_df.index, cgen_33, color='k', marker='o', label='$CGEN_{33}$', linewidth=1.0, zorder=1, markersize=8, linestyle=':', fillstyle='none')
        ax.plot(self.results_df.index, cgen_12, color='k', marker='s', label='$CGEN_{12}$', linewidth=1.0, zorder=0, markersize=4, linestyle='--')
        ax.plot(self.results_df.index, cgen_23, color='k', marker='s', label='$CGEN_{23}$', linewidth=1.0, zorder=1, markersize=8, linestyle='--', fillstyle='none')
        ax.plot(self.results_df.index, cgen_13, color='k', marker='d', label='$CGEN_{13}$', linewidth=1.0, zorder=0, markersize=6, linestyle='-.')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Generalized Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.nf_max + 1, 1))
        ax.set_yticks(np.arange(0, 40, 5))
        ax.set_xlim(0.8, self.nf_max + 0.2)
        ax.set_ylim(0, 40)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper right', fontsize=9, markerscale=1.0, frameon=True)        
        plt.tight_layout()


    def plot_generalized_capacitance(self):
        """ Gera o gráfico de convergência a partir dos resultados armazenados, replicando a figura de referência. """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        self.nf_max = self.results_df.index.max()

        # Extrai os componentes da matriz de capacitância 'C0' e converte para pF/m.
        # Por convenção, os elementos Cgen_ij (i != j) são negativos. 
        try:
            cgen_11 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: +x[0, 0] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            cgen_22 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: +x[1, 1] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            cgen_33 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: +x[2, 2] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            # cgen_12 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: -x[0, 1] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            # cgen_13 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: -x[0, 2] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            # cgen_23 = self.results_df['CGEN (RIBBON.FOR)'].apply(lambda x: -x[1, 2] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
        except (TypeError, IndexError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C0 (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados com o estilo da figura de referência
        ax.plot(self.results_df.index, cgen_11, color='k', marker='o', label='$CGEN_{11}$', linewidth=1.0, zorder=0, markersize=3, linestyle=':')
        ax.plot(self.results_df.index, cgen_22, color='k', marker='o', label='$CGEN_{22}$', linewidth=1.0, zorder=0, markersize=3, linestyle='-')
        ax.plot(self.results_df.index, cgen_33, color='k', marker='o', label='$CGEN_{33}$', linewidth=1.0, zorder=1, markersize=8, linestyle=':', fillstyle='none')
        # ax.plot(self.results_df.index, cgen_12, color='k', marker='s', label='$CGEN_{12}$', linewidth=1.0, zorder=0, markersize=4, linestyle='--')
        # ax.plot(self.results_df.index, cgen_23, color='k', marker='s', label='$CGEN_{23}$', linewidth=1.0, zorder=1, markersize=8, linestyle='--', fillstyle='none')
        # ax.plot(self.results_df.index, cgen_13, color='k', marker='*', label='$CGEN_{13}$', linewidth=1.0, zorder=0, markersize=8, linestyle='-.')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Generalized Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.nf_max + 1, 1))
        ax.set_yticks(np.arange(0, 40, 5))
        ax.set_xlim(0.8, self.nf_max + 0.2)
        ax.set_ylim(0, 40)
        ax.set_title('')
        ax.grid(False)
        ax.legend(loc='upper right', fontsize=9, markerscale=1.0, frameon=True)        
        plt.tight_layout()


    def plot_free_space_capacitance(self):
        """
        Gera o gráfico de convergência a partir dos resultados armazenados,
        replicando a figura de referência.
        """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        self.nf_max = self.results_df.index.max()

        # Extrai os componentes da matriz de capacitância 'C0' e converte para pF/m.
        # Por convenção, os elementos C0ij (i != j) são negativos. 
        try:
            ribbon_c011 = self.results_df['C0 (RIBBON.FOR)'].apply(lambda x: +x[0, 0] * self.c_factor if isinstance(x, np.ndarray) else np.nan)
            bare_c011   = self.results_df['C0 (BARE_WIRE.PY)'].apply(lambda x: x * self.c_factor if isinstance(x, float) else np.nan)
            
        except (TypeError, IndexError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C0 (RIBBON.FOR)' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))

        # Plotagem dos dados com o estilo da figura de referência
        ax.plot(self.results_df.index, ribbon_c011, color='k', marker='o', linestyle=':', markersize=3, label='$C0_{11}$ (RIBBON.FOR)')
        ax.plot(self.results_df.index, bare_c011,   color='k', marker='o', linestyle=':', markersize=8, label='$C0_{11}$ (BARE_WIRE.PY)', fillstyle='none')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Bare-Wire Capacitance (pF/m)', fontsize=12)
        # ax.set_xticks(np.arange(1, self.nf_max + 1, 1))
        # ax.set_yticks(np.arange(10, 24, 2))
        # ax.set_xlim(0.8, self.nf_max + 0.2)
        # ax.set_ylim(10, 24)
        ax.set_title('')
        ax.grid(False)
        ax.legend()        
        plt.tight_layout()