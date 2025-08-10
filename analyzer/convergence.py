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
        self.sum_max = SUM_MAX
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
            
            # Fortran Instance
            self._prepare_fortran_runner(temp_mtl)
            self.runner_silent.run_fortran(self.fortran_base_params)

            # MoM Instance
            bare_wire_mtl = copy.deepcopy(self.mtl_copy)
            bare_wire_mtl['type'] = 'bare_wires'
            for key in bare_wire_mtl.keys():
                if isinstance(key, int):
                    bare_wire_mtl[key]['fourier_order'] = k
                    bare_wire_mtl[key]['sheath'] = None

            mom_bare_wires = MulticonductorBareWireSystems(bare_wire_mtl)
            mom_bare_wires.run_simulation()

            # Coleta de resultados
            results.append({
                'k': k,
                'L (RIBBON.FOR)':       self.runner_silent.L_matrix if self.fortran_base_params is not None else np.nan,
                'C (RIBBON.FOR)':       self.runner_silent.C_matrix if self.runner_silent.C_matrix is not None else np.nan,
                'C0 (RIBBON.FOR)':      self.runner_silent.C0_matrix if self.runner_silent.C0_matrix is not None else np.nan,
                'CGEN (RIBBON.FOR)':    self.runner_silent.CGEN_matrix if self.runner_silent.CGEN_matrix is not None else np.nan,
                'C0 (BARE_WIRE.PY)':    mom_bare_wires.C_maxwellian if mom_bare_wires.C_maxwellian is not None else np.nan,
                'CGEN (BARE_WIRE.PY)':  mom_bare_wires.C_generalized if mom_bare_wires.C_generalized is not None else np.nan,
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


    def plot_bifilar_generalized_capacitance(self):
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
            c11_pymom  = self._extract_matrix_element('CGEN (BARE_WIRE.PY)', row=0, col=0)
            c22_pymom  = self._extract_matrix_element('CGEN (BARE_WIRE.PY)', row=1, col=1)
            c12_pymom  = self._extract_matrix_element('CGEN (BARE_WIRE.PY)', row=0, col=1)
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
        """ Gera o gráfico de convergência da capacitância generalizada, adaptando-se ao número de condutores do sistema. """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Determina o número de condutores (N) a partir do dicionário mtl
        # Conta quantas chaves no dicionário são inteiros (0, 1, 2, ...)
        N = len([key for key in self.mtl_copy.keys() if isinstance(key, int)])

        # Dicionários para armazenar as séries de dados extraídas
        c_ribbon = {}
        c_pymom = {}

        try:
            # Extração de dados de forma programática
            for i in range(N):
                # Elementos da diagonal principal (auto-capacitâncias C_ii)
                c_ribbon[f'c_{i}{i}'] = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=i, col=i)
                c_pymom[f'c_{i}{i}'] = self._extract_matrix_element('CGEN (BARE_WIRE.PY)', row=i, col=i)
                
                # Elementos fora da diagonal (capacitâncias mútuas C_ij)
                for j in range(i + 1, N):
                    c_ribbon[f'c_{i}{j}'] = self._extract_matrix_element('CGEN (RIBBON.FOR)', row=i, col=j)
                    c_pymom[f'c_{i}{j}'] = self._extract_matrix_element('CGEN (BARE_WIRE.PY)', row=i, col=j)

        except (TypeError, IndexError, ValueError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'CGEN' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5)) # Aumentar a figura para acomodar mais legendas

        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1

        # 1. Encontra o valor máximo de NF para os dados do Fortran
        max_nf_fortran = fortran_nf_axis.max()

        # 2. Cria uma máscara para filtrar os dados do Python que excedem esse limite
        mask_py = python_nf_axis <= max_nf_fortran
        
        # Marcadores e estilos para diferenciar as plots
        markers = ['o', 's', '^', 'd', 'v', '<', '>']
        linestyles = [':', '-.', '--', '-']

        # Plotagem programática dos dados
        legend_items_count = 0
        for i in range(N):
            # Plot da diagonal principal
            ax.plot(fortran_nf_axis, c_ribbon[f'c_{i}{i}'], color='k', marker=markers[i % len(markers)], 
                    linestyle=linestyles[0], label=f'$C_{{{i}{i}}} (.FOR)$')
            ax.plot(python_nf_axis[mask_py], c_pymom[f'c_{i}{i}'][mask_py], color='gray', marker=markers[i % len(markers)], 
                    linestyle=linestyles[1], label=f'$C_{{{i}{i}}} (.PY)$', fillstyle='none')
            legend_items_count += 2

            # Plot dos elementos fora da diagonal
            for j in range(i + 1, N):
                ax.plot(fortran_nf_axis, c_ribbon[f'c_{i}{j}'], color='k', marker=markers[(i+j) % len(markers)], 
                        linestyle=linestyles[2], label=f'$C_{{{i}{j}}} (.FOR)$')
                ax.plot(python_nf_axis[mask_py], c_pymom[f'c_{i}{j}'][mask_py], color='gray', marker=markers[(i+j) % len(markers)], 
                        linestyle=linestyles[3], label=f'$C_{{{i}{j}}} (.PY)$', fillstyle='none')
                legend_items_count += 2

        # Configurações do gráfico
        ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=12)
        ax.set_ylabel('Generalized Capacitance Matrix (pF/m)', fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # 4. Define os limites e os ticks do eixo X com base no NF máximo do Fortran
        ax.set_xlim(0.8, max_nf_fortran + 0.2)
        ax.set_xticks(np.arange(1, max_nf_fortran + 1, 1))
        ax.legend(loc='upper center', fontsize=9, ncol=legend_items_count, bbox_to_anchor=(0.5, 1.1), fancybox=True)
        plt.tight_layout()


    def plot_free_space_capacitance_matrix(self):
        """ Gera o gráfico de convergência da capacitância do espaço livro, adaptando-se ao número de condutores do sistema. """
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return

        # Garante que nf_max seja consistente com os dados
        assert self.results_df.index.max() == self.sum_max - 1, \
            f"A simulação não rodou até o valor máximo esperado de k={self.sum_max - 1}"

        # Determina o número de condutores (N) a partir do dicionário mtl
        # Conta quantas chaves no dicionário são inteiros (0, 1, 2, ...)
        N = len([key for key in self.mtl_copy.keys() if isinstance(key, int)])

        # Dicionários para armazenar as séries de dados extraídas
        c_ribbon = {}
        c_pymom = {}

        try:
            # Extração de dados de forma programática
            for i in range(N-1):
                # Elementos da diagonal principal (auto-capacitâncias C_ii)
                c_ribbon[f'c_{i}{i}'] = self._extract_matrix_element('C0 (RIBBON.FOR)', row=i, col=i)
                c_pymom[f'c_{i}{i}'] = self._extract_matrix_element('C0 (BARE_WIRE.PY)', row=i, col=i)

                # Elementos fora da diagonal (capacitâncias mútuas C_ij)
                for j in range(i + 1, N-1):
                    c_ribbon[f'c_{i}{j}'] = self._extract_matrix_element('C0 (RIBBON.FOR)', row=i, col=j)
                    c_pymom[f'c_{i}{j}'] = self._extract_matrix_element('C0 (BARE_WIRE.PY)', row=i, col=j)

        except (TypeError, IndexError, ValueError) as e:
            print(f"Erro ao extrair elementos da matriz de capacitância: {e}")
            print("Verifique se as simulações foram executadas e se a matriz 'C' foi populada.")
            return

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5)) # Aumentar a figura para acomodar mais legendas

        fortran_nf_axis = self.results_df.index + 1
        python_nf_axis = 2 * self.results_df.index + 1

        # 1. Encontra o valor máximo de NF para os dados do Fortran
        max_nf_fortran = fortran_nf_axis.max()

        # 2. Cria uma máscara para filtrar os dados do Python que excedem esse limite
        mask_py = python_nf_axis <= max_nf_fortran
        
        # Marcadores e estilos para diferenciar as plots
        markers = ['o', 's', '^', 'd', 'v', '<', '>']
        linestyles = [':', '-.', '--', '-']

        # Plotagem programática dos dados
        legend_items_count = 0
        for i in range(N-1):
            # Plot da diagonal principal
            ax.plot(fortran_nf_axis, c_ribbon[f'c_{i}{i}'], color='k', marker=markers[i % len(markers)], 
                    linestyle=linestyles[0], label=f'$C_{{{i+1}{i+1}}} (.FOR)$')
            ax.plot(python_nf_axis[mask_py], c_pymom[f'c_{i}{i}'][mask_py], color='gray', marker=markers[i % len(markers)], 
                    linestyle=linestyles[1], label=f'$C_{{{i+1}{i+1}}} (.PY)$', fillstyle='none')
            legend_items_count += 2

            # Plot dos elementos fora da diagonal
            for j in range(i + 1, N-1):
                ax.plot(fortran_nf_axis, c_ribbon[f'c_{i}{j}'], color='k', marker=markers[(i+j) % len(markers)], 
                        linestyle=linestyles[2], label=f'$C_{{{i+1}{j+1}}} (.FOR)$')
                ax.plot(python_nf_axis[mask_py], c_pymom[f'c_{i}{j}'][mask_py], color='gray', marker=markers[(i+j) % len(markers)], 
                        linestyle=linestyles[3], label=f'$C_{{{i+1}{j+1}}} (.PY)$', fillstyle='none')
                legend_items_count += 2

        # Configurações do gráfico
        ax.set_xlabel('Number of Fourier Coefficients (NF)', fontsize=12)
        ax.set_ylabel('Capacitance Matrix (pF/m)', fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # 4. Define os limites e os ticks do eixo X com base no NF máximo do Fortran
        ax.set_xlim(0.8, max_nf_fortran + 0.2)
        ax.set_xticks(np.arange(1, max_nf_fortran + 1, 1))
        ax.legend(loc='upper center', fontsize=9, ncol=legend_items_count, bbox_to_anchor=(0.5, 1.1), fancybox=True)
        plt.tight_layout()