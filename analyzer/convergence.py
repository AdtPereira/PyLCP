import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any
from scipy.constants import epsilon_0

from mtl_data.mtl import MulticonductorTransmissionLine as MTL
from mtl_paul.py_fortran import FortranRunner


class ConvergenceAnalyzer(MTL):
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
        super().__init__(mtl)
        self.project_root = project_root
        self.nf_max = NF_MAX
        self.results_df = None
        self.c_factor = 1e12 # Fator de conversão para pF/m

        # Extrai parâmetros e prepara o executor do Fortran
        self._prepare_fortran_runner()
        self._calculate_analytical_solution()

    def _prepare_fortran_runner(self):
        """Prepara os parâmetros e o executor para a simulação Fortran."""
        center_0 = np.array(self.mtl[0]['center_point'])
        center_1 = np.array(self.mtl[1]['center_point'])
        
        self.fortran_base_params = {
            'N':    sum([key for key in self.mtl.keys() if isinstance(key, int)]),
            'NF':   self.NF,
            'IREF': self.idx_ref,
            'RW':   self.mtl[self.idx_ref]['radius'][1],
            'TD':   self.mtl[self.idx_ref]['sheath']['thickness'],
            'ER':   self.mtl[self.idx_ref]['sheath']['relative_permittivity'],
            'S':    np.linalg.norm(center_0 - center_1)
        }
        
        fortran_exe_path = self.project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.fortran_runner = FortranRunner(exe_path=str(fortran_exe_path), silent=False)
        self.fortran_runner_silent = FortranRunner(exe_path=str(fortran_exe_path), silent=True)
        self.D_over_R = self.fortran_base_params['S'] / self.fortran_base_params['RW']

    def _calculate_analytical_solution(self):
        """Calcula a solução analítica para fios nus como referência."""
        R = self.fortran_base_params['RW']
        D = self.fortran_base_params['S']
        self.analytical_bare_capacitance = (np.pi * epsilon_0) / np.arccosh(D / (2 * R))

    def run_simulations(self):
        """Executa o laço de convergência para ambas as simulações e armazena os resultados."""
        print(f"\nIniciando estudo de convergência até k={self.nf_max}...")
        
        results = []
        for NF in range(1, self.nf_max + 1):
            print(f"  Executando para k = {NF}...")
            
            temp_mtl = copy.deepcopy(self.mtl)
            for key in self.mtl.keys():
                temp_mtl[key]['fourier_order'] = NF
                temp_mtl[key]['sheath']['fourier_order'] = NF
            
            # Simulação Fortran
            temp_fortran_params = self.fortran_base_params.copy()
            temp_fortran_params['NF'] = NF
            self.fortran_runner_silent.run_fortran(temp_fortran_params)
            
            # Coleta de resultados
            results.append({
                'NF': NF,
                'C (RIBBON.FOR)': self.fortran_runner_silent.C_matrix if self.fortran_runner_silent.C_matrix is not None else np.nan,
                'C0 (RIBBON.FOR)': self.fortran_runner_silent.C0_matrix if self.fortran_runner_silent.C0_matrix is not None else np.nan,
                'CGEN (RIBBON.FOR)': self.fortran_runner_silent.CGEN_matrix if self.fortran_runner_silent.CGEN_matrix is not None else np.nan,
            })

        self.results_df = pd.DataFrame(results).set_index('NF')

    def run_single_simulation(self):
        """Executa uma simulação única para um valor específico de k."""
        
        self.fortran_runner.run_fortran(self.fortran_base_params)
        # py_mom = TwoCoatedWireSystem(self.mtl)
        # py_mom.run_simulation()

        dist1, dist2 = 59, 8
        print("\n" + "="*dist1)
        print("="*dist2 + " THREE-WIRE RIBBON CABLE SYSTEM SIMULATION " + "="*dist2)
        print("="*dist1)
        
        print("\n")
        print("="*dist2 + "                 RIBBON.FOR                " + "="*dist2)
        if self.fortran_runner.L_matrix is not None:
            print("\nMatriz de Indutância (L):")
            print(self.fortran_runner.L_matrix)
            
            print("\nMatriz de Capacitância (C):")
            print(self.fortran_runner.C_matrix)

            print("\nMatriz de Capacitância no Vácuo (C0):")
            print(self.fortran_runner.C0_matrix)

            print("\nMatriz de Capacitância Generalizada (CGEN):")
            print(self.fortran_runner.CGEN_matrix)
        else:
            print("\nNenhum resultado foi analisado. Verifique os logs de erro.")

        # print("\n")
        # print("="*dist2 + "           TWOCOATEDWIRESYSTEM.PY          " + "="*dist2)
        # if py_mom.C_generalized is not None:
        #     print("\nMatriz de Capacitância Generalizada (CGEN):")
        #     print(py_mom.C_generalized)
        # else:
        #     print("\nNenhum resultado foi analisado. Verifique os logs de erro.")

    def plot_results(self):
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
        fig, ax = plt.subplots(figsize=(8, 6))

        # Plotagem dos dados com o estilo da figura de referência
        ax.plot(self.results_df.index, c11, color='k', marker='o', linestyle='-', label='C11')
        ax.plot(self.results_df.index, c22, color='k', marker='o', linestyle=':', label='C22', markerfacecolor='white', markeredgecolor='k')
        ax.plot(self.results_df.index, c12, color='k', marker='s', linestyle='-.', label='C12')

        # Configuração dos eixos para corresponder à imagem de referência
        ax.set_xlabel('Number of Fourier Coefficients', fontsize=12)
        ax.set_ylabel('Capacitance (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(1, self.nf_max + 1, 1))
        ax.set_xlim(0.8, self.nf_max + 0.2)
        ax.set_ylim(15, 40)
        
        # Remove título e grade para um visual mais limpo
        ax.set_title('')
        ax.grid(False)
        ax.legend(fontsize=11)        
        plt.tight_layout()
