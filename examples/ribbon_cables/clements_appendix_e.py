import sys
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.constants import epsilon_0
from typing import Dict, Any

# --- ETAPA 1: CONFIGURAÇÃO DE CAMINHOS E IMPORTAÇÕES ---
def setup_project_paths() -> Path:
    """
    Determina a raiz do projeto e adiciona os diretórios necessários ao path do Python.
    
    Returns:
        Path: O objeto Path para a raiz do projeto.
    """
    try:
        script_dir = Path(__file__).resolve().parent
        project_root = script_dir.parents[1]
        
        # Adiciona os diretórios dos módulos ao path do sistema
        sys.path.append(str(project_root))
        sys.path.append(str(project_root / 'mom'))
        
        print("Caminhos do projeto configurados com sucesso.")
        return project_root
    except IndexError:
        raise FileNotFoundError(
            "Não foi possível encontrar a raiz do projeto. "
            "Certifique-se de que o script está em 'examples/coated_wires'."
        )

# Executa a configuração de caminhos
project_root = setup_project_paths()

try:
    from mom.coated_wire_systems import TwoCoatedWireSystem
    from mtl_paul.py_fortran import FortranRunner
    from mtl_data.models import BIFILAR_COATED_WIRE_S4000 as MTL
    print("Módulos e modelo de dados importados com sucesso.")
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

class ConvergenceAnalyzer:
    """
    Encapsula a lógica para executar e analisar o estudo de convergência
    de capacitância, comparando MoM Python e Fortran.
    """
    C_FACTOR = 1e12  # Fator de conversão para pF/m

    def __init__(self, mtl_config: Dict[str, Any], NF_MAX: int = 10):
        """
        Inicializa o analisador de convergência.

        Args:
            mtl_config (Dict[str, Any]): Dicionário com a configuração do modelo MTL.
            nf_max (int): Número máximo de coeficientes/ordem harmônica para testar.
        """
        self.mtl_config = mtl_config
        self.nf_max = NF_MAX
        self.k_range = range(1, NF_MAX + 1)
        self.results_df = None

        # Extrai parâmetros e prepara o executor do Fortran
        self._prepare_fortran_runner()
        self._calculate_analytical_solution()

    def _prepare_fortran_runner(self):
        """Prepara os parâmetros e o executor para a simulação Fortran."""
        center_0 = np.array(self.mtl_config[0]['center_point'])
        center_1 = np.array(self.mtl_config[1]['center_point'])
        
        self.fortran_base_params = {
            'N': 2,
            'IREF': 1,
            'RW': self.mtl_config[0]['radius'][1],
            'TD': self.mtl_config[0]['sheath']['thickness'],
            'ER': self.mtl_config[0]['sheath']['relative_permittivity'],
            'S': np.linalg.norm(center_0 - center_1)
        }
        
        fortran_exe_path = project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
        self.fortran_runner = FortranRunner(exe_path=str(fortran_exe_path))
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

        for k in self.k_range:
            print(f"  Executando para k = {k}...")
            
            # Simulação MoM Python
            temp_mtl = copy.deepcopy(self.mtl_config)
            for key in [0, 1]:
                temp_mtl[key]['fourier_order'] = k
                temp_mtl[key]['sheath']['fourier_order'] = k
            
            sim_python = TwoCoatedWireSystem(temp_mtl)
            sim_python.run_simulation()
            
            # Simulação Fortran
            temp_fortran_params = self.fortran_base_params.copy()
            temp_fortran_params['NF'] = k
            self.fortran_runner.run_fortran(temp_fortran_params)
            
            # Coleta de resultados
            results.append({
                'NF': k,
                'C (TwoCoatedWireSystem)': sim_python.C_maxwellian * self.C_FACTOR,
                'C (RIBBON.FOR)': self.fortran_runner.C_matrix[0, 0] * self.C_FACTOR if self.fortran_runner.C_matrix is not None else np.nan,
                'C0 (RIBBON.FOR)': self.fortran_runner.C0_matrix[0, 0] * self.C_FACTOR if self.fortran_runner.C0_matrix is not None else np.nan
            })

        self.results_df = pd.DataFrame(results).set_index('NF')

    def plot_results(self):
        """Gera o gráfico de convergência a partir dos resultados armazenados."""
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return
            
        PLOT_CONFIG = {
            'C (TwoCoatedWireSystem)': {'color': 'blue', 'marker': 'o', 'linestyle': '-', 'label': 'MoM Python C (Fios Revestidos)'},
            'C (RIBBON.FOR)': {'color': 'green', 'marker': 'x', 'linestyle': ':', 'label': 'Fortran C (Fios Revestidos)'},
            'C0 (RIBBON.FOR)': {'color': 'purple', 'marker': 's', 'linestyle': '--', 'label': 'Fortran C0 (Fios Nus)', 'fillstyle': 'none'},
        }

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plota cada curva a partir do DataFrame
        for col, style in PLOT_CONFIG.items():
            ax.plot(self.results_df.index, self.results_df[col], **style)

        # Plota a linha da solução analítica
        ax.axhline(y=self.analytical_bare_capacitance * self.C_FACTOR, color='red', linestyle='--', 
                   label=f'Analítica (Fios Nus) = {self.analytical_bare_capacitance*self.C_FACTOR:.4f} pF/m')
        
        ax.set_title(f'Estudo de Convergência da Capacitância para D/R = {self.D_over_R:.2f}', fontsize=14)
        ax.set_xlabel('Ordem Harmônica (k) / Coeficientes de Fourier (NF)', fontsize=12)
        ax.set_ylabel('Capacitância (pF/m)', fontsize=12)
        ax.set_xticks(np.arange(0, self.nf_max + 2, 2))
        ax.set_xlim(0, self.nf_max + 1)
        ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.7)
        ax.legend(fontsize=10)
        plt.tight_layout()

    def display_results_table(self):
        """Exibe a tabela de resultados formatada no terminal."""
        if self.results_df is None:
            print("Execute as simulações primeiro com 'run_study()'.")
            return
            
        pd.options.display.float_format = '{:.4f}'.format
        DIST1 = 75
        DIST2 = 10
        print("\n\n" + "="*DIST1)
        print(" " * DIST2 + "TWO DIELECTRIC-COATED WIRE SYSTEM CONVERGENCE")
        print(" " * DIST2 + 
                f"D = {self.fortran_base_params['S']:.3f} m, "
                f"RW = {self.fortran_base_params['RW']:.3f} m, "
                f"TD = {self.fortran_base_params['TD']:.3f} m, "
                f"ER = {self.fortran_base_params['ER']:.1f}")
        print(" " * DIST2 + f"Exact Capacitance (Bare Wires) = {self.analytical_bare_capacitance * self.C_FACTOR:.4f} pF/m")
        print(" " * DIST2 + "Results in [pF/m]")
        print("="*DIST1)
        print(self.results_df)
        print("="*DIST1)

if __name__ == '__main__':
    analyzer = ConvergenceAnalyzer(mtl_config=MTL, NF_MAX=15)
    analyzer.run_simulations()
    analyzer.plot_results()
    analyzer.display_results_table()
    plt.show()
