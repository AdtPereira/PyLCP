import re
import subprocess
import numpy as np
from pathlib import Path

class FortranRunner:
    """
    Uma classe para controlar a execução de um programa Fortran legado,
    gerenciando seus arquivos de entrada e saída.
    """
    def __init__(self, exe_path: str, input_filename: str = 'RIBBON.IN', output_filename: str = 'PUL.DAT', silent: bool = True):
        """
        Inicializa o executor do Fortran.

        Args:
            exe_path (str): O caminho para o arquivo executável do Fortran.
            input_filename (str): O nome do arquivo de entrada que o executável espera.
            output_filename (str): O nome do arquivo de saída que o executável gera.
        """
        self.exe_path = Path(exe_path).resolve()
        self.input_filename = self.exe_path.parent / input_filename
        self.output_filename = self.exe_path.parent / output_filename
        self.silent = silent

        if not self.exe_path.is_file():
            raise FileNotFoundError(f"O arquivo executável não foi encontrado em: {self.exe_path}")

        self.L_matrix = None
        self.C_matrix = None
        self.C0_matrix = None

    def _write_input_file(self, params: dict):
        """
        Escreve os parâmetros de simulação no arquivo de entrada do Fortran.
        """
        if not self.silent:
            print(f"Escrevendo arquivo de entrada em: {self.input_filename}")

        param_order = ['N', 'NF', 'IREF', 'RW', 'TD', 'ER', 'S']
        
        with open(self.input_filename, 'w') as f:
            for key in param_order:
                if key not in params:
                    raise ValueError(f"Parâmetro obrigatório '{key}' não encontrado no dicionário de parâmetros.")
                f.write(f"{params[key]}\n")

    def _run_executable(self):
        """
        Executa o programa Fortran compilado.
        """
        if not self.silent:
            print(f"Executando: {self.exe_path}...")
        try:
            result = subprocess.run(
                [str(self.exe_path)],
                cwd=self.exe_path.parent,
                check=True,
                capture_output=True,
                text=True,
                timeout=30
            )
            if not self.silent:
                print("Execução do Fortran concluída com sucesso.")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Erro durante a execução do programa Fortran.")
            print(f"Código de saída: {e.returncode}")
            print(f"Saída de erro (stderr):\n{e.stderr}")
            return False
        except Exception as e:
            print(f"Ocorreu um erro inesperado ao tentar executar o Fortran: {e}")
            return False

    def _parse_output_file(self):
        """
        Lê e analisa o arquivo de saída para extrair as matrizes de resultado.
        """
        if not self.silent:
            print(f"Analisando arquivo de saída: {self.output_filename}")

        if not self.output_filename.is_file():
            print(f"Erro: Arquivo de saída '{self.output_filename}' não foi gerado.")
            return

        results = {'L': {}, 'C': {}, 'C0': {}}
        max_index = 0

        with open(self.output_filename, 'r') as f:
            for line in f:
                line_stripped = line.strip()
                if not line_stripped:
                    continue

                match = re.match(r"(\d+)\s+(\d+)\s+([0-9.E+-]+)\s+=\s*([LC0]+)\(", line_stripped)
                if match:
                    i, j, value, key = match.groups()
                    i, j, value = int(i), int(j), float(value)
                    max_index = max(max_index, i, j)
                    if i > j: i, j = j, i
                    results[key][(i, j)] = value
        
        matrix_size = max_index
        if matrix_size > 0:
            self.L_matrix = np.zeros((matrix_size, matrix_size))
            self.C_matrix = np.zeros((matrix_size, matrix_size))
            self.C0_matrix = np.zeros((matrix_size, matrix_size))
            
            for (i, j), val in results.get('L', {}).items():
                self.L_matrix[i-1, j-1] = self.L_matrix[j-1, i-1] = val
            for (i, j), val in results.get('C', {}).items():
                self.C_matrix[i-1, j-1] = self.C_matrix[j-1, i-1] = val
            for (i, j), val in results.get('C0', {}).items():
                self.C0_matrix[i-1, j-1] = self.C0_matrix[j-1, i-1] = val
        
        if not self.silent:
            print("Análise do arquivo de saída concluída.")

    def run_simulation(self, params: dict):
        """
        Orquestra o processo completo: escreve o input, executa e analisa o output.
        """
        try:
            self._write_input_file(params)
            if self._run_executable():
                self._parse_output_file()
            else:
                if not self.silent:
                    print("Simulação abortada devido a erro na execução.")
        except Exception as e:
            print(f"Ocorreu um erro inesperado durante a simulação: {e}")


# Exemplo de uso da classe FortranRunner
# if __name__ == '__main__':
#     import os
#     os.system('cls' if os.name == 'nt' else 'clear')
#     simulation_params = {
#         'N': 2,
#         'NF': 10,
#         'IREF': 1,
#         'RW': 1.000E-02,
#         'TD': 1.000E-02,
#         'ER': 3.0,
#         'S': 4.000E-02
#     }

#     print("Iniciando simulação controlada por Python a partir da raiz do projeto...")
    
#     script_dir = Path(__file__).resolve().parent
#     project_root = script_dir.parent
#     fortran_exe_path = project_root / 'mtl_paul' / 'RIBBON' / 'RIBBON.EXE'
    
#     runner = FortranRunner(exe_path=fortran_exe_path)
#     runner.run_simulation(simulation_params)

#     print("\n" + "="*40)
#     print("      Resultados da Simulação Fortran")
#     print("="*40)
#     if runner.L_matrix is not None:
#         print("\nMatriz de Indutância (L):")
#         print(runner.L_matrix)
        
#         print("\nMatriz de Capacitância (C):")
#         print(runner.C_matrix)
        
#         print("\nMatriz de Capacitância no Vácuo (C0):")
#         print(runner.C0_matrix)
#     else:
#         print("\nNenhum resultado foi analisado. Verifique os logs de erro.")