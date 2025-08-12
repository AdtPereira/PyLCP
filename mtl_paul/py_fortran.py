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

        # Inicializa todos os atributos da matriz para refletir a nomenclatura final
        self.A_matrix = None
        self.B_matrix = None
        self.C_matrix = None
        self.D_matrix = None
        self.CGEN0_matrix = None
        self.IND_matrix = None
        self.CAP_matrix = None
        self.CAP0_matrix = None
        self.CGEN_matrix = None
        self.NF = None

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

        raw_values = {}
        max_indices = {}
        
        pattern = re.compile(r"(\d+)\s+(\d+)\s+([0-9.E+-]+)\s+=\s*([A-Z0-9]+)\(")

        with open(self.output_filename, 'r') as f:
            for line in f:
                match = pattern.match(line.strip())
                if match:
                    i_str, j_str, val_str, key = match.groups()
                    i, j, value = int(i_str), int(j_str), float(val_str)
                    
                    if key not in raw_values:
                        raw_values[key] = {}
                        max_indices[key] = 0

                    raw_values[key][(i, j)] = value
                    max_indices[key] = max(max_indices[key], i, j)
        
        # Define quais matrizes são simétricas com base nos nomes no arquivo PUL.DAT
        symmetric_keys = ['IND', 'CAP', 'CAP0', 'CGEN', 'CGEN0']

        for key, values in raw_values.items():
            size = max_indices.get(key, 0)
            if size > 0:
                matrix = np.zeros((size, size))
                is_symmetric = key in symmetric_keys
                
                for (i, j), val in values.items():
                    matrix[i-1, j-1] = val
                    if is_symmetric:
                        matrix[j-1, i-1] = val
                
                setattr(self, f"{key}_matrix", matrix)

        if not self.silent:
            print("Análise do arquivo de saída concluída.")

    def run_fortran(self, params: dict):
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