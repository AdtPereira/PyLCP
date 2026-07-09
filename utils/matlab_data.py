import sys
import os
import numpy as np
from pathlib import Path

import scipy.io as sio


class MatlabDataReader:
    """
    Uma classe dedicada para ler e analisar todos os arquivos .mat do MATLAB
    do diretório 'Results' de um caso específico, de forma análoga ao
    ComsolDataReader (utils/comsol_data.py).

    Ela descobre e carrega automaticamente todos os arquivos .mat, retornando-os
    em um dicionário estruturado.
    """

    def __init__(self, script_file_path: str, autoShow: bool = True):
        """
        Inicializa o leitor identificando o diretório 'Results' alvo.

        Args:
            script_file_path (str): O __file__ do script que está chamando
                                     (usado para localizar o diretório do caso).
            autoShow (bool): Se True, exibe um resumo dos dados carregados.
        """
        project_root = Path(script_file_path).resolve().parents[2]
        if str(project_root) not in sys.path:
            sys.path.insert(0, str(project_root))

        case_name = os.path.splitext(os.path.basename(script_file_path))[0]

        self.project_root = project_root
        self.case_name = case_name
        self.results_path = project_root / 'testData' / case_name / 'Results'

        print(f"Project root configured at: {project_root}")
        print(f"\nInstanciando MatlabDataReader para o caso '{case_name}' ---")

        if not self.results_path.is_dir():
            print(f"  Aviso: diretório Results não encontrado para o caso '{case_name}'. Dados MATLAB ignorados.")

        self.data = self.load_all_results()

        if autoShow and self.data:
            self.show_summary()

    def load_all_results(self) -> dict:
        """
        Verifica o diretório 'Results', carrega todos os arquivos .mat e os retorna
        como um dicionário.

        Returns:
            Um dicionário onde as chaves são os nomes dos arquivos (sem a extensão
            .mat). O valor é o array NumPy da variável, caso o arquivo contenha uma
            única variável de dados (caso comum), ou um dicionário {nome: array}
            caso contenha mais de uma.
        """
        if not self.results_path.is_dir():
            return {}

        mat_files = list(self.results_path.glob('*.mat'))
        data = {}

        if not mat_files:
            print(f"Aviso: Nenhum arquivo .mat encontrado em {self.results_path}")
            return {}

        print(f"Encontrado(s) {len(mat_files)} arquivo(s) .mat no diretório Results do caso '{self.case_name}'.")

        for file_path in mat_files:
            file_stem = file_path.stem
            try:
                data[file_stem] = self._parse_single_file(file_path)
            except Exception as e:
                print(f"Erro ao analisar o arquivo {file_path.name}: {e}")

        return data

    def _parse_single_file(self, file_path: Path):
        """
        Carrega um único arquivo .mat (formato MATLAB level 5, via scipy.io.loadmat)
        e retorna suas variáveis, ignorando as chaves de metadados do MATLAB
        ('__header__', '__version__', '__globals__').

        Matrizes 3D quadradas exportadas no formato MATLAB (N, N, num_freq) são
        reordenadas para (num_freq, N, N), a convenção usada pelas demais matrizes
        PUL do pyLCP (ver, por ex., InternalPerUnitParameters).
        """
        print(f"  -> Carregando e analisando: {file_path.name}...")

        raw = sio.loadmat(file_path)
        variables = {}
        for name, value in raw.items():
            if name.startswith('__'):
                continue
            if isinstance(value, np.ndarray) and value.ndim == 3 and value.shape[0] == value.shape[1]:
                value = np.moveaxis(value, -1, 0)
            variables[name] = value

        if len(variables) == 1:
            return next(iter(variables.values()))
        return variables

    def show_summary(self):
        """
        Exibe um resumo de todos os dados carregados (nomes de variáveis, shapes
        e dtypes).
        """
        if not self.data:
            print("Nenhum dado carregado para exibir o resumo. Execute 'load_all_results()' primeiro.")
            return

        print(f"\n--- Resumo dos Dados Carregados para o Caso '{self.case_name}' ---")
        for name, value in self.data.items():
            print(f"\n==================================================")
            print(f"  Arquivo: '{name}.mat'")
            print(f"==================================================")

            variables = value if isinstance(value, dict) else {name: value}
            for var_name, arr in variables.items():
                shape = getattr(arr, 'shape', None)
                dtype = getattr(arr, 'dtype', type(arr))
                print(f"  {var_name}: shape={shape}, dtype={dtype}")
