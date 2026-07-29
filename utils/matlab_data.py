import sys
import os
import numpy as np
from pathlib import Path
from typing import Sequence

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

    @staticmethod
    def _reorder_conductor_matrix(matrix, conductor_order: Sequence[int]):
        """
        O MATLAB exporta os condutores agrupados por tipo: [core_A, core_B, core_C,
        sheath_A, sheath_B, sheath_C]. O pyLCP monta suas matrizes (interna e
        quasi-TEM) agrupadas por cabo: [core_A, sheath_A, core_B, sheath_B,
        core_C, sheath_C] (ver np.kron(np.identity(N), Zij) em
        InternalPerUnitParameters.matrices()). Sem essa reordenação, M[p,q]
        (pyLCP) e Z[p,q] (MATLAB) apontam para pares de condutores fisicamente
        diferentes para os mesmos índices (p, q).
        """
        if matrix is None:
            return None
        return matrix[:, conductor_order, :][:, :, conductor_order]

    def get_scc_scenario_data(self, prefix: str, conductor_order: Sequence[int]) -> dict:
        """
        Monta o dicionário de dados de referência do MATLAB (no formato usado por
        pul_data['matlab']) para um caso SCC (núcleo + bainha), a partir dos
        arquivos exportados com o prefixo `prefix` (ex.: 'andreata' ->
        'andreata_frequency_range', 'andreata_series_impedance_matrix', ...).

        `conductor_order` reordena os condutores da convenção MATLAB (agrupada
        por tipo) para a convenção pyLCP (agrupada por cabo); ver
        `_reorder_conductor_matrix`. As duas matrizes de retorno pelo solo (Zg,
        Pg) seguem o mesmo layout tipo-agrupado das demais -- núcleo e bainha
        do mesmo cabo têm entradas redundantes (idênticas), já que o retorno
        pelo solo só depende da posição do cabo, não de qual condutor dentro
        dele -- por isso usam a mesma `_reorder_conductor_matrix`, sem recorte
        especial.

        Retorna 'frequencies' como None se o arquivo correspondente não for
        encontrado -- o fallback (ex.: para pul_data['frequencies']) fica a
        cargo do chamador.
        """
        frequencies = self.data.get(f'{prefix}_frequency_range')

        internal_z = self._reorder_conductor_matrix(self.data.get(f'{prefix}_internal_impedance_matrix'), conductor_order)
        series_z = self._reorder_conductor_matrix(self.data.get(f'{prefix}_series_impedance_matrix'), conductor_order)
        shunt_y = self._reorder_conductor_matrix(self.data.get(f'{prefix}_shunt_admittance_matrix'), conductor_order)
        internal_y = self._reorder_conductor_matrix(self.data.get(f'{prefix}_internal_admittance_matrix'), conductor_order)
        earth_return_pg = self._reorder_conductor_matrix(self.data.get(f'{prefix}_earth_return_potential_coefficient_matrix'), conductor_order)
        earth_return_zg = self._reorder_conductor_matrix(self.data.get(f'{prefix}_earth_return_impedance_matrix'), conductor_order)

        return {
            'frequencies': frequencies.flatten() if frequencies is not None else None,
            'scenarios': {
                'measured': {
                    'internal_impedance_matrix': internal_z,
                    'series_impedance_matrix': series_z,
                    'shunt_admittance_matrix': shunt_y,
                    'internal_admittance_matrix': internal_y,
                    'earth_return_potential_coefficient_matrix': earth_return_pg,
                    'earth_return_impedance_matrix': earth_return_zg,
                },
            },
        }
