import re
import sys
import os
import re
import numpy as np
import pandas as pd
from pathlib import Path

import scipy.constants as sc
from scipy.linalg import lu_factor, lu_solve
from utils.case_utils import *
from analytical_forms.single_core_cable import InternalPerUnitParameters, _expand_by_block_sizes

class MergedComsolDataReader:
    """
    A dedicated class to read and parse COMSOL text files containing
    two or more datasets concatenated in sequence.

    It automatically detects headers, splits the data where the first column
    resets to zero, and provides clean pandas DataFrames as output.

    Attributes:
        df1 (pd.DataFrame): DataFrame containing the first dataset.
        df2 (pd.DataFrame): DataFrame containing the second dataset.
        column_names (list): A list of the cleaned column names found.
    """
    def __init__(self, file_path: str):
        """
        Initializes the reader and processes the file.

        Args:
            file_path (str): The full path to the .txt data file.
        """
        self.file_path = Path(file_path)
        self.df1 = pd.DataFrame()
        self.df2 = pd.DataFrame()
        self.column_names = []        
        self._parse_file()

    def _parse_file(self):
        """
        The core private method that handles the entire file parsing logic.
        """
        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {self.file_path}")

        header_lines_raw = []
        data_lines = []
        
        # This logic correctly separates the header from the data.
        with open(self.file_path, 'r', encoding='utf-8') as f:
            is_header_section = True
            for line in f:
                stripped_line = line.strip()
                if not stripped_line:
                    continue
                
                is_data = False
                try:
                    # A data line typically starts with a number.
                    float(stripped_line.split()[0])
                    is_data = True
                except (ValueError, IndexError):
                    is_data = False
                
                if is_header_section and not is_data:
                    header_lines_raw.append(stripped_line)
                else:
                    # Once we find the first data line, everything after it is data.
                    is_header_section = False
                    data_lines.append(stripped_line)

        # 1. Filter out metadata to isolate potential column name lines.
        metadata_keywords = ['Model:', 'Version:', 'Date:', 'Table:', 'Dimension:', 'Nodes:', 'Expressions:', 'Description:']
        column_name_candidates = []
        for line in header_lines_raw:
            clean_line = line.replace('%', '').strip()
            if clean_line and not any(keyword in line for keyword in metadata_keywords):
                column_name_candidates.append(clean_line)

        if not column_name_candidates:
            raise ValueError("Could not parse valid column headers from the file.")

        # 2. Join all candidate lines and then split them by multiple spaces.
        # This robustly handles headers on the same line OR across multiple lines.
        full_header_string = " ".join(column_name_candidates)
        raw_column_names = re.split(r'\s{2,}', full_header_string.strip())

        # 3. Find the split point between the two datasets.
        split_index = -1
        for i, line in enumerate(data_lines):
            # Check if the line starts with '0' followed by a space or tab.
            if i > 0 and (line.lstrip().startswith('0 ') or line.lstrip().startswith('0\t')):
                split_index = i
                break
        
        if split_index == -1:
            raise ValueError("Could not find a data split point in the file.")

        data_lines_1 = data_lines[:split_index]
        data_lines_2 = data_lines[split_index:]

        # 4. Create both DataFrames.
        self.df1 = self._create_dataframe(data_lines_1, raw_column_names)
        self.df2 = self._create_dataframe(data_lines_2, raw_column_names)
        
        if not self.df1.empty:
            self.column_names = self.df1.columns.tolist()

    def _create_dataframe(self, lines: list, cols: list) -> pd.DataFrame:
        """
        Helper method to convert raw data lines into a clean DataFrame.
        """
        num_cols = len(cols)
        # Join all lines and then split to handle numbers that might be broken by newlines.
        all_values_str = " ".join(lines).split()
        
        if not all_values_str:
             return pd.DataFrame()
        
        # Ensure the total number of data points is a multiple of the number of columns.
        if len(all_values_str) % num_cols != 0:
            raise ValueError(
                f"Data length mismatch: {len(all_values_str)} values "
                f"is not divisible by {num_cols} columns."
            )
            
        data_array = np.array(all_values_str, dtype=float).reshape(-1, num_cols)
        df = pd.DataFrame(data_array, columns=cols)
        
        # Clean up the column names for easier use.
        clean_names = {}
        for col in df.columns:
            # Remove units (e.g., '(nC/m^2)').
            new_name = re.sub(r'\s*\([^)]+\)', '', col)
            # Convert to lowercase and remove leading/trailing whitespace.
            new_name = new_name.lower().strip()
            # Replace spaces with underscores.
            new_name = re.sub(r'\s+', '_', new_name) 
            clean_names[col] = new_name
        df.rename(columns=clean_names, inplace=True)
        return df

class ComsolDataReader:
    """
    Uma classe dedicada para ler e analisar todos os arquivos .txt do COMSOL
    do diretório 'Results' de um caso específico.

    Ela descobre, lê e analisa automaticamente todos os arquivos .txt, retornando-os
    em um dicionário estruturado.
    """

    def __init__(self, script_file_path: str, autoShow: bool = True):
        """
        Inicializa o leitor identificando o diretório 'Results' alvo.

        Args:
            project_root (Path): O diretório raiz do projeto pyLCP.
            case_name (str): O nome do caso de teste específico (ex: 'coated_bifilar_s40').
        """
        project_root = Path(script_file_path).resolve().parents[2]
        if str(project_root) not in sys.path:
            sys.path.insert(0, str(project_root))
            
        case_name = os.path.splitext(os.path.basename(script_file_path))[0]        
        
        self.project_root = project_root
        self.case_name = case_name
        self.results_path = project_root / 'testData' / case_name / 'Results'
        self.data = self.load_all_results()

        if autoShow and self.data:
            self.show_summary()

        print(f"Project root configured at: {project_root}")
        print(f"\nInstanciando ComsolDataReader para o caso '{case_name}' ---")
        
        if not self.results_path.is_dir():
            print(f"  Aviso: diretório Results não encontrado para o caso '{case_name}'. Dados COMSOL ignorados.")
        
    def load_all_results(self) -> dict[str, pd.DataFrame]:
        """
        Verifica o diretório 'Results', carrega todos os arquivos .txt e os retorna
        como um dicionário de DataFrames.

        Returns:
            Um dicionário onde as chaves são os nomes dos arquivos (sem a extensão .txt)
            e os valores são os DataFrames do pandas analisados.
        """
        txt_files = list(self.results_path.glob('*.txt'))
        data = {}
        
        if not txt_files:
            print(f"Aviso: Nenhum arquivo .txt encontrado em {self.results_path}")
            return {}

        print(f"Encontrado(s) {len(txt_files)} arquivo(s) .txt no diretório Results do caso '{self.case_name}'.")

        for file_path in txt_files:
            file_stem = file_path.stem
            try:
                data[file_stem] = self._parse_single_file(file_path)
            except Exception as e:
                print(f"Erro ao analisar o arquivo {file_path.name}: {e}")

        return data

    def _parse_single_file(self, file_path: Path) -> pd.DataFrame:
        """
        Método privado para ler e analisar um único arquivo .txt do COMSOL.
        Contém a lógica de análise principal.
        """
        print(f"  -> Carregando e analisando: {file_path.name}...")
        
        header_lines, data_lines = [], []
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('%'):
                    if not any(keyword in line for keyword in ['Model:', 'Version:', 'Date:', 'Table:']):
                        header_lines.append(line)
                elif line.strip():
                    data_lines.append(line.strip())

        full_header_str = ' '.join([h.replace('%', '').strip() for h in header_lines])
        parts = re.split(r'(\([^)]+\))', full_header_str)
        
        raw_names = []
        i = 0
        while i < len(parts) - 1:
            var_name = parts[i].strip()
            unit = parts[i+1].strip()
            if var_name:
                raw_names.append(f"{var_name} {unit}")
            i += 2
        num_cols = len(raw_names)

        # Make raw column names unique before DataFrame creation so that
        # duplicate exports (e.g. two excitations in the same file) don't
        # cause df[col] to return a DataFrame instead of a Series.
        seen_raw = {}
        unique_raw = []
        for name in raw_names:
            if name in seen_raw:
                seen_raw[name] += 1
                unique_raw.append(f"{name}#{seen_raw[name]}")
            else:
                seen_raw[name] = 0
                unique_raw.append(name)

        all_values_str = " ".join(data_lines).split()
        if not all_values_str: raise ValueError("Nenhum dado encontrado no arquivo.")
        if len(all_values_str) % num_cols != 0:
            raise ValueError(f"Incompatibilidade de dados: {len(all_values_str)} valores não é múltiplo de {num_cols} colunas.")

        data_array = np.array(all_values_str).reshape(-1, num_cols)
        df = pd.DataFrame(data_array, columns=unique_raw)

        for col in df.columns:
            if df[col].astype(str).str.contains('i').any():
                df[col] = df[col].astype(str).str.replace('i', 'j', regex=False).apply(complex)
            else:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        # Clean names: strip units and special chars, then re-deduplicate with _N suffix.
        seen_clean = {}
        final_names = []
        for col in df.columns:
            base = re.sub(r'#\d+$', '', col)  # remove disambiguation marker
            clean = re.sub(r'[^a-z0-9_]+', '_', re.sub(r'\s*\([^)]+\)', '', base.lower())).strip('_')
            if clean in seen_clean:
                seen_clean[clean] += 1
                final_names.append(f"{clean}_{seen_clean[clean]}")
            else:
                seen_clean[clean] = 0
                final_names.append(clean)
        df.columns = final_names

        return df

    def show_summary(self, head_rows: int = 5):
        """
        Exibe um resumo de todos os DataFrames carregados, mostrando o head e info de cada um.

        Args:
            head_rows (int): O número de linhas a serem exibidas do cabeçalho de cada DataFrame.
        """
        if not self.data:
            print("Nenhum dado carregado para exibir o resumo. Execute 'load_all_results()' primeiro.")
            return

        print(f"\n--- Resumo dos Dados Carregados para o Caso '{self.case_name}' ---")
        for name, df in self.data.items():
            print(f"\n==================================================")
            print(f"  Arquivo: '{name}.txt'")
            print(f"==================================================")
            
            print(f"\n--- Head ---")
            print(df.head(head_rows))
            
            print(f"\n--- Info ---")
            df.info()
            print("\n")

class ComsolPostProcessor:
    """
    Classe para processar dados COMSOL específicos para linhas de transmissão
    com configuração flat de 3 cabos (A, B, C) e construir matrizes de impedância
    e admitância de retorno à terra.
    """
    def __init__(self, script_file_path: str, autoShow: bool = True):
        self.cmsl_reader = ComsolDataReader(script_file_path, autoShow=autoShow)
        
    def get_general_parameters(self, cmsl_file_name: str = None) -> dict:
        """
        Retorna os parâmetros gerais extraídos dos dados COMSOL.
        """
        data = self.cmsl_reader.data.get(cmsl_file_name, None)
        if data is None:
            print(f"  Aviso: arquivo COMSOL '{cmsl_file_name}.txt' não encontrado. Dados COMSOL ignorados.")
            return None
        freq = np.asarray(data['freq'])

        return {
            'frequencies': freq,
            'angular_frequencies': 2 * np.pi * freq,
        }

    def load_scc_earth_return_and_internal_scenarios(
        self, pul_data: dict, internal_mtl_model, internal_form: str = 'approximation',
        earth_return_key: str = 'cmsl_ground_return_impedance',
    ) -> None:
        """
        Popula pul_data['comsol'] in-place com os dois blocos de dados COMSOL
        que os casos SCC (andreata_case1/2/3) sempre carregam juntos:

        1. Retorno à terra por cenário (chaves já presentes em
           pul_data['comsol']['scenarios'], ex. 'rho_g_100_epsr1_1_mf'): para
           cada uma, monta earth_return_parameters + quasi_tem_matrices a
           partir de internal_mtl_model (usado só para os parâmetros
           internos -- o retorno à terra em si vem do COMSOL).
        2. Impedância interna combinada (núcleo + blindagem), lida de
           'cmsl_internal_impedance_matrix.txt': data1() = excitação pelo
           núcleo (self do núcleo + mútua), data2() = excitação pela
           blindagem (self da blindagem) -- ver
           get_scc_internal_impedance_matrix_combined(). Não depende do
           solo/retorno à terra, por isso é carregada incondicionalmente,
           mesmo que o bloco 1 acima não tenha dados.

        NOTA: os dois blocos compartilham pul_data['comsol']['frequencies']
        -- hoje inofensivo nos três casos porque só um dos dois arquivos
        costuma existir por vez; se algum caso vier a ter os dois
        simultaneamente, vão precisar da mesma grade de frequência.
        """
        print("Construindo matrizes COMSOL...")
        cmsl_params = self.get_general_parameters(earth_return_key)
        if cmsl_params is not None:
            pul = InternalPerUnitParameters(internal_mtl_model, cmsl_params['frequencies'])
            internal_matrices = pul.matrices(internal_form=internal_form)
            pul_data['comsol'].update(cmsl_params)
            pul_data['comsol']['internal_matrices'] = internal_matrices

            for key, value in pul_data['comsol']['scenarios'].items():
                print(f"  -> Processando COMSOL para: {key}")
                earth_return = self.get_earth_return_parameters(key)
                quasi_tem = self.get_quasi_tem_approx_matrices(internal_matrices, earth_return)
                value['earth_return_parameters'] = earth_return
                value['quasi_tem_matrices'] = quasi_tem
        else:
            print("  Aviso: Processamento COMSOL ignorado (dados não disponíveis).")
            pul_data['comsol'] = {}

        print("Carregando dados COMSOL de impedância interna...")
        scc_internal_cmsl = self.get_scc_internal_impedance_matrix_combined()
        if scc_internal_cmsl is not None:
            pul_data['comsol'].setdefault('scenarios', {})
            pul_data['comsol']['frequencies'] = scc_internal_cmsl['frequencies']
            pul_data['comsol']['scenarios'].update(scc_internal_cmsl['scenarios'])
        else:
            print("  Aviso: dados COMSOL de impedância interna não disponíveis.")

    def get_bifilar_data(self, excitation_type: str = 'average'):
        N = 2
        general_data = self.get_general_parameters('cmsl_series_impedance_matrix')
        if general_data is None:
            return None
        data = self.cmsl_reader.data['cmsl_series_impedance_matrix']
        freq = general_data['frequencies']
        Zp = np.zeros((len(freq), N, N), dtype=complex)
        Zs = np.zeros((len(freq), N-1, N-1), dtype=complex)

        Zp[:, 0, 0] = data['v11']
        Zp[:, 1, 1] = data['v22'] 

        if excitation_type == 'conductor_1':
            Zp[:, 0, 1] = data['v12']
        elif excitation_type == 'conductor_2':
            Zp[:, 0, 1] = data['v21']
        elif excitation_type == 'average':
            Zp[:, 0, 1] = 0.5 * (data['v21'] + data['v12'])

        Zp[:, 1, 0] = Zp[:, 0, 1]

        Zs[:, 0, 0] = data['single_coil_voltage']

        return {
            'partial_impedance_matrix': Zp,
            'series_impedance_matrix': Zs,
        }        

    def get_coaxial_series_impedance_matrix(self):
        """
        Constrói a matriz de impedância [n, 3, 3] simétrica para uma configuração
        flat de 3 cabos (A, B, C), a partir dos dados do COMSOL.

        A montagem assume uma configuração simétrica:
        - Z_AA = Z_BB = Z_CC (self, de 'vcoil_1')
        - Z_AB = Z_BA = Z_BC = Z_CB (mutual adjacente, de 'vcoil_2')
        - Z_AC = Z_CA (mutual externa, de 'vcoil_3')

        :param base_key: A chave base do cenário COMSOL 
                        (ex: 'rho_g_100_epsr1_1_mf').
        :param model: O modelo de linha de transmissão multiconductor.
        :return: Uma tupla (Z0, freq), onde Z0 é a matriz [n, 3, 3] e 
                freq é o vetor de frequências [n]. Retorna (None, None) se 
                os dados não forem encontrados.
        """
        N = 2
        general_data = self.get_general_parameters('cmsl_coaxial_cable_impedance')
        data = self.cmsl_reader.data['cmsl_coaxial_cable_impedance']
        freq = general_data['frequencies']
        Zs = np.zeros((len(freq), N-1, N-1), dtype=complex)

        # Zs: series impedance matrix
        Zs[:, 0, 0] = data['coil_impedance']

        return Zs
    
    def get_coaxial_cable_parameters(self):
        """
        Constrói a matriz de impedância [n, 3, 3] simétrica para uma configuração
        flat de 3 cabos (A, B, C), a partir dos dados do COMSOL.

        A montagem assume uma configuração simétrica:
        - Z_AA = Z_BB = Z_CC (self, de 'vcoil_1')
        - Z_AB = Z_BA = Z_BC = Z_CB (mutual adjacente, de 'vcoil_2')
        - Z_AC = Z_CA (mutual externa, de 'vcoil_3')

        :param base_key: A chave base do cenário COMSOL 
                        (ex: 'rho_g_100_epsr1_1_mf').
        :param model: O modelo de linha de transmissão multiconductor.
        :return: Uma tupla (Z0, freq), onde Z0 é a matriz [n, 3, 3] e 
                freq é o vetor de frequências [n]. Retorna (None, None) se 
                os dados não forem encontrados.
        """
        general_data = self.get_general_parameters('cmsl_coaxial_cable_impedance')
        data = self.cmsl_reader.data['cmsl_coaxial_cable_impedance']
        jw = 1j * general_data['angular_frequencies']

        # z11: internal impedance of core outer surface
        z11 = data['r11'] + jw * data['l11']

        # z12: core outer insulator impedance
        z12 = data['r12'] + jw * data['l12']

        # z2i: internal impedance of sheath inner surface
        z2i = data['r2i'] + jw * data['l2i']

        return {
            'Zcs': data['coil_impedance'],
            'z11': z11,
            'z12': z12,
            'z2i': z2i,
        }
    
    def get_scc_internal_impedance_elements(self, excitation_type: str = 'core'):
        """
        Constrói a matriz de impedância [n, 3, 3] simétrica para uma configuração
        flat de 3 cabos (A, B, C), a partir dos dados do COMSOL.

        A montagem assume uma configuração simétrica:
        - Z_AA = Z_BB = Z_CC (self, de 'vcoil_1')
        - Z_AB = Z_BA = Z_BC = Z_CB (mutual adjacente, de 'vcoil_2')
        - Z_AC = Z_CA (mutual externa, de 'vcoil_3')

        :param base_key: A chave base do cenário COMSOL 
                        (ex: 'rho_g_100_epsr1_1_mf').
        :param model: O modelo de linha de transmissão multiconductor.
        :return: Uma tupla (Z0, freq), onde Z0 é a matriz [n, 3, 3] e 
                freq é o vetor de frequências [n]. Retorna (None, None) se 
                os dados não forem encontrados.
        """
        N = 2
        general_data = self.get_general_parameters('cmsl_series_impedance_core_excitation')
        if general_data is None:
            return None
        if 'cmsl_series_impedance_sheath_excitation' not in self.cmsl_reader.data:
            print("  Aviso: arquivo COMSOL 'cmsl_series_impedance_sheath_excitation.txt' não encontrado. Dados COMSOL ignorados.")
            return None
        core = self.cmsl_reader.data['cmsl_series_impedance_core_excitation']
        sheath = self.cmsl_reader.data['cmsl_series_impedance_sheath_excitation']
        freq = general_data['frequencies']
        jw = 1j * 2 * np.pi * freq

        # z11: internal impedance of core outer surface
        z11_cr = core['r11'] + jw * core['l11']
        z11_sh = sheath['r2i'] + jw * sheath['l2i']

        # z12: core insulator impedance
        z12_cr = core['r12'] + jw * core['l12']
        z12_sh = sheath['r12'] + jw * sheath['l12']

        # z2i: internal impedance of sheath inner surface
        z2i_cr = core['r2i'] + jw * core['l2i']
        z2i_sh = sheath['r11'] + jw * sheath['l11']

        # z13: sheath insulator impedance
        z13_cr, z13_sh = 0, 0
        if 'l13' in core and 'l13' in sheath:
            z13_cr = jw * core['l13']
            z13_sh = jw * sheath['l13']

        # z14: air gap impedance
        z14_cr, z14_sh = 0, 0
        if 'l14' in core and 'l14' in sheath:   
            z14_cr = jw * core['l14']
            z14_sh = jw * sheath['l14']

        # z15: HDPE tube impedance
        z15_cr, z15_sh = 0, 0
        if 'l15' in core and 'l15' in sheath:
            z15_cr = jw * core['l15']
            z15_sh = jw * sheath['l15']

        Zi_energy = np.zeros((len(freq), N, N), dtype=complex)
        Zi_js = np.zeros((len(freq), N, N), dtype=complex)
        Z11 = z11_cr + z12_cr + z2i_cr + z13_cr + z14_cr + z15_cr
        Z22 = z11_sh + z12_sh + z2i_sh + z13_sh + z14_sh + z15_sh
        
        Zi_energy[:, 0, 0] = Z11
        Zi_energy[:, 1, 1] = Z22

        Zi_js[:, 0, 0] = core['core_voltage']
        Zi_js[:, 1, 1] = sheath['sheath_voltage']

        if excitation_type == 'core':
            mutual_js = core['sheath_voltage']
            mutual_energy = 0.5 * z2i_cr
            Zi_js[:, 0, 1] = core['sheath_voltage']
            Zi_js[:, 1, 0] = core['sheath_voltage']
            Zi_energy[:, 0, 1] = mutual_energy
            Zi_energy[:, 1, 0] = mutual_energy
        
        elif excitation_type == 'sheath':
            mutual_js = sheath['core_voltage']
            mutual_energy = z2i_sh
            Zi_js[:, 0, 1] = sheath['core_voltage']
            Zi_js[:, 1, 0] = sheath['core_voltage']
            Zi_energy[:, 0, 1] = mutual_energy
            Zi_energy[:, 1, 0] = mutual_energy
        
        else:
            raise ValueError("excitation_type deve ser 'core' ou 'sheath'.")
        
        return {
            'js_method': Zi_js,
            'energy_method': Zi_energy,
            'self_core_js': core['core_voltage'],
            'self_sheath_js': sheath['sheath_voltage'],
            'mutual_js': mutual_js,
            'self_core_energy': Z11,
            'self_sheath_energy': Z22,
            'mutual_energy': mutual_energy
        }

    def get_internal_impedance_elements(self) -> dict:
        """
        Retorna a impedância interna medida (r11 + jwL11) do arquivo
        'cmsl_internal_impedance.txt' para um condutor sólido único.
        """
        general_data = self.get_general_parameters('cmsl_internal_impedance')
        if general_data is None:
            return None
        data = self.cmsl_reader.data['cmsl_internal_impedance']
        jw = 1j * general_data['angular_frequencies']

        return {
            'frequencies': general_data['frequencies'],
            'r11': data['r11'].to_numpy(),
            'l11': data['l11'].to_numpy(),
            'Zi_measured': data['r11'].to_numpy() + jw * data['l11'].to_numpy(),
        }

    def get_bare_and_hollow_wire_internal_impedance(self) -> dict:
        """
        Retorna a impedância interna medida via método Js (tensão da bobina de
        excitação sob corrente unitária) do arquivo
        'cmsl_bare_and_hollow_wire_internal_impedance.txt', para os condutores
        sólido (bare wire) e oco (hollow/tubular) de single_deConti.
        """
        general_data = self.get_general_parameters('cmsl_bare_and_hollow_wire_internal_impedance')
        if general_data is None:
            return None
        data = self.cmsl_reader.data['cmsl_bare_and_hollow_wire_internal_impedance']

        return {
            'frequencies': general_data['frequencies'],
            'Zi_measured': data['solid_conductor_coil_voltage'].to_numpy(),
            'Zi_hollow': data['hollow_conductor_coil_voltage'].to_numpy(),
        }

    def get_scc_internal_impedance_matrix(self) -> dict:
        """
        Returns the SCC internal impedance matrix [Zi] (core + sheath, 2x2)
        measured via the Js method (coil voltage under core/sheath excitation),
        from 'internal_impedance_matrix_core_excitation.txt' and
        'internal_impedance_matrix_sheath_excitation.txt'.

        Returned in the standard 'scenarios' shape expected by
        BasePlotter._get_data_from_source(source='comsol'), under a single
        synthetic scenario key 'measured' (there is only one measured matrix,
        not per-analytical-scenario data).
        """
        N = 2
        general_data = self.get_general_parameters('internal_impedance_matrix_core_excitation')
        if general_data is None:
            return None
        if 'internal_impedance_matrix_sheath_excitation' not in self.cmsl_reader.data:
            print("  Warning: file 'internal_impedance_matrix_sheath_excitation.txt' not found. COMSOL data ignored.")
            return None

        core = self.cmsl_reader.data['internal_impedance_matrix_core_excitation']
        sheath = self.cmsl_reader.data['internal_impedance_matrix_sheath_excitation']
        freq = general_data['frequencies']

        Zi = np.zeros((len(freq), N, N), dtype=complex)
        Zi[:, 0, 0] = core['core_voltage']
        Zi[:, 1, 1] = sheath['sheath_voltage']
        Zi[:, 0, 1] = core['sheath_voltage']
        Zi[:, 1, 0] = core['sheath_voltage']

        return {
            'frequencies': freq,
            'scenarios': {
                'measured': {'impedance_matrix': Zi},
            },
        }

    def get_scc_internal_impedance_matrix_combined(self) -> dict:
        """
        Returns the SCC internal impedance matrix [Zi] (core + sheath, 2x2)
        measured via the Js method, from a single combined file
        'cmsl_internal_impedance_matrix.txt' holding both excitations in one
        table: data1(...) = core excitation (mf.VCoil_core_i0: core coil
        voltage, self-term; mf.VCoil_sheath_0: sheath coil voltage, mutual
        term, sheath held at zero current) and data2(...) = sheath excitation
        (mf.VCoil_core_0: core coil voltage, mutual term, core held at zero
        current; mf.VCoil_sheath_i0: sheath coil voltage, self-term).

        The four column names each contain their own inner parenthesis (e.g.
        'data1(mf.VCoil_core_i0)'), which the generic parser's unit-stripping
        regex collapses down to positional placeholders 'data1'/'data1_1'/
        'data2'/'data2_1' — read here by column order (as they appear in the
        file), not by semantic name. The mutual term is taken from the core-
        excitation reading ('data1_1'), matching the convention used by
        get_scc_internal_impedance_matrix() for the two-file format.

        Returned in the standard 'scenarios' shape expected by
        BasePlotter._get_data_from_source(source='comsol').
        """
        N = 2
        general_data = self.get_general_parameters('cmsl_internal_impedance_matrix')
        if general_data is None:
            return None

        data = self.cmsl_reader.data['cmsl_internal_impedance_matrix']
        freq = general_data['frequencies']

        Zi = np.zeros((len(freq), N, N), dtype=complex)
        Zi[:, 0, 0] = data['data1']    # data1(mf.VCoil_core_i0): core excitation, core voltage (self)
        Zi[:, 0, 1] = data['data1_1']  # data1(mf.VCoil_sheath_0): core excitation, sheath voltage (mutual)
        Zi[:, 1, 0] = data['data1_1']
        Zi[:, 1, 1] = data['data2_1']  # data2(mf.VCoil_sheath_i0): sheath excitation, sheath voltage (self)

        return {
            'frequencies': freq,
            'scenarios': {
                'measured': {'impedance_matrix': Zi},
            },
        }

    def get_shunt_capacitance_elements(self) -> dict:
        """
        Returns C_11 (core) and C_22 (sheath) self-capacitances per unit length
        from the COMSOL shunt parameters file, computed via the energy method.
        Both values are frequency-independent; the first row is used.
        """
        df = self.cmsl_reader.data['cmsl_shunt_params']
        return {
            'c11': float(df['ccc_energy'].iloc[0]),
            'c22': float(df['css_energy'].iloc[0]),
        }

    def get_earth_return_parameters(self, base_key: str):
        """
        Constrói a matriz de impedância [n, 3, 3] simétrica para uma configuração
        flat de 3 cabos (A, B, C), a partir dos dados do COMSOL.

        A montagem assume uma configuração simétrica:
        - Z_AA = Z_BB = Z_CC (self, de 'vcoil_1')
        - Z_AB = Z_BA = Z_BC = Z_CB (mutual adjacente, de 'vcoil_2')
        - Z_AC = Z_CA (mutual externa, de 'vcoil_3')

        :param base_key: A chave base do cenário COMSOL 
                        (ex: 'rho_g_100_epsr1_1_mf').
        :param model: O modelo de linha de transmissão multiconductor.
        :return: Uma tupla (Z0, freq), onde Z0 é a matriz [n, 3, 3] e 
                freq é o vetor de frequências [n]. Retorna (None, None) se 
                os dados não forem encontrados.
        """
        N = 3
        general_data = self.get_general_parameters('cmsl_ground_return_impedance')
        data = self.cmsl_reader.data['cmsl_ground_return_impedance']

        freq = general_data['frequencies']
        jw = 1j * general_data['angular_frequencies']

        Zg = np.zeros((len(freq), N, N), dtype=complex)
        Yg = np.zeros_like(Zg, dtype=complex)
        Pg = np.zeros_like(Zg, dtype=complex)

        # 3. Parsear a base_key para extrair rho_g e eps_r
        match = re.search(r"rho_g_(\d+)_epsr1_(\d+)_mf", base_key)
        
        if not match:
            print(f"Erro: Não foi possível extrair os parâmetros (rho_g, eps_r) da chave '{base_key}'.")
            return None, None
            
        sigma1 = 1.0 / float(match.group(1))    # Condutividade do solo (S/m)
        epsr1 = float(match.group(2))           # Permissividade Relativa do solo

        # Squared Ground propagation constant 
        gamma_earth = np.sqrt(jw * sc.mu_0 * (sigma1 + jw * epsr1 * sc.epsilon_0))

        # 3. Definir as chaves de dados com base no mapeamento fornecido
        key_self = f"{base_key}_vcoil_1" # Z_AA, Z_BB, Z_CC
        key_adj  = f"{base_key}_vcoil_2" # Z_AB, Z_BC
        key_ext  = f"{base_key}_vcoil_3" # Z_AC
        required_keys = [key_self, key_adj, key_ext]
        
        # 4. Verificar se todas as chaves de dados necessárias existem
        if not all(key in data for key in required_keys):
            print(f"Erro: Faltando uma ou mais chaves para a base_key '{base_key}' nos dados COMSOL.")
            print(f"Chaves necessárias: {required_keys}")
            print(f"Chaves disponíveis: {list(data.keys())}")
            return None, None

        # Diagonal (Self-impedances)
        Zg[:, 0, 0] = data[key_self]  # Z_AA
        Zg[:, 1, 1] = data[key_self]  # Z_BB
        Zg[:, 2, 2] = data[key_self]  # Z_CC

        # Termos adjacentes (A-B e B-C)
        Zg[:, 0, 1] = data[key_adj]   # Z_AB
        Zg[:, 1, 0] = data[key_adj]   # Z_BA
        Zg[:, 1, 2] = data[key_adj]   # Z_BC
        Zg[:, 2, 1] = data[key_adj]   # Z_CB
        
        # Termos externos (A-C)
        Zg[:, 0, 2] = data[key_ext]   # Z_AC
        Zg[:, 2, 0] = data[key_ext]   # Z_CA

        # Earth-Return Admittance based on Vance (1978) formulation
        for i in range(len(freq)):
            gamma_earth2 = gamma_earth[i]**2 * np.eye(N)
            Yg[i, :, :] = gamma_earth2 @ np.linalg.inv(Zg[i, :, :])
            Pg[i, :, :] = jw[i] * np.linalg.inv(Yg[i, :, :])

        return {
            'impedance_matrix': Zg,
            'potential_coefficient': Pg,
            'admittance_matrix': Yg,
            'gamma_earth': gamma_earth,
        }
    
    def get_quasi_tem_approx_matrices(self, internal_matrices, earth_return_params):
        """
        Assembles the final PUL matrices for a vector of frequencies.
        """
        general_data = self.get_general_parameters('cmsl_ground_return_impedance')
        freq = general_data['frequencies']
        jw = 1j * general_data['angular_frequencies']

        # Conductors per cable, in the same order as the ground-return (N x N)
        # matrices below. Uniform ([M]*N) for homogeneous MTL types; for
        # 'scc-flat-ecc' the SCC and ECC cables have different M, so the
        # ground-return coupling (which only depends on cable position, not on
        # which/how many conductors that cable has) must be tiled block-by-block
        # with each cable's own size rather than a single uniform M. Mirrors
        # PerUnitParameters.quasi_tem_approx_matrices (analytical_forms/single_core_cable.py).
        block_sizes = internal_matrices['block_sizes']

        z0_jk = earth_return_params['impedance_matrix']         # Shape (num_freq, N, N)
        pg_jk = earth_return_params['potential_coefficient']    # Shape (num_freq, N, N)
        yg_jk = earth_return_params['admittance_matrix']        # Shape (num_freq, N, N)

        Zi = internal_matrices['impedance_matrix']              # Shape (num_freq, sum(M), sum(M))
        Pi = internal_matrices['potential_coefficient_matrix']  # Shape (sum(M), sum(M))
        Yi = internal_matrices['shunt_admittance_matrix']       # Shape (num_freq, sum(M), sum(M))

        Zg = _expand_by_block_sizes(z0_jk, block_sizes)
        Pg = _expand_by_block_sizes(pg_jk, block_sizes)
        Yg = np.zeros_like(Zi, dtype=complex)

        # Series impedance is a simple element-wise addition
        Zs = Zi + Zg

        # Shunt Admittance Matrix calculation
        # Pi is 2D, Pe is 3D. Use broadcasting to add them.
        Psh = Pi[np.newaxis, :, :] + Pg
        
        # The linear solve must be looped over the frequency axis
        Ysh = np.zeros_like(Psh, dtype=complex)
        for i in range(len(freq)):
            Ysh[i, :, :] = jw[i] * np.linalg.inv(Psh[i, :, :])

        return {
            'earth_return_impedance_matrix': Zg,
            'earth_return_potential_coefficient': Pg,
            'earth_return_admittance_matrix': Yg,
            'potential_coefficient': Psh,
            'series_impedance_matrix': Zs,
            'shunt_admittance_matrix': Ysh,
        }