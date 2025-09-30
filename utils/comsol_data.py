import re
import numpy as np
import pandas as pd
from pathlib import Path

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

    def __init__(self, project_root: Path, case_name: str):
        """
        Inicializa o leitor identificando o diretório 'Results' alvo.

        Args:
            project_root (Path): O diretório raiz do projeto pyLCP.
            case_name (str): O nome do caso de teste específico (ex: 'coated_bifilar_s40').
        """
        self.project_root = project_root
        self.case_name = case_name
        self.results_path = project_root / 'testData' / case_name / 'Results'
        self.data = {}

        if not self.results_path.is_dir():
            raise FileNotFoundError(
                f"O diretório Results não foi encontrado para o caso '{case_name}' em: {self.results_path}"
            )
        
    def load_all_results(self) -> dict[str, pd.DataFrame]:
        """
        Verifica o diretório 'Results', carrega todos os arquivos .txt e os retorna
        como um dicionário de DataFrames.

        Returns:
            Um dicionário onde as chaves são os nomes dos arquivos (sem a extensão .txt)
            e os valores são os DataFrames do pandas analisados.
        """
        txt_files = list(self.results_path.glob('*.txt'))
        if not txt_files:
            print(f"Aviso: Nenhum arquivo .txt encontrado em {self.results_path}")
            return {}

        print(f"Encontrado(s) {len(txt_files)} arquivo(s) .txt no diretório Results do caso '{self.case_name}'.")

        for file_path in txt_files:
            file_stem = file_path.stem
            try:
                self.data[file_stem] = self._parse_single_file(file_path)
            except Exception as e:
                print(f"Erro ao analisar o arquivo {file_path.name}: {e}")

        return self.data

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
        
        column_names = []
        i = 0
        while i < len(parts) - 1:
            var_name = parts[i].strip()
            unit = parts[i+1].strip()
            if var_name:
                column_names.append(f"{var_name} {unit}")
            i += 2
        num_cols = len(column_names)

        all_values_str = " ".join(data_lines).split()
        if not all_values_str: raise ValueError("Nenhum dado encontrado no arquivo.")
        if len(all_values_str) % num_cols != 0:
            raise ValueError(f"Incompatibilidade de dados: {len(all_values_str)} valores não é múltiplo de {num_cols} colunas.")
        
        data_array = np.array(all_values_str).reshape(-1, num_cols)
        df = pd.DataFrame(data_array, columns=column_names)

        for col in df.columns:
            if df[col].astype(str).str.contains('i').any():
                df[col] = df[col].str.replace('i', 'j', regex=False).apply(complex)
            else:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        clean_names = {col: re.sub(r'[^a-z0-9_]+', '_', re.sub(r'\s*\([^)]+\)', '', col.lower())).strip('_') for col in df.columns}
        df.rename(columns=clean_names, inplace=True)
        
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
