# root_folder/model/case_utils.py
""" It is a Python script that contains the main function and a class called User. """

import re
import os
import json
import numpy as np
import pandas as pd
import scipy.special as ss
from pathlib import Path
from typing import Optional, Tuple

UNITS_DATA = {
    'meter':      {'scale': 1,       'label': 'm'},
    'centimeter': {'scale': 100,     'label': 'cm'},
    'millimeter': {'scale': 1000,    'label': 'mm'},
    'mil':        {'scale': 39370.1, 'label': 'mil'},
}

def complex_formatter(x):
    """ Format complex numbers in scientific notation."""
    if x == 0:
        return " 0 "
    else:
        return f"{x.real:.3e}{x.imag:+.3e}j"

def set_numpy_print_options():
    """ Set the numpy print options to use the custom formatter."""
    # Define a custom formatter for the numpy print options
    my_formatter = {
        'float_kind': lambda x: f"{x:.2e}",
        'complex_kind': complex_formatter
    }
    np.set_printoptions(formatter=my_formatter)

def matrix_to_string(matrix: np.ndarray) -> str:
    """
    Formats a NumPy matrix into a multi-line string for display.
    Example:
     [ 1.2345e+01  -2.3456e-02]
     [ 3.4567e-03   4.5678e+04]
    """
    lines = []
    num_rows = matrix.shape[0]
    # Format each number with fixed width for alignment
    s_rows = [[f"{val:11.4e}" for val in row] for row in matrix]
    
    for i, row in enumerate(s_rows):
        s_row = "  ".join(row)
        # Add brackets to the first and last lines
        # prefix = " [" if i == 0 else "  "
        # suffix = "]" if i == num_rows - 1 else ""
        # lines.append(f"{prefix}{s_row}{suffix}")
        lines.append(s_row)        
    return "\n".join(lines)

def format_scientific_notation(value, precision=1):
    """
    Formats a number into scientific notation using LaTeX style (e.g., 5.8 x 10^7).
    """
    e_notation = f'{value:.{precision}e}'
    base, exponent = e_notation.split('e')
    exponent_int = int(exponent)
    return fr'{base} \times 10^{{{exponent_int}}}'

def matrix_viewer(matrix: np.ndarray, title: str, columns_name: Optional[list] = None) -> pd.DataFrame:
    """
    Formata e exibe uma matriz NumPy como um DataFrame do pandas com um título.

    Args:
        matriz_numpy (np.ndarray): A matriz de entrada do NumPy.
        titulo (str): O título a ser exibido acima da matriz formatada.
        colunas (Optional[list]): Uma lista opcional de nomes para as colunas.

    Returns:
        pd.DataFrame: A matriz formatada como um DataFrame do pandas.
    """
    pd.options.display.float_format = '{:.5e}'.format
    df = pd.DataFrame(matrix, columns=columns_name)
    
    print("\n")
    print(f"{title}:")
    print("-" * (len(title) + 1))
    if columns_name is None:
        print(df.to_string(header=False, index=False))
    else:
        print(df)

def verify_kelvin_functions(q = 1.5):
    """
    This script verifies the output of scipy.special.kelvin() by comparing it
    against the fundamental mathematical definitions that relate Kelvin functions
    to the modified Bessel functions I_0 and I_1 with a complex argument.
    
    This version correctly unpacks the complex tuple returned by ss.kelvin()
    as per the official SciPy documentation:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kelvin.html
    """
    
    print(f"--- Verification for q = {q} ---")
    
    # --- 1. Call ss.kelvin() and correctly extract real values ---
    print("\n[Reference] Unpacking values from ss.kelvin(q) according to SciPy docs:")
    
    # ss.kelvin() returns a tuple of 4 complex numbers: (Be, Ke, Bep, Kep)
    Be_complex, Ke_complex, Bep_complex, Kep_complex = ss.kelvin(q)
    
    # The real-valued functions are the real/imaginary parts of the complex results.
    # ber and bei come from the first element (Be).
    ber_ref = Be_complex.real
    bei_ref = Be_complex.imag
    # ber' and bei' come from the third element (Bep).
    ber_p_ref = Bep_complex.real
    bei_p_ref = Bep_complex.imag
    
    print(f"ber(q)  = Re[kelvin(q)[0]] = {ber_ref:>10.6f}")
    print(f"bei(q)  = Im[kelvin(q)[0]] = {bei_ref:>10.6f}")
    print(f"ber'(q) = Re[kelvin(q)[2]] = {ber_p_ref:>10.6f}")
    print(f"bei'(q) = Im[kelvin(q)[2]] = {bei_p_ref:>10.6f}")
    
    # --- 2. Verification using the Bessel function I_0 ---
    # According to the definition: I_0(q*sqrt(j)) = ber(q) + j*bei(q)
    complex_arg = q * (1j**0.5)
    i0_complex = ss.iv(0, complex_arg)
    
    ber_from_bessel = i0_complex.real
    bei_from_bessel = i0_complex.imag
    
    print("\n[Test 1] Values derived from I_0(q * sqrt(j)):")
    print(f"Re[I_0] = {ber_from_bessel:>10.6f} -> Matches ber(q)? {np.isclose(ber_ref, ber_from_bessel)}")
    print(f"Im[I_0] = {bei_from_bessel:>10.6f} -> Matches bei(q)? {np.isclose(bei_ref, bei_from_bessel)}")

    # --- 3. Verification of the derivatives using Bessel function I_1 ---
    # According to the definition: sqrt(j)*I_1(q*sqrt(j)) = ber'(q) + j*bei'(q)
    i1_complex_term = (1j**0.5) * ss.iv(1, complex_arg)

    ber_p_from_bessel = i1_complex_term.real
    bei_p_from_bessel = i1_complex_term.imag

    print("\n[Test 2] Values derived from sqrt(j) * I_1(q * sqrt(j)):")
    print(f"Re[...] = {ber_p_from_bessel:>10.6f} -> Matches ber'(q)? {np.isclose(ber_p_ref, ber_p_from_bessel)}")
    print(f"Im[...] = {bei_p_from_bessel:>10.6f} -> Matches bei'(q)? {np.isclose(bei_p_ref, bei_p_from_bessel)}")

def save_figure_multiformat(fig, results_dir, base_filename, formats=['png', 'pdf']):
    """
    Saves a figure to multiple file formats in the case's results directory.

    Args:
        fig (matplotlib.figure.Figure): The figure object to save.
        base_filename (str): The base name for the output file, without extension.
        formats (list, optional): A list of file extensions to save as. 
                                    Defaults to ['png', 'svg'].
    """
    # Loop through each format and save the figure
    for fmt in formats:
        # Create the full filename with the current format's extension
        full_filename = f"{base_filename}.{fmt}"
        
        # Create the full path to the results directory
        file_path = os.path.join(results_dir, full_filename)
        
        # Prepare keyword arguments for savefig to handle format-specific options
        save_kwargs = {'bbox_inches': 'tight'}
        if fmt == 'png':
            save_kwargs['dpi'] = 300
        
        # Save the figure using the specific options
        fig.savefig(file_path, **save_kwargs)
        
        print(f"Plot saved to {file_path}")

def load_json_parameters(script_file_path, show_content=False):
    """
    Dynamically loads parameters from a '.json' file.
    It assumes the '.in.json' file has the same base name as the
    calling script and is located in the same directory.

    Args:
        script_file_path (str): The __file__ attribute from the calling script.
        show_content (bool): If True, prints the content of the loaded
                             dictionary to the console. Defaults to False.

    Returns:
        dict: A dictionary with the parameters loaded from the JSON file.
    """
    script_path = Path(script_file_path)
    script_dir = script_path.resolve().parent
    case_name = script_path.stem
    model_file_name = f"{case_name}.json"
    model_path = script_dir / model_file_name
    
    if not model_path.exists():
        raise FileNotFoundError(f"The parameter file could not be found at: {model_path}")

    with open(model_path, 'r') as f:
        parameters = json.load(f)
    
    print(f"Successfully loaded parameters from: {model_path}")
    
    # --- NOVO BLOCO DE CÓDIGO ---
    # Se o parâmetro show_content for True, imprime o dicionário
    if show_content:
        print(f"\n--- Content of {model_file_name} ---")
        # Usa json.dumps para uma impressão "pretty print"
        print(json.dumps(parameters, indent=2))
        print("--------------------------------------\n")
        
    return parameters

def load_comsol_results(script_path: str, comsol_tag: str = '') -> pd.DataFrame:
    """
    Loads COMSOL data by deriving the case name and filename from the script path.

    This acts as a convenient wrapper for `read_comsol_results`, making the
    main script calls cleaner and more consistent with `load_json_parameters`.

    Args:
        script_path (str): The path of the calling script (typically __file__).
        comsol_tag (str): The suffix to append to the case name to form the
                          .txt filename. Defaults to 'cmsl_1'.

    Returns:
        pd.DataFrame: A pandas DataFrame with the COMSOL simulation data.
    """
    case_name = Path(script_path).stem
    file_name = f"{case_name}{comsol_tag}.txt"
    print(f"Loading COMSOL data {file_name}")
    return read_comsol_results(case_name, file_name)

def read_comsol_results(case_name: str, file_name: str, base_dir: str = "testData") -> pd.DataFrame:
    """
    Reads data from a COMSOL-exported .txt file, dynamically parsing
    column names and their units from the header. This robust version
    handles single or multi-line headers and automatically cleans column names.
    """
    file_path = Path(base_dir) / case_name / "Results" / file_name

    if not file_path.exists():
        raise FileNotFoundError(f"The specified file was not found: {file_path}")

    header_lines = []
    data_lines = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('%'):
                # Skip metadata lines
                if not any(keyword in line for keyword in ['Model:', 'Version:', 'Date:', 'Table:']):
                    header_lines.append(line)
            elif line.strip():
                data_lines.append(line.strip())

    # --- 1. Robust Header Parsing with a More General Regex ---
    full_header_str = ' '.join([h.replace('%', '').strip() for h in header_lines])
    
    # This regex is more flexible and finds any text within parentheses,
    # such as (Hz), (H/m), or (Ω).
    parts = re.split(r'(\([^)]+\))', full_header_str)
    
    column_names = []
    i = 0
    # Group the variable name with its unit
    # Example: ['freq ', '(Hz)', ' r11 ', '(Ω/m)'] -> ['freq (Hz)', 'r11 (Ω/m)']
    while i < len(parts) - 1:
        var_name = parts[i].strip()
        unit = parts[i+1].strip()
        if var_name:
            column_names.append(f"{var_name} {unit}")
        i += 2
            
    num_cols = len(column_names)

    # --- 2. Data Processing and Reshaping ---
    all_values_str = " ".join(data_lines).split()
    
    if len(all_values_str) == 0:
        raise ValueError("No data found in the file.")
    
    if len(all_values_str) % num_cols != 0:
        raise ValueError(
            f"Data mismatch: Total values ({len(all_values_str)}) "
            f"is not a multiple of the number of columns ({num_cols})."
        )
    
    data_array = np.array(all_values_str).reshape(-1, num_cols)

    # --- 3. DataFrame Creation and Generalized Data Type Conversion ---
    df = pd.DataFrame(data_array, columns=column_names)

    for col in df.columns:
        # If any value in the column contains 'i', treat it as complex.
        # This is more robust than checking units.
        if df[col].astype(str).str.contains('i').any():
            df[col] = df[col].str.replace('i', 'j', regex=False).apply(complex)
        else:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # --- 4. Extensible Column Name Cleaning ---
    clean_names = {}
    for col in df.columns:
        new_name = col.lower()
        # Remove units in parentheses and clean up the name
        new_name = re.sub(r'\s*\([^)]+\)', '', new_name) # Remove "(unit)"
        new_name = re.sub(r'[^a-z0-9_]+', '_', new_name) # Replace special chars with underscore
        new_name = new_name.strip('_') # Clean leading/trailing underscores
        clean_names[col] = new_name
    df.rename(columns=clean_names, inplace=True)
    
    return df

# def read_and_split_comsol_results(file_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
#     """
#     Reads data from a COMSOL .txt file containing two concatenated datasets
#     and splits them into two separate DataFrames. The split point is detected
#     when the first column's value (e.g., arc length) resets to zero.

#     Args:
#         file_path (str): The full path to the .txt data file.

#     Returns:
#         Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two pandas
#                                            DataFrames, one for each dataset.
#     """
#     path_obj = Path(file_path)
#     if not path_obj.exists():
#         raise FileNotFoundError(f"The specified file was not found: {path_obj}")

#     header_lines = []
#     data_lines = []
#     with open(path_obj, 'r', encoding='utf-8') as f:
#         for line in f:
#             if line.startswith('%'):
#                 if not any(keyword in line for keyword in ['Model:', 'Version:', 'Date:', 'Table:']):
#                     header_lines.append(line)
#             elif line.strip():
#                 data_lines.append(line.strip())

#     # --- 1. Parse Header (similar to the original function) ---
#     full_header_str = ' '.join([h.replace('%', '').strip() for h in header_lines])
#     parts = re.split(r'(\([^)]+\))', full_header_str)
    
#     column_names = []
#     i = 0
#     while i < len(parts) - 1:
#         var_name = parts[i].strip()
#         unit = parts[i+1].strip()
#         if var_name:
#             column_names.append(f"{var_name} {unit}")
#         i += 2
    
#     if not column_names:
#         raise ValueError("Could not parse column headers from the file.")

#     # --- 2. Find Split Point and Separate Data ---
#     split_index = -1
#     for i, line in enumerate(data_lines):
#         # Find the second occurrence of the start of a dataset (arc length = 0)
#         if i > 0 and line.lstrip().startswith('0 '):
#             split_index = i
#             break
            
#     if split_index == -1:
#         raise ValueError("Could not find a split point in the data file.")

#     data_lines_1 = data_lines[:split_index]
#     data_lines_2 = data_lines[split_index:]

#     # --- 3. Create a helper function to process each data subset ---
#     def create_dataframe(lines, cols):
#         num_cols = len(cols)
#         all_values_str = " ".join(lines).split()
        
#         if len(all_values_str) % num_cols != 0:
#             raise ValueError("Data mismatch during processing.")
            
#         data_array = np.array(all_values_str, dtype=float).reshape(-1, num_cols)
#         df = pd.DataFrame(data_array, columns=cols)
        
#         # Clean column names
#         clean_names = {}
#         for col in df.columns:
#             new_name = col.lower()
#             new_name = re.sub(r'\s*\([^)]+\)', '', new_name)
#             new_name = re.sub(r'[^a-z0-9_]+', '_', new_name)
#             new_name = new_name.strip('_')
#             clean_names[col] = new_name
#         df.rename(columns=clean_names, inplace=True)
#         return df

#     # --- 4. Create both DataFrames ---
#     df1 = create_dataframe(data_lines_1, column_names)
#     df2 = create_dataframe(data_lines_2, column_names)
    
#     return df1, df2

# def load_split_comsol_results(script_path: str, comsol_tag: str = '') -> pd.DataFrame:
#     """
#     Loads COMSOL data by deriving the case name and filename from the script path.

#     This acts as a convenient wrapper for `read_comsol_results`, making the
#     main script calls cleaner and more consistent with `load_json_parameters`.

#     Args:
#         script_path (str): The path of the calling script (typically __file__).
#         comsol_tag (str): The suffix to append to the case name to form the
#                           .txt filename. Defaults to 'cmsl_1'.

#     Returns:
#         pd.DataFrame: A pandas DataFrame with the COMSOL simulation data.
#     """
#     case_name = Path(script_path).stem
#     file_name = f"{case_name}{comsol_tag}.txt"
#     print(f"Loading COMSOL data {file_name}")
#     return read_split_comsol_results(case_name, file_name)

# def read_split_comsol_results(case_name: str, file_name: str, base_dir: str = "testData") -> pd.DataFrame:
#     """
#     Reads data from a COMSOL-exported .txt file, dynamically parsing
#     column names and their units from the header. This robust version
#     handles single or multi-line headers and automatically cleans column names.
#     """
#     file_path = Path(base_dir) / case_name / "Results" / file_name

#     if not file_path.exists():
#         raise FileNotFoundError(f"The specified file was not found: {file_path}")

#     header_lines = []
#     data_lines = []
#     with open(file_path, 'r', encoding='utf-8') as f:
#         for line in f:
#             line_content = line.strip()
#             if line_content.startswith('%'):
#                 # Ignore metadata comments, but keep header lines
#                 if not any(keyword in line for keyword in ['Model:', 'Version:', 'Date:', 'Table:']):
#                     # Check if the line seems to be a real header and not just '%'
#                     if len(line_content) > 1:
#                         header_lines.append(line_content)
#             elif line_content:
#                 data_lines.append(line_content)

#     # --- BLOCO CORRIGIDO ---
#     # Process each header line individually instead of joining them.
#     # This correctly handles headers spread across multiple lines.
#     column_names = [h.replace('%', '').strip() for h in header_lines]
#     if not column_names:
#         raise ValueError("Could not parse column headers from the file.")
#     num_cols = len(column_names)
    
#     # --- 2. Data Processing and Reshaping ---
#     all_values_str = " ".join(data_lines).split()
    
#     if len(all_values_str) == 0:
#         raise ValueError("No data found in the file.")
    
#     if len(all_values_str) % num_cols != 0:
#         raise ValueError(
#             f"Data mismatch: Total values ({len(all_values_str)}) "
#             f"is not a multiple of the number of columns ({num_cols})."
#         )
    
#     # Use float type for initial conversion
#     data_array = np.array(all_values_str, dtype=float).reshape(-1, num_cols)

#     # --- 3. DataFrame Creation ---
#     df = pd.DataFrame(data_array, columns=column_names)

#     # --- 4. Extensible Column Name Cleaning ---
#     clean_names = {}
#     for col in df.columns:
#         new_name = col.lower()
#         new_name = re.sub(r'\s*\([^)]+\)', '', new_name) # Remove "(unit)"
#         new_name = re.sub(r'[^a-z0-9_]+', '_', new_name) # Replace special chars with underscore
#         new_name = new_name.strip('_') # Clean leading/trailing underscores
#         clean_names[col] = new_name
#     df.rename(columns=clean_names, inplace=True)
    
#     return df

# def read_and_split_comsol_results(file_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
#     """
#     Reads data from a COMSOL .txt file containing two concatenated datasets
#     and splits them into two separate DataFrames. The split point is detected
#     when the first column's value (e.g., arc length) resets to zero.
#     """
#     path_obj = Path(file_path)
#     if not path_obj.exists():
#         raise FileNotFoundError(f"The specified file was not found: {path_obj}")

#     header_lines = []
#     data_lines = []
#     with open(path_obj, 'r', encoding='utf-8') as f:
#         for line in f:
#             line_content = line.strip()
#             if line_content.startswith('%'):
#                 if not any(keyword in line for keyword in ['Model:', 'Version:', 'Date:', 'Table:']):
#                     if len(line_content) > 1:
#                         header_lines.append(line_content)
#             elif line_content:
#                 data_lines.append(line_content)

#     # --- BLOCO CORRIGIDO ---
#     # Process each header line individually.
#     column_names = [h.replace('%', '').strip() for h in header_lines]
#     if not column_names:
#         raise ValueError("Could not parse column headers from the file.")

#     # --- 2. Find Split Point and Separate Data ---
#     split_index = -1
#     for i, line in enumerate(data_lines):
#         # Find the second occurrence of a dataset start (arc length = 0)
#         # We check for a line that starts with '0' and is not the very first line.
#         if i > 0 and line.lstrip().startswith('0 '):
#             split_index = i
#             break
            
#     if split_index == -1:
#         raise ValueError("Could not find a split point in the data file.")

#     data_lines_1 = data_lines[:split_index]
#     data_lines_2 = data_lines[split_index:]

#     # --- 3. Helper function to create each DataFrame ---
#     def create_dataframe(lines, cols):
#         num_cols = len(cols)
#         all_values_str = " ".join(lines).split()
        
#         if not all_values_str:
#              return pd.DataFrame(columns=cols) # Return empty DataFrame if no data
        
#         if len(all_values_str) % num_cols != 0:
#             raise ValueError("Data mismatch during subset processing.")
            
#         data_array = np.array(all_values_str, dtype=float).reshape(-1, num_cols)
#         df = pd.DataFrame(data_array, columns=cols)
        
#         # Clean column names
#         clean_names = {}
#         for col in df.columns:
#             new_name = col.lower()
#             new_name = re.sub(r'\s*\([^)]+\)', '', new_name)
#             new_name = re.sub(r'[^a-z0-9_]+', '_', new_name)
#             new_name = new_name.strip('_')
#             clean_names[col] = new_name
#         df.rename(columns=clean_names, inplace=True)
#         return df

#     # --- 4. Create both DataFrames ---
#     df1 = create_dataframe(data_lines_1, column_names)
#     df2 = create_dataframe(data_lines_2, column_names)
    
#     return df1, df2

# def load_split_comsol_results(script_path: str, comsol_tag: str = ''):
#     """
#     Convenience wrapper for read_and_split_comsol_results.
#     """
#     case_name = Path(script_path).stem
#     file_name = f"{case_name}{comsol_tag}.txt"
#     file_path = Path("testData") / case_name / "Results" / file_name
#     print(f"Loading and splitting COMSOL data from {file_name}")
#     return read_and_split_comsol_results(str(file_path))

# def read_and_split_comsol_results(file_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
#     """
#     Reads a COMSOL .txt file with two concatenated datasets and splits them
#     into two separate DataFrames.
#     """
#     path_obj = Path(file_path)
#     if not path_obj.exists():
#         raise FileNotFoundError(f"The specified file was not found: {path_obj}")

#     header_lines = []
#     data_lines = []
#     with open(path_obj, 'r', encoding='utf-8') as f:
#         for line in f:
#             line_content = line.strip()
#             if line_content.startswith('%'):
#                 # --- CORRECTED LOGIC ---
#                 # A true header line is assumed to contain units in parentheses.
#                 if '(' in line_content and ')' in line_content:
#                     header_lines.append(line_content)
#             elif line_content:
#                 data_lines.append(line_content)

#     column_names = [h.replace('%', '').strip() for h in header_lines]
#     if not column_names:
#         raise ValueError("Could not parse valid column headers from the file.")

#     split_index = -1
#     for i, line in enumerate(data_lines):
#         if i > 0 and line.lstrip().startswith('0 '):
#             split_index = i
#             break
            
#     if split_index == -1:
#         raise ValueError("Could not find a split point in the data file.")

#     data_lines_1 = data_lines[:split_index]
#     data_lines_2 = data_lines[split_index:]

#     def create_dataframe(lines, cols):
#         num_cols = len(cols)
#         all_values_str = " ".join(lines).split()
#         if not all_values_str: return pd.DataFrame(columns=cols)
#         if len(all_values_str) % num_cols != 0: raise ValueError("Data mismatch during subset processing.")
        
#         data_array = np.array(all_values_str, dtype=float).reshape(-1, num_cols)
#         df = pd.DataFrame(data_array, columns=cols)
        
#         clean_names = {}
#         for col in df.columns:
#             new_name = col.lower().replace('(','').replace(')','').replace(' ', '_')
#             new_name = re.sub(r'[^a-z0-9_]+', '_', new_name).strip('_')
#             clean_names[col] = new_name
#         df.rename(columns=clean_names, inplace=True)
#         return df

#     df1 = create_dataframe(data_lines_1, column_names)
#     df2 = create_dataframe(data_lines_2, column_names)
    
#     return df1, df2

# def read_and_split_comsol_results(file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
#     """
#     Reads data from a COMSOL .txt file containing two concatenated datasets
#     and splits them into two separate DataFrames. This version uses robust
#     header detection and column name cleaning.
#     """
#     path_obj = Path(file_path)
#     if not path_obj.exists():
#         raise FileNotFoundError(f"The specified file was not found: {path_obj}")

#     header_lines_raw = []
#     data_lines = []
#     with open(path_obj, 'r', encoding='utf-8') as f:
#         for line in f:
#             stripped_line = line.strip()
#             if not stripped_line:
#                 continue
#             if stripped_line.startswith('%'):
#                 header_lines_raw.append(stripped_line)
#             else:
#                 data_lines.append(stripped_line)

#     # --- 1. Robust Header Identification ---
#     # Filter out common COMSOL metadata to isolate only column name lines.
#     metadata_keywords = ['Model:', 'Version:', 'Date:', 'Table:', 'Dimension:', 'Nodes:', 'Expressions:', 'Description:']
#     column_names = []
#     for line in header_lines_raw:
#         clean_line = line.replace('%', '').strip()
#         if clean_line and not any(keyword in line for keyword in metadata_keywords):
#             column_names.append(clean_line)

#     if not column_names:
#         raise ValueError("Could not parse valid column headers from the file.")

#     # --- 2. Find Split Point and Separate Data ---
#     split_index = -1
#     for i, line in enumerate(data_lines):
#         # Find the second occurrence of a dataset start (arc length = 0)
#         if i > 0 and (line.lstrip().startswith('0 ') or line.lstrip().startswith('0\t')):
#             split_index = i
#             break
            
#     if split_index == -1:
#         # If no split is found, assume one dataset and raise a warning or handle as needed.
#         # For this specific problem, we expect a split.
#         raise ValueError("Could not find a split point in the data file.")

#     data_lines_1 = data_lines[:split_index]
#     data_lines_2 = data_lines[split_index:]

#     # --- 3. Helper function to create each DataFrame ---
#     def create_dataframe(lines, cols):
#         num_cols = len(cols)
#         # Join lines and then split by any whitespace
#         all_values_str = " ".join(lines).split()
        
#         if not all_values_str:
#              return pd.DataFrame(columns=cols)
        
#         if len(all_values_str) % num_cols != 0:
#             raise ValueError(f"Data mismatch: {len(all_values_str)} values for {num_cols} columns.")
            
#         data_array = np.array(all_values_str, dtype=float).reshape(-1, num_cols)
#         df = pd.DataFrame(data_array, columns=cols)
        
#         # --- 4. Robust Column Name Cleaning ---
#         clean_names = {}
#         for col in df.columns:
#             # First, remove unit in parentheses, if it exists
#             new_name = re.sub(r'\s*\([^)]+\)', '', col)
#             # Then, convert to lower case, strip whitespace, and replace internal whitespace with a single underscore
#             new_name = new_name.lower().strip()
#             new_name = re.sub(r'\s+', '_', new_name)
#             clean_names[col] = new_name
#         df.rename(columns=clean_names, inplace=True)
#         return df

#     # --- 5. Create both DataFrames ---
#     df1 = create_dataframe(data_lines_1, column_names)
#     df2 = create_dataframe(data_lines_2, column_names)
    
#     return df1, df2
