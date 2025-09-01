""" It is a Python script that contains the main function and a class called User. """

import numpy as np
import pandas as pd
from typing import Optional

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
