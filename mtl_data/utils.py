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

# Set the print options to use the custom formatter
def set_numpy_print_options():
    """ Set the numpy print options to use the custom formatter."""
    # Define a custom formatter for the numpy print options
    my_formatter = {
        'float_kind': lambda x: f"{x:.2e}",
        'complex_kind': complex_formatter
    }
    np.set_printoptions(formatter=my_formatter)


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
