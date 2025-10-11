# root_folder/model/case_utils.py
""" It is a Python script that contains the main function and a class called User. """

import os
import numpy as np
import pandas as pd
import scipy.special as ss
from typing import Optional

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
