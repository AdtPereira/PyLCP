import os
import sys
import time
import copy
import pandas as pd
from tabulate import tabulate
from pathlib import Path
import matplotlib.pyplot as plt

# RAIZ DO PROJETO E DIRETÓRIOS
try:
    os.system('cls' if os.name == 'nt' else 'clear')
    script_dir = Path(__file__).resolve().parent
    print(f"Script directory: {script_dir}")
    project_root = script_dir.parents[1]
    print(f"Project root: {project_root}")
    sys.path.append(str(project_root))
    print("Caminhos do projeto configurados com sucesso.")
except IndexError:
    raise FileNotFoundError(
        "Não foi possível encontrar a raiz do projeto. "
        "Certifique-se de que o script está em 'examples/coated_wires'."
    )

# IMPORTAÇÕES DOS MÓDULOS E MODELO DE DADOS
try:
    from mtl_main.models_wires import SINGLE_WIRE_R468 as MODEL
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from mtl_main.utils import *
    from analytical_forms.overhead_lines import PerUnitParameters
    print("Módulos e modelo de dados importados com sucesso.")
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

def log_pul_parameters(pul_parameters):
    """
    Generates and prints a terminal log report of the calculated PUL parameters
    using a pandas DataFrame.

    Args:
        pul_parameters (dict): The dictionary containing the calculated parameters,
                               with frequencies as keys.
    """
    report_data = []
    dim_scale = 1e3  # Scale factor for m to km

    # Iterate through the frequencies and their corresponding data
    for freq, params in sorted(pul_parameters.items()):
        row = {
            'Frequency (Hz)': f'{freq:,.0f}' # Format frequency with commas
        }
        
        param_mapping = {
            'Series Imp. (Ohm/km)': 'zs',
            'Shunt Adm. (S/km)': 'ysh',
            'Earth-Return (Ohm/km)': 'zg'
        }

        for name_in_report, key_in_dict in param_mapping.items():
            row[name_in_report] = params[key_in_dict][0, 0] * dim_scale            
        report_data.append(row)

    df = pd.DataFrame(report_data)
    
    # Define a lambda function to format complex numbers into strings
    complex_formatter = lambda c: f"{c.real:.3e} {'-' if c.imag < 0 else '+'} {abs(c.imag):.3e}j"
    complex_cols = df.select_dtypes(include='complex128').columns
    df[complex_cols] = df[complex_cols].map(complex_formatter)
    # 'psql' is a clean, readable table format. Others include 'grid', 'fancy_grid', etc.
    table = tabulate(df, headers='keys', tablefmt='psql', stralign="center")

    # --- Print the Report ---
    print("Per-Unit-Length (PUL) Parameters Report")
    print(table)

if __name__ == "__main__":
    """ Main function to orchestrate the analysis, calculation, and visualization. """
    st = time.time()

    frequency = {
        'Analytically': [60, 100e3, 1e6],
        'Numerically': [60, 100e3, 1e6],
    }

    pul_parameters = {}
    mtl = MulticonductorTransmissionLine(MODEL)
    for f in frequency['Analytically']:
        pul_parameters[f] = PerUnitParameters(mtl, f).pul_quasi_tem(form='deri')

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")
    log_pul_parameters(pul_parameters)
    MTLRepresentation(mtl, units='millimeter').ground_return_systems()
    plt.show()