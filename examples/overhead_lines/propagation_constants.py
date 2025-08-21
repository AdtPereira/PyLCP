""" 
Este script executa simulações de impedância de linhas de transmissão coaxiais
usando tanto uma abordagem analítica (formulação de Patel) quanto uma abordagem numérica
(Método dos Momentos - MoM-SO). Ele gera gráficos comparativos dos resultados
obtidos por ambas as metodologias.

Este arquivo é parte do projeto PyLCP, que é um pacote Python para análise de linhas de transmissão.

REFERENCES:
[1] PATEL, Utkarsh R. A Surface Admittance Approach For Fast Calculation of the 
    Series Impedance of Cables Including Skin, Proximity, and Ground Return Effects.
    2014. University of Toronto, Graduate Department of The Edward S. Rogers Sr. 
    Department of Electrical & Computer Engineering. 

[2] U. R. Patel, B. Gustavsen and P. Triverio, "An Equivalent Surface Current Approach
    for the Computation of the Series Impedance of Power Cables with Inclusion of Skin
    and Proximity Effects," in IEEE Transactions on Power Delivery, vol. 28, no. 4, pp.
    2474-2482, Oct. 2013, doi: 10.1109/TPWRD.2013.2267098.

[3] U. R. Patel, B. Gustavsen and P. Triverio, "Application of the MoM-SO Method for 
    Accurate Impedance Calculation of Single-Core Cables Enclosed by a Conducting Pipe," 
    Proc. International Conference on Power Systems Transients (IPST 2013), Vancouver, 
    Canada July 18-20, 2013. https://www.ipstconf.org/Proc_IPST2013.php

[4] A. Ametani, "A General Formulation of Impedance and Admittance of Cables," in IEEE
    Transactions on Power Apparatus and Systems, vol. PAS-99, no. 3, pp. 902-910, May
    1980, doi: 10.1109/TPAS.1980.319718.

[5] A. Ametani, "Wave Propagation Characteristics of Cables," in IEEE Transactions on
    Power Apparatus and Systems, vol. PAS-99, no. 2, pp. 499-505, March 1980, 
    doi: 10.1109/TPAS.1980.319685.
"""
import os
import sys
import time
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# RAIZ DO PROJETO E DIRETÓRIOS
os.system('cls' if os.name == 'nt' else 'clear')
try:
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
    from mtl_main.models_wires import SINGLE_OHTL_CONTI as MODEL
    from mtl_main.source import MulticonductorTransmissionLine
    from mtl_main.graphics import MTLRepresentation
    from analytical_formulation.overhead_lines import PerUnitParameters
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)


def plot_gamma(freq, pul_dict, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """
    plt.figure(figsize=(8, 5))

    quasitem = [1E3 * item['γ'][p, q] for item in pul_dict['quasitem']]
    quasitem_log = [1E3 * item['γ'][p, q] for item in pul_dict['quasitem_log']]
    nakagawa = [1E3 * item['γ'][p, q] for item in pul_dict['nakagawa']]
    sunde = [1E3 * item['γ'][p, q] for item in pul_dict['sunde']]
    carson = [1E3 * item['γ'][p, q] for item in pul_dict['carson']]

    # Plot the results
    plt.plot(freq['Analytically'], np.real(quasitem),       color='black',  linestyle='-',  label='Quasi-TEM (Integral Eq.)')
    plt.plot(freq['Analytically'], np.real(quasitem_log),   color='blue',   linestyle='--', label='Quasi-TEM (Approx. Log.)')
    plt.plot(freq['Analytically'], np.real(nakagawa),       color='red',    linestyle='--', label='Wise (1948) - Nakagawa (1981)')
    plt.plot(freq['Analytically'], np.real(sunde),          color='green',  linestyle='--', label='Sunde (1968)', )
    plt.plot(freq['Analytically'], np.real(carson),         color='black',  linestyle=':',  label='Carson (1926)')

    plt.xscale('log')
    plt.xlim(1E3, 1E10)
    plt.ylim(0, 4)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'Attenuation Constant, $\alpha$ (Np/km)')
    plt.title('Attenuation Constant of the single overhead line\n'
              r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},'
              r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
    plt.grid(True)


if __name__ == "__main__":
    """ Main function to perform the calculations and display the results."""
    st = time.time()

    # Calculate the series impedance for each frequency
    frequency = {
        'Analytically': np.logspace(0, 10, num=200),
        'Numerically': np.logspace(0, 7, num=30)
    }

    # Dictionary to hold the series impedance calculations
    pul_parameters = {
        'quasitem': [],
        'quasitem_log': [],
        'nakagawa': [],
        'sunde': [],
        'carson': [],
    }

    # The geometric model is constant, so we create the object once for efficiency.
    mtl_model = MulticonductorTransmissionLine(MODEL)

    # Analytical Formulation
    for f in frequency['Analytically']:
        pul = PerUnitParameters(mtl_model, f, sigma_1=1/200, er_1=5)

        # Internal and external impedance matrices
        pul_parameters['quasitem'].append(pul.pul_extended_theory())
        pul_parameters['quasitem_log'].append(pul.pul_extended_theory(type_form='quasitem_log'))
        pul_parameters['nakagawa'].append(pul.pul_extended_theory(type_form='nakagawa'))
        pul_parameters['sunde'].append(pul.pul_extended_theory(type_form='sunde'))
        pul_parameters['carson'].append(pul.pul_extended_theory(type_form='carson'))

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")
    plot_gamma(frequency, pul_parameters, p=0, q=0)
    MTLRepresentation(mtl_model, units='meter').isolated_wires()
    plt.show()
