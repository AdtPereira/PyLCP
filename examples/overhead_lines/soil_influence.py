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
import copy
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
    from mtl_main.models_wires import SINGLE_OHTL_XUE as MODEL
    from mtl_main.source import MulticonductorTransmissionLine
    from mtl_main.graphics import MTLRepresentation
    from analytical_formulation.overhead_lines import PerUnitParameters
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)


def plot_gamma(freq, pul, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """

    plt.figure(figsize=(8, 5))
    f = np.array(freq['Analytically'])
    w = 2 * np.pi * f

    gmma_1 = np.array([item['γ'][p, q] for item in pul['rho_1u']['nakagawa']])
    gmma_2 = np.array([item['γ'][p, q] for item in pul['rho_200']['nakagawa']])
    gmma_3 = np.array([item['γ'][p, q] for item in pul['rho_5k']['nakagawa']])
    gmma_4 = np.array([item['γ'][p, q] for item in pul['rho_1u']['quasi_tem']])
    gmma_5 = np.array([item['γ'][p, q] for item in pul['rho_200']['quasi_tem']])
    gmma_6 = np.array([item['γ'][p, q] for item in pul['rho_5k']['quasi_tem']])

    # Plot the results
    plt.plot(freq['Analytically'], np.real(1E3 * gmma_4), color='grey', linestyle='-') # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(1E3 * gmma_5), color='grey', linestyle='-') # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(1E3 * gmma_6), color='grey', linestyle='-') # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(1E3 * gmma_1), label=r'$\rho_g=1\;\mu\Omega m$', color='blue', linestyle='--') # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(1E3 * gmma_2), label=r'$\rho_g=200\;\Omega m$', color='green', linestyle='--') # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(1E3 * gmma_3), label=r'$\rho_g=5000\;\Omega m$', color='red', linestyle='--') # pylint: disable=line-too-long

    # Optional: Additional plotting configurations like labels, grid, etc.
    plt.xscale('log')
    plt.xlim(1E0, 1E10)
    plt.ylim(0, 1.6)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'Attenuation Constant, $\alpha$ (Np/km)')
    plt.title('Attenuation Constant of the single overhead line\n'
              r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m}, \epsilon_{rg}=5$ [3]')
    plt.grid(False)

    # Plot the results
    plt.figure(figsize=(8, 5))
    plt.plot(f, w / np.imag(1E6 * gmma_1), label=r'$\rho_g=1\;\mu\Omega m$', color='black', linestyle='--')
    plt.plot(f, w / np.imag(1E6 * gmma_2), label=r'$\rho_g=200\;\Omega m$', color='black', linestyle='-')
    plt.plot(f, w / np.imag(1E6 * gmma_3), label=r'$\rho_g=5000\;\Omega m$', color='red', linestyle='-')

    # Optional: Additional plotting configurations like labels, grid, etc.
    plt.xscale('log')
    plt.xlim(1E0, 1E10)
    plt.ylim(50, 350)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'Phase Velocity\; $(m/\mu s)$')
    plt.title('Phase Velocity of the single overhead line (Log. Approx. Quasi-TEM)\n'
              r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m}, \epsilon_{rg}=5$ [3]')
    plt.grid(False)


if __name__ == "__main__":
    """ Main function to perform the calculations and display the results."""
    st = time.time()

    # Calculate the series impedance for each frequency
    frequency = {
        'Analytically': np.logspace(0, 10, num=400),
        'Numerically': np.logspace(0, 7, num=30)
    }

    # Dictionary to hold the series impedance calculations
    pul_parameters = {
        'rho_1u': {'nakagawa': [], 'quasi_tem': []},
        'rho_200': {'nakagawa': [], 'quasi_tem': []},
        'rho_5k': {'nakagawa': [], 'quasi_tem': []},
    }

    # The geometric model is constant, so we create the object once for efficiency.
    model_a = copy.deepcopy(MODEL)
    model_a[0]['relative_permittivity'] = 5
    model_a[0]['conductivity'] = 1E6
    mtl_model_a = MulticonductorTransmissionLine(model_a)

    model_b = copy.deepcopy(MODEL)
    model_b[0]['relative_permittivity'] = 5
    model_b[0]['conductivity'] = 1/2E2
    mtl_model_b = MulticonductorTransmissionLine(model_b)

    model_c = copy.deepcopy(MODEL)
    model_c[0]['relative_permittivity'] = 5
    model_c[0]['conductivity'] = 1/5E3
    mtl_model_c = MulticonductorTransmissionLine(model_c)

    # Analytical Formulation
    for f in frequency['Analytically']:
        pul1 = PerUnitParameters(mtl_model_a, f)
        pul2 = PerUnitParameters(mtl_model_b, f)
        pul3 = PerUnitParameters(mtl_model_c, f)

        pul_parameters['rho_1u']['nakagawa'].append(pul1.pul_extended_theory(type_form='nakagawa'))
        pul_parameters['rho_1u']['quasi_tem'].append(pul1.pul_extended_theory())
        pul_parameters['rho_200']['nakagawa'].append(pul2.pul_extended_theory(type_form='nakagawa'))
        pul_parameters['rho_200']['quasi_tem'].append(pul2.pul_extended_theory())
        pul_parameters['rho_5k']['nakagawa'].append(pul3.pul_extended_theory(type_form='nakagawa'))
        pul_parameters['rho_5k']['quasi_tem'].append(pul3.pul_extended_theory())

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")
    plot_gamma(frequency, pul_parameters, p=0, q=0)
    MTLRepresentation(mtl_model_a, units='meter').ground_return_systems()
    plt.show()
    