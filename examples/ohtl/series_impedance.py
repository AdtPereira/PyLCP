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
    from mtl_data.wire_models import SINGLE_OHTL_XUE as MTL
    from mtl_data.graphics import MTLRepresentation
    from ohtl.pul_parameters import PerUnitParameters
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)


def plot_series_impedance(freq, pul_dict, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """
    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    # fig.suptitle('')

    f = np.array(freq['Analytically'])
    z_a = np.array([item['Zs'][p, q] for item in pul_dict['a']])
    z_b = np.array([item['Zs'][p, q] for item in pul_dict['b']])
    z_c = np.array([item['Zs'][p, q] for item in pul_dict['c']])
    l_a = np.imag(z_a) / (2 * np.pi * f)  # H/m
    l_b = np.imag(z_b) / (2 * np.pi * f)  # H/m
    l_c = np.imag(z_c) / (2 * np.pi * f)  # H/m

    # Extracting the impedance elements from zi_matrix
    ax1.plot(f, 1E3 * np.real(z_a), label=r'$\rho_e = 100 \;\Omega m, \epsilon_r=1$', color='black', linestyle='-')
    ax1.plot(f, 1E3 * np.real(z_b), label=r'$\rho_e = 100 \;\Omega m, \epsilon_r=20$', color='black', linestyle='--')
    ax1.plot(f, 1E3 * np.real(z_c), label=r'$\rho_e = 2000 \;\Omega m, \epsilon_r=1$', color='black', linestyle='-.')

    # Additional plotting configurations
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(1E3, 1E9)
    ax1.set_ylim(1E0, 1E5)
    ax1.legend()
    ax1.set_xlabel('Frequency (Hz)')
    ax1.set_ylabel(r'$R_s \, (\Omega/km)$')
    ax1.grid(False)
    ax1.set_title('P.u.l. series resistance of the single overhead line\n'
                 r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m}, \rho = 1.68 \times 10^{-8} \, \mathrm{\Omega m}$ [1]')

    # Extracting the impedance elements from zi_matrix
    ax2.plot(f, 1E6 * l_a, label=r'$\rho_e = 100 \;\Omega m, \epsilon_r=1$', color='black', linestyle='-')
    ax2.plot(f, 1E6 * l_b, label=r'$\rho_e = 100 \;\Omega m, \epsilon_r=20$', color='black', linestyle='--')
    ax2.plot(f, 1E6 * l_c, label=r'$\rho_e = 2000 \;\Omega m, \epsilon_r=1$', color='black', linestyle='-.')

    ax2.set_xscale('log')
    ax2.set_xlim(1E3, 1E9)
    ax2.set_ylim(1, 2.5)
    ax2.legend()
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel(r'$L_s \, (mH/km)$')
    ax2.grid(False)
    ax2.set_title('P.u.l. series inductance of the single overhead line\n'
                 r'$r_1 = 0.01 \, \mathrm{m}, h_1 = 10 \, \mathrm{m}, \rho = 1.68 \times 10^{-8} \, \mathrm{\Omega m}$ [1]')
    plt.tight_layout(rect=[0, 0, 1, 0.96])


if __name__ == "__main__":
    """ Main function to perform the calculations and display the results."""
    st = time.time()

    # Calculate the series impedance for each frequency
    frequency = {'Analytically': np.logspace(3, 9, num=200), 'Numerically': np.logspace(0, 7, num=30)}

    # Dictionary to hold the series impedance calculations
    pul_dict = {'a': [], 'b': [], 'c': []}

    # Analytical Formulation
    for f in frequency['Analytically']:
        pul_a = PerUnitParameters(MTL, f, sigma_1=1/1E2, er_1=1)
        pul_b = PerUnitParameters(MTL, f, sigma_1=1/1E2, er_1=20)
        pul_c = PerUnitParameters(MTL, f, sigma_1=1/2E3, er_1=1)

        # Internal and external impedance matrices
        pul_dict['a'].append(pul_a.pul_extended_theory())
        pul_dict['b'].append(pul_b.pul_extended_theory())
        pul_dict['c'].append(pul_c.pul_extended_theory())

    # End the timer
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")
    plot_series_impedance(frequency, pul_dict, p=0, q=0)
    MTLRepresentation(MTL, units='meter').bared_and_coated_wires()
    plt.show()
    