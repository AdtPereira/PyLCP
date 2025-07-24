"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

REFERENCES:
[1] 

"""

import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from ohtl.pul_parameters import PerUnitParameters
from data.models import MTL_MODELS
from mom_so.mtl_graphics import MTLRepresentation as graph

# Multiconductor Transmission Line choices
MTL = MTL_MODELS['overhead']['xue']
MTL = MTL_MODELS['overhead']['deConti']


def plot_gamma(freq, pul_dict, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """

    quasitem = [1E3 * item['γ'][p, q] for item in pul_dict['quasitem']] # pylint: disable=line-too-long # Np/km
    quasitem_log = [1E3 * item['γ'][p, q] for item in pul_dict['quasitem_log']] # pylint: disable=line-too-long # Np/km
    nakagawa = [1E3 * item['γ'][p, q] for item in pul_dict['nakagawa']] # pylint: disable=line-too-long # Np/km
    sunde = [1E3 * item['γ'][p, q] for item in pul_dict['sunde']] # pylint: disable=line-too-long # Np/km
    carson = [1E3 * item['γ'][p, q] for item in pul_dict['carson']] # pylint: disable=line-too-long # Np/km

    # Plot the results
    plt.plot(freq['Analytically'], np.real(quasitem), label='Quasi-TEM (Integral Eq.)',
             color='black', linestyle='-')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(quasitem_log), label='Quasi-TEM (Approx. Log.)',
             color='blue', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(nakagawa), label='Wise (1948) - Nakagawa (1981)',
             color='red', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(sunde), label='Sunde (1968)',
             color='green', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(carson), label='Carson (1926)',
             color='black', linestyle=':')  # pylint: disable=line-too-long

    # Optional: Additional plotting configurations like labels, grid, etc.
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
    plt.show()


def main():
    """ Main function to perform the calculations and display the results."""

    # Clears the console screen and starts the timer
    os.system('cls' if os.name == 'nt' else 'clear')
    print("Calculations were started ... ...")
    start_time = time.time()

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

    # Analytical Formulation
    for f in frequency['Analytically']:
        pul = PerUnitParameters(MTL, f, sigma_1=1/200, er_1=5)

        # Internal and external impedance matrices
        pul_parameters['quasitem'].append(pul.pul_extended_theory()) # pylint: disable=line-too-long
        pul_parameters['quasitem_log'].append(pul.pul_extended_theory(type_form='quasitem_log')) # pylint: disable=line-too-long
        pul_parameters['nakagawa'].append(pul.pul_extended_theory(type_form='nakagawa')) # pylint: disable=line-too-long
        pul_parameters['sunde'].append(pul.pul_extended_theory(type_form='sunde')) # pylint: disable=line-too-long
        pul_parameters['carson'].append(pul.pul_extended_theory(type_form='carson')) # pylint: disable=line-too-long

    # End the timer
    elapsed_time = time.time() - start_time
    print(f"End of the routine! Time spent on simulation: {elapsed_time:.2f} seconds.\n")  # pylint: disable=line-too-long

    # Display the geometry of the transmission line
    graph(MTL).wires_and_cables()

    # Plot the series resistance as a function of frequency
    plot_gamma(frequency, pul_parameters, p=0, q=0)


if __name__ == "__main__":
    main()
