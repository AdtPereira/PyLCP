"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

REFERENCES:
[1] 

"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
from ohtl.pul_parameters import PerUnitParameters
from data.systems import MULTICONDUCTOR_TRANSMISSION_LINE
from data.multiconductor import GraphicRepresentation as graph

# Multiconductor Transmission Line choices
# 0 - A single overhead conductor, h1 = 10 m and r1 = 1 cm
MTL = MULTICONDUCTOR_TRANSMISSION_LINE[0]


def plot_gamma(freq, pul, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """

    # gamma = []
    # for z, y in zip(Z, Y):
    #     if z.shape != y.shape:
    #         raise ValueError(f"The dimensions of the matrices do not match: {z.shape} and {y.shape}") # pylint: disable=line-too-long
    #     gamma.append(np.sqrt(z @ y))

    # gamma = [np.sqrt(z @ y) if z.shape == y.shape else None for z, y in zip(Z, Y)]

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
    plt.show()

    # Plot the results
    plt.plot(f, w / np.imag(1E6 * gmma_1),
             label=r'$\rho_g=1\;\mu\Omega m$', color='black', linestyle='--')
    plt.plot(f, w / np.imag(1E6 * gmma_2),
             label=r'$\rho_g=200\;\Omega m$', color='black', linestyle='-')
    plt.plot(f, w / np.imag(1E6 * gmma_3),
             label=r'$\rho_g=5000\;\Omega m$', color='red', linestyle='-')

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
    plt.show()


def main():
    """ Main function to perform the calculations and display the results."""

    # Clears the console screen and starts the timer
    os.system('cls' if os.name == 'nt' else 'clear')
    print("Calculations were started ... ...")
    start_time = time.time()

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

    # Analytical Formulation
    for f in frequency['Analytically']:
        pul_1 = PerUnitParameters(MTL, f, sigma_e=1/1E-6, er_e=5)
        pul_2 = PerUnitParameters(MTL, f, sigma_e=1/2E+2, er_e=5)
        pul_3 = PerUnitParameters(MTL, f, sigma_e=1/5E+3, er_e=5)

        # Internal and external impedance matrices
        pul_parameters['rho_1u']['nakagawa'].append(pul_1.pul_extended_theory(type_form='nakagawa')) # pylint: disable=line-too-long
        pul_parameters['rho_200']['nakagawa'].append(pul_2.pul_extended_theory(type_form='nakagawa')) # pylint: disable=line-too-long
        pul_parameters['rho_5k']['nakagawa'].append(pul_3.pul_extended_theory(type_form='nakagawa')) # pylint: disable=line-too-long
        pul_parameters['rho_1u']['quasi_tem'].append(pul_1.pul_extended_theory()) # pylint: disable=line-too-long
        pul_parameters['rho_200']['quasi_tem'].append(pul_2.pul_extended_theory()) # pylint: disable=line-too-long
        pul_parameters['rho_5k']['quasi_tem'].append(pul_3.pul_extended_theory()) # pylint: disable=line-too-long

    # End the timer
    elapsed_time = time.time() - start_time
    print(f"End of the routine! Time spent on simulation: {
          elapsed_time:.2f} seconds.\n")

    # Display the geometry of the transmission line
    graph(MTL, sigma_1=1/1E-6, er_1=5).wires_and_cables(line_type='overhead')

    # Plot the series resistance as a function of frequency
    plot_gamma(frequency, pul_parameters, p=0, q=0)


if __name__ == "__main__":
    main()
