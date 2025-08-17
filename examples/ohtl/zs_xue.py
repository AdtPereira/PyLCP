import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from ohtl.pul_parameters import PerUnitParameters
from mtl_data.models import MTL_MODELS
from mtl_data.graphics import MTLRepresentation as graph
#from mom_so import green, patel

# Multiconductor Transmission Line choices
MTL = MTL_MODELS['overhead']['xue']


def plot_impedance(freq, pul_dict, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """

    # Frequency array
    f = np.array(freq['Analytically'])

    # Extracting the impedance elements from zi_matrix
    z_a = np.array([item['Zs'][p, q] for item in pul_dict['a']])
    z_b = np.array([item['Zs'][p, q] for item in pul_dict['b']])
    z_c = np.array([item['Zs'][p, q] for item in pul_dict['c']])
    l_a = np.imag(z_a) / (2 * np.pi * f)  # H/m
    l_b = np.imag(z_b) / (2 * np.pi * f)  # H/m
    l_c = np.imag(z_c) / (2 * np.pi * f)  # H/m

    plt.plot(f, 1E3 * np.real(z_a), label=r'$\rho_e = 100 \;\Omega m, \epsilon_r=1$', color='black', linestyle='-')
    plt.plot(f, 1E3 * np.real(z_b), label=r'$\rho_e = 100 \;\Omega m, \epsilon_r=20$', color='black', linestyle='--')
    plt.plot(f, 1E3 * np.real(z_c), label=r'$\rho_e = 2000 \;\Omega m, \epsilon_r=1$', color='black', linestyle='-.')

    # Additional plotting configurations
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1E3, 1E9)
    plt.ylim(1E0, 1E5)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'$R_s \, (\Omega/km)$')
    plt.grid(True)

    # Adjust the layout of the plots
    plt.suptitle('P.u.l. series resistance of the single overhead line\n'
                 r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m}, \rho = 1.68 \times 10^{-8} \, \mathrm{\Omega m}$ [1]')
    plt.tight_layout()
    plt.show()

    # Extracting the impedance elements from zi_matrix
    plt.plot(f, 1E6 * l_a, label=r'$\rho_e = 100 \;\Omega m, \epsilon_r=1$', color='black', linestyle='-')
    plt.plot(f, 1E6 * l_b, label=r'$\rho_e = 100 \;\Omega m, \epsilon_r=20$', color='black', linestyle='--')
    plt.plot(f, 1E6 * l_c, label=r'$\rho_e = 2000 \;\Omega m, \epsilon_r=1$', color='black', linestyle='-.')

    # Additional plotting configurations
    plt.xscale('log')
    plt.xlim(1E3, 1E9)
    plt.ylim(1, 2.5)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'$L_s \, (mH/km)$')
    plt.grid(True)

    # Adjust the layout of the plots
    plt.suptitle('P.u.l. series inductance of the single overhead line\n'
                 r'$r_1 = 0.01 \, \mathrm{m}, h_1 = 10 \, \mathrm{m}, \rho = 1.68 \times 10^{-8} \, \mathrm{\Omega m}$ [1]') # pylint: disable=line-too-long
    plt.tight_layout()
    plt.show()


def main():
    """ Main function to perform the calculations and display the results."""

    # Clears the console screen and starts the timer
    os.system('cls' if os.name == 'nt' else 'clear')
    print("Calculations were started ... ...")
    start_time = time.time()

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
    elapsed_time = time.time() - start_time
    print(f"End of the routine! Time spent on simulation: {elapsed_time:.2f} seconds.\n")

    # Display the geometry of the transmission line
    graph(MTL).wires_and_cables()

    # Plot the series resistance as a function of frequency
    plot_impedance(frequency, pul_dict, p=0, q=0)


if __name__ == "__main__":
    main()
