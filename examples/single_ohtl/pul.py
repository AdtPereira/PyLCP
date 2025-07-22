import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# Adiciona a raiz do projeto ao PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..')))

from data.systems import MTL_MODELS
from data.graph import GraphicRepresentation as graph

from ohtl.pul_parameters import PerUnitParameters
from mom_so import green, patel

# Multiconductor Transmission Line choices
MTL = MTL_MODELS['overhead']['xue']
MTL = MTL_MODELS['overhead']['deConti']

def plot_zi(freq, z_dict, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """

    # Extracting the impedance elements from zi_matrix
    exact = np.array([item[p][q] for item in z_dict['zi_exact']])
    approx = np.array([item[p][q] for item in z_dict['zi_approx']])
    ri_cc = np.array([item[p][q] for item in z_dict['ri_cc']])
    zi_hf = np.array([item[p][q] for item in z_dict['zi_hf']])
    zi_mom = np.array([item[p][q] for item in z_dict['zi_momso']])

    plt.plot(freq['Analytically'], np.real(1E3 * exact), label=fr'$R_{{{ p+1}}}$ (Exact-Form)', color='black', linestyle='-')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(1E3 * approx), label=fr'$R_{{{p+1}}}$ (Closed-Form)', color=( 1, 0, 0, 0.5), linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(1E3 * ri_cc), label=fr'$R_{{cc({p+1})}}$', color=(0, 1, 0, 0.5), linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.real(1E3 * zi_hf), label=fr'$Z_{{hf({p+1})}}$', color=(0, 0, 1, 0.5), linestyle='--')  # pylint: disable=line-too-long
    plt.scatter(freq['Numerically'], np.real(1E3 * zi_mom), label=fr'$R_{{{p+1}}}$ (MoM-SO [1])', color='blue', marker='x', s=40)  # pylint: disable=line-too-long

    # Optional: Additional plotting configurations like labels, grid, etc.
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1E1, 1E6)
    plt.ylim(1E-2, 1E2)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Series Resistance p.u.l. (Ω/km)')
    plt.title('P.u.l. series resistance of the single overhead line\n'
              r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},'
              r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
    plt.grid(True)
    plt.show()


def plot_zs(freq, z_dict, zg_dict, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """

    # Extracting the impedance elements from zi_matrix
    zi = np.array([item[p, q] for item in z_dict['zi_exact']])
    zg = np.array([item[p, q] for item in zg_dict['quasi_tem']])
    ze = np.array([item[p, q] for item in z_dict['ze']])
    zs = zi + ze + zg

    plt.subplot(1, 2, 1)
    plt.plot(freq['Analytically'], np.absolute(1E3 * zs), label='$Z_s$', color='black', linestyle='-')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.absolute(1E3 * zi), label='$Z_i$', color='green', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.absolute(1E3 * ze), label='$Z_e$', color='red', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.absolute(1E3 * zg), label='$Z_g$', color='blue', linestyle='--')  # pylint: disable=line-too-long

    # Additional plotting configurations
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1E0, 1E7)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('$|Z_s|$ (Ω/km)')
    plt.grid(True)

    # Extracting the impedance elements from zi_matrix
    plt.subplot(1, 2, 2)
    plt.plot(freq['Analytically'], np.angle(zs, deg=True), label='$Z_s$', color='black', linestyle='-')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.angle(zi, deg=True), label='$Z_i$', color='green', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.angle(ze, deg=True), label='$Z_e$', color='red', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.angle(zg, deg=True), label='$Z_g$', color='blue', linestyle='--')  # pylint: disable=line-too-long

    # Additional plotting configurations
    plt.xscale('log')
    plt.xlim(1E0, 1E7)
    plt.ylim(0, 100)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase Angle (deg)')
    plt.grid(True)

    # Adjust the layout of the plots
    plt.suptitle('P.u.l. series impedance of the single overhead line\n'
                 r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},'
                 r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
    plt.tight_layout()
    plt.show()


def plot_zg(freq, zg_dict, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """

    # Extracting the impedance elements from zg dictionary
    quasitem = [item[p, q] for item in zg_dict['quasi_tem']]
    quasitem_log = [item[p, q] for item in zg_dict['quasitem_log']]
    sunde = [item[p, q] for item in zg_dict['sunde']]
    sunde_log = [item[p, q] for item in zg_dict['sunde_log']]
    carson = [item[p, q] for item in zg_dict['carson']]
    deri = [item[p, q] for item in zg_dict['deri']]

    plt.subplot(1, 2, 1)
    # Integral Equations
    plt.plot(freq['Analytically'], np.absolute(quasitem), label='Integral Equations',
             color='black', linestyle='-')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.absolute(sunde),
             color='black', linestyle='-')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.absolute(carson),
             color='black', linestyle='-')  # pylint: disable=line-too-long

    # Logarithmic Approximation
    plt.plot(freq['Analytically'], np.absolute(quasitem_log), label='Quasi-TEM (Approx. Log.)',
             color='green', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.absolute(sunde_log), label='Sunde (1968) (Approx. Log.)',
             color='red', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.absolute(deri), label='A. Deri (1981)',
             color='blue', linestyle='--')  # pylint: disable=line-too-long

    # Additional plotting configurations
    plt.xscale('log')
    plt.xlim(1E3, 1E7)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('$|Z_g|$ (Ω/m)')
    plt.grid(True)

    plt.subplot(1, 2, 2)
    # Integral Equations
    plt.plot(freq['Analytically'], np.angle(quasitem, deg=True), label='Integral Equations', color='black', linestyle='-')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.angle(sunde, deg=True),
             color='black', linestyle='-')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.angle(carson, deg=True),
             color='black', linestyle='-')  # pylint: disable=line-too-long

    # Logarithmic Approximation
    plt.plot(freq['Analytically'], np.angle(quasitem_log, deg=True), label='Quasi-TEM (Approx. Log.)', color='green', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.angle(sunde_log, deg=True), label='Sunde (1968) (Approx. Log.)', color='red', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.angle(deri, deg=True), label='A. Deri (1981)', color='blue', linestyle='--')  # pylint: disable=line-too-long

    # Additional plotting configurations
    plt.xscale('log')
    plt.xlim(1E3, 1E7)
    plt.ylim(0, 90)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase Angle (deg)')
    plt.grid(True)

    # Adjust the layout of the plots
    plt.suptitle('P.u.l. earth-return impedance of the single overhead line from classical theory\n'
                 r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},'
                 r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
    plt.tight_layout()
    plt.show()


def plot_ys(freq, ys_dict, p, q):
    """
    This function plots the series resistance as a function of frequency.

    Parameters:
    freq (array): Frequency array.
    zi_matrix (list of matrices): Matrix containing impedance values.
    p (int): Row index in the impedance matrix.
    q (int): Column index in the impedance matrix.
    """

    # Extracting the impedance elements from zi_matrix
    ye = np.array([item[p, q] for item in ys_dict['ye']])
    yg = np.array([item[p, q] for item in ys_dict['yg_approx']])

    plt.subplot(1, 2, 1)
    plt.plot(freq['Analytically'], np.absolute(1E3 * ye), label='$Y_e$',
             color='red', linestyle='--')  # pylint: disable=line-too-long
    plt.plot(freq['Analytically'], np.absolute(1E3 * yg), label='$Y_g$',
             color='blue', linestyle='--')  # pylint: disable=line-too-long

    # Additional plotting configurations
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1E0, 1E7)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('$|Y|$ (S/km)')
    plt.grid(True)

    # Extracting the impedance elements from zi_matrix
    plt.subplot(1, 2, 2)
    plt.plot(freq['Analytically'], np.angle(ye, deg=True),
             label='$Y_e$', color='red', linestyle='--')
    plt.plot(freq['Analytically'], np.angle(yg, deg=True),
             label='$Y_g$', color='blue', linestyle='--')

    # Additional plotting configurations
    plt.xscale('log')
    plt.xlim(1E0, 1E7)
    plt.ylim(0, 100)
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase Angle (deg)')
    plt.grid(True)

    # Adjust the layout of the plots
    plt.suptitle('P.u.l. longitudinal admittance of the single overhead line\n'
                 r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},'
                 r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
    plt.tight_layout()
    plt.show()


def main():
    """ Main function to perform the calculations and display the results."""

    # Clears the console screen and starts the timer
    os.system('cls' if os.name == 'nt' else 'clear')
    print("Calculations were started ... ...")
    start_time = time.time()

    # Calculate the series impedance for each frequency
    frequency = {
        'Analytically': np.logspace(0, 7, num=200),
        'Numerically': np.logspace(0, 7, num=30)
    }

    # Dictionary to hold the earth return impedance calculations
    zg_results = {
        'quasi_tem': [],
        'quasitem_log': [],
        'sunde': [],
        'sunde_log': [],
        'carson': [],
        'deri': []
    }

    # Dictionary to hold the series impedance calculations
    z_results = {
        'zi_exact': [],
        'zi_approx': [],
        'ri_cc': [],
        'zi_hf': [],
        'zi_momso': [],
        'ze': []
    }

    # Analytical Formulation
    for f in frequency['Analytically']:
        pul = PerUnitParameters(MTL, f, sigma_1=1/200, er_1=5)

        # Internal and external impedance matrices
        z_results['zi_exact'].append(pul.internal_impedance(type_form='bessel')[0])
        z_results['zi_approx'].append(pul.internal_impedance()[0])
        z_results['ri_cc'].append(pul.internal_impedance()[1])
        z_results['zi_hf'].append(pul.internal_impedance()[2])
        z_results['ze'].append(pul.external_impedance())

        # Store the results in the dictionary
        zg_results['quasi_tem'].append(pul.earth_return_impedance())
        zg_results['sunde'].append(pul.earth_return_impedance(type_form='sunde'))
        zg_results['carson'].append(pul.earth_return_impedance(type_form='carson'))
        zg_results['quasitem_log'].append(pul.earth_return_impedance(type_form='quasitem_log'))
        zg_results['sunde_log'].append(pul.earth_return_impedance(type_form='approx_log'))
        zg_results['deri'].append(pul.earth_return_impedance(type_form='deri'))

    # MoM-SO routine
    green_list = ['Analytically', 'Numerically']
    green_matrix = green.QuasiStatic(MTL).g_tanaka(mode=green_list[0])
    for f in frequency['Numerically']:
        mom_so = patel.HomogeneousLosslessMedium(MTL, f)
        z_results['zi_momso'].append(mom_so.z_partial(green_matrix))

    # End the timer
    elapsed_time = time.time() - start_time
    print(f"End of the routine! Time spent on simulation: {elapsed_time:.2f} seconds.\n")

    # Display the geometry of the transmission line
    graph(MTL).wires_and_cables(line_type='overhead')

    # Plot the series resistance as a function of frequency
    plot_zi(frequency, z_results, p=0, q=0)
    plot_zg(frequency, zg_results, p=0, q=0)
    plot_zs(frequency, z_results, zg_results, p=0, q=0)


if __name__ == "__main__":
    main()
