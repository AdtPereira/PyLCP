""" This module contains the plotting functions for the two-wire line system """

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .patel import HomogeneousLosslessMedium, LosslessPostProcessing
from .green import QuasiStatic
from .geometry import FreeSpace, AuxiliaryGeometry


class TwoWire(FreeSpace, AuxiliaryGeometry):
    """ This class contains the plotting functions for the system. """

    def __init__(self, mtl_dict, f, f_mom, green_evaluation='Analytically'):
        super().__init__(mtl_dict)
        self.f = f
        self.f_mom = f_mom
        self.green = green_evaluation
        self.analytical_results = None
        self.numerical_results = None
        self._perform_calculations(mtl_dict)

        # Distance between the conductors for graphical representation
        # self.d = self.distance_matrices(self.mtl)[0][0][1]
        self.d = 0

    def _perform_calculations(self, mtl_dict):
        """
        Performs analytical and numerical formulation calculations and stores the results.
        """

        # self.analytical_results = self._analytical_formulation(mtl_dict)
        self.numerical_results = self._numerical_formulation(mtl_dict)

    def _analytical_formulation(self, mtl_dict):
        """
        Performs analytical formulation for the given MTL and frequencies.
        """

        resistance_hf = []
        external_inductance = []
        series_impedance = []
        for f in self.f:
            parameters = analytic.TwoWire(mtl_dict, f).pul_parameters()
            resistance_hf.append(parameters[0][0, 1])
            external_inductance.append(parameters[1][0, 1])
            series_impedance.append(parameters[2][0, 1])

        return resistance_hf, external_inductance, series_impedance

    def _numerical_formulation(self, mtl_dict):
        """ Calculates the series impedance given frequencies using numerical formulation. """

        series_impedance = []

        # Green's matrix
        green = QuasiStatic(mtl_dict).g_tanaka(mode=self.green)

        # Print the Green's matrix
        # print("Greens' Matrix: \n", green)

        # Post-processing parameters
        post_processing = LosslessPostProcessing(mtl_dict)

        # Calculate the series impedance for each frequency
        for f in self.f_mom:
            mom_so_patel = HomogeneousLosslessMedium(mtl_dict, f)
            z_partial = mom_so_patel.z_partial(green)
            zs = post_processing.z_total(z_partial)
            series_impedance.append(zs[0][0])
        return series_impedance

    def plot_series_resistance(self):
        """ This function plots the series resistance as a function of frequency. """

        # Analytical Series Resistance
        plt.plot(self.f, np.real(self.analytical_results[2]),
                 label='Analytical (no proximity)', color='black', linestyle='-')

        # Asymptotic Series Resistance
        plt.plot(self.f, self.analytical_results[0],
                 label='Analytical (high-freq)', color='red', linestyle='--')

        # MoM Series Resistance
        plt.scatter(self.f_mom, np.real(self.numerical_results),
                    label='MoM-SO [1]', color='blue', marker='x', s=40)

        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1, 1E6)
        plt.ylim(1E-5, 2E-2)
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Resistance p.u.l. (Ω/m)')
        plt.title('Figure 2.4: P.u.l. resistance of the two-wire line of Sec. 2.6.1\n'
                  f'for D = {self.d} m, a = 0.01 m, and σ = 5.8E7 S/m [1]')
        plt.grid(False)
        plt.show()

    def plot_series_inductance(self):
        """ This function plots the series inductance as a function of frequency. """

        # Analytical Series Inductance
        plt.plot(self.f, 1E6 * np.imag(self.analytical_results[2]) / (2 * np.pi * self.f),
                 label='Analytical (no proximity)', color='black', linestyle='-')

        # Asymptotic Series Inductance
        plt.plot(self.f, 1E6 * np.array(self.analytical_results[1]),
                 label='Analytical (high-freq)', color='red', linestyle='--')

        # MoM Series Inductance
        plt.scatter(self.f_mom, 1E6 * np.imag(self.numerical_results) / (2 * np.pi * self.f_mom),
                    label='MoM-SO [1]', color='blue', marker='x', s=40)

        plt.xscale('log')
        plt.legend()
        plt.xlim(1, 1E6)
        if self.d == 0.1:
            plt.ylim(0.90, 1.06)
        elif self.d == 0.025:
            plt.ylim(0.25, 0.55)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Inductance p.u.l. (uH/m)')
        plt.title('Figure 2.5: P.u.l. inductance of the two-wire line of Sec. 2.6.1\n'
                  f'for D = {self.d} m, a = 0.01 m, and σ = 5.8E7 S/m [1]')
        plt.grid(False)
        plt.show()


class CoaxialCable(FreeSpace):
    """ 
    This class contains the plotting functions for the system.

    Args:
        mtl (str): The MTL (Multi-Terminal Line) object.
        f_analytic (list): A list of frequencies for analytical formulation.
        f_mom (list): A list of frequencies for numerical formulation.

    Attributes:
        d (float): Distance between the conductors [m].
    """

    def __init__(self, mtl_dict, f, f_mom, green_evaluation='Analytically'):
        super().__init__(mtl_dict)
        self.f = f
        self.f_mom = f_mom
        self.green = green_evaluation
        self.analytical = None
        self.numerical = None
        self.comsol = None
        self._perform_calculations(mtl_dict)
        self._comsol_data()

    def _perform_calculations(self, mtl_dict):
        """ Performs analytical and numerical formulation calculations and stores the results. """

        self.analytical = self._analytical_formulation(mtl_dict)
        self.numerical = self._numerical_formulation(mtl_dict)

    def _analytical_formulation(self, mtl_dict):
        """ Performs analytical formulation for the given MTL and frequencies. """
        series_impedance = []
        ametani_impedances = []

        for f in self.f:
            # Patel's Formulation
            z = analytic.SingleCoreCable(mtl_dict, f).pul_parameters()

            # Ametani's Formulation
            zz11, zz12, zz22 = analytic.Ametani(
                mtl_dict, f).impedance_two_layered_conductor()

            series_impedance.append(z)
            ametani_impedances.append([zz11, zz12, zz22])

        return series_impedance, ametani_impedances

    def _numerical_formulation(self, mtl_dict):
        """ Calculates the series impedance given frequencies using numerical formulation. """

        series_impedance = []

        # Green's matrix
        green = QuasiStatic(mtl_dict).g_tanaka(mode=self.green)

        # Post-processing parameters
        post_processing = mom_so.HomogeneousLosslessMediumPostProcessing(
            mtl_dict)

        # Calculate the series impedance for each frequency
        for f in self.f_mom:
            mom_so_patel = mom_so.HomogeneousLosslessMedium(mtl_dict, f)
            z_partial = mom_so_patel.z_partial(green)
            zs = post_processing.z_matrix(z_partial)
            series_impedance.append(zs[0][0])
        return series_impedance

    def _comsol_data(self):
        # Read the data from the file
        file_path = 'C:\\Users\\adilt\\OneDrive\\01 ACADEMIA\\06 MODELOS\\7.MoM-SO\\data'
        resistance = pd.read_csv(
            file_path+'\\comsol_resistance_coax.txt', sep=r'\s+', comment='%')
        inductance = pd.read_csv(
            file_path+r'\\comsol_inductance_coax.txt', sep=r'\s+', comment='%')

        # Rename the columns
        resistance.columns = [
            'freq (Hz)', 'Analytic (DC)', 'Analytic (HF)', 'COMSOL (mf/ec)']
        inductance.columns = [
            'freq (Hz)', 'Analytic (DC)', 'Analytic (HF)', 'COMSOL (mf/ec)']

        self.comsol = [
            resistance['freq (Hz)'], resistance['COMSOL (mf/ec)'], inductance['COMSOL (mf/ec)']]

    def plot_series_resistance(self):
        """
        This function plots the series resistance as a function of frequency.
        """

        # Analytical Series Resistance
        plt.plot(self.f, np.real(self.analytical[0]),
                 label='Analytic', color='black', linestyle='-')

        plt.plot(self.comsol[0], self.comsol[1],
                 label='COMSOL', color='red', marker='o', linestyle='None', markersize=2)

        # MoM Series Resistance
        plt.scatter(self.f_mom, np.real(self.numerical),
                    label='MoM-SO [1]', color='blue', marker='x', s=55)

        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1E0, 1E6)
        plt.ylim(1E-5, 1E-2)
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Resistance p.u.l. (Ω/m)')
        plt.title('Figure 2.6: P.u.l. resistance of a coaxial cable of Sec. 2.6.2\n'
                  'a = 22 mm, b = 39.5 mm, c = 44 mm [1]')
        plt.grid(False)
        plt.show()

    def plot_series_resistance_ametani(self):
        """This function plots the series resistance as a function of frequency."""
        zz11 = np.array([data[0] for data in self.analytical[1]])
        zz12 = np.array([data[1] for data in self.analytical[1]])
        zz22 = np.array([data[2] for data in self.analytical[1]])

        # Analytical Series Resistance
        plt.plot(self.f, np.real(self.analytical[0]),
                 label='Analytic', color='black', linestyle='-')

        plt.plot(self.f, np.real(zz11 - 2*zz12 + zz22),
                 label=r'Re($Z_{11}$ - $2*Z_{12}$ + $Z_{22}$)/$\omega$',
                 color='red', linestyle='--')

        plt.plot(self.f, np.real(zz11),
                 label=r'Re($Z_{11}$)/$\omega$',
                 color='black', linestyle=':')

        plt.plot(self.f, np.real(zz12),
                 label=r'Re($Z_{12}$)/$\omega$',
                 color='blue', linestyle=':')

        plt.plot(self.f, np.real(zz22),
                 label=r'Re($Z_{22}$)/$\omega$',
                 color='green', linestyle=':')

        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1E0, 1E6)
        plt.ylim(1E-5, 1E-2)
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Resistance p.u.l. (Ω/m)')
        plt.title('Figure 2.6: P.u.l. resistance of a coaxial cable of Sec. 2.6.2\n'
                  'a = 22 mm, b = 39.5 mm, c = 44 mm [1]')
        plt.grid(False)
        plt.show()

    def plot_series_inductance(self):
        """
        This function plots the series inductance as a function of frequency.
        """

        # Analytical Series Inductance
        plt.plot(self.f, 1E6 * np.imag(self.analytical[0]) / (2 * np.pi * self.f),
                 label='Analytic', color='black', linestyle='-')

        plt.plot(self.comsol[0], self.comsol[2],
                 label='COMSOL', color='red', marker='o', linestyle='None', markersize=2)

        # MoM Series Inductance
        plt.scatter(self.f_mom, 1E6 * np.imag(self.numerical) / (2 * np.pi * self.f_mom),
                    label='MoM-SO [1]', color='blue', marker='x', s=55)

        plt.xscale('log')
        plt.legend()
        plt.xlim(1E0, 1E6)
        plt.ylim(0.11, 0.19)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Inductance p.u.l. (uH/m)')
        plt.title('Figure 2.6: P.u.l. inductance of a coaxial cable of Sec. 2.6.2\n'
                  'a = 22 mm, b = 39.5 mm, c = 44 mm [1]')
        plt.grid(False)
        plt.show()

    def plot_series_inductance_ametani(self):
        """
        This function plots the series inductance as a function of frequency.
        """
        zz11 = np.array([data[0] for data in self.analytical[1]])
        zz12 = np.array([data[1] for data in self.analytical[1]])
        zz22 = np.array([data[2] for data in self.analytical[1]])

        # Analytical Series Inductance
        plt.plot(self.f, 1E6 * np.imag(self.analytical[0]) / (2 * np.pi * self.f),
                 label='Analytic', color='black', linestyle='-')

        plt.plot(self.f, 1E6 * np.imag(zz11 - 2*zz12 + zz22) / (2 * np.pi * self.f),
                 label=r'Im($Z_{11}$ - $2*Z_{12}$ + $Z_{22}$)/$\omega$',
                 color='red', linestyle='--')

        plt.plot(self.f, 1E6 * np.imag(zz11) / (2 * np.pi * self.f),
                 label=r'($Z_{11}$)/$\omega$', color='black', linestyle=':')

        plt.plot(self.f, 1E6 * np.imag(zz12) / (2 * np.pi * self.f),
                 label=r'($Z_{12}$)/$\omega$', color='blue', linestyle=':')

        plt.plot(self.f, 1E6 * np.imag(zz22) / (2 * np.pi * self.f),
                 label=r'Im($Z_{22}$)/$\omega$', color='green', linestyle=':')

        plt.xscale('log')
        plt.legend()
        plt.xlim(1E0, 1E6)
        plt.ylim(0, 0.19)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Inductance p.u.l. (uH/m)')
        plt.title('Figure 2.6: P.u.l. inductance of a coaxial cable of Sec. 2.6.2\n'
                  'a = 22 mm, b = 39.5 mm, c = 44 mm [1]')
        plt.grid(False)
        plt.show()


class EnclosureGIB(FreeSpace, AuxiliaryGeometry):
    """ This class contains the plotting functions for the system. """

    def __init__(self, mtl_dict, f, f_mom, green_evaluation='Analytically'):
        super().__init__(mtl_dict)
        self.f = f
        self.f_mom = f_mom
        self.green = green_evaluation
        self.analytical_results = None
        self.numerical_results = None
        self._perform_calculations(mtl_dict)

        # Distance between the conductors for graphical representation
        self.d = self.distance_matrices(self.mtl)[0][0][1]

    def _perform_calculations(self, mtl_dict):
        """ Performs analytical and numerical formulation calculations and stores the results. """

        # self.analytical_results = self._analytical_formulation(mtl_dict)
        self.numerical_results = self._numerical_formulation(mtl_dict)

    def _analytical_formulation(self, mtl_dict):
        """ Performs analytical formulation for the given MTL and frequencies. """

        resistance_hf = []
        external_inductance = []
        series_impedance = []
        for f in self.f:
            parameters = analytic.TwoWire(mtl_dict, f).pul_parameters()
            resistance_hf.append(parameters[0][0, 1])
            external_inductance.append(parameters[1][0, 1])
            series_impedance.append(parameters[2][0, 1])

        return resistance_hf, external_inductance, series_impedance

    def _numerical_formulation(self, mtl_dict):
        """ Calculates the series impedance given frequencies using numerical formulation. """

        series_impedance = []

        # Green's matrix
        green = QuasiStatic(mtl_dict).g_tanaka(mode=self.green)

        # Print the Green's matrix
        # print("Greens' Matrix: \n", green)

        # Post-processing parameters
        post_processing = mom_so.HomogeneousLosslessMediumPostProcessing(
            mtl_dict)

        # Calculate the series impedance for each frequency
        for f in self.f_mom:
            mom_so_patel = mom_so.HomogeneousLosslessMedium(mtl_dict, f)
            z_partial = mom_so_patel.z_partial(green)
            zs = post_processing.z_matrix(z_partial)
            series_impedance.append(zs[0][0])
        return series_impedance

    def plot_series_resistance(self):
        """ This function plots the series resistance as a function of frequency. """

        # Analytical Series Resistance
        # plt.plot(self.f, np.real(self.analytical_results[2]),
        #          label='Analytical (no proximity)', color='black', linestyle='-')

        # Asymptotic Series Resistance
        # plt.plot(self.f, self.analytical_results[0],
        #          label='Analytical (high-freq)', color='red', linestyle='--')

        # MoM Series Resistance
        plt.scatter(self.f_mom, np.real(self.numerical_results),
                    label='MoM-SO [1]', color='blue', marker='x', s=40)

        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1, 1E6)
        # plt.ylim(1E-5, 2E-2)
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Resistance p.u.l. (Ω/m)')
        plt.title('Figure 2.4: P.u.l. resistance of the two-wire line of Sec. 2.6.1\n'
                  f'for D = {self.d} m, a = 0.01 m, and σ = 5.8E7 S/m [1]')
        plt.grid(False)
        plt.show()

    def plot_series_inductance(self):
        """ This function plots the series inductance as a function of frequency. """

        # Analytical Series Inductance
        # plt.plot(self.f, 1E6 * np.imag(self.analytical_results[2]) / (2 * np.pi * self.f),
        #          label='Analytical (no proximity)', color='black', linestyle='-')

        # Asymptotic Series Inductance
        # plt.plot(self.f, 1E6 * np.array(self.analytical_results[1]),
        #          label='Analytical (high-freq)', color='red', linestyle='--')

        # MoM Series Inductance
        plt.scatter(self.f_mom, 1E6 * np.imag(self.numerical_results) / (2 * np.pi * self.f_mom),
                    label='MoM-SO [1]', color='blue', marker='x', s=40)

        plt.xscale('log')
        plt.legend()
        plt.xlim(1, 1E6)
        # if self.d == 0.1:
        #     plt.ylim(0.90, 1.06)
        # elif self.d == 0.025:
        #     plt.ylim(0.25, 0.55)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Inductance p.u.l. (uH/m)')
        plt.title('Figure 2.5: P.u.l. inductance of the two-wire line of Sec. 2.6.1\n'
                  f'for D = {self.d} m, a = 0.01 m, and σ = 5.8E7 S/m [1]')
        plt.grid(False)
        plt.show()
