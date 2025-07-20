""" This module contains the plotting functions for the two-wire line system. """

import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
import mom_so_analytical_formulation as analytic
import mom_so_green_matrices as green_matrices
#from mom_so_patel import MultilayeredLossyMedium as lossy_medium
from mom_so_geometry import UndergroundSystem, AuxiliaryGeometry


class BuriedSingleCoreCable(UndergroundSystem, AuxiliaryGeometry):
    """ This class performs the calculations for a single-core cable buried underground. """

    def __init__(self, mtl_dict, f, f_mom, green_evaluation='Analytically'):
        super().__init__(mtl_dict)
        self.f = f
        self.f_mom = f_mom
        self.green = green_evaluation
        self.analytical = None
        self.mom_numerical = None
        self.mom_analytical = None
        self.comsol = None
        self._perform_calculations(mtl_dict)
        # self._comsol_data()

    def _perform_calculations(self, mtl_dict):
        """ Performs analytical and numerical formulation calculations and stores the results. """

        # self.analytical = self._analytical_formulation(mtl_dict)
        # self.mom_numerical = self._mom_numerical_formulation(mtl_dict)
        self.mom_analytical = self._mom_analytical_formulation(mtl_dict)

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

    def _mom_analytical_formulation(self, mtl_dict):
        """ Calculates the series impedance for the given frequencies. """

        series_impedance = []

        # Post-processing parameters
        # post_processing = mom_so.PostProcessingParameters(mtl_dict)

        # Function to calculate conductor-to-conductor distances
        def conductors2conductors_dqp(green, p, q):
            return AuxiliaryGeometry().distance_vector_dqp(
                green.conductor_surfaces, p, green.conductor_surfaces, q)

        # Function to calculate conductor-to-hole distances
        def conductors2hole_dqp(green, p, q):
            return AuxiliaryGeometry().distance_vector_dqp(
                green.hole_surfaces, p, green.conductor_surfaces, q)

        # Calculate the series impedance for each frequency
        for f in self.f_mom:

            # Green's matrices as a function of frequency
            green = green_matrices.FullWaveAnalytically(mtl_dict, f)
            green2 = green_matrices.FullWaveNumerically(mtl_dict, f)

            print("\nGEOMETRY:")
            print("Principal Dimension N = SUM (2*Np+1): ", green.N)
            print("Principal Dimension Nhat = SUM (2*N_hat+1): ", green.Nhat)
            print("Surface type [1-conductor, 0-hole]:", green.surfaces_type)

            print("\nCONDUCTORS:")
            print(green.conductor_surfaces)

            print("Distance vector dpp:",
                  conductors2conductors_dqp(green, p=0, q=0))
            if len(green.conductor_surfaces) == 2:
                print(conductors2conductors_dqp(green, p=0, q=1))

            print("\nHOLES:")
            print(green.hole_surfaces)

            print("Distance vector dp-hole:",
                  conductors2hole_dqp(green, p=0, q=0))

            print("\nFREQUENCY DEPENDENT PARAMETERS:")
            print("Frequency [Hz]: ", f)
            print("Wave number, k_hat [rad/m]: ", green.khat[0])

            print("\n[GREEN'S MATRICES]:")
            print("\n[Gc_hat]:")
            print(green.ghatc())

            print("\n[Gc_hat] [Numerically]:")
            print(green2.ghatc())

            print("\n[G0_hat]:")
            print(green.g0_matrices()[0])

            print("\n[G0_til]:")
            print(green.g0_matrices()[1])

            print("\n[H_hat]:")
            print(green.hhat())

            print("\n[D1]:")
            print(green.d_matrices()[0])

            print("\n[D2]:")
            print(green.d_matrices()[1])

            print("\n[T]:")
            print(green.t_matrix())

            print("\n[CONDUCTORS SURFACE ADMITTANCE OPERATOR]:")

            print("\n[Ys]:")
            print(green.ys())

            print("\n[Yhat_s]:")
            print(green.yhats())

            print("\n Psi Matrix:")
            print(green.psi())

            print("\n[IMPEDANCE MATRIX]:")
            print("\n[U]:")
            print(np.array([[1]]))

            print("\n[Z]:")
            z = green.z_partial(green.psi())
            series_impedance.append(z)
            print(z)

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
        plt.scatter(self.f_mom, np.real(self.mom_analytical),
                    label='MoM-SO [1]', color='blue', marker='x', s=40)

        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1, 1E6)
        #plt.ylim(1E-5, 2E-2)
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Resistance p.u.l. (Ω/m)')
        # plt.title('Figure 2.4: P.u.l. resistance of the two-wire line of Sec. 2.6.1\n'
        #           f'for D = {self.d} m, a = 0.01 m, and σ = 5.8E7 S/m [1]')
        plt.grid(False)
        plt.show()