# import numpy as np
# from scipy.special import jv, jvp #, iv, kv
# from scipy.constants import mu_0
# import matplotlib.pyplot as plt

# from mtl_data.mtl import MulticonductorTransmissionLine


# class SingleCoreCable(FreeSpace):
#     """ This class contains the analytical formulation of the system. """

#     def __init__(self, mtl_dict, frequency):
#         """
#         Initialize the AnalyticalFormulation class.

#         Parameters:
#         conductor (list): List of dictionaries containing the properties of the conductors.
#         frequency (float): The frequency of the system.
#         """
#         super().__init__(mtl_dict)

#         # Conductivity of the conductors [np.array]
#         self.sigma = np.array([c['conductivity'] for c in self.mtl])

#         # Angular frequency, rad/s [float]
#         self.jw = 1j * 2 * np.pi * frequency

#         # Propagation Constant [np.array]
#         self.gamma = np.sqrt(self.jw * mu_0 * self.sigma)

#         # Coaxial Cable radii
#         core = [c for c in self.mtl if c['conductor_name'] == 'core']
#         sheath = [c for c in self.mtl if c['conductor_name'] == 'sheath']

#         self.a = core[0]['radius'][1] if core else None
#         self.b, self.c = (sheath[0]['radius'] if sheath else (None, None))

#     def external_inductance(self):
#         """
#         Calculate the external inductance for a lossless coaxial
#         cable, L'.
#         # Equation 2.70 [1]
#         """
#         return mu_0 / (2 * np.pi) * np.log(self.b / self.a)

#     def internal_impedance(self):
#         """
#         Calculate the internal impedance of the inner conductor
#         Za (omega).
#         # Equation 2.71 [1]
#         """
#         # Propagation Constant of the inner conductor
#         gama_a = self.gamma[0] * self.a

#         # Intrinsic Impedance of the inner conductor
#         eta = np.sqrt(self.jw * mu_0 / self.sigma)[0]

#         # Internal Impedance of the inner conductor
#         za = eta / (2 * np.pi * self.a) * iv(0, gama_a) / iv(1, gama_a)

#         return za

#     def external_impedance(self):
#         """
#         Calculate the external impedance of the inner conductor
#         Zb (omega).
#         # Equation 2.72 [1]
#         """
#         # Propagation Constant of the inner conductor
#         gama_b = self.gamma[0] * self.b
#         gama_c = self.gamma[0] * self.c

#         # Intrinsic Impedance of the inner conductor
#         eta = np.sqrt(self.jw * mu_0 / self.sigma)[0]

#         numerator = iv(0, gama_b) * kv(1, gama_c) + (
#             kv(0, gama_b) * iv(1, gama_c))

#         denominator = iv(1, gama_c) * kv(1, gama_b) - (
#             iv(1, gama_b) * kv(1, gama_c))

#         # External Impedance of the inner conductor
#         zb = eta / (2 * np.pi * self.b) * numerator / denominator

#         return zb

#     def pul_parameters(self):
#         """
#         This function calculates the series resistance of the system using
#         the high frequency approximation.

#         Returns:
#         tuple: A tuple containing the high frequency resistance, external inductance,
#         and matrix impedance.
#         """
#         # Matrix Impedance, z (Ω/m)
#         # Equation (2.69) [1]
#         l_ext = self.external_inductance()
#         za = self.internal_impedance()
#         zb = self.external_impedance()
#         zs = self.jw * l_ext + za + zb

#         return zs


# class Ametani(FreeSpace):
#     """ This class contains the analytical formulation of the system. """

#     def __init__(self, mtl_dict, frequency):
#         """
#         Initialize the AnalyticalFormulation class.

#         Parameters:
#         conductor (list): List of dictionaries containing the properties of the conductors.
#         frequency (float): The frequency of the system.
#         """
#         super().__init__(mtl_dict)

#         # Conductivity of the conductors [np.array]
#         self.sigma = np.array([c['conductivity'] for c in self.mtl])

#         # Angular frequency, rad/s [float]
#         self.jw = 1j * 2 * np.pi * frequency

#         # Permeability of the medium [np.array]
#         self.mu = np.array([mu_0 * c['relative_permeability']
#                            for c in self.mtl])

#         # Coaxial Cable radii
#         self.a = None
#         self.b = None
#         self.b_prime = None
#         self.c = None

#         # Call the function to configure the parameters
#         self._scc_from_data()

#     def _scc_from_data(self):
#         """
#         This function configures the parameters of the coaxial cable.
#         """

#         for conductor in self.mtl:
#             if conductor['conductor_name'] == 'core':
#                 self.a = conductor['radius'][0]
#                 self.b = conductor['radius'][1]

#             elif conductor['conductor_name'] == 'sheath':
#                 self.b_prime = conductor['radius'][0]
#                 self.c = conductor['radius'][1]

#     def parameter_m(self, mu, sigma):
#         """
#         Calculate the parameter m for the two-layered conductor.
#         """
#         return np.sqrt(self.jw * mu * sigma)

#     def impedance_two_layered_conductor(self):
#         """
#         Calculate the impedance of a two-layered conductor.
#         """

#         # Intermediate variables
#         m1 = self.parameter_m(self.mu[0], self.sigma[0])
#         m2 = self.parameter_m(self.mu[1], self.sigma[1])
#         x1 = m1 * self.a
#         x2 = m1 * self.b
#         x3 = m2 * self.b_prime
#         x4 = m2 * self.c

#         # Intermediate variables A, B, E, F, R
#         aa = kv(1, x1) * iv(0, x2) + iv(1, x1) * kv(0, x2)
#         bb = kv(1, x1) * iv(1, x2) - iv(1, x1) * kv(1, x2)
#         ee = iv(0, x3) * kv(1, x4) + kv(0, x3) * iv(1, x4)
#         ff = kv(1, x3) * iv(1, x4) - iv(1, x3) * kv(1, x4)
#         rr = kv(1, x3) * iv(0, x4) + iv(1, x3) * kv(0, x4)

#         # Calculating z10, z2i, z2m, z20, z12
#         rho1 = 1 / self.sigma[0]
#         rho2 = 1 / self.sigma[1]

#         # Solid conductor case
#         if x1 == 0:
#             z10 = (m1 * rho1 / (2 * np.pi * self.b)) * iv(0, x2) / iv(1, x2)
#         else:
#             z10 = (m1 * rho1 / (2 * np.pi * self.b)) * aa / bb

#         z2i = (m2 * rho2 / (2 * np.pi * self.b_prime)) * ee / ff
#         z2m = rho2 / (2 * np.pi * self.b_prime * self.c * ff)
#         z20 = (m2 * rho2 / (2 * np.pi * self.c)) * rr / ff
#         z12 = self.jw * (mu_0 / 2 / np.pi) * np.log(self.b_prime / self.b)

#         # Calculating Z11, Z12, Z22
#         zz22 = z20
#         zz12 = zz22 - z2m
#         zz11 = z10 + z12 + z2i - z2m + zz12

#         return zz11, zz12, zz22

# class CoaxialCable(FreeSpace):
#     """ 
#     This class contains the plotting functions for the system.

#     Args:
#         mtl (str): The MTL (Multi-Terminal Line) object.
#         f_analytic (list): A list of frequencies for analytical formulation.
#         f_mom (list): A list of frequencies for numerical formulation.

#     Attributes:
#         d (float): Distance between the conductors [m].
#     """

#     def __init__(self, mtl_dict, f, f_mom, green_evaluation='Analytically'):
#         super().__init__(mtl_dict)
#         self.f = f
#         self.f_mom = f_mom
#         self.green = green_evaluation
#         self.analytical = None
#         self.numerical = None
#         self.comsol = None
#         self._perform_calculations(mtl_dict)
#         self._comsol_data()

#     def _perform_calculations(self, mtl_dict):
#         """ Performs analytical and numerical formulation calculations and stores the results. """

#         self.analytical = self._analytical_formulation(mtl_dict)
#         self.numerical = self._numerical_formulation(mtl_dict)

#     def _analytical_formulation(self, mtl_dict):
#         """ Performs analytical formulation for the given MTL and frequencies. """
#         series_impedance = []
#         ametani_impedances = []

#         for f in self.f:
#             # Patel's Formulation
#             z = analytic.SingleCoreCable(mtl_dict, f).pul_parameters()

#             # Ametani's Formulation
#             zz11, zz12, zz22 = analytic.Ametani(
#                 mtl_dict, f).impedance_two_layered_conductor()

#             series_impedance.append(z)
#             ametani_impedances.append([zz11, zz12, zz22])

#         return series_impedance, ametani_impedances

#     def _numerical_formulation(self, mtl_dict):
#         """ Calculates the series impedance given frequencies using numerical formulation. """

#         series_impedance = []

#         # Green's matrix
#         green = QuasiStatic(mtl_dict).g_tanaka(green_evaluation=self.green)

#         # Post-processing parameters
#         post_processing = mom_so.HomogeneousLosslessMediumPostProcessing(mtl_dict)

#         # Calculate the series impedance for each frequency
#         for f in self.f_mom:
#             mom_so_patel = mom_so.HomogeneousLosslessMedium(mtl_dict, f)
#             z_partial = mom_so_patel.z_partial(green)
#             zs = post_processing.z_matrix(z_partial)
#             series_impedance.append(zs[0][0])
#         return series_impedance

#     def _comsol_data(self):
#         # Read the data from the file
#         file_path = 'C:\\Users\\adilt\\OneDrive\\01 ACADEMIA\\06 MODELOS\\7.MoM-SO\\data'
#         resistance = pd.read_csv(
#             file_path+'\\comsol_resistance_coax.txt', sep=r'\s+', comment='%')
#         inductance = pd.read_csv(
#             file_path+r'\\comsol_inductance_coax.txt', sep=r'\s+', comment='%')

#         # Rename the columns
#         resistance.columns = [
#             'freq (Hz)', 'Analytic (DC)', 'Analytic (HF)', 'COMSOL (mf/ec)']
#         inductance.columns = [
#             'freq (Hz)', 'Analytic (DC)', 'Analytic (HF)', 'COMSOL (mf/ec)']

#         self.comsol = [
#             resistance['freq (Hz)'], resistance['COMSOL (mf/ec)'], inductance['COMSOL (mf/ec)']]

#     def plot_series_resistance(self):
#         """
#         This function plots the series resistance as a function of frequency.
#         """

#         # Analytical Series Resistance
#         plt.plot(self.f, np.real(self.analytical[0]),
#                  label='Analytic', color='black', linestyle='-')

#         plt.plot(self.comsol[0], self.comsol[1],
#                  label='COMSOL', color='red', marker='o', linestyle='None', markersize=2)

#         # MoM Series Resistance
#         plt.scatter(self.f_mom, np.real(self.numerical),
#                     label='MoM-SO [1]', color='blue', marker='x', s=55)

#         plt.xscale('log')
#         plt.yscale('log')
#         plt.xlim(1E0, 1E6)
#         plt.ylim(1E-5, 1E-2)
#         plt.legend()
#         plt.xlabel('Frequency (Hz)')
#         plt.ylabel('Series Resistance p.u.l. (Ω/m)')
#         plt.title('Figure 2.6: P.u.l. resistance of a coaxial cable of Sec. 2.6.2\n'
#                   'a = 22 mm, b = 39.5 mm, c = 44 mm [1]')
#         plt.grid(False)
#         plt.show()

#     def plot_series_resistance_ametani(self):
#         """This function plots the series resistance as a function of frequency."""
#         zz11 = np.array([data[0] for data in self.analytical[1]])
#         zz12 = np.array([data[1] for data in self.analytical[1]])
#         zz22 = np.array([data[2] for data in self.analytical[1]])

#         # Analytical Series Resistance
#         plt.plot(self.f, np.real(self.analytical[0]),
#                  label='Analytic', color='black', linestyle='-')

#         plt.plot(self.f, np.real(zz11 - 2*zz12 + zz22),
#                  label=r'Re($Z_{11}$ - $2*Z_{12}$ + $Z_{22}$)/$\omega$',
#                  color='red', linestyle='--')

#         plt.plot(self.f, np.real(zz11),
#                  label=r'Re($Z_{11}$)/$\omega$',
#                  color='black', linestyle=':')

#         plt.plot(self.f, np.real(zz12),
#                  label=r'Re($Z_{12}$)/$\omega$',
#                  color='blue', linestyle=':')

#         plt.plot(self.f, np.real(zz22),
#                  label=r'Re($Z_{22}$)/$\omega$',
#                  color='green', linestyle=':')

#         plt.xscale('log')
#         plt.yscale('log')
#         plt.xlim(1E0, 1E6)
#         plt.ylim(1E-5, 1E-2)
#         plt.legend()
#         plt.xlabel('Frequency (Hz)')
#         plt.ylabel('Series Resistance p.u.l. (Ω/m)')
#         plt.title('Figure 2.6: P.u.l. resistance of a coaxial cable of Sec. 2.6.2\n'
#                   'a = 22 mm, b = 39.5 mm, c = 44 mm [1]')
#         plt.grid(False)
#         plt.show()

#     def plot_series_inductance(self):
#         """
#         This function plots the series inductance as a function of frequency.
#         """

#         # Analytical Series Inductance
#         plt.plot(self.f, 1E6 * np.imag(self.analytical[0]) / (2 * np.pi * self.f),
#                  label='Analytic', color='black', linestyle='-')

#         plt.plot(self.comsol[0], self.comsol[2],
#                  label='COMSOL', color='red', marker='o', linestyle='None', markersize=2)

#         # MoM Series Inductance
#         plt.scatter(self.f_mom, 1E6 * np.imag(self.numerical) / (2 * np.pi * self.f_mom),
#                     label='MoM-SO [1]', color='blue', marker='x', s=55)

#         plt.xscale('log')
#         plt.legend()
#         plt.xlim(1E0, 1E6)
#         plt.ylim(0.11, 0.19)
#         plt.xlabel('Frequency (Hz)')
#         plt.ylabel('Series Inductance p.u.l. (uH/m)')
#         plt.title('Figure 2.6: P.u.l. inductance of a coaxial cable of Sec. 2.6.2\n'
#                   'a = 22 mm, b = 39.5 mm, c = 44 mm [1]')
#         plt.grid(False)
#         plt.show()

#     def plot_series_inductance_ametani(self):
#         """
#         This function plots the series inductance as a function of frequency.
#         """
#         zz11 = np.array([data[0] for data in self.analytical[1]])
#         zz12 = np.array([data[1] for data in self.analytical[1]])
#         zz22 = np.array([data[2] for data in self.analytical[1]])

#         # Analytical Series Inductance
#         plt.plot(self.f, 1E6 * np.imag(self.analytical[0]) / (2 * np.pi * self.f),
#                  label='Analytic', color='black', linestyle='-')

#         plt.plot(self.f, 1E6 * np.imag(zz11 - 2*zz12 + zz22) / (2 * np.pi * self.f),
#                  label=r'Im($Z_{11}$ - $2*Z_{12}$ + $Z_{22}$)/$\omega$',
#                  color='red', linestyle='--')

#         plt.plot(self.f, 1E6 * np.imag(zz11) / (2 * np.pi * self.f),
#                  label=r'($Z_{11}$)/$\omega$', color='black', linestyle=':')

#         plt.plot(self.f, 1E6 * np.imag(zz12) / (2 * np.pi * self.f),
#                  label=r'($Z_{12}$)/$\omega$', color='blue', linestyle=':')

#         plt.plot(self.f, 1E6 * np.imag(zz22) / (2 * np.pi * self.f),
#                  label=r'Im($Z_{22}$)/$\omega$', color='green', linestyle=':')

#         plt.xscale('log')
#         plt.legend()
#         plt.xlim(1E0, 1E6)
#         plt.ylim(0, 0.19)
#         plt.xlabel('Frequency (Hz)')
#         plt.ylabel('Series Inductance p.u.l. (uH/m)')
#         plt.title('Figure 2.6: P.u.l. inductance of a coaxial cable of Sec. 2.6.2\n'
#                   'a = 22 mm, b = 39.5 mm, c = 44 mm [1]')
#         plt.grid(False)
#         plt.show()


# class EnclosureGIB(FreeSpace, AuxiliaryGeometry):
    # """ This class contains the plotting functions for the system. """

    # def __init__(self, mtl_dict, f, f_mom, green_evaluation='Analytically'):
    #     super().__init__(mtl_dict)
    #     self.f = f
    #     self.f_mom = f_mom
    #     self.green = green_evaluation
    #     self.analytical_results = None
    #     self.numerical_results = None
    #     self._perform_calculations(mtl_dict)

    #     # Distance between the conductors for graphical representation
    #     self.d = self.distance_matrices(self.mtl)[0][0][1]

    # def _perform_calculations(self, mtl_dict):
    #     """ Performs analytical and numerical formulation calculations and stores the results. """

    #     #self.analytical_results = self._analytical_formulation(mtl_dict)
    #     self.numerical_results = self._numerical_formulation(mtl_dict)

    # def _analytical_formulation(self, mtl_dict):
    #     """ Performs analytical formulation for the given MTL and frequencies. """

    #     resistance_hf = []
    #     external_inductance = []
    #     series_impedance = []
    #     for f in self.f:
    #         parameters = analytic.TwoWire(mtl_dict, f).pul_parameters()
    #         resistance_hf.append(parameters[0][0, 1])
    #         external_inductance.append(parameters[1][0, 1])
    #         series_impedance.append(parameters[2][0, 1])

    #     return resistance_hf, external_inductance, series_impedance

    # def _numerical_formulation(self, mtl_dict):
    #     """ Calculates the series impedance given frequencies using numerical formulation. """

    #     series_impedance = []

    #     # Green's matrix
    #     green = QuasiStatic(mtl_dict).g_tanaka(green_evaluation=self.green)

    #     # Print the Green's matrix
    #     # print("Greens' Matrix: \n", green)

    #     # Post-processing parameters
    #     post_processing = mom_so.HomogeneousLosslessMediumPostProcessing(mtl_dict)

    #     # Calculate the series impedance for each frequency
    #     for f in self.f_mom:
    #         mom_so_patel = mom_so.HomogeneousLosslessMedium(mtl_dict, f)
    #         z_partial = mom_so_patel.z_partial(green)
    #         zs = post_processing.z_matrix(z_partial)
    #         series_impedance.append(zs[0][0])
    #     return series_impedance

    # def plot_series_resistance(self):
    #     """ This function plots the series resistance as a function of frequency. """

    #     # Analytical Series Resistance
    #     # plt.plot(self.f, np.real(self.analytical_results[2]),
    #     #          label='Analytical (no proximity)', color='black', linestyle='-')

    #     # Asymptotic Series Resistance
    #     # plt.plot(self.f, self.analytical_results[0],
    #     #          label='Analytical (high-freq)', color='red', linestyle='--')

    #     # MoM Series Resistance
    #     plt.scatter(self.f_mom, np.real(self.numerical_results),
    #                 label='MoM-SO [1]', color='blue', marker='x', s=40)

    #     plt.xscale('log')
    #     plt.yscale('log')
    #     plt.xlim(1, 1E6)
    #     # plt.ylim(1E-5, 2E-2)
    #     plt.legend()
    #     plt.xlabel('Frequency (Hz)')
    #     plt.ylabel('Series Resistance p.u.l. (Ω/m)')
    #     plt.title('Figure 2.4: P.u.l. resistance of the two-wire line of Sec. 2.6.1\n'
    #               f'for D = {self.d} m, a = 0.01 m, and σ = 5.8E7 S/m [1]')
    #     plt.grid(False)
    #     plt.show()

    # def plot_series_inductance(self):
    #     """ This function plots the series inductance as a function of frequency. """

    #     # Analytical Series Inductance
    #     # plt.plot(self.f, 1E6 * np.imag(self.analytical_results[2]) / (2 * np.pi * self.f),
    #     #          label='Analytical (no proximity)', color='black', linestyle='-')

    #     # Asymptotic Series Inductance
    #     # plt.plot(self.f, 1E6 * np.array(self.analytical_results[1]),
    #     #          label='Analytical (high-freq)', color='red', linestyle='--')

    #     # MoM Series Inductance
    #     plt.scatter(self.f_mom, 1E6 * np.imag(self.numerical_results) / (2 * np.pi * self.f_mom),
    #                 label='MoM-SO [1]', color='blue', marker='x', s=40)

    #     plt.xscale('log')
    #     plt.legend()
    #     plt.xlim(1, 1E6)
    #     # if self.d == 0.1:
    #     #     plt.ylim(0.90, 1.06)
    #     # elif self.d == 0.025:
    #     #     plt.ylim(0.25, 0.55)
    #     plt.xlabel('Frequency (Hz)')
    #     plt.ylabel('Series Inductance p.u.l. (uH/m)')
    #     plt.title('Figure 2.5: P.u.l. inductance of the two-wire line of Sec. 2.6.1\n'
    #               f'for D = {self.d} m, a = 0.01 m, and σ = 5.8E7 S/m [1]')
    #     plt.grid(False)
    #     plt.show()

