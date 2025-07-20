"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

REFERENCES:
[1] 

"""


import numpy as np
from scipy.special import jv, jvp #, iv, kv
from scipy.constants import mu_0
import matplotlib.pyplot as plt
from data.multiconductor import MulticonductorTransmissionLine as Multiconductor


class Bifilar(Multiconductor):
    """ This class contains the analytical formulation of the system. """

    def __init__(self, mtl):
        super().__init__(mtl)

        # Kelvin Functions
        self.kelvin_exp = np.exp(1j * 3 * np.pi / 4)

    def ber(self, xi):
        """ Kelvin ber(xi) function """
        return np.real(jv(0, xi * self.kelvin_exp))

    def bei(self, xi):
        """ Kelvin bei(xi) function """
        return np.imag(jv(0, xi * self.kelvin_exp))

    def ber_prime(self, xi):
        """ Kelvin ber'(xi) derivative function """
        return np.real(self.kelvin_exp * jvp(0, xi * self.kelvin_exp, 1))

    def bei_prime(self, xi):
        """ Kelvin bei'(xi) derivative function """
        return np.imag(self.kelvin_exp * jvp(0, xi * self.kelvin_exp, 1))

    def series_impedance(self, f):
        """
        This function calculates the series resistance of the system using 
        the high frequency approximation.

        Returns:
        tuple: A tuple containing the high frequency resistance, external inductance, 
        and matrix impedance.
        """
        # Angular frequency, rad/s [float]
        w = 2 * np.pi * f

        # Skin Depth [np.array]
        delta = np.sqrt(1 / (w / 2 * mu_0 * self.sigma))

        # Surface Resistance [np.array]
        Rs = 1 / (self.sigma * delta) # pylint: disable=invalid-name

        # Outer Radii of the conductors [np.array]
        ap = np.array([cp['radius'][1] for cp in self.mtl])

        # Matrix Distance [np.array]
        D = self.D[0]  # pylint: disable=invalid-name

        # High Frequency Resistance and External Inductance [np.array]
        N = len(self.mtl)-1  # pylint: disable=invalid-name
        Rhf = np.zeros((N, N))  # pylint: disable=invalid-name
        Lext = np.zeros_like(Rhf)  # pylint: disable=invalid-name
        Zi = np.zeros_like(Rhf, dtype=complex)  # pylint: disable=invalid-name

        # Constant Term and Bessel argument
        Xi = np.sqrt(2) * ap / delta # pylint: disable=invalid-name
        constant_term = 1 / (np.sqrt(2) * np.pi * ap * self.sigma * delta)

        for p in range(N):
            # 1st solution: High Frequency Approximation
            # These formulas account for proximity effect
            # only at high frequencies.These formulas
            # account for proximity effect only at high
            # frequencies.

            # Common fraction term
            D_2a = D[p][p+1] / 2 / ap[p] # pylint: disable=invalid-name

            # Surface resistance
            Rs_pia = Rs[p] / np.pi / ap[p] # pylint: disable=invalid-name

            # High Frequency Resistance (Ω/m)
            # Equation (2.64) [1]
            Rhf[p] = Rs_pia * D_2a / np.sqrt(D_2a ** 2 - 1)

            # External Inductance (H/m)
            # Equation (2.65) [1]
            Lext[p] = mu_0 / np.pi * np.arccosh(D_2a)

            # 2nd solution: Internal Impedance Matrix, z_int (Ω/m)
            # These formulas captures skin effect, but proximity
            # effect is neglected
            # Equation (2.67) [1]
            ber_bei = self.ber(Xi[p]) + 1j * self.bei(Xi[p])
            beip_berp = self.bei_prime(Xi[p]) - 1j * self.ber_prime(Xi[p])

            # Internal Impedance Matrix, Zi (Ω/m)
            Zi[p] = constant_term[p] * ber_bei / beip_berp

        # Matrix Impedance, Zs (Ω/m)
        # Equation (2.68) [1]
        Zs = 2 * Zi + 1j * w * Lext  # pylint: disable=invalid-name

        return Zs, Rhf, Lext

    def plot_series_resistance(self, f, z, rhf, data):
        """
        This function plots the series resistance as a function of frequency.

        Parameters:
        freq (array): Frequency array.
        zi_matrix (list of matrices): Matrix containing impedance values.
        p (int): Row index in the impedance matrix.
        q (int): Column index in the impedance matrix.
        """

        # Extracting the data from the list_data
        p = data[0]
        Np = data[1] # pylint: disable=invalid-name
        D = data[2][0,1] # pylint: disable=invalid-name

        # Extracting the impedance elements
        zs = np.array([item[p] for item in z[0]])
        rhf = np.array([item[p] for item in rhf])

        # Asymptotic Series Resistance
        plt.plot(f[0], 1E3 * rhf, label='Asymptotic',
                color='red', linestyle='--')

        # Closed-Form Approximation Series Resistance
        plt.plot(f[0], 1E3 * np.real(zs), label='Skin Effect Only',
                color='black', linestyle='-')

        # MoM Series Resistance
        plt.scatter(f[1], 1E3 * np.real(z[1]),
                    label=fr'MoM-SO: $Np=Nq={Np}$ [1]',
                    color='blue', marker='x', s=50)

        # Optional: Additional plotting configurations like labels, grid, etc.
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1E0, 1E7)
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Resistance p.u.l. (Ω/km)')
        plt.title('Figure 2.4: P.u.l. series resistance, $R_{int}$, of the bifilar overhead line\n'
                r'$r_1=r_2=0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},' fr'D={D}\,m,'
                r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
        plt.grid(True)
        plt.show()

    def plot_series_inductance(self, f, z, lext_hf, data):
        """ This function plots the series inductance as a function of frequency. """

        # Extracting the data from the list_data
        p = data[0]
        Np = data[1]  # pylint: disable=invalid-name
        D = data[2][0, 1]  # pylint: disable=invalid-name

        # Extracting the impedance elements
        zs = np.array([item[p] for item in z[0]])
        lext = np.array([np.imag(zs) / (2 * np.pi * f) for zs, f in zip(zs, f[0])])
        lext_hf = np.array([item[p] for item in lext_hf])
        lext_mom = np.array([np.imag(z) / (2 * np.pi * f)
                            for z, f in zip(z[1], f[1])])

        # Closed-Form Approximation Series Inductance
        plt.plot(f[0], 1E6 * lext, label='Analytical (Skin Effect Only)',
                color='black', linestyle='-')

        # Asymptotic Series Inductance
        plt.plot(f[0], 1E6 * lext_hf, label='Analytical (Asymptotic)',
                color='red', linestyle='--')

        # MoM Series Inductance
        plt.scatter(f[1], 1E6 * lext_mom, label=fr'MoM-SO ($Np=Nq={Np}$) [1]',
                    color='blue', marker='x', s=40)

        plt.xscale('log')
        plt.legend()
        plt.xlim(1E0, 1E7)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Inductance p.u.l. (mH/km)')
        plt.title('Figure 2.5: P.u.l. inductance of the bifilar overhead line\n'
                r'$r_1=r_2=0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},' fr'D={D}\,m,'
                r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
        plt.grid(False)
        plt.show()


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
