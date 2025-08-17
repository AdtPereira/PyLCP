"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

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

[4] PAUL, Clayton R. Analysis of multiconductor transmission lines. 2. ed. Hoboken,
    N.J.: John Wiley & Sons, Inc., c2008.

"""

import numpy as np
from scipy.special import iv, kv
from scipy.constants import mu_0
from mtl_data.mtl import MulticonductorTransmissionLine


class CoaxialCable(MulticonductorTransmissionLine):
    """ This class contains the analytical formulation of the system. """

    def __init__(self, mtl, frequency):
        """
        Initialize the AnalyticalFormulation class.

        Parameters:
        conductor (list): List of dictionaries containing the properties of the conductors.
        frequency (float): The frequency of the system.
        """
        super().__init__(mtl)

        # Angular frequency, rad/s [float]
        self.jw = 1j * 2 * np.pi * frequency

        # Propagation Constant [np.array]
        self.gamma = np.sqrt(self.jw * mu_0 * self.sigma)

        # Coaxial Cable radii
        for key, conductor in self.mtl.items():
            if isinstance(key, int):  # Ensures the key is an integer
                if conductor['conductor_name'] == 'core':
                    self.a = conductor['radius'][1]
                elif conductor['conductor_name'] == 'sheath':
                    self.b, self.c = conductor['radius']

    # Equation 2.70 [1]
    # Equation 4.51 [4]
    def external_inductance(self):
        """ Calculate the external inductance for a lossless coaxial cable, L'. """
        return mu_0 / (2 * np.pi) * np.log(self.b / self.a)

    # Equation 2.71 [1]
    def internal_impedance(self):
        """ Calculate the internal impedance of the inner conductor Za (omega). """
        # Propagation Constant of the inner conductor
        gama_a = self.gamma[0] * self.a

        # Intrinsic Impedance of the inner conductor
        eta = np.sqrt(self.jw * mu_0 / self.sigma)[0]

        # Internal Impedance of the inner conductor
        za = eta / (2 * np.pi * self.a) * iv(0, gama_a) / iv(1, gama_a)

        return za

    # Equation 2.72 [1]
    def external_impedance(self):
        """ Calculate the external impedance of the inner conductor Zb (omega). """
        # Propagation Constant of the inner conductor
        gama_b = self.gamma[0] * self.b
        gama_c = self.gamma[0] * self.c

        # Intrinsic Impedance of the inner conductor
        eta = np.sqrt(self.jw * mu_0 / self.sigma)[0]

        numerator = iv(0, gama_b) * kv(1, gama_c) + (kv(0, gama_b) * iv(1, gama_c))
        denominator = iv(1, gama_c) * kv(1, gama_b) - (iv(1, gama_b) * kv(1, gama_c))

        # External Impedance of the inner conductor
        zb = eta / (2 * np.pi * self.b) * numerator / denominator

        return zb

    # Equation (2.69) [1]
    def pul_parameters(self):
        """
        This function calculates the series resistance of the system using 
        the high frequency approximation.

        Returns:
            tuple: A tuple containing the high frequency resistance, external inductance, 
            and matrix impedance.
        """
        # Matrix Impedance, z (Ω/m)
        l_ext = self.external_inductance()
        za = self.internal_impedance()
        zb = self.external_impedance()
        zs = self.jw * l_ext + za + zb

        return zs


class Ametani(MulticonductorTransmissionLine):
    """ This class contains the analytical formulation of the system. """

    def __init__(self, mtl, frequency):
        """
        Initialize the AnalyticalFormulation class.

        Parameters:
        conductor (list): List of dictionaries containing the properties of the conductors.
        frequency (float): The frequency of the system.
        """
        super().__init__(mtl)

        # Angular frequency, rad/s [float]
        self.jw = 1j * 2 * np.pi * frequency

        # Call the function to configure the parameters
        self._scc_from_data()

    def _scc_from_data(self):
        """ This function configures the parameters of the coaxial cable.
        """

        for key, conductor in self.mtl.items():
            if isinstance(key, int):  # Ensures the key is an integer
                if conductor['conductor_name'] == 'core':
                    self.a, self.b = conductor['radius']
                elif conductor['conductor_name'] == 'sheath':
                    self.b_prime, self.c = conductor['radius']

    def parameter_m(self, mu, sigma):
        """
        Calculate the parameter m for the two-layered conductor.
        """
        return np.sqrt(self.jw * mu * sigma)

    def impedance_two_layered_conductor(self):
        """ 
        Calculate the impedance of a two-layered conductor.
        """

        # Intermediate variables
        m1 = self.parameter_m(self.mu[0], self.sigma[0])
        m2 = self.parameter_m(self.mu[1], self.sigma[1])
        x1 = m1 * self.a
        x2 = m1 * self.b
        x3 = m2 * self.b_prime
        x4 = m2 * self.c

        # Intermediate variables A, B, E, F, R
        aa = kv(1, x1) * iv(0, x2) + iv(1, x1) * kv(0, x2)
        bb = kv(1, x1) * iv(1, x2) - iv(1, x1) * kv(1, x2)
        ee = iv(0, x3) * kv(1, x4) + kv(0, x3) * iv(1, x4)
        ff = kv(1, x3) * iv(1, x4) - iv(1, x3) * kv(1, x4)
        rr = kv(1, x3) * iv(0, x4) + iv(1, x3) * kv(0, x4)

        # Calculating z10, z2i, z2m, z20, z12
        rho1 = 1 / self.sigma[0]
        rho2 = 1 / self.sigma[1]

        # Solid conductor case
        if x1 == 0:
            z10 = (m1 * rho1 / (2 * np.pi * self.b)) * iv(0, x2) / iv(1, x2)
        else:
            z10 = (m1 * rho1 / (2 * np.pi * self.b)) * aa / bb

        z2i = (m2 * rho2 / (2 * np.pi * self.b_prime)) * ee / ff
        z2m = rho2 / (2 * np.pi * self.b_prime * self.c * ff)
        z20 = (m2 * rho2 / (2 * np.pi * self.c)) * rr / ff
        z12 = self.jw * (mu_0 / 2 / np.pi) * np.log(self.b_prime / self.b)

        # Calculating Z11, Z12, Z22
        zz22 = z20
        zz12 = zz22 - z2m
        zz11 = z10 + z12 + z2i - z2m + zz12

        return zz11, zz12, zz22
