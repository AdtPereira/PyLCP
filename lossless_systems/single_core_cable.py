"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

Below is a high-level overview of the script components:

Imports and Global Variables:

Required libraries and global variables are imported and defined.

BIFILAR_TL: A list containing properties of the two conductors.
Classes:

Geometry: Handles basic geometry calculations, such as distance matrices between conductor centers.
ParametersWithFrequency: Extends Geometry to include frequency-dependent parameters.
GreensMatrices: Uses the geometry to compute Green's matrices, which are essential for the method 
of moments.
MoMSuperficialOperator: Implements the method of moments, calculating matrices like U, Ys, G, and Z.
AnalyticalFormulation: Provides analytical formulations for high-frequency resistance, external 
inductance, and impedance.
Plotter: Handles plotting of series resistance and inductance against frequency.

Functions:

clear_screen: Clears the console screen.
main: The main function orchestrates the scattering calculations and plotting. It performs the 
following steps:
Clears the screen.
Initializes objects for the method of moments and analytical formulations.
Computes series resistance, external inductance, and impedance over a range of frequencies.
Plots the results using the Plotter class.
Detailed Class and Function Explanations
Geometry
__init__: Initializes the geometry of the system based on conductor properties.
distance_matrices: Calculates matrices for distances and angles between conductor centers.
ParametersWithFrequency
__init__: Extends the Geometry class to include frequency-dependent parameters such as 
conductivity, permeability, and permittivity.
ynp_operator: Calculates the surface admittance operator for a conductor.
GreensMatrices
__init__: Initializes Green's matrices using the geometry of the system.
dissertation and ieee_paper: Calculate Green's functions using different methods.
sub_matrices: Generates Green's sub-matrices.
MoMSuperficialOperator
__init__: Extends ParametersWithFrequency to initialize the method of moments parameters.
matrix_u: Constructs matrix U.
matrix_ys: Constructs matrix Ys.
matrix_g and matrix_g_ieee: Constructs matrix G using different methods.
matrix_z: Computes the impedance matrix Z.
AnalyticalFormulation
__init__: Initializes the analytical formulation based on the conductor properties and frequency.
pul_parameters: Calculates high-frequency resistance, external inductance, and impedance.
Plotter
__init__: Initializes the plotting class with frequency and impedance data.
series_resistance: Plots series resistance against frequency.
series_inductance: Plots series inductance against frequency.

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

"""

import numpy as np
from scipy.special import iv, kv
from scipy.constants import mu_0

from mom_so.mtl import MulticonductorTransmissionLine


class SingleCoreCable(MulticonductorTransmissionLine):
    """ This class contains the analytical formulation of the system. """

    def __init__(self, mtl, frequency):
        """
        Initialize the AnalyticalFormulation class.

        Parameters:
        conductor (list): List of dictionaries containing the properties of the conductors.
        frequency (float): The frequency of the system.
        """
        super().__init__(mtl)

        # Conductivity of the conductors [np.array]
        self.sigma = np.array([c['conductivity'] for c in self.mtl])

        # Angular frequency, rad/s [float]
        self.jw = 1j * 2 * np.pi * frequency

        # Propagation Constant [np.array]
        self.gamma = np.sqrt(self.jw * mu_0 * self.sigma)

        # Coaxial Cable radii
        core = [c for c in self.mtl if c['conductor_name'] == 'core']
        sheath = [c for c in self.mtl if c['conductor_name'] == 'sheath']

        self.a = core[0]['radius'][1] if core else None
        self.b, self.c = (sheath[0]['radius'] if sheath else (None, None))

    # Equation 2.70 [1]
    # Equation 4.50 [2]
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
        """ 
        This function configures the parameters of the coaxial cable.
        """

        for conductor in self.mtl:
            if conductor['conductor_name'] == 'core':
                self.a = conductor['radius'][0]
                self.b = conductor['radius'][1]

            elif conductor['conductor_name'] == 'sheath':
                self.b_prime = conductor['radius'][0]
                self.c = conductor['radius'][1]

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
