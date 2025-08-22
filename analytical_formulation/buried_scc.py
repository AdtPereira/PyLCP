"""
REFERENCES:
[1] XUE, Haoyan. General Formulation and Accurate Evaluation of Earth-Return Parameters
    for Overhead / Underground Cables. PhD thesis, Department of Electrical Engineering,
    École Polytechnique de Montréal, Université de Montréal, August 2018. 

[2] PAUL, Clayton R. Analysis of multiconductor transmission lines. 2. ed. Hoboken,
    N.J.: John Wiley & Sons, Inc., c2008.
"""

import numpy as np
from scipy.special import kv
import scipy.constants as sc
from mtl_main.source import MulticonductorTransmissionLine

class PerUnitParameters:    
    """ This class calculates PUL parameters using an MTL geometry model. """

    def __init__(self, model: MulticonductorTransmissionLine, f: float):
        # MTL Geometry Model
        self.mtl = model

        # Soil Relative Permittivity
        self.er_1 = model.mtl_ref[0]['relative_permittivity']

        # Soil conductivity (S/m)
        self.sigma_1 = model.mtl_ref[0]['conductivity']

        # Soil Relative Permeability
        self.mur_1 = model.mtl_ref[0]['relative_permeability']

        # External Conductance (S/m)
        self.ge = model.mtl_ref[0]['external_conductance']

        # Soil resistivity (ohm.m)
        self.rho_1 = 1 / self.sigma_1

        # Angular frequency (rad/s)
        self.jw = 1j * 2 * np.pi * f

        # Constant vacuum terms
        self.jw_mu0_2pi = self.jw * sc.mu_0 / (2 * np.pi)
        self.jw_2pi_e0 = self.jw * 2 * np.pi * sc.epsilon_0

        # Air wave number - Equation (2.15) [1]
        self.k_air2 = - self.jw * sc.mu_0 * self.jw * sc.epsilon_0

        # Earth wave number - Equation (2.15) [1]
        self.k_earth2 = - self.jw * self.mur_1 * sc.mu_0 * (self.sigma_1 + self.jw * self.er_1 * sc.epsilon_0)

    # 3.3.2.2 Earth-return impedance and admittance formulas based on quasi-TEM assumption [1]
