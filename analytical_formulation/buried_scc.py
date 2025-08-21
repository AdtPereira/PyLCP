import numpy as np
from scipy.special import kv
import scipy.constants as sc
from mtl_main.source import MulticonductorTransmissionLine

class PerUnitParameters:    
    """ This class calculates PUL parameters using an MTL geometry model. """

    def __init__(self, model: MulticonductorTransmissionLine, f: float, sigma_1: float, er_1: float = 1, mur_1: float = 1, ge: float = 0):
        # MTL Geometry Model
        self.mtl = model

        # Soil Relative Permittivity [np.array]
        self.er_1 = er_1

        # Soil conductivity (S/m) [np.array]
        self.sigma_1 = sigma_1

        # Soil resistivity (ohm.m)    
        self.rho_1 = 1 / sigma_1

        # Soil Relative Permeability [np.array]
        self.mur_1 = mur_1

        # External Conductance (S/m) [np.array]
        self.ge = ge

        # Angular frequency (rad/s)
        self.jw = 1j * 2 * np.pi * f            

            