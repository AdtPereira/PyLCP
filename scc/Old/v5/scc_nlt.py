import numpy as np
from scipy.constants import mu_0 as MUO, epsilon_0 as EO
from scc_functions import MatrixOperation as mo, TrigonometricOperation as tg, MatrixDistances as md

class MonoNetworkTopology():

    def __init__(self, Z, Y):
        self.gama = np.sqrt(Z*Y)
        self.Yc = np.sqrt(Y/Z)   