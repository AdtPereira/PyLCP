import numpy as np
from scipy.constants import mu_0 as MUO, epsilon_0 as EO
from scipy.special import kv  # Modified second-order Bessel function
from scc_user_functions import MatrixDistances

class GroundReturnImpedance:

    def __init__(self, jw, scc, syst):        
        # Pre-calculated constants
        self.jwu = jw * MUO
        self.y0 = jw * np.sqrt(MUO * EO)
        self.y1 = np.sqrt(jw * syst.mu1 * (syst.sgm1 + jw * syst.eps1))
        
        # System attributes
        self.nph = syst.nph
        self.h = scc.pos_y   
        self.x = scc.pos_x 
        self.d, self.D = MatrixDistances.Distance(scc, syst)
        
    def DeConti(self):
        # Compact Internal variables
        y1 = self.y1
        y0 = self.y0
        
        # Ground Return Impedance
        Zg = np.zeros((self.nph, self.nph), dtype='complex_')

        for m in range(self.nph):
            for n in range(self.nph):
                r = np.sqrt((self.x[m] - self.x[n])**2 + (self.h[m] - self.h[n])**2)
                
                Zg[m][n] = self.jwu / (2 * np.pi) * (
                    kv(0, y1 * self.d[m][n]) +
                    (y1 - y0) / (y0 + y1) * np.exp(-y1 * (self.h[m] + self.h[n])) * (2 / (4 + (y1 * r)**2))
                )
        
        return Zg
