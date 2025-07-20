import numpy as np
from scipy.constants import mu_0 as MUO, epsilon_0 as EO
from scc_user_functions import MatrixOper as mo
from scc_user_functions import Trigonometric as tg
from scc_user_functions import MatrixDistances as md

# Modified second-order Bessel function
from scipy.special import kv

class InternalParameters:

    def __init__(self, jw):
        self.jw = jw

    def CoreImpedanceAprox(self, r_out, rho):
        
        # Internal parameters
        m = np.sqrt(self.jw * MUO / rho)
        coth_val = tg.coth(0.777 * m * r_out)
        
        # Core Impedance
        return (rho*m/(2*np.pi*r_out) * coth_val + 0.356*rho/(np.pi*(r_out ** 2)))
  
    def InternalImpedanceAprox(self, r_in, r_out, rho):

        # Internal parameters
        m = np.sqrt(self.jw * MUO / rho)
        tck = r_out - r_in
        coth_val = tg.coth(m * tck)
        Rio = r_in + r_out

        # Internal surface Conductor (core - Sheath) Impedance 
        Zxi = rho*m / (2*np.pi * r_in)  * coth_val - rho/(2*np.pi * r_in * Rio)
        
        # Outer surface Conductor (core - Sheath) Impedance  
        Zxo = rho*m / (2*np.pi * r_out) * coth_val + rho/(2*np.pi * r_out * Rio)
        
        # Mutual Conductor (core - Sheath) Impedance 
        Zxm = rho*m / (np.pi * Rio) * tg.csch(m * tck)

        return Zxi, Zxo, Zxm

    def InsulationImpedance(self, r_out, r_in):        
        return (self.jw * MUO / (2*np.pi) * np.log(r_out / r_in))   

class SeriesImpedanceMatrix(InternalParameters):
    
    def __init__(self, jw):
        super().__init__(jw)          

    def LoopImpedance(self, syst, Zg):
        if syst.ncc == 1: 
            Zc = self.CoreImpedanceAprox(syst.ra, syst.rho_c)
            Zeins = self.InsulationImpedance(syst.rf, syst.re)
            zL = Zc + Zeins + Zg
            return zL
        
        else:
            return 0

    def MatrixTransformation(self, syst, zL, Zg):
        n = syst.ncc
        nz = np.zeros((n, n), dtype=complex)
        a = np.eye(n) + np.diag(-1 * np.ones(n - 1), -1)
        
        if syst.nph == 1:
            A = a
            ZL = zL

        elif syst.nph == 2:
            A = np.hstack((np.vstack((a, nz)), np.vstack((nz, a))))
            zL12 = np.zeros((n, n), dtype=complex)
            zL12[-1][-1] = Zg[0][1]
            ZL = np.hstack((np.vstack((zL, zL12)), np.vstack((zL12, zL))))

        elif syst.nph == 3:
            A = np.hstack((np.vstack((a, nz, nz)), np.vstack((nz, a, nz)), np.vstack((nz, nz, a))))
            zL12 = np.zeros((n, n), dtype=complex)
            zL12[-1][-1] = Zg[0][1]
            zL13 = np.zeros((n, n), dtype=complex)
            zL13[-1][-1] = Zg[0][2]
            zL23 = np.zeros((n, n), dtype=complex)
            zL23[-1][-1] = Zg[1][2]
            ZL = np.hstack((np.vstack((zL, zL12, zL13)), np.vstack((zL12, zL, zL23)), np.vstack((zL13, zL23, zL))))
        return A, ZL
    
    def SeriesImpedance(self, A, ZL):        
        Zs = mo.RInv(mo.LInv(A.T, ZL), A) 
        return Zs
    
class GroundReturnImpedance:

    def __init__(self, jw, syst):        
        # Pre-calculated constants
        self.jwu = jw * MUO
        self.y0 = jw * np.sqrt(MUO * EO)
        self.y1 = np.sqrt(jw * syst.mu1 * (syst.sgm1 + jw * syst.eps1))
        
        # System attributes
        self.d, self.D = md.Distance(syst)
        
    def DeConti(self, syst):
        # Compact Internal variables
        y1 = self.y1
        y0 = self.y0
        k = syst.nph
        x = syst.pos_x
        h = syst.pos_y
        
        # Ground Return Impedance
        Zg = np.zeros((k, k), dtype='complex_')

        for m in range(k):
            for n in range(k):
                # Horizontal separation between the cable centers
                r = np.sqrt((x[m] - x[n])**2 + (h[m] - h[n])**2)
                
                Zg[m][n] = self.jwu / (2 * np.pi) * (
                    kv(0, y1 * self.d[m][n]) +
                    (y1 - y0) / (y0 + y1) * np.exp(-y1 * (h[m] + h[n])) * (2 / (4 + (y1 * r)**2))
                )
        
        return Zg    