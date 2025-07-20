import numpy as np
from scipy.constants import mu_0 as MUO
from scc_user_functions import MatrixDef as md
from scc_user_functions import Trigonometric as tg

class InternalImpedanceAprox:

    def __init__(self, jw, scc):

        # Core Resistivity
        rho_cor = scc.rho_cor        
        a = scc.a

        # Core impedance
        if rho_cor != 0:
            m = np.sqrt(jw * MUO / rho_cor)
            coth_val = tg.coth(0.777 * m * a)
            self.Zc = rho_cor * m / (2 * np.pi * a) * coth_val + 0.356 * rho_cor / (np.pi * (a ** 2))
        else:
            self.Zc = 0

        # Sheath Resistivity
        rho_sth = scc.rho_sth
        b, c = scc.b, scc.c

        # Sheath impedance
        if rho_sth != 0: 
            m = np.sqrt(jw * MUO / rho_sth)
            tck = c - b
            coth_val = tg.coth(m * tck)
            self.Zsi = rho_sth * m / (2 * np.pi * b) * coth_val - rho_sth / (2 * np.pi * b * (b + c))
            self.Zso = rho_sth * m / (2 * np.pi * c) * coth_val + rho_sth / (2 * np.pi * c * (b + c))
            self.Zsm = rho_sth * m / (np.pi * (b + c)) * tg.csch(m * tck)
        else:
            self.Zsi = self.Zso = self.Zsm = 0

        # Armor Resistivity
        rho_arm = scc.rho_arm
        d, e = scc.d, scc.e

        # Armor Impedance
        if rho_arm != 0:
            m = np.sqrt(jw * MUO / rho_arm)
            tck = e - d
            coth_val = tg.coth(m * tck)
            self.Zai = rho_arm * m / (2 * np.pi * d) * coth_val - rho_arm / (2 * np.pi * d * (d + e))
            self.Zao = rho_arm * m / (2 * np.pi * e) * coth_val + rho_arm / (2 * np.pi * e * (d + e))
            self.Zam = rho_arm * m / (np.pi * (d + e)) * tg.csch(m * tck)
        else:
            self.Zai = self.Zao = self.Zam = 0
     
        # Primary insulation impedance between core and sheath
        self.Zpins = jw * MUO / (2 * np.pi) * np.log(b / a)
        
        # Secondary insulation impedance between sheath and armor
        self.Zsins = jw * MUO / (2 * np.pi) * np.log(d / c)

        # External insulation impedance
        self.Zeins = jw * MUO / (2 * np.pi) * np.log(scc.f / e)

class SeriesImpedanceMatrix:

    def LoopImpedance(self, zI, Zg, syst):
        if syst.ncc == 1:
            return zI.Zc + zI.Zeins + Zg[0][0]
        else:
            return 0

    def MatrixTransformation(self, zI, Zg, syst):
        n = syst.ncc
        nz = np.zeros((n, n), dtype=complex)
        a = np.eye(n) + np.diag(-1 * np.ones(n - 1), -1)
        zL = self.LoopImpedance(zI, Zg, syst)

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
    
    def SeriesImpedance(self, zI, Zg, syst):
        A, ZL = self.MatrixTransformation(zI, Zg, syst)
        Zs = md.RInv(md.LInv(A.T, ZL), A) 
        return Zs