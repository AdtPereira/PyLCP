import numpy as np
from scipy.special import kv
from scipy.constants import mu_0 as MUO, epsilon_0 as EO
from scc.scc_functions import MatrixOperation as mo, TrigonometricOperation as trig, MatrixDistances as md
from mtl_main.source import MulticonductorTransmissionLine

class InternalParameters:

    def __init__(self, jw):
        self.jw = jw

    def CoreImpedanceAprox(self, r_out, rho):
        
        # Internal parameters
        m = np.sqrt(self.jw * MUO / rho)
        coth_val = trig.coth(0.777 * m * r_out)
        
        # Core Impedance
        return (rho*m/(2*np.pi*r_out) * coth_val + 0.356*rho/(np.pi*(r_out ** 2)))
  
    def InternalImpedanceAprox(self, r_in, r_out, rho):

        # Internal parameters
        m = np.sqrt(self.jw * MUO / rho)
        tck = r_out - r_in
        coth_val = trig.coth(m * tck)
        Rio = r_in + r_out

        # Internal surface Conductor (core - Sheath) Impedance 
        Zxi = rho*m / (2*np.pi * r_in)  * coth_val - rho/(2*np.pi * r_in * Rio)
        
        # Outer surface Conductor (core - Sheath) Impedance  
        Zxo = rho*m / (2*np.pi * r_out) * coth_val + rho/(2*np.pi * r_out * Rio)
        
        # Mutual Conductor (core - Sheath) Impedance 
        Zxm = rho*m / (np.pi * Rio) * trig.csch(m * tck)

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
    
class GroundReturnAdmittance(GroundReturnImpedance):

    def __init__(self, jw, syst):  
        super().__init__(jw, syst) 
        self.P1 = jw/2/np.pi/(syst.sgm1 + jw*syst.eps1)  
        self.jw = jw   

    def LowFrequency(self, syst): 
        if syst.ncc == 1:
            Ye = np.array([[self.jw * EO * 2 * np.pi * syst.er_ei / np.log(syst.rf/syst.re)]])            

        elif syst.ncc == 2:
            Y1 = self.jw * EO * 2 * np.pi * syst.er_pi / np.log(syst.rb/syst.ra)
            Y2 = self.jw * EO * 2 * np.pi * syst.er_ei / np.log(syst.rf/syst.re)
            Ye = np.array([[ Y1,   -Y1], \
                          [-Y1, Y1+Y2]])  

        elif syst.ncc == 3:  
            Y1 = self.jw * EO * 2 * np.pi * syst.er_pi / np.log(syst.rb/syst.ra)
            Y2 = self.jw * EO * 2 * np.pi * syst.er_si / np.log(syst.rd/syst.rc)
            Y3 = self.jw * EO * 2 * np.pi * syst.er_ei / np.log(syst.rf/syst.re)
            Ye = np.array([[ Y1,   -Y1,     0], \
                          [-Y1, Y1+Y2,   -Y2], \
                          [  0,   -Y2, Y2+Y3]])
        return Ye

    def DeConti(self, syst):        
        # Compact Internal variables
        y1 = self.y1
        y0 = self.y0
        k = syst.nph

        # Ground-Return Coefficients
        Pg = np.zeros((k, k), dtype='complex_')
        U = np.eye(k)

        for m in range(k):
            for n in range(k):
                               
                Pg[m][n] = self.P1 * (kv(0, y1 * self.d[m][n]) +
                    (y1**2 - y0**2) / (y1**2 + y0**2) * (kv(0, y1 * self.D[m][n]))
                )

        Yg = self.jw * mo.LInv(Pg, U)

        return Pg, Yg 

class PerUnitParameters:    
    def __init__(self, model: MulticonductorTransmissionLine, s: float):
        # MTL Geometry Model
        self.mtl = model

        self.s = s

        # scc model
        self.ra = Model['ra']['value']
        self.rb = Model['rb']['value']
        self.rc = Model['rc']['value']
        self.rd = Model['rd']['value']
        self.re = Model['re']['value']
        self.rf = Model['rf']['value']  
        self.rho_c = Model['rho_cor']['value']
        self.rho_s = Model['rho_sth']['value']
        self.rho_a = Model['rho_arm']['value']
        self.er_ei = Model['er_eins']['value']
        self.er_pi = Model['er_pins']['value']
        self.er_si = Model['er_sins']['value']
        self.pos_x = Model['pos_x']['value']
        self.pos_y = Model['pos_y']['value']        

        # Soil resistivity (ohm.m)    
        self.rho1 = rhog                
        
        # Soil relative permittivity
        self.eps_r1 = erg                

        # Soil permeability (H/m)
        self.mu1 = MUO                  

        # Soil Conductivity (S/m)
        self.sgm1 = 1/rhog         

        # Soil permittivity (F/m)
        self.eps1 = EO * erg    

        # Number of internal conductors (core - sheath - armor)
        if self.rho_c != 0 and self.rho_s == 0 and self.rho_a == 0:
            self.ncc = 1

        elif self.rho_c != 0 and self.rho_s != 0 and self.rho_a == 0:
            self.ncc = 2

        elif self.rho_c != 0 and self.rho_s != 0 and self.rho_a != 0:
            self.ncc = 3

        else:
            self.ncc = 0

        # Number of single cables (single-phase - bi-phase - three-phase)
        self.nph = len(self.pos_y)

        # Number of total conductor (multi-phase transmission line)
        self.nc = self.ncc*self.nph


    def SeriesImpedance(self):
        # Ground Return Impedance
        Zg = GroundReturnImpedance(self.s, self.syst).DeConti(self.syst)
        
    #     # Loop Impedance Matrix (NODA,2008)
    #     zL = SeriesImpedanceMatrix(self.s).LoopImpedance(self.syst, Zg[0][0])
    #     A, ZL = SeriesImpedanceMatrix(self.s).MatrixTransformation(self.syst, zL, Zg)
        
    #     # Series Impedance Matrix
    #     self.Zs = SeriesImpedanceMatrix(self.s).SeriesImpedance(A, ZL)
            
    # def ShuntAdmittance(self):
    #     self.Ye = GroundReturnAdmittance(self.s, self.syst).LowFrequency(self.syst)
            
    # def PropagationFunction(self):
    #     self.gama = np.sqrt(self.Zs * self.Ye)  
    #     self.Yc = np.sqrt(self.Ye / self.Zs)

    # def QuadripoleParameters(self, RS, LX):
    #     self.Ykk = self.Yc * trig.coth(self.gama*LX)
    #     self.Ykm = -self.Yc * trig.csch(self.gama*LX)
    #     L1 = np.hstack((self.Ykk + 1/RS, self.Ykm))
    #     L2 = np.hstack((self.Ykm       , self.Ykk))
    #     return np.vstack((L1, L2))