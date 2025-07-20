from numpy import *
from scipy.constants import mu_0 as MUO, epsilon_0 as EO

class PRY_M01_1SCC_1C:

    def __init__(self):

        self.a = 0.012975               # Geometric Parameters
        self.b = 0.012975
        self.c = 0.012975
        self.d = 0.012975
        self.e = 0.012975
        self.f = 0.039315

        self.rho_cor = 1.934753E-8      # Core Resistivity (ohm.m)
        self.rho_sth = 0 #2.2E-7        # Sheath Resistivity (ohm.m)
        self.rho_arm = 0 #1.8E-7        # Armor Resistivity (ohm.m)

        self.er_eins = 2.963538 #2.3   # External Relative permittivity Insulation
        self.er_pins = 0               # Primary Relative permittivity Insulation
        self.er_sins = 0 #2.5          # Secondary Relative permittivity Insulation            
    
        self.pos_y = [1.7, 1.7, 1.7]    # cable Depth buried (positive numbers)
        self.pos_x = [0.0, 0.3, 0.6]    # Horizontal separation between the cable centers

class DC_M02_3SCC_1C:

    def __init__(self):

        self.a = 0.0190                 # Geometric Parameters
        self.b = 0.0425
        self.c = 0.0425
        self.d = 0.0425
        self.e = 0.0425
        self.f = 0.0425

        self.rho_cor = 1.7E-8           # Core Resistivity (ohm.m)
        self.rho_sth = 0                # Sheath Resistivity (ohm.m)
        self.rho_arm = 0                # Armor Resistivity (ohm.m)

        self.er_eins = 3.5             # External Relative permittivity Insulation
        self.er_pins = 0               # Primary Relative permittivity Insulation
        self.er_sins = 0               # Secondary Relative permittivity Insulation            

        self.pos_y = [1.00, 1.00, 1.00] # cable Depth buried (positive numbers)
        self.pos_x = [0.00, 0.35, 0.70] # Horizontal separation between the cable centers

class SystemType:

    def __init__(self, SystId, SCC, rhog, erg):
        
        # Ground resistivity (ohm.m)    
        self.syst_id = SystId   

        # Ground resistivity (ohm.m)    
        self.rho1 = rhog                
        
        # Ground relative permittivity
        self.eps_r1 = erg                

        # Ground permeability (H/m)
        self.mu1 = MUO                  

        # Ground Conductivity (S/m)
        self.sgm1 = 1/self.rho1         

        # Ground permittivity (F/m)
        self.eps1 = EO * self.eps_r1     

        # Number of internal conductors (core - sheath - armor)
        if SCC.rho_cor != 0 and SCC.rho_sth == 0 and SCC.rho_arm == 0:
            self.ncc = 1

        elif SCC.rho_cor != 0 and SCC.rho_sth != 0 and SCC.rho_arm == 0:
            self.ncc = 2

        else:
            self.ncc = 3

        # Number of single cables (single-phase - bi-phase - three-phase)
        self.nph = len(SCC.pos_y)

        # Number of total conductor (multi-phase transmission line)
        self.nc = self.ncc*self.nph