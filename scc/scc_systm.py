from scipy.constants import mu_0 as MUO, epsilon_0 as EO

class SystemType:

    def __init__(self, Model, rhog, erg, Syst_id):
        
        # scc model
        self.id = Syst_id
        self.Name = Model['name']['value']
        self.Description = Model['name']['description']
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