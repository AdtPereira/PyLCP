import datetime
import numpy as np
from scipy.constants import mu_0 as MUO, epsilon_0 as EO

class Result:
    @staticmethod
    def show_results(syst, Gri, Zg, Zi, A, ZL, Zs):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"{syst.Name}_{timestamp}.log"
        
        with open(filename, 'w') as Log_File:
            Log_File.write(                               
                # Writing system information
                f"-----------------------------------------------------------\n"
                f"Universidade Federal de Minas Gerais (UFMG)\n"
                f"Programa de Pós-Graduação em Engenharia Elétrica (PPGEE)\n\n"                

                f"Routine for calculating the Electrical Parameters of single\n" 
                f"core underground cables (SCC) \n"
                f"v3\n\n"
                                
                f"System Cable: {syst.Name}\n"
                f"Author: Adilton Pereira (c)\n"
                f"Belo Horizonte, {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                f"Brazil\n"
                f"-----------------------------------------------------------\n\n"
                                
                # Writing physical constants
                f"UNIVERSAL PHYSICAL CONSTANTS\n"
                f"  Vacuum magnetic permeability. mu_0 = {MUO:.4e} H/m\n"
                f"  Vacuum electrical permittivity. epsilon_0 = {EO:.4e} F/m\n\n"                
            )

            # Writing single core cables parameters
            Log_File.write("SINGLE CORE CABLES (SCC) PARAMETERS\n")
            for param, value in vars(syst).items():
                Log_File.write(f" {param} = {value}\n")
            Log_File.write("------------------------------------------------------------\n\n")
                      
            # Writing SCC internal impedance
            Log_File.write("SCC INTERNAL IMPEDANCE\n")            
            for param, value in Zi.items():
                if value != 0:
                    Log_File.write(f" {param} = {value:.4e} ohm/m\n")
            
            Log_File.write(
                # Writing matrix distances
                "\n"
                "ANGULAR FREQUENCY [Hz]\n"
                f"  s = {Gri.jwu/MUO} \n"
                f"  f = {Gri.jwu/MUO/1j/2/np.pi} \n"

                # Writing matrix distances
                "\nGLOBAL MATRIX DISTANCES [m]\n"
                f"  \nd = \n{Gri.d} \n"
                f"  \nD = \n{Gri.D} \n\n"
                
                # Writing physical parameters
                "PHYSICAL PARAMETERS\n"
                f"  Air Propagation constant. y0 = {Gri.y0:.4e}\n"
                f"  Ground Propagation constant. y1 = {Gri.y1:.4e}\n\n"
                
                # Writing De Conti et al. (2023) Closed-Form Expressions
                "GROUND RETURN IMPEDANCE MATRIX [ohm/m]\n"
                "De Conti et al. (2023) Closed-Form Expressions\n\n"
                f" zG = \n {Zg}\n\n"
                
                # Writing series impedance matrix (NODA, 2008)
                "SERIES IMPEDANCE MATRIX (NODA, 2008)\n"
                f" Transformation Matrix\n A = \n {A}\n\n"
                f" Loop impedance [ohm/m]\n ZL = \n {ZL}\n\n"
                f" Series Impedance [ohm/m]\n Zs = \n {Zs}\n\n"
            )