import datetime
import numpy as np
from scipy.constants import mu_0 as MUO, epsilon_0 as EO

class LogResult:
    def Parameters(syst, Gri, Zg, Zi, zL, A, ZL, Zs,  Ye, Pg, Yg, file_path):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"{file_path}{syst.Name}_{timestamp}.log"
        
        with open(filename, 'w') as Log_File:
            Log_File.write(                               
                # Writing system information
                f"==================================================================================\n"
                f"Universidade Federal de Minas Gerais (UFMG)\n"
                f"Programa de Pós-Graduação em Engenharia Elétrica (PPGEE)\n\n"                

                f"Routine for calculating the Electrical Parameters of single core underground\n" 
                f"cables (SCC) arrangements.\n"
                f"v5\n\n"
                                
                f"Author: Adilton Pereira (c)\n"
                f"Belo Horizonte, {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                f"Brazil"

                f"\n\n"
                f"References:\n"
                f"[1] A. De Conti, N. Duarte and R. Alipio, 'Closed-Form Expressions for the \n" 
                f"    Calculation of the Ground-Return Impedance and Admittance of Underground\n"
                f"    Cables,' in IEEE Transactions on Power Delivery, vol. 38, no. 4,\n"
                f"    pp. 2891-2900, Aug. 2023, doi: 10.1109/TPWRD.2023.3264614.\n"
                f"==================================================================================\n\n"

                f"SYSTEM CABLE NAME: {syst.Name}\n"
                f"{syst.Description}\n\n"

                # Writing physical constants
                f"UNIVERSAL PHYSICAL CONSTANTS\n"
                f"  Vacuum magnetic permeability. mu_0 = {MUO:.4e} H/m\n"
                f"  Vacuum electrical permittivity. epsilon_0 = {EO:.4e} F/m\n\n"                
            )

            # Writing single core cables parameters
            Log_File.write("SINGLE CORE CABLES (SCC) AND SYSTEM PARAMETERS\n")
            for param, value in vars(syst).items():
                Log_File.write(f" {param} = {value}\n")
            Log_File.write("\n")
                      
            # Writing SCC internal impedance
            Log_File.write("SCC INTERNAL IMPEDANCE\n")            
            for param, value in Zi.items():
                if value != 0:
                    Log_File.write(f" {param} = {value:.4e} ohm/m\n")
            
            Log_File.write(
                # LOOP IMPEDANCE PARAMETER [ohm/m] (NODA, 2008)
                f" Phase-Loop Impedance. zL = {zL} \n"

                # Writing matrix distances
                "\n"
                "ANGULAR FREQUENCY\n"
                f"  jw = {Gri.jwu/MUO:.2e} rad/s\n"
                f"  f = {Gri.jwu/MUO/1j/2/np.pi:.2e} [Hz]\n"

                # Writing matrix distances
                "\nGLOBAL MATRIX DISTANCES [m]"
                f"  \nd = \n{Gri.d} \n"
                f"  \nD = \n{Gri.D} \n\n"
                
                # Writing physical parameters
                "PHYSICAL PARAMETERS\n"
                f"  Air Propagation constant. y0 = {Gri.y0:.4e}\n"
                f"  Ground Propagation constant. y1 = {Gri.y1:.4e}\n\n"
                
                # Writing De Conti et al. (2023) Closed-Form Expressions
                "GROUND RETURN IMPEDANCE MATRIX [1] [ohm/m]\n"
                f" zG = \n {Zg}\n\n"
                
                # Writing series impedance matrix (NODA, 2008)
                "SERIES IMPEDANCE MATRIX (NODA, 2008)\n"
                f" Transformation Matrix\n A = \n {A}\n\n"
                f" Loop impedance Matrix [ohm/m]\n ZL = \n {ZL}\n\n"
                f" Series Impedance Matrix [ohm/m]\n Zs = \n {Zs}\n\n"

                # Writing low frequency admittance 
                "LOW FREQUENCY ADMITTANCE MATRIX\n"
                f" Ye = \n {Ye}\n\n"                

                # Writing ground-return potential coefficients and admittance matrices
                "GROUND-RETURN ADMITTANCE MATRIX\n"
                "Ground-return potential coefficients\n"
                f" Pg = \n {Pg}\n\n"
                "Ground-return admittance\n"
                f" Yg = \n {Yg}\n\n"
            )
    
    def Transient(syst, s, Zs, Ye, file_path):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"{file_path}{syst.Name}_{timestamp}.log"
        
        with open(filename, 'w') as Log_File:
            ## Basic Information
            Log_File.write(                               
                # Writing system information
                f"==================================================================================\n"
                f"Universidade Federal de Minas Gerais (UFMG)\n"
                f"Programa de Pós-Graduação em Engenharia Elétrica (PPGEE)\n\n"                

                f"Transient in Transmission Lines of\n" 
                f"single core underground cables (SCC) arrangements.\n"
                f"v5\n\n"
                                
                f"Author: Adilton Pereira (c)\n"
                f"Belo Horizonte, {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                f"Brazil"

                f"\n\n"
                f"References:\n"
                f"[1] A. De Conti, N. Duarte and R. Alipio, 'Closed-Form Expressions for the \n" 
                f"    Calculation of the Ground-Return Impedance and Admittance of Underground\n"
                f"    Cables,' in IEEE Transactions on Power Delivery, vol. 38, no. 4,\n"
                f"    pp. 2891-2900, Aug. 2023, doi: 10.1109/TPWRD.2023.3264614.\n"
                f"==================================================================================\n\n"

                f"SYSTEM CABLE NAME: {syst.Name}\n"
                f"{syst.Description}\n\n"

                # Writing physical constants
                f"UNIVERSAL PHYSICAL CONSTANTS\n"
                f"  Vacuum magnetic permeability. mu_0 = {MUO:.4e} H/m\n"
                f"  Vacuum electrical permittivity. epsilon_0 = {EO:.4e} F/m\n\n"                
            )

            ## Per-Unit Transmission Line Matrices          
            Log_File.write(
                # FREQUENCY
                f" s = {s} \n"

                # SERIES IMPEDANCE MATRIX [ohm/m]
                f" Zs \n {Zs} \n"

                # SHUNT ADMITTANCE MATRIX [ohm/m]
                f" Ye \n {Ye} \n"                
            )