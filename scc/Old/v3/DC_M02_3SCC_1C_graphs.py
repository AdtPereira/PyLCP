import os
import numpy as np
import scipy.io as sio

from scc_data import scc_models_list
from scc_systm import SystemType
from scc_impedance import GroundReturnImpedance as gri
from scc_user_graph import DC_M02_3SCC_1C_gph1 as gph1, DC_M02_3SCC_1C_gph2 as gph2

def clear_screen():
    """Clears the console screen."""
    os.system('cls' if os.name == 'nt' else 'clear')

def load_matlab_data(file_path):
    """Loads MATLAB data from the .mat file."""
    mat_data = sio.loadmat(file_path)
    
    return mat_data['f']['A1'][0][0][0], mat_data['Z31']['A5'][0][0][0], \
           mat_data['Z31']['B5'][0][0][0], mat_data['Z31']['C5'][0][0][0]

def calculate_Zg_impedance(s_values, system_types):
    """Calculates ground return impedance for different system types."""
    
    # Define an empty dictionary "Impedance_data"
    Zg_data = {syst.id: [] for syst in system_types}
    
    for s in s_values:
        for syst in system_types:     
            Zg = gri(s, syst).DeConti(syst)
            Zg_data[syst.id].append(Zg[2][0])
    
    return Zg_data
    
def main():
    clear_screen()
    
    # Define system types and soil parameters
    scc = scc_models_list[1]
    system_types = [SystemType(Model=scc, rhog=100, erg=10, Syst_id='#1'),
                    SystemType(Model=scc, rhog=1000, erg=10, Syst_id='#2'),
                    SystemType(Model=scc, rhog=10000, erg=10, Syst_id='#3')]

    # Generate logarithmic frequency range
    dec_start = 4
    dec_end = 7
    num_points = 300
    f = np.logspace(dec_start, dec_end, num=num_points)
    s = 1j * (2 * np.pi * f)

    # Calculate ground return impedance list for all systems
    impedance_py = calculate_Zg_impedance(s, system_types)

    # Load MATLAB data
    impedance_m = load_matlab_data(r'C:\Users\adilt\OneDrive\01 ACADEMIA\05 PPGEE\EEE003 MODELOS\1.DECONTI\MatLabData.mat')

    # Plot impedance data
    gph1(f, impedance_py, impedance_m)
    gph2(f, impedance_py, impedance_m)

if __name__ == "__main__":
    main()