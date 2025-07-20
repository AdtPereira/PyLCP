import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

from scc_data import DC_M02_3SCC_1C as model, SystemType
from scc_impedance_ground_return import GroundReturnImpedance

def clear_screen():
    """Clears the console screen."""
    os.system('cls' if os.name == 'nt' else 'clear')

def load_matlab_data(file_path):
    """Loads MATLAB data from the .mat file."""
    mat_data = sio.loadmat(file_path)
    return mat_data['f']['A1'][0][0][0], mat_data['Z31']['A5'][0][0][0], \
           mat_data['Z31']['B5'][0][0][0], mat_data['Z31']['C5'][0][0][0]

def calculate_impedance(s_values, scc, system_types):
    """Calculates ground return impedance for different system types."""
    impedance_data = {system.syst_id: [] for system in system_types}
    for s in s_values:
        for system_type in system_types:
            impedance = GroundReturnImpedance(s, scc, system_type).DeConti()
            impedance_data[system_type.syst_id].append(impedance[2][0])
    return impedance_data

def plot_impedance_abs(frequencies, impedance_data, matlab_data=None):
    """Plots impedance data."""
    plt.figure(figsize=(8, 6))
    for system_id, impedance_values in impedance_data.items():
        plt.semilogx(frequencies, np.abs(impedance_values), label=f'System {system_id} Python Data')
    if matlab_data:
        mat_f, mat_Z31A, mat_Z31B, mat_Z31C = matlab_data
        plt.semilogx(mat_f, np.abs(mat_Z31A), color='black', linestyle='--', label='MATLAB Data')
        plt.semilogx(mat_f, np.abs(mat_Z31B), color='black', linestyle='--')
        plt.semilogx(mat_f, np.abs(mat_Z31C), color='black', linestyle='--')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('|Zg| (Ω/m)')
    plt.legend(loc='upper left')
    plt.xlim(1E4, 1E7)
    plt.ylim(0, 25)
    plt.grid(False)
    plt.show()

def plot_impedance_angle(frequencies, impedance_data, matlab_data=None):
    """Plots impedance data."""
    plt.figure(figsize=(8, 6))
    for system_id, impedance_values in impedance_data.items():
        plt.semilogx(frequencies, np.angle(impedance_values)*180/np.pi, label=f'System {system_id} Python Data')
    if matlab_data:
        mat_f, mat_Z31A, mat_Z31B, mat_Z31C = matlab_data
        plt.semilogx(mat_f, np.angle(mat_Z31A)*180/np.pi, color='black', linestyle='--', label='MATLAB Data')
        plt.semilogx(mat_f, np.angle(mat_Z31B)*180/np.pi, color='black', linestyle='--')
        plt.semilogx(mat_f, np.angle(mat_Z31C)*180/np.pi, color='black', linestyle='--')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Angle of Zg (Degrees)')
    plt.legend(loc='lower left')
    plt.xlim(1E4, 1E7)
    plt.ylim(20, 90)
    plt.grid(False)
    plt.show()
    
def main():
    clear_screen()

    # Define single core cable model
    scc = model()

    # Define system types and soil parameters
    system_types = [SystemType(r'$\rho_g=100\; \Omega.m$', scc, 1E2, 10),
                    SystemType(r'$\rho_g=1000\; \Omega.m$', scc, 1E3, 10), 
                    SystemType(r'$\rho_g=10000\; \Omega.m$', scc, 1E4, 10)]

    # Generate logarithmic frequency range
    decade_start = 4
    decade_end = 7
    num_points = 300
    f = np.logspace(decade_start, decade_end, num=num_points)
    s = 1j * (2 * np.pi * f)

    # Calculate ground return impedance
    impedance_data = calculate_impedance(s, scc, system_types)

    # Load MATLAB data
    matlab_data = load_matlab_data(r'C:\Users\adilt\OneDrive\01 ACADEMIA\05 PPGEE\EEE003 MODELOS\1.DECONTI\MatLabData.mat')

    # Plot impedance data
    plot_impedance_abs(f, impedance_data, matlab_data)
    plot_impedance_angle(f, impedance_data, matlab_data)

if __name__ == "__main__":
    main()