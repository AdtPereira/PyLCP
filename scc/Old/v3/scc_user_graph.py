import matplotlib.pyplot as plt
import numpy as np

def DC_M02_3SCC_1C_gph1(frequencies, impedance_data, matlab_data=None):
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

def DC_M02_3SCC_1C_gph2(frequencies, impedance_data, matlab_data=None):
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