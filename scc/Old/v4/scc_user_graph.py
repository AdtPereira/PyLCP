import matplotlib.pyplot as plt
import numpy as np
import os

def plot_data(py_f, py_data, m_data=None, impedance=True, abs=True, ylabel='', ylim=None, yscale=1, legend='upper left', graph_path=None, file_name=None):
    """
    Plot and Save a Graph at especific directory path

    Parameters:
    - frequencies: Array-like, frequencies for x-axis.
    - py_data: Dictionary, Python data for plotting.
    - m_data: Tuple, optional, MATLAB data for comparison.
    - impedance: bool, optional, whether plotting impedance data (True) or potential data (False).
    - abs: bool, optional, whether plotting absolute values (True) or angles (False) of the data.
    - ylabel: str, optional, label for the y-axis.
    - ylim: tuple, optional, limits for the y-axis.
    - yscale: float, optional, scaling factor for the y-axis.
    - legend: str, optional, location of the legend.

    Returns:
    - None
    """
    plt.figure(figsize=(8, 6))
    
    for system_id, values in py_data.items():
        if abs:
            y_values = np.abs(values) * yscale
        else:
            y_values = np.angle(values) * 180 / np.pi        
        plt.semilogx(py_f, y_values, label=f'System {system_id} Python Data')

    if m_data:
        mat_f, *mat_values = m_data
        if impedance:
            if abs:
                plt.semilogx(mat_f, np.abs(mat_values[0]), color='black', linestyle='--', label='MATLAB Data')
                plt.semilogx(mat_f, np.abs(mat_values[1]), color='black', linestyle='--')
                plt.semilogx(mat_f, np.abs(mat_values[2]), color='black', linestyle='--')
            else:
                plt.semilogx(mat_f, np.angle(mat_values[0]) * 180 / np.pi, color='black', linestyle='--', label='MATLAB Data')
                plt.semilogx(mat_f, np.angle(mat_values[1]) * 180 / np.pi, color='black', linestyle='--')
                plt.semilogx(mat_f, np.angle(mat_values[2]) * 180 / np.pi, color='black', linestyle='--')
        else:
            if abs:
                plt.semilogx(mat_f, np.abs(mat_values[3]) * yscale, color='black', linestyle='--', label='MATLAB Data')
                plt.semilogx(mat_f, np.abs(mat_values[4]) * yscale, color='black', linestyle='--')
                plt.semilogx(mat_f, np.abs(mat_values[5]) * yscale, color='black', linestyle='--')
            else:
                plt.semilogx(mat_f, np.angle(mat_values[3]) * 180 / np.pi, color='black', linestyle='--', label='MATLAB Data')
                plt.semilogx(mat_f, np.angle(mat_values[4]) * 180 / np.pi, color='black', linestyle='--')
                plt.semilogx(mat_f, np.angle(mat_values[5]) * 180 / np.pi, color='black', linestyle='--')

    plt.xlabel('Frequency (Hz)')
    plt.ylabel(ylabel)
    plt.legend(loc=legend)
    plt.xlim(1E4, 1E7)
    plt.ylim(ylim)
    plt.grid(False)

    # Save the plot if graph_path is provided
    if graph_path:
        os.makedirs(graph_path, exist_ok=True)  # Create directory if it doesn't exist
        plt.savefig(file_name)
        print(f"Graph saved at {graph_path}")

    plt.show()

def Graph1(frequencies, py_data, m_data=None, syst_name=None, graph_path=None):
    """Plots impedance data."""
    ylabel = '|Zg| (Ω/m)' 
    ylim = [0, 40] if syst_name == 'DC_M04_3SCC_1C' else [0, 25]
    file_name = graph_path + syst_name + '_ImpedanceAbs.svg'    
    plot_data(frequencies, py_data, m_data, ylabel=ylabel, ylim=ylim, graph_path=graph_path, file_name=file_name)

def Graph2(frequencies, py_data, m_data=None, syst_name=None, graph_path=None):
    """Plots impedance data."""
    ylabel = 'Angle of Zg (Degrees)'   
    ylim = [55, 85] if syst_name == 'DC_M04_3SCC_1C' else [20, 90]
    legend = 'lower left'
    file_name = graph_path + syst_name + '_ImpedanceAngle.svg'
    plot_data(frequencies, py_data, m_data, abs=False, ylabel=ylabel, ylim=ylim, legend=legend, graph_path=graph_path, file_name=file_name)    

def Graph3(frequencies, py_data, m_data=None, syst_name=None, graph_path=None):
    """Plots potential data."""
    ylabel = r'$|P_g| \times 10^9 (\Omega ms^{-1})$'
    yscale = 1E-9
    ylim = [0, 15] if syst_name == 'DC_M04_3SCC_1C' else [0, 14] if syst_name == 'DC_M02_3SCC_1C' else [0, 12]
    file_name = graph_path + syst_name + '_PotentialAbs.svg'
    plot_data(frequencies, py_data, m_data, impedance=False, ylabel=ylabel, ylim=ylim, yscale=yscale, graph_path=graph_path, file_name=file_name)

def Graph4(frequencies, py_data, m_data=None, syst_name=None, graph_path=None):
    """Plots potential data."""
    ylabel = 'Angle of Pg (Degrees)'
    ylim = [-50, 90] if syst_name == 'DC_M04_3SCC_1C' else [-90, 90]
    legend = 'lower left'
    file_name = graph_path + syst_name + '_PotentialAngle.svg'
    plot_data(frequencies, py_data, m_data, impedance=False, abs=False, ylabel=ylabel, ylim=ylim, legend=legend, graph_path=graph_path, file_name=file_name)