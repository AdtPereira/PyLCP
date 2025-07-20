import datetime
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0 as MUO, epsilon_0 as EO

class ParameterAsLog:
    
    def Parameters(syst, Gri, Zg, Zi, zL, A, ZL, Zs, Ye, Pg, Yg, file_path):
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

                f"Transient Analysis in single core underground cables (SCC)\n" 
                f"Transmission Lines.\n"
                f"v6\n\n"
                                
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

class ParameterAsGraph:

    def Matlab_Python_plot_data(py_f, py_data, m_data=None, impedance=True, abs=True, ylabel='', ylim=None, yscale=1, legend='upper left', graph_path=None, file_name=None):
        """
        Plot and Save a Graph at specific directory path

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

    def Matlab_Python_Graph1(frequencies, py_data, m_data=None, syst_name=None, graph_path=None):
        """Plots impedance data."""
        ylabel = '|Zg| (Ω/m)' 
        ylim = [0, 40] if syst_name == 'DC_M04_3SCC_1C' else [0, 25]
        file_name = graph_path + syst_name + '_ImpedanceAbs.svg'    
        ParameterAsGraph.Matlab_Python_plot_data(frequencies, py_data, m_data, ylabel=ylabel, ylim=ylim, graph_path=graph_path, file_name=file_name)

    def Matlab_Python_Graph2(frequencies, py_data, m_data=None, syst_name=None, graph_path=None):
        """Plots impedance data."""
        ylabel = 'Angle of Zg (Degrees)'   
        ylim = [55, 85] if syst_name == 'DC_M04_3SCC_1C' else [20, 90]
        legend = 'lower left'
        file_name = graph_path + syst_name + '_ImpedanceAngle.svg'
        ParameterAsGraph.Matlab_Python_plot_data(frequencies, py_data, m_data, abs=False, ylabel=ylabel, ylim=ylim, legend=legend, graph_path=graph_path, file_name=file_name)    

    def Matlab_Python_Graph3(frequencies, py_data, m_data=None, syst_name=None, graph_path=None):
        """Plots potential data."""
        ylabel = r'$|P_g| \times 10^9 (\Omega ms^{-1})$'
        yscale = 1E-9
        ylim = [0, 15] if syst_name == 'DC_M04_3SCC_1C' else [0, 14] if syst_name == 'DC_M02_3SCC_1C' else [0, 12]
        file_name = graph_path + syst_name + '_PotentialAbs.svg'
        ParameterAsGraph.Matlab_Python_plot_data(frequencies, py_data, m_data, impedance=False, ylabel=ylabel, ylim=ylim, yscale=yscale, graph_path=graph_path, file_name=file_name)

    def Matlab_Python_Graph4(frequencies, py_data, m_data=None, syst_name=None, graph_path=None):
        """Plots potential data."""
        ylabel = 'Angle of Pg (Degrees)'
        ylim = [-50, 90] if syst_name == 'DC_M04_3SCC_1C' else [-90, 90]
        legend = 'lower left'
        file_name = graph_path + syst_name + '_PotentialAngle.svg'
        ParameterAsGraph.Matlab_Python_plot_data(frequencies, py_data, m_data, impedance=False, abs=False, ylabel=ylabel, ylim=ylim, legend=legend, graph_path=graph_path, file_name=file_name)            