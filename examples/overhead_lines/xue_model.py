"""
REFERENCES:
[1] XUE, Haoyan. General Formulation and Accurate Evaluation of Earth-Return Parameters
    for Overhead / Underground Cables. PhD thesis, Department of Electrical Engineering,
    École Polytechnique de Montréal, Université de Montréal, August 2018.
"""
import os
import sys
import copy
import time
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# --- PROJECT ROOT AND DIRECTORIES CONFIGURATION ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    script_dir = Path(__file__).resolve().parent
    print(f"Script directory: {script_dir}")
    project_root = script_dir.parents[1]
    print(f"Project root: {project_root}")
    sys.path.append(str(project_root))
    print("Project paths configured successfully.")
except IndexError:
    raise FileNotFoundError(
        "Could not find the project root. "
        "Make sure the script is in 'examples/coated_wires'."
    )

# --- MODULES AND DATA MODEL IMPORTS ---
try:
    from mtl_main.models_wires import SINGLE_OHTL_XUE as MODEL
    from mtl_main.source import MulticonductorTransmissionLine
    from mtl_main.graphics import MTLRepresentation
    from analytical_formulation.overhead_lines import PerUnitParameters
    print("Modules and data model imported successfully.") 
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

class XueModelPlotter:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, freq, pul, p=0, q=0):
        self.freq_data = freq
        self.pul_data = pul
        self.p = p
        self.q = q
        
        self.f = self.freq_data['Analytically']
        
        # --- Centralized Plotting Configuration ---
        base_title_info = (
            r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m}, '
            r'\rho = 1.68 \times 10^{-8} \, \mathrm{\Omega m}$ [1]'
        )
        self.plot_configs = {
            'fig42': {
                'resistance_title': f'Figure 4.2: P.u.l. series resistance\n{base_title_info}',
                'inductance_title': f'Figure 4.2: P.u.l. series inductance\n{base_title_info}',
                'series_to_plot': [
                    {'key': 'p100_er1', 'label': r'$\rho_e = 100 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-'},
                    {'key': 'p100_er20', 'label': r'$\rho_e = 100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--'},
                    {'key': 'p2000_er1', 'label': r'$\rho_e = 2000 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-.'}
                ]
            },
            'fig43': {
                'resistance_title': f'Figure 4.3: P.u.l. series resistance\n{base_title_info}',
                'inductance_title': f'Figure 4.3: P.u.l. series inductance\n{base_title_info}',
                'series_to_plot': [
                    {'key': 'p100_er1', 'label': r'$\rho_e = 100 \;\Omega m, \epsilon_r=1$ (Extended)', 'color': 'black', 'linestyle': '-'},
                    {'key': 'p100_er20', 'label': r'$\rho_e = 100 \;\Omega m, \epsilon_r=20$ (Extended)', 'color': 'black', 'linestyle': '-.'},
                    {'key': 'p100_er1_classical', 'label': r'$\rho_e = 100 \;\Omega m, \epsilon_r=1$ (Classical)', 'color': 'red', 'linestyle': '--'}
                ]
            }
        }

    def _plot_resistance_inductance_subplot(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

        for series in config['series_to_plot']:
            # Extract and process data
            zs_raw = np.array([item['zs'][self.p, self.q] for item in self.pul_data[series['key']]])
            rs = np.real(zs_raw) * 1e3  # Resistance in Ohm/km
            ls = np.imag(zs_raw) / (2 * np.pi * self.f) * 1e6  # Inductance in mH/km
            
            # Plotting with specified style
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, rs, **style)
            ax2.plot(self.f, ls, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1E3, 1E9)
        ax1.set_ylim(1E0, 1E5)
        ax1.legend()
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_s \, (\Omega/km)$')
        ax1.grid(False)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1E3, 1E9)
        ax2.set_ylim(1, 2.5)
        ax2.legend()
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_s \, (mH/km)$')
        ax2.grid(False)
        ax2.set_title(config['inductance_title'])

        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def plot_fig42(self):
        """Plots the data corresponding to Figure 4.2 from the reference."""
        self._plot_resistance_inductance_subplot('fig42')

    def plot_fig43(self):
        """Plots the data corresponding to Figure 4.3 from the reference."""
        self._plot_resistance_inductance_subplot('fig43')

if __name__ == "__main__":
    st = time.time()
    frequency = {'Analytically': np.logspace(3, 9, num=200)}

    # --- Data-Driven Calculation Setup ---
    # Define models for different physical scenarios
    mtl_model_a = MulticonductorTransmissionLine(MODEL)
    
    model_b_data = copy.deepcopy(MODEL)
    model_b_data[0]['relative_permittivity'] = 20
    mtl_model_b = MulticonductorTransmissionLine(model_b_data)

    model_c_data = copy.deepcopy(MODEL)
    model_c_data[0]['conductivity'] = 0.0005 # rho = 2000 Ohm.m
    mtl_model_c = MulticonductorTransmissionLine(model_c_data)

    # Define the calculation scenarios
    scenarios = {
        'p100_er1': {'mtl': mtl_model_a, 'form': 'nakagawa'},
        'p100_er20': {'mtl': mtl_model_b, 'form': 'nakagawa'},
        'p2000_er1': {'mtl': mtl_model_c, 'form': 'nakagawa'},
        'p100_er1_classical': {'mtl': mtl_model_a, 'form': 'carson'}
    }
    
    # Initialize dictionary to hold results
    pul = {key: [] for key in scenarios}

    # --- Refactored Calculation Loop ---
    for f in frequency['Analytically']:
        for key, params in scenarios.items():
            pul_instance = PerUnitParameters(params['mtl'], f)
            result = pul_instance.pul_extended_theory(form=params['form'])
            pul[key].append(result)

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")

    # --- Plotting ---
    plotter = XueModelPlotter(frequency, pul)
    plotter.plot_fig42()
    plotter.plot_fig43()

    MTLRepresentation(mtl_model_a, units='millimeter').ground_return_systems()
    plt.show()