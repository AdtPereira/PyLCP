# C:\git\PyLCP\testData\scc_deConti_flat\scc_deConti_flat.py

import sys
import os
from pathlib import Path
import time
import numpy as np
import copy
import matplotlib.pyplot as plt

try:
    os.system('cls' if os.name == 'nt' else 'clear')
    project_root = Path(__file__).resolve().parents[2]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root)) 
    print(f"Project root configured at: {project_root}")
except IndexError:
    raise RuntimeError("Could not find project root. Ensure the directory structure is correct.")

try:
    from utils.case_utils import load_model 
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_formulation.scc import PerUnitParameters
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

class ModelPlotter:
    """
    (O conteúdo desta classe permanece inalterado)
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, mtl_model, freq, pul):
        self.freq_data = freq
        self.pul_data = pul
        self.f = self.freq_data['Analytically']
        self.w = 2 * np.pi * self.f
        self.ro = mtl_model.surfaces[0]['radius']
        self.h1 = mtl_model.surfaces[0]['center_point'][1]

        self.comparison_series = [
            {'key': 'p100',   'label': r'$\rho_e=100 \;\Omega m$',  'color': 'black', 'linestyle': '-', 'marker': 'o'},
            {'key': 'p1000',  'label': r'$\rho_e=1000 \;\Omega m$', 'color': 'black', 'linestyle': '-', 'marker': 's'},
            {'key': 'p10000', 'label': r'$\rho_e=10000 \;\Omega m$','color': 'black', 'linestyle': '-', 'marker': '^'},
            {'key': 'p100_deConti', 'label': '', 'color': 'red', 'linestyle': '--', 'marker': None},
            {'key': 'p1000_deConti', 'label': '', 'color': 'red', 'linestyle': '--', 'marker': None},
            {'key': 'p10000_deConti', 'label': '', 'color': 'red', 'linestyle': '--', 'marker': None},
        ]

        self.plot_configs = {
            'impedance': {
                'suptitle': r"Flat arrangement's mutual ground-return impedance for $\varepsilon_{r1} = 10$",
                'norm_title': 'Absolute Impedance',
                'angle_title': 'Angle of Impedance',
                'p': 2,
                'q': 0,
                'x_lim': {'norm': (1E4, 1E7), 'angle': (1E4, 1E7)},
                'y_lim': {'norm': (0, 25), 'angle': (20, 90)},
                'series_to_plot': self.comparison_series
            },
            'potential': {
                'suptitle': r"Flat arrangement's mutual ground-return potential coefficients for $\varepsilon_{r1} = 10$",
                'norm_title': 'Absolute Value',
                'angle_title': 'Angle of Potential Coefficient',
                'p': 2,
                'q': 0,
                'x_lim': {'norm': (1E4, 1E7), 'angle': (1E4, 1E7)},
                'y_lim': {'norm': (0, 14), 'angle': (-90, 90)},
                'series_to_plot': self.comparison_series
            },            
        }

    def _plot_impedance_subplots(self, config_key):
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zg = np.array([item['earth-return_impedance_matrix'][p, q] for item in self.pul_data[series['key']]])
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'],
                     'marker': series['marker'], 'markersize': 4, 'linewidth': 1.0}
            ax1.plot(self.f, np.abs(zg), **style)
            ax2.plot(self.f, np.angle(zg, deg=True), **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['norm'])
        ax1.set_ylim(config['y_lim']['norm'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$|Zg_{{{p+1}{q+1}}}| \, (\Omega/m)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['norm_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['angle'])
        ax2.set_ylim(config['y_lim']['angle'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'Angle of $Zg_{{{p+1}{q+1}}}$ (Degrees)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['angle_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _plot_potential_subplots(self, config_key):
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zg = np.array([item['earth-return_potential_coefficient'][p, q] for item in self.pul_data[series['key']]])
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'],
                     'marker': series['marker'], 'markersize': 4, 'linewidth': 1.0}
            ax1.plot(self.f, np.abs(zg) * 1e-9, **style)
            ax2.plot(self.f, np.angle(zg, deg=True), **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['norm'])
        ax1.set_ylim(config['y_lim']['norm'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$|Pg_{{{p+1}{q+1}}}| \times 10^9 \, (\Omega m s^{{-1}})$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['norm_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['angle'])
        ax2.set_ylim(config['y_lim']['angle'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'Angle of $Pg_{{{p+1}{q+1}}}$ (Degrees)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['angle_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def plot_impedance_comparison(self):
        """Plots the data corresponding to Figure 3 from the reference."""
        self._plot_impedance_subplots('impedance')

    def plot_potential_comparison(self):
        """Plots the data corresponding to Figure 6 from the reference."""
        self._plot_potential_subplots('potential')

def main():
    """
    Main function to run the simulation and plotting.
    """
    st = time.time()    
    MODEL = load_model(__file__)
    frequency = {'Analytically': np.logspace(4, 7, num=80)}

    mtl_model_a = MulticonductorTransmissionLine(MODEL)
    
    model_b = copy.deepcopy(MODEL)
    model_b[0]['conductivity'] = 0.001
    mtl_model_b = MulticonductorTransmissionLine(model_b)
    
    model_c = copy.deepcopy(MODEL)
    model_c[0]['conductivity'] = 0.0001
    mtl_model_c = MulticonductorTransmissionLine(model_c)

    scenarios = {
        'p100':  {'mtl': mtl_model_a, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p1000': {'mtl': mtl_model_b, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p10000':{'mtl': mtl_model_c, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p100_deConti':  {'mtl': mtl_model_a, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        'p1000_deConti': {'mtl': mtl_model_b, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        'p10000_deConti':{'mtl': mtl_model_c, 'zg_form': 'deconti', 'yg_form': 'deconti'},
    }
    
    pul_parameters = {key: [] for key in scenarios}
    for f in frequency['Analytically']:
        print(f"Calculating for {f:.1e} Hz ...")
        for key, params in scenarios.items():
            pul = PerUnitParameters(params['mtl'], f)
            pul_parameters[key].append(
                pul.ground_return_parameters(zg_form=params['zg_form'], yg_form=params['yg_form'])
            )

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = ModelPlotter(mtl_model_a, frequency, pul_parameters)
    plotter.plot_impedance_comparison()
    plotter.plot_potential_comparison()    
    MTLRepresentation(mtl_model_a, units='millimeter').ground_return_systems()
    plt.show()

if __name__ == "__main__":
    main()