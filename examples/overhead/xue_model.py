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
import scipy.constants as sc
import matplotlib.pyplot as plt
from pathlib import Path

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
    from mtl_main.models_wires import SINGLE_WIRE_R10_XUE as MODEL
    from mtl_main.source import MulticonductorTransmissionLine
    from mtl_main.graphics import MTLRepresentation
    from analytical_forms.overhead_lines import PerUnitParameters
    print("Modules and data model imported successfully.") 
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

class XueModelPlotter:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, mtl_model, freq, pul, p=0, q=0):
        self.freq_data = freq
        self.pul_data = pul
        self.p = p
        self.q = q
        self.f = self.freq_data['Analytically']
        self.w = 2 * np.pi * self.f
        self.ro = mtl_model.surfaces[0]['radius']
        self.h1 = mtl_model.surfaces[0]['center_point'][1]

        self.nakagawa_series = [
            {'key': 'p100',         'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20',    'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--'},
            {'key': 'p2000',        'label': r'$\rho_e=2000 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-.'}
        ]

        self.nakagawa_carson_series = [
            {'key': 'p100',         'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$ (Nakagawa)',   'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20',    'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$ (Nakagawa)',  'color': 'black', 'linestyle': '-.'},
            {'key': 'p100_carson',  'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$ (Carson)',     'color': 'red',   'linestyle': '--'}
        ]

        self.attenuation_constant_series = [
            {'key': 'p100',         'label': r'$\rho_e=100 \;\Omega m$ (Nakagawa)', 'color': 'black',   'linestyle': '-'},
            {'key': 'p100_carson',  'label': r'$\rho_e=100 \;\Omega m$ (Carson)',   'color': 'red',     'linestyle': '-'},
            {'key': 'p2000',        'label': r'$\rho_e=2000 \;\Omega m$ (Nakagawa)','color': 'black',   'linestyle': '--'},
            {'key': 'p2000_carson', 'label': r'$\rho_e=2000 \;\Omega m$ (Carson)',  'color': 'red',     'linestyle': '--'},
        ]

        self.plot_configs = {
            'fig42': {
                'suptitle': 'Figure 4.2: P.u.l. series impedance with Nakagawa formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (1, 2.5)},
                'series_to_plot': self.nakagawa_series
            },
            'fig43': {
                'suptitle': 'Figure 4.3: P.u.l. series impedance comparison [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'y_lim': {'resistance': (1E0, 1E4), 'inductance': (1.4, 2.2)},
                'series_to_plot': self.nakagawa_carson_series
            },
            'fig45': {
                'suptitle': 'Figure 4.5: P.u.l. shunt admittance with Nakagawa formulation [1]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'y_lim': {'conductance': (-0.3, 0.1), 'capacitance': (6.6, 7.4)},
                'series_to_plot': self.nakagawa_series
            },
            'fig46': {
                'suptitle': 'Figure 4.6: P.u.l. shunt admittance comparison [1]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'y_lim': {'conductance': (-0.06, 0.02), 'capacitance': (7.22, 7.32)},
                'series_to_plot': self.nakagawa_carson_series
            },
            'fig47': {
                'suptitle': 'Figure 4.7: Propagation constant with Nakagawa formulation [1]',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'y_lim': {'attenuation': (1E-3, 1E1), 'phase_velocity': (0.8, 1.1)},
                'series_to_plot': self.nakagawa_series
            },
            'fig48': {
                'suptitle': r'Figure 4.8: Propagation constant comparison for $\varepsilon_r = 1$ [1]',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'y_lim': {'attenuation': (1E-3, 1E2), 'phase_velocity': (0.8, 1.1)},
                'series_to_plot': self.attenuation_constant_series
            },
        }

    def _plot_impedance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zs_raw = np.array([item['zs'][self.p, self.q] for item in self.pul_data[series['key']]])
            rs = np.real(zs_raw) * 1e3  # Resistance in Ohm/km
            ls = np.imag(zs_raw) / (2 * np.pi * self.f) * 1e6  # Inductance in mH/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, rs, **style)
            ax2.plot(self.f, ls, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1E3, 1E9)
        ax1.set_ylim(config['y_lim']['resistance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_s \, (\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1E3, 1E9)
        ax2.set_ylim(config['y_lim']['inductance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_s \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _plot_admittance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            ysh_raw = np.array([item['ysh'][self.p, self.q] for item in self.pul_data[series['key']]])
            cond = np.real(ysh_raw) * 1e3  # Conductance in S/km
            cap = np.imag(ysh_raw) / (2 * np.pi * self.f) * 1e12  # Capacitance in nF/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, cond, **style)
            ax2.plot(self.f, cap, **style)

        # Configure left subplot (Conductance)
        ax1.set_xscale('log')
        ax1.set_xlim(1E3, 1E9)
        ax1.set_ylim(config['y_lim']['conductance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$G \, (S/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['conductance_title'])

        # Configure right subplot (Capacitance)
        ax2.set_xscale('log')
        ax2.set_xlim(1E3, 1E9)
        ax2.set_ylim(config['y_lim']['capacitance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$C \, (nF/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['capacitance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _plot_propagation_constant_subplots(self, config_key):
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            gamma_v = np.array([item['gamma_v'][self.p, self.q] for item in self.pul_data[series['key']]])
            alfa = gamma_v.real * 1e3
            phase_vel = self.w / gamma_v.imag / sc.c
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, alfa, **style)
            ax2.plot(self.f, phase_vel, **style)

        # Configure left subplot (attenuation)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1E3, 1E9)
        ax1.set_ylim(config['y_lim']['attenuation'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel( r'Attenuation Constant, $\alpha_{\nu}$ (Np/km)')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['attenuation_title'])

        # Configure right subplot (phase velocity)
        ax2.set_xscale('log')
        ax2.set_xlim(1E3, 1E9)
        ax2.set_ylim(config['y_lim']['phase_velocity'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Phase Velocity, $c_{\nu}/c_0$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['phase_velocity_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def plot_fig42(self):
        """Plots the data corresponding to Figure 4.2 from the reference."""
        self._plot_impedance_subplots('fig42')

    def plot_fig43(self):
        """Plots the data corresponding to Figure 4.3 from the reference."""
        self._plot_impedance_subplots('fig43')

    def plot_fig45(self):
        """Plots the data corresponding to Figure 4.5 from the reference."""
        self._plot_admittance_subplots('fig45')

    def plot_fig46(self):
        """Plots the data corresponding to Figure 4.6 from the reference."""
        self._plot_admittance_subplots('fig46')

    def plot_fig47(self):
        """Plots the data corresponding to Figure 4.7 from the reference."""
        self._plot_propagation_constant_subplots('fig47')

    def plot_fig48(self):
        """Plots the data corresponding to Figure 4.8 from the reference."""
        self._plot_propagation_constant_subplots('fig48')

if __name__ == "__main__":
    st = time.time()
    frequency = {'Analytically': np.logspace(3, 9, num=100)}

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
        'p100':         {'mtl': mtl_model_a, 'form': 'nakagawa'},
        'p100_er20':    {'mtl': mtl_model_b, 'form': 'nakagawa'},
        'p2000':        {'mtl': mtl_model_c, 'form': 'nakagawa'},
        'p100_carson':  {'mtl': mtl_model_a, 'form': 'carson'},
        'p2000_carson': {'mtl': mtl_model_c, 'form': 'carson'}
    }
    
    pul_parameters = {key: [] for key in scenarios}
    for f in frequency['Analytically']:
        for key, params in scenarios.items():
            pul_parameters[key].append(
                PerUnitParameters(params['mtl'], f).pul_quasi_tem(form=params['form'])
            )

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")
    plotter = XueModelPlotter(mtl_model_a, frequency, pul_parameters)
    plotter.plot_fig42()
    plotter.plot_fig43()
    plotter.plot_fig45()
    plotter.plot_fig46()
    plotter.plot_fig47()
    plotter.plot_fig48()
    MTLRepresentation(mtl_model_a, units='millimeter').ground_return_systems()
    plt.show()