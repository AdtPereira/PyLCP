import os
import sys
import time
import copy
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import scipy.constants as sc

class LimaModels:
    """
    A highly refactored class to handle plotting for the Lima model results.
    It uses a configuration-driven approach to minimize code duplication,
    with specific settings for each plot type.
    """
    def __init__(self, pul_parameters):
        self.pul_data = pul_parameters
        self.f = pul_parameters['frequencies']
        self.w = 2 * np.pi * self.f

        self.figsize = (12, 5)
        
        suptitle_suffix = (
            'of the single overhead line\n'
            r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m}, \sigma = 6.496 \times 10^7\,\mathrm{S/m}$ [2]'
        )

        self.comparison_forms = [
            {'key': 'quasi_tem', 'color': 'black', 'linestyle': '-', 'label': 'Quasi-TEM (Integral Eq.)'},
            {'key': 'quasi_tem_log', 'color': 'blue', 'linestyle': '--', 'label': 'Quasi-TEM (Approx. Log.)'},
            {'key': 'nakagawa', 'color': 'red', 'linestyle': '-', 'label': 'Wise (1948) & Nakagawa (1981)'},
            {'key': 'sunde', 'color': 'green', 'linestyle': '-', 'label': 'Sunde (1968)'},
            {'key': 'carson', 'color': 'black', 'linestyle': ':', 'label': 'Carson (1926)'}
        ]

        self.nakagawa_form = [
            {'key': 'nakagawa_1u', 'color': 'black', 'linestyle': '--', 'label': r'$\rho_g = 1 \,    \mu \Omega m$'},
            {'key': 'nakagawa',    'color': 'black', 'linestyle': '-',  'label': r'$\rho_g = 200 \,  \Omega m$'},
            {'key': 'nakagawa_5k', 'color': 'red',   'linestyle': ':',  'label': r'$\rho_g = 5000 \, \Omega m$'},
        ]

        self.quasi_tem_form = [
            {'key': 'quasi_tem_1u', 'color': 'black', 'linestyle': '--', 'label': r'$\rho_g = 1 \,    \mu \Omega m$'},
            {'key': 'quasi_tem',    'color': 'black', 'linestyle': '-',  'label': r'$\rho_g = 200 \,  \Omega m$'},
            {'key': 'quasi_tem_5k', 'color': 'red',   'linestyle': ':',  'label': r'$\rho_g = 5000 \, \Omega m$'},
        ]
        
        self.plot_configs = {
            'propagation_constant_forms': {
                'suptitle': f'Propagation constant {suptitle_suffix}',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1E2, 1E10), 'phase_velocity': (1E0, 1E10)},
                'y_lim': {'attenuation': (0, 4), 'phase_velocity': (0.5, 1.1)},
                'series_to_plot': self.comparison_forms
            },
            'propagation_constant_nakagawa': {
                'suptitle': f'Nakagawa Form. (1981) Propagation constant {suptitle_suffix}',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1E2, 1E10), 'phase_velocity': (1E0, 1E10)},
                'y_lim': {'attenuation': (0, 1.6), 'phase_velocity': (0.5, 1.1)},
                'series_to_plot': self.nakagawa_form
            },
            'propagation_constant_quasi_tem': {
                'suptitle': f'Quasi-TEM (Integral Form.) Propagation constant {suptitle_suffix}',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1E2, 1E10), 'phase_velocity': (1E0, 1E10)},
                'y_lim': {'attenuation': (0, 1.6), 'phase_velocity': (0.5, 1.1)},
                'series_to_plot': self.quasi_tem_form
            },
            'characteristic_impedance_matrix': {
                'suptitle': f'Characteristic Impedance Matrix {suptitle_suffix}',
                'norm_title': 'Absolute Impedance',
                'angle_title': 'Angle of Impedance',
                'p': 0, 'q': 0,
                'x_lim': {'norm': (1E0, 1E10), 'angle': (1E0, 1E10)},
                'y_lim': {'norm': (300, 1100), 'angle': (-40, 5)},
                'series_to_plot': self.comparison_forms
            },
            'characteristic_impedance_nakagawa': {
                'suptitle': f'Nakagawa Form. (1981) Characteristic Impedance Matrix {suptitle_suffix}',
                'norm_title': 'Absolute Impedance',
                'angle_title': 'Angle of Impedance',
                'p': 0, 'q': 0,
                'x_lim': {'norm': (1E0, 1E10), 'angle': (1E0, 1E10)},
                'y_lim': {'norm': (300, 1100), 'angle': (-40, 5)},
                'series_to_plot': self.nakagawa_form
            },
            'characteristic_impedance_quasi_tem': {
                'suptitle': f'Quasi-TEM (Integral Form.) Characteristic Impedance Matrix {suptitle_suffix}',
                'norm_title': 'Absolute Impedance',
                'angle_title': 'Angle of Impedance',
                'p': 0, 'q': 0,
                'x_lim': {'norm': (1E0, 1E10), 'angle': (1E0, 1E10)},
                'y_lim': {'norm': (300, 1100), 'angle': (-40, 5)},
                'series_to_plot': self.quasi_tem_form
            },
        }

    def _propagation_constant_subplots(self, config_key):
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            gamma_v_3d = self.pul_data[series['key']]['propagation_voltage_matrix']
            gamma_v = gamma_v_3d[:, p, q]   
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, gamma_v.real * 1e3, **style)
            ax2.plot(self.f, self.w / gamma_v.imag / sc.c, **style)

        # Configure left subplot (attenuation)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['attenuation'])
        ax1.set_ylim(config['y_lim']['attenuation'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel( r'Attenuation Constant, $\alpha_{\nu}$ (Np/km)')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['attenuation_title'])

        # Configure right subplot (phase velocity)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['phase_velocity'])
        ax2.set_ylim(config['y_lim']['phase_velocity'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Phase Velocity, $c_{\nu}/c_0$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['phase_velocity_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _characteristic_impedance_subplots(self, config_key):
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zc_3d = self.pul_data[series['key']]['characteristic_impedance_matrix']
            zc = zc_3d[:, p, q]
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': 1.0}
            ax1.plot(self.f, np.abs(zc), **style)
            ax2.plot(self.f, np.angle(zc, deg=True), **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['norm'])
        ax1.set_ylim(config['y_lim']['norm'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$|Zc_{{{p+1}{q+1}}}| \, (\Omega/m)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['norm_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['angle'])
        ax2.set_ylim(config['y_lim']['angle'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'Angle of $Zc_{{{p+1}{q+1}}}$ (Degrees)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['angle_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def propagation_constant_forms(self):
        self._propagation_constant_subplots('propagation_constant_forms')

    def propagation_constant_nakagawa(self):
        self._propagation_constant_subplots('propagation_constant_nakagawa')

    def propagation_constant_quasi_tem(self):
        self._propagation_constant_subplots('propagation_constant_quasi_tem')

    def characteristic_impedance_matrix(self):
        self._characteristic_impedance_subplots('characteristic_impedance_matrix')

    def characteristic_impedance_nakagawa(self):
        self._characteristic_impedance_subplots('characteristic_impedance_nakagawa')

    def characteristic_impedance_quasi_tem(self):
        self._characteristic_impedance_subplots('characteristic_impedance_quasi_tem')
