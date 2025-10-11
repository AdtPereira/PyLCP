import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

from pathlib import Path
from utils.case_utils import *

class XueModels:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    This version is adapted for the vectorized data structure.
    """
    def __init__(self, file_path: str, pul_data: dict, autoSave: bool = True):
        """
        Initializes the XueModels class with the provided per-unit-length data.
        :param file_path: Path to the current file.
        :param pul_data: Dictionary containing per-unit-length data and frequencies.
        """
        self.script_path = Path(file_path)
        self.autoSave = autoSave
        self.pul_data = pul_data
        self.f = pul_data['frequencies']
        self.w = 2 * np.pi * self.f
        self.cmsl = self.pul_data['comsol'] if 'comsol' in self.pul_data else None
        
        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.script_path.stem, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

        self.figsize = (12, 5)

        comsol_layout_template = [
            {
                'base_key': 'rho_g_100_epsr1_1_mf',
                'label': r'COMSOL ($\rho_e=100, \epsilon_r=1$)',
                'marker': 'o', 's': 30, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 10
            },
            {
                'base_key': 'rho_g_100_epsr1_20_mf',
                'label': r'COMSOL ($\rho_e=100, \epsilon_r=20$)',
                'marker': 'o', 's': 10, 'facecolors': 'black', 'edgecolors': 'none', 'zorder': 10
            },
            {
                'base_key': 'rho_g_500_epsr1_1_mf',
                'label': r'COMSOL ($\rho_e=500, \epsilon_r=1$)',
                'marker': '*', 's': 12, 'color': 'black', 'zorder': 10
            }
        ]

        def generate_comsol_series(template, coil_suffix):
            """
            Gera uma lista de dicionários para plotagem, combinando um template de layout
            com um sufixo específico para a chave de dados.
            """
            series_list = []
            for item_template in template:
                new_item = item_template.copy()
                new_item['data_key'] = f"{new_item.pop('base_key')}_{coil_suffix}"
                series_list.append(new_item)
            return series_list
        
        self.overhead_line = [
            {'key': 'p100',         'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20',    'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--'},
            {'key': 'p2000',        'label': r'$\rho_e=2000 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-.'}
        ]

        self.sc_cables = [
            {
                'key': 'p100',
                'type': {
                    'series_term': {'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                    'ground-return_term': {'label': '', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0}
            }},
            {
                'key': 'p100_er20',
                'type': {
                    'series_term': {'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--', 'linewidth': 1.5},
                    'ground-return_term': {'label': '', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0}
            }},
            {
                'key': 'p500',
                'type': {
                    'series_term': {'label': r'$\rho_e=500 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-.', 'linewidth': 1.5},
                    'ground-return_term': {'label': 'Ground-Return Term', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0}
            }},
        ]        

        self.xue_series_norm = [
            {'key': 'p100', 'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-'},
            {'key': 'p100_deconti', 'label': 'De Conti Approx.', 'color': 'red', 'linestyle': ':'},
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
                'suptitle': 'Figure 4.2: P.u.l. series impedance with Nakagawa formulation [Xue, 2018]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 0, 'q': 0,
                'x_lim': {'resistance': (1E3, 1E9), 'inductance': (1E3, 1E9)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (1, 2.5)},
                'series_to_plot': self.overhead_line
            },
            'fig43': {
                'suptitle': 'Figure 4.3: P.u.l. series impedance comparison [Xue, 2018]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 0, 'q': 0,
                'x_lim': {'resistance': (1E3, 1E9), 'inductance': (1E3, 1E9)},
                'y_lim': {'resistance': (1E0, 1E4), 'inductance': (1.4, 2.2)},
                'series_to_plot': self.nakagawa_carson_series
            },
            'fig45': {
                'suptitle': 'Figure 4.5: P.u.l. shunt admittance with Nakagawa formulation [Xue, 2018]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 0, 'q': 0,
                'x_lim': {'conductance': (1E3, 1E9), 'capacitance': (1E3, 1E9)},
                'y_lim': {'conductance': (-0.3, 0.1), 'capacitance': (6.6, 7.4)},
                'series_to_plot': self.overhead_line
            },
            'fig46': {
                'suptitle': 'Figure 4.6: P.u.l. shunt admittance comparison [Xue, 2018]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 0, 'q': 0,
                'x_lim': {'conductance': (1E3, 1E9), 'capacitance': (1E3, 1E9)},
                'y_lim': {'conductance': (-0.06, 0.02), 'capacitance': (7.22, 7.32)},
                'series_to_plot': self.nakagawa_carson_series
            },
            'fig47': {
                'suptitle': 'Figure 4.7: Propagation constant with Nakagawa formulation [Xue, 2018]',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1E3, 1E9), 'phase_velocity': (1E3, 1E9)},
                'y_lim': {'attenuation': (1E-3, 1E1), 'phase_velocity': (0.8, 1.1)},
                'series_to_plot': self.overhead_line
            },
            'fig48': {
                'suptitle': r'Figure 4.8: Propagation constant comparison for $\varepsilon_r = 1$ [Xue, 2018]',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1E3, 1E9), 'phase_velocity': (1E3, 1E9)},
                'y_lim': {'attenuation': (1E-3, 1E2), 'phase_velocity': (0.8, 1.1)},
                'series_to_plot': self.attenuation_constant_series
            },
            'fig419': {
                'suptitle': 'Figure 4.19: P.u.l. Self-impedance of phase - a sheath with Magalhães/Xue formulation [Xue, 2018]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, 'q': 1,
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E1, 1E5), 'inductance': (0.5, 2.0)},
                'y_ticks': {'resistance': np.arange(1E0, 1E5, 1E1), 'inductance': np.arange(0.5, 2.1, 0.5)},
                'series_to_plot': self.sc_cables,
                'comsol_series_to_plot': generate_comsol_series(comsol_layout_template, 'vcoil_1')
            },
            'fig421': {
                'suptitle': 'Figure 4.21: P.u.l. Mutual-impedance between phase - a and phase - b sheaths with Magalhães/Xue formulation [Xue, 2018]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, 'q': 3,
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E1, 1E5), 'inductance': (0.0, 1.5)},
                'y_ticks': {'resistance': np.arange(1E1, 1E5, 1E1), 'inductance': np.arange(0.0, 1.6, 0.5)},
                'series_to_plot': self.sc_cables,
                'comsol_series_to_plot': generate_comsol_series(comsol_layout_template, 'vcoil_2')
            },
            'fig421a': {
                'suptitle': 'Figure 4.21a: P.u.l. Mutual-impedance between phase - a and phase - c sheaths with Magalhães/Xue formulation [Xue, 2018]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, 'q': 5,
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E1, 1E5), 'inductance': (0.0, 1.5)},
                'y_ticks': {'resistance': np.arange(1E1, 1E5, 1E1), 'inductance': np.arange(0.0, 1.6, 0.5)},
                'series_to_plot': self.sc_cables,
                'comsol_series_to_plot': generate_comsol_series(comsol_layout_template, 'vcoil_3')
            },
            'fig423': {
                'suptitle': 'Figure 4.23: P.u.l. Self-admittance of phase - a sheath with Magalhães/Xue formulation [Xue, 2018]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 1, 'q': 1,
                'x_lim': {'conductance': (1E3, 1E7), 'capacitance': (1E3, 1E7)},
                'y_lim': {'conductance': (0.0, 20), 'capacitance': (0.0, 3.0)},
                'y_ticks': {'conductance': np.arange(0.0, 21, 5), 'capacitance': np.arange(0.0, 3.1, 1.0)},
                'series_to_plot': self.sc_cables
            },
            'fig425': {
                'suptitle': 'Figure 4.25: P.u.l. Mutual-admittance between phase - a and phase - b sheaths with Magalhães/Xue formulation [Xue, 2018]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 1, 'q': 3,
                'x_lim': {'conductance': (1E3, 1E7), 'capacitance': (1E3, 1E7)},
                'y_lim': {'conductance': (-8.0, 2.0), 'capacitance': (-0.8, 0.2)},
                'y_ticks': {'conductance': np.arange(-8.0, 2.1, 2), 'capacitance': np.arange(-0.8, 0.3, 0.2)},
                'series_to_plot': self.sc_cables
            },
            'fig425a': {
                'suptitle': 'Figure 4.25a: P.u.l. Mutual-admittance between phase - a and phase - c sheaths with Magalhães/Xue formulation [Xue, 2018]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 1, 'q': 5,
                'x_lim': {'conductance': (1E3, 1E7), 'capacitance': (1E3, 1E7)},
                'y_lim': {'conductance': (-4.0, 2.0), 'capacitance': (-0.6, 0.2)},
                'y_ticks': {'conductance': np.arange(-4.0, 2.1, 2), 'capacitance': np.arange(-0.6, 0.3, 0.2)},
                'series_to_plot': self.sc_cables
            },
        }
    
    def _scc_impedance_subplots(self, graph_key, graph_form):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)
        
        if graph_form == 'resistance_and_inductance':

            # Analytical plotting loop
            for series in config['series_to_plot']:
                zs = self.pul_data[series['key']]['series_impedance_matrix']
                zg = self.pul_data[series['key']]['earth-return_impedance_matrix']
                
                ax1.plot(self.f, np.real(zs[:, p, q]) * 1e3, **series['type']['series_term'])
                ax2.plot(self.f, np.imag(zs[:, p, q]) / self.w * 1e6, **series['type']['series_term'])
                
                ax1.plot(self.f, np.real(zg[:, p, q]) * 1e3, **series['type']['ground-return_term'])
                ax2.plot(self.f, np.imag(zg[:, p, q]) / self.w * 1e6, **series['type']['ground-return_term'])

            # COMSOL plotting loop
            if self.cmsl is not None and 'comsol_series_to_plot' in config:
                comsol_data_source = self.cmsl['cmsl_ground_return_impedance']
                f_cmsl = comsol_data_source['freq']
                w_cmsl = 2 * np.pi * f_cmsl
                
                for series in config['comsol_series_to_plot']:
                    zg = comsol_data_source[series['data_key']]                    
                    style = {k: v for k, v in series.items() if k != 'data_key'}

                    ax1.scatter(f_cmsl, np.real(zg) * 1e3, **style)
                    ax2.scatter(f_cmsl, np.imag(zg) / w_cmsl * 1e6, **style)

            # Configure left subplot (Resistance)
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlim(config['x_lim']['resistance'])
            ax1.set_ylim(config['y_lim']['resistance'])
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$Rs_{{{p+1}{q+1}}} \, (\Omega/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['resistance_title'])
            
            ax2.set_xscale('log')
            ax2.set_xlim(config['x_lim']['inductance'])
            ax2.set_ylim(config['y_lim']['inductance'])
            ax2.set_yticks(config['y_ticks']['inductance'])
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$Ls_{{{p+1}{q+1}}} \, (mH/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['inductance_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

            if self.autoSave:
                save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}_{graph_form}')

    def _scc_admittance_subplots(self, config_key, graph_form):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)
        
        if graph_form == 'conductance_and_capacitance':
            
            # Analytical plotting loop
            for series in config['series_to_plot']:
                ysh = self.pul_data[series['key']]['shunt_admittance_matrix']
                
                ax1.plot(self.f, np.real(ysh[:, p, q]) * 1e3, **series['type']['series_term'])
                ax2.plot(self.f, np.imag(ysh[:, p, q]) / self.w * 1e9, **series['type']['series_term'])

            ax1.set_xscale('log')
            ax1.set_xlim(config['x_lim']['conductance'])
            ax1.set_ylim(config['y_lim']['conductance'])
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$G_{{{p+1}{q+1}}} \, (S/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['conductance_title'])
            
            ax2.set_xscale('log')
            ax2.set_xlim(config['x_lim']['capacitance'])
            ax2.set_ylim(config['y_lim']['capacitance'])
            ax2.set_yticks(config['y_ticks']['capacitance'])
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$C_{{{p+1}{q+1}}} \, (\mu F/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['capacitance_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

            if self.autoSave:
                save_figure_multiformat(fig, self.results_dir, base_filename=f'{config_key}_{graph_form}')

    def _overhead_impedance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zs_3d = self.pul_data[series['key']]['series_impedance_matrix']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, zs_3d[:, p, q].real * 1e3, **style)
            ax2.plot(self.f, zs_3d[:, p, q].imag / self.w * 1e6, **style)

        if self.cmsl is not None:
            cmsl = self.cmsl['cmsl_ground_return_impedance_h5']  
            zs = cmsl['coil_impedance']
            f, w = cmsl['freq'], 2 * np.pi * cmsl['freq']
            ax1.scatter(f, np.real(zs) * 1e3, label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)
            ax2.scatter(f, np.imag(zs) / w * 1e6, label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(config['x_lim']['resistance'])
        # ax1.set_ylim(config['y_lim']['resistance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_s \, (\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['inductance'])
        # ax2.set_ylim(config['y_lim']['inductance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_s \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _overhead_admittance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            ysh_3d = self.pul_data[series['key']]['shunt_admittance_matrix']
            ysh = ysh_3d[:, p, q]    
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, ysh.real * 1e3 , **style)
            ax2.plot(self.f, ysh.imag / self.w * 1e12, **style)

        # Configure left subplot (Conductance)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['conductance'])
        ax1.set_ylim(config['y_lim']['conductance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$G \, (S/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['conductance_title'])

        # Configure right subplot (Capacitance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['capacitance'])
        ax2.set_ylim(config['y_lim']['capacitance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$C \, (nF/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['capacitance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _overhead_propagation_constant_subplots(self, config_key):
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
        ax1.set_yscale('log')
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

    def plot_fig42(self):
        """Plots the data corresponding to Figure 4.2 from the reference."""
        self._overhead_impedance_subplots('fig42')

    def plot_fig43(self):
        """Plots the data corresponding to Figure 4.3 from the reference."""
        self._overhead_impedance_subplots('fig43')

    def plot_fig45(self):
        """Plots the data corresponding to Figure 4.5 from the reference."""
        self._overhead_admittance_subplots('fig45')

    def plot_fig46(self):
        """Plots the data corresponding to Figure 4.6 from the reference."""
        self._overhead_admittance_subplots('fig46')

    def plot_fig47(self):
        """Plots the data corresponding to Figure 4.7 from the reference."""
        self._overhead_propagation_constant_subplots('fig47')

    def plot_fig48(self):
        """Plots the data corresponding to Figure 4.8 from the reference."""
        self._overhead_propagation_constant_subplots('fig48')

    def plot_fig419(self, graph_form='resistance_and_inductance'):
        """Plots the data corresponding to Figure 4.19 from the reference."""
        self._scc_impedance_subplots('fig419', graph_form)

    def plot_fig421(self, graph_form='resistance_and_inductance'):
        """Plots the data corresponding to Figure 4.21 from the reference."""
        self._scc_impedance_subplots('fig421', graph_form)

    def plot_fig421a(self, graph_form='resistance_and_inductance'):
        """Plots the data corresponding to Figure 4.21a from the reference."""
        self._scc_impedance_subplots('fig421a', graph_form)

    def plot_fig423(self, graph_form='conductance_and_capacitance'):
        """Plots the data corresponding to Figure 4.23 from the reference."""
        self._scc_admittance_subplots('fig423', graph_form)

    def plot_fig425(self, graph_form='conductance_and_capacitance'):
        """Plots the data corresponding to Figure 4.25 from the reference."""
        self._scc_admittance_subplots('fig425', graph_form)

    def plot_fig425a(self, graph_form='conductance_and_capacitance'):
        """Plots the data corresponding to Figure 4.25a from the reference."""
        self._scc_admittance_subplots('fig425a', graph_form)