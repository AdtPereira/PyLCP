import os
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from pathlib import Path
from utils.case_utils import *

import os
import re
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
        A highly refactored class to handle plotting for the Xue model results.
        It uses a configuration-driven approach to generate complex subplot figures.
        This version is adapted for the vectorized data structure.
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
        
        self.xlim = tuple(pul_data['frequencies'][[0, -1]])
        self.figsize = (12, 5)

        comsol_layout_template = [
            {
                'base_key': 'rho_g_100_epsr1_1_mf',
                'label': r'COMSOL ($\rho_e=100, \epsilon_r=1$)',
                'marker': 'o', 's': 35, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 10
            },
            {
                'base_key': 'rho_g_100_epsr1_20_mf',
                'label': r'COMSOL ($\rho_e=100, \epsilon_r=20$)',
                'marker': 'o', 's': 15, 'facecolors': 'black', 'edgecolors': 'none', 'zorder': 10
            },
            {
                'base_key': 'rho_g_500_epsr1_1_mf',
                'label': r'COMSOL ($\rho_e=500, \epsilon_r=1$)',
                'marker': 's', 's': 15, 'color': 'black', 'zorder': 10
            }
        ]

        self.overhead_line = [
            {'key': 'p100',         'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20',    'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--'},
            {'key': 'p2000',        'label': r'$\rho_e=2000 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-.'}
        ]
        
        self.xue_series = [
            {
                'key': 'p100_er1',
                'type': {
                    'main': {'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                }
            },
            {
                'key': 'p100_er20',
                'type': {
                    'main': {'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--', 'linewidth': 1.5},
                }
            },
            {
                'key': 'p500_er1',
                'type': {
                    'main': {'label': r'$\rho_e=500 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-.', 'linewidth': 1.5},
                }
            },
        ] 

        self.vance_series = [
            {
                'key': 'p100_er1_vance',
                'type': {
                    'main': {'label': 'Vance Form. (1978)', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0},
                }
            },
            {
                'key': 'p100_er20_vance',
                'type': {
                    'main': {'label': '', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0},
                }
            },
            {
                'key': 'p500_er1_vance',
                'type': {
                    'main': {'label': '', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0},
                }
            },
        ] 
        
        self.deconti_series = [
            {
                'key': 'p100_er1_deconti',
                'type': {
                    'main': {'label': 'De Conti et al. Form. (2023)', 'color': 'red', 'linestyle': ':', 'linewidth': 1.0},
                }
            },
            {
                'key': 'p100_er20_deconti',
                'type': {
                    'main': {'label': '', 'color': 'red', 'linestyle': ':', 'linewidth': 1.0},
                }
            },
            {
                'key': 'p500_er1_deconti',
                'type': {
                    'main': {'label': '', 'color': 'red', 'linestyle': ':', 'linewidth': 1.0},
                }
            },
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

        self.plot_configuration = {
            'fig42': {
                'suptitle': 'Figure 4.2: P.u.l. series impedance with Nakagawa formulation [Xue, 2018]',
                'p': 0, 'q': 0,
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (1, 2.5)},
                'series_to_plot': self.overhead_line
            },
            'fig43': {
                'suptitle': 'Figure 4.3: P.u.l. series impedance comparison [Xue, 2018]',
                'p': 0, 'q': 0,
                'y_lim': {'resistance': (1E0, 1E4), 'inductance': (1.4, 2.2)},
                'series_to_plot': self.nakagawa_carson_series
            },
            'fig45': {
                'suptitle': 'Figure 4.5: P.u.l. shunt admittance with Nakagawa formulation [Xue, 2018]',
                'p': 0, 'q': 0,
                'y_lim': {'conductance': (-0.3, 0.1), 'capacitance': (6.6, 7.4)},
                'series_to_plot': self.overhead_line
            },
            'fig46': {
                'suptitle': 'Figure 4.6: P.u.l. shunt admittance comparison [Xue, 2018]',
                'p': 0, 'q': 0,
                'y_lim': {'conductance': (-0.06, 0.02), 'capacitance': (7.22, 7.32)},
                'series_to_plot': self.nakagawa_carson_series
            },
            'fig47': {
                'suptitle': 'Figure 4.7: Propagation constant with Nakagawa formulation [Xue, 2018]',
                'p': 0, 'q': 0,
                'y_lim': {'attenuation': (1E-3, 1E1), 'phase_velocity': (0.8, 1.1)},
                'series_to_plot': self.overhead_line
            },
            'fig48': {
                'suptitle': r'Figure 4.8: Propagation constant comparison for $\varepsilon_r = 1$ [Xue, 2018]',
                'p': 0, 'q': 0,
                'y_lim': {'attenuation': (1E-3, 1E2), 'phase_velocity': (0.8, 1.1)},
                'series_to_plot': self.attenuation_constant_series
            },
            'fig419': {
                'suptitle': 'Figure 4.19: P.u.l. Self-impedance of phase - a sheath [Xue, 2018]',
                'p': 1, 'q': 1,
                'y_lim': {'resistance': (1E1, 1E5), 'inductance': (0.5, 2.0)},
                'y_ticks': {'resistance': np.arange(1E0, 1E5, 1E1), 'inductance': np.arange(0.5, 2.1, 0.5)},
                'series_to_plot': self.xue_series + self.deconti_series,
                'comsol_series_to_plot': comsol_layout_template
            },
            'fig421a': {
                'suptitle': 'Figure 4.21a: P.u.l. Mutual-impedance between phase - a and phase - b sheaths [Xue, 2018]',
                'p': 1, 'q': 3,
                'y_lim': {'resistance': (1E1, 1E5), 'inductance': (0.0, 1.5)},
                'y_ticks': {'resistance': np.arange(1E1, 1E5, 1E1), 'inductance': np.arange(0.0, 1.6, 0.5)},
                'series_to_plot': self.xue_series + self.deconti_series,
                'comsol_series_to_plot': comsol_layout_template
            },
            'fig421b': {
                'suptitle': 'Figure 4.21b: P.u.l. Mutual-impedance between phase - a and phase - c sheaths [Xue, 2018]',
                'p': 1, 'q': 5,
                'y_lim': {'resistance': (1E1, 1E5), 'inductance': (0.0, 1.5)},
                'y_ticks': {'resistance': np.arange(1E1, 1E5, 1E1), 'inductance': np.arange(0.0, 1.6, 0.5)},
                'series_to_plot': self.xue_series + self.deconti_series,
                'comsol_series_to_plot': comsol_layout_template
            },
            'fig423': {
                'suptitle': 'Figure 4.23: P.u.l. Self-admittance of phase - a sheath [Xue, 2018]',
                'p': 1, 'q': 1,
                'y_lim': {'conductance': (0.0, 20), 'capacitance': (0.0, 3.0)},
                'y_ticks': {'conductance': np.arange(0, 21, 5), 'capacitance': np.arange(0, 3.1, 0.5)},
                'series_to_plot': self.xue_series + self.vance_series + self.deconti_series,
            },
            'fig425a': {
                'suptitle': 'Figure 4.25: P.u.l. Mutual-admittance between phase - a and phase - b sheaths [Xue, 2018]',
                'p': 1, 'q': 3,
                'y_lim': {'conductance': (-8.0, 2.0), 'capacitance': (-0.8, 0.2)},
                'y_ticks': {'conductance': np.arange(-8.0, 2.1, 2), 'capacitance': np.arange(-0.8, 0.3, 0.2)},
                'series_to_plot': self.xue_series + self.vance_series + self.deconti_series,
            },
            'fig425b': {
                'suptitle': 'Figure 4.25a: P.u.l. Mutual-admittance between phase - a and phase - c sheaths [Xue, 2018]',
                'p': 1, 'q': 5,
                'y_lim': {'conductance': (-4.0, 2.0), 'capacitance': (-0.6, 0.2)},
                'y_ticks': {'conductance': np.arange(-4.0, 2.1, 2), 'capacitance': np.arange(-0.6, 0.3, 0.2)},
                'series_to_plot': self.xue_series + self.vance_series + self.deconti_series,
            },
            'earth_return_impedance_self': {
                'suptitle': 'P.u.l. Self earth-return impedance of phase - a [Xue, 2018]',
                'p': 0, 'q': 0,
                'series_to_plot': self.xue_series,
                'comsol_series_to_plot': comsol_layout_template
            },
            'earth_return_impedance_mutual_ab': {
                'suptitle': 'P.u.l. Mutual earth-return impedance between phase - a and phase - b [Xue, 2018]',
                'p': 0, 'q': 1,
                'series_to_plot': self.xue_series,
                'comsol_series_to_plot': comsol_layout_template
            },
            'earth_return_impedance_mutual_ac': {
                'suptitle': 'P.u.l. Mutual earth-return impedance between phase - a and phase - c [Xue, 2018]',
                'p': 0, 'q': 2,
                'series_to_plot': self.xue_series,
                'comsol_series_to_plot': comsol_layout_template
            },
            'earth_return_admittance_self': {
                'suptitle': 'P.u.l. Self earth-return admittance of phase - a [Xue, 2018]',
                'p': 0, 'q': 0,
                'series_to_plot': self.xue_series + self.vance_series,
                'comsol_series_to_plot': comsol_layout_template,
            },
            'earth_return_admittance_mutual_ab': {
                'suptitle': 'P.u.l. Mutual earth-return admittance between phase - a and phase - b [Xue, 2018]',
                'p': 0, 'q': 1,
                'series_to_plot': self.xue_series + self.vance_series,
                'comsol_series_to_plot': comsol_layout_template,
            },
            'earth_return_admittance_mutual_ac': {
                'suptitle': 'P.u.l. Mutual earth-return admittance between phase - a and phase - c [Xue, 2018]',
                'p': 0, 'q': 2,
                'series_to_plot': self.xue_series + self.vance_series,
                'comsol_series_to_plot': comsol_layout_template,
            },
        }

    def _scc_series_impedance_subplots(self, graph_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configuration[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)
        
        for series in config['series_to_plot']:
            zs = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['series_impedance_matrix']
            ax1.plot(self.f, np.real(zs[:, p, q]) * 1e3, **series['type']['main'])
            ax2.plot(self.f, np.imag(zs[:, p, q]) / self.w * 1e6, **series['type']['main'])
            
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(self.xlim)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$Rs_{{{p+1}{q+1}}} \, (\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$Ls_{{{p+1}{q+1}}} \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'series_impedance_{graph_key}')

    def _scc_shunt_admittance_subplots(self, graph_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configuration[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)
        
        for series in config['series_to_plot']:
            Ysh = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['shunt_admittance_matrix']
            ax1.plot(self.f, np.real(Ysh[:, p, q]) * 1e3, **series['type']['main'])
            ax2.plot(self.f, np.imag(Ysh[:, p, q]) / self.w * 1e9, **series['type']['main'])

        ax1.set_xscale('log')
        ax1.set_xlim(self.xlim)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$G_{{{p+1}{q+1}}} \, (S/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$C_{{{p+1}{q+1}}} \, (\mu F/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'shunt_admittance_{graph_key}')

    def _scc_earth_return_impedance_subplots(self, graph_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configuration[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)        
        
        for series in config['series_to_plot']:
            zg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['impedance_matrix']
            ax1.plot(self.f, np.real(zg[:, p, q]) * 1e3, **series['type']['main'])
            ax2.plot(self.f, np.imag(zg[:, p, q]) / self.w * 1e6, **series['type']['main'])

        # --- Bloco COMSOL MODIFICADO ---
        # Lendo dados pré-processados de self.cmsl (pul_data['comsol'])
        if self.cmsl is not None and 'comsol_series_to_plot' in config:
            for series in config['comsol_series_to_plot']:
                # Obter a chave de dados, ex: 'p100_er1'
                comsol_key = series.get('base_key')
                assert comsol_key is not None, "base_key must be provided in comsol_series_to_plot"

                # Acessar os dados COMSOL pré-processados
                if comsol_key in self.cmsl:
                    comsol_data = self.cmsl[comsol_key]
                    freq = comsol_data['freq']
                    cmsl_zg = comsol_data['Zg']
                    
                    # Obter o estilo de plotagem, removendo a chave de dados
                    style = {k: v for k, v in series.items() if k != 'base_key'}
                    
                    # Plotar os dados
                    ax1.scatter(freq, np.real(cmsl_zg[:, p, q]) * 1e3, **style)
                    ax2.scatter(freq, np.imag(cmsl_zg[:, p, q]) / (2 * np.pi * freq) * 1e6, **style)
                else:
                    print(f"Aviso: Chave de dados COMSOL '{comsol_key}' não encontrada em pul_data['comsol'].")
            
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(self.xlim)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$Rg_{{{p+1}{q+1}}} \, (\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$Lg_{{{p+1}{q+1}}} \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')

    def _scc_earth_return_admittance_subplots(self, graph_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configuration[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)        
        
        for series in config['series_to_plot']:
            Yg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['admittance_matrix']
            ax1.plot(self.f, np.real(Yg[:, p, q]) * 1e3, **series['type']['main'])
            ax2.plot(self.f, np.imag(Yg[:, p, q]) / self.w * 1e6, **series['type']['main'])

        # --- Bloco COMSOL MODIFICADO ---
        # Lendo dados pré-processados de self.cmsl (pul_data['comsol'])
        if self.cmsl is not None and 'comsol_series_to_plot' in config:
            for series in config['comsol_series_to_plot']:
                # Obter a chave de dados, ex: 'p100_er1'
                comsol_key = series.get('base_key')
                assert comsol_key is not None, "base_key must be provided in comsol_series_to_plot"

                # Acessar os dados COMSOL pré-processados
                if comsol_key in self.cmsl:
                    comsol_data = self.cmsl[comsol_key]
                    freq = comsol_data['freq']
                    cmsl_yg = comsol_data['Yg'] # <-- Mudança aqui para Yg
                    
                    # Obter o estilo de plotagem, removendo a chave de dados
                    style = {k: v for k, v in series.items() if k != 'base_key'}
                    
                    # Plotar os dados
                    ax1.scatter(freq, np.real(cmsl_yg[:, p, q]) * 1e3, **style)
                    ax2.scatter(freq, np.imag(cmsl_yg[:, p, q]) / (2 * np.pi * freq) * 1e6, **style)
                else:
                    print(f"Aviso: Chave de dados COMSOL '{comsol_key}' não encontrada em pul_data['comsol'].")

        ax1.set_xscale('log')
        ax1.set_xlim(self.xlim)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$Gg_{{{p+1}{q+1}}} \, (S/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$Cg_{{{p+1}{q+1}}} \, (mF/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')

    def _overhead_series_impedance_subplots(self, graph_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configuration[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zs = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['series_impedance_matrix']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, zs[:, p, q].real * 1e3, **style)
            ax2.plot(self.f, zs[:, p, q].imag / self.w * 1e6, **style)

        if self.cmsl is not None:
            cmsl = self.cmsl['cmsl_ground_return_impedance_h5']  
            zs = cmsl['coil_impedance']
            f, w = cmsl['freq'], 2 * np.pi * cmsl['freq']
            ax1.scatter(f, np.real(zs) * 1e3, label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)
            ax2.scatter(f, np.imag(zs) / w * 1e6, label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(self.xlim)
        # ax1.set_ylim(config['y_lim']['resistance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_s \, (\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        # ax2.set_ylim(config['y_lim']['inductance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_s \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _overhead_shunt_admittance_subplots(self, graph_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configuration[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            ysh = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['shunt_admittance_matrix']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, np.real(ysh[:, p, q]) * 1e3 , **style)
            ax2.plot(self.f, np.imag(ysh[:, p, q]) / self.w * 1e12, **style)

        # Configure left subplot (Conductance)
        ax1.set_xscale('log')
        ax1.set_xlim(self.xlim)
        # ax1.set_ylim(config['y_lim']['conductance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$G \, (S/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['conductance_title'])

        # Configure right subplot (Capacitance)
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        # ax2.set_ylim(config['y_lim']['capacitance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$C \, (nF/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['capacitance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _overhead_propagation_constant_subplots(self, graph_key):
        config = self.plot_configuration[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            gamma_v = self.pul_data['scenarios'][series['key']]['propagation_voltage_matrix']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, np.real(gamma_v[:, p, q]) * 1e3, **style)
            ax2.plot(self.f, self.w / np.imag(gamma_v[:, p, q]) / sc.c, **style)

        # Configure left subplot (attenuation)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(self.xlim)
        # ax1.set_ylim(config['y_lim']['attenuation'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel( r'Attenuation Constant, $\alpha_{\nu}$ (Np/km)')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['attenuation_title'])

        # Configure right subplot (phase velocity)
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        # ax2.set_ylim(config['y_lim']['phase_velocity'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Phase Velocity, $c_{\nu}/c_0$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['phase_velocity_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def overhead_series_impedance_matrix(self, graph_key_list):
        for key in graph_key_list:
            self._overhead_series_impedance_subplots(key)

    def overhead_shunt_admittance_matrix(self, graph_key_list):
        for key in graph_key_list:
            self._overhead_shunt_admittance_subplots(key)

    def overhead_propagation_constant(self, graph_key_list):
        for key in graph_key_list:
            self._overhead_propagation_constant_subplots(key)

    def scc_series_impedance_matrix(self, graph_key_list):
        for key in graph_key_list:
            self._scc_series_impedance_subplots(key)

    def scc_shunt_admittance_matrix(self, graph_key_list):
        for key in graph_key_list:
            self._scc_shunt_admittance_subplots(key)

    def scc_earth_return_impedance_matrix(self):
        for key in ['earth_return_impedance_self', 'earth_return_impedance_mutual_ab', 'earth_return_impedance_mutual_ac']:
            self._scc_earth_return_impedance_subplots(key)

    def scc_earth_return_admittance_matrix(self):
        for key in ['earth_return_admittance_self', 'earth_return_admittance_mutual_ab', 'earth_return_admittance_mutual_ac']:
            self._scc_earth_return_admittance_subplots(key)
