import os
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from utils.case_utils import *

class SingleCoreCableModels:
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
        
        self.rho_1 = 1000
        self.epsr_1 = 1.0
        self.figsize = (12, 5)
        self.xlim = tuple(pul_data['frequencies'][[0, -1]])

        def _generate_comsol_series(template, coil_suffix):
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

        self.impedance_series_to_plot = [
            {
                'key': 'p100_xue',
                'type': {
                    'self-core': {'label': 'Self-Core (Integral Form.)', 'color': 'red', 'linestyle': '-', 'linewidth': 2.0},
                    'self-sheath': {'label': 'Self-Sheath (Integral Form.)', 'color': 'blue', 'linestyle': '-', 'linewidth': 2.0},
                    'mutual-core-sheath': {'label': 'Mutual Core-Sheath (Integral Form.)', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 2.0}
                }
            },
            {
                'key': 'p100_deconti',
                'type': {
                    'self-core': {'label': 'De Conti Approx. Form.', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
                    'self-sheath': {'label': '', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
                    'mutual-core-sheath': {'label': '', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
                }
            },
        ] 

        self.admittance_series_to_plot = [
            {
                'key': 'p100_xue',
                'type': {
                    'self-core': {'label': 'Self-Core (Integral Form.)', 'color': 'red', 'linestyle': '-', 'linewidth': 2.0},
                    'self-sheath': {'label': 'Self-Sheath (Integral Form.)', 'color': 'blue', 'linestyle': '-', 'linewidth': 2.0},
                    'mutual-core-sheath': {'label': 'Mutual Core-Sheath (Integral Form.)', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 2.0}
                }
            },
            {
                'key': 'p100_vance',
                'type': {
                    'self-core': {'label': 'Vance Approx. Form.', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
                    'self-sheath': {'label': '', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
                    'mutual-core-sheath': {'label': '', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
                }
            },
        ]

        self.scc_series_params = [
            {
                'key': 'p100_xue',
                'type': {
                    'series_term': {'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
                    'shunt_term': {'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
                    'earth_return_term': {'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0}
                }
            },
            {
                'key': 'p100_deconti',
                'type': {
                    'series_term': {'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
                    'shunt_term': {'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
                    'earth_return_term': {'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0}
                }
            },
        ]    

        self.scc_shunt_params = [
            {
                'key': 'p100_xue',
                'type': {
                    'series_term': {'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
                    'shunt_term': {'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
                    'earth_return_term': {'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0}
                }
            },
            {
                'key': 'p100_deconti',
                'type': {
                    'series_term': {'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
                    'shunt_term': {'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
                    'earth_return_term': {'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0}
                }
            },
            {
                'key': 'p100_vance',
                'type': {
                    'series_term': {'label': 'Vance Approx. Form.', 'color': 'blue', 'linestyle': '-.', 'linewidth': 1.0},
                    'shunt_term': {'label': 'Vance Approx. Form.', 'color': 'blue', 'linestyle': '-.', 'linewidth': 1.0},
                    'earth_return_term': {'label': 'Vance Approx. Form.', 'color': 'blue', 'linestyle': '-.', 'linewidth': 1.0}
                }
            },
        ] 

        self.scc_params = [
            {
                'key': 'p100_xue',
                'type': {
                    'series_term': {
                        'label': 'Series Impedance', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0, 'zorder': 1},
                    'shunt_term': {
                        'label': 'Shunt Admittance', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0, 'zorder': 1},
                    'potential_term': {
                        'label': 'Equivalent Potential Coeff.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0, 'zorder': 1},
                    'internal_term': {
                        'label': 'Internal', 'color': 'blue', 'linestyle': '--', 'linewidth': 1.5, 'zorder': 2},
                    'earth_return_term': {
                        'label': 'Earth-return (Integral Form.)', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.5, 'zorder': 3},
                    'series_composition': {
                        'label': 'Series Composition', 'color': 'red', 'linestyle': ':', 'linewidth': 1.5, 'zorder': 4},
                }
            },
        ]  

        self.xue_plot_configs = {
            'series_impedance_matrix': {
                'suptitle': r'P.u.l. series impedance matrix of single core coaxial cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'series_to_plot': self.impedance_series_to_plot,
            },
            'shunt_admittance_matrix': {
                'suptitle': r'P.u.l. shunt admittance matrix of single core coaxial cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'series_to_plot': self.admittance_series_to_plot,
            },
            'series_impedance_composition_core': {
                'suptitle': r'P.u.l. core self-impedance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'resistance_title': r'P.u.l. series resistance',
                'inductance_title': r'P.u.l. series inductance',
                'series_to_plot': self.scc_params,
                'comsol_series_to_plot': _generate_comsol_series(comsol_layout_template, 'vcoil_1')
            },
            'series_impedance_composition_sheath': {
                'suptitle': r'P.u.l. sheath self-impedance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'resistance_title': r'P.u.l. series resistance',
                'inductance_title': r'P.u.l. series inductance',
                'series_to_plot': self.scc_params,
                'comsol_series_to_plot': _generate_comsol_series(comsol_layout_template, 'vcoil_2')
            },
            'series_impedance_composition_core_sheath': {
                'suptitle': r'P.u.l. core-sheath mutual impedance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'resistance_title': r'P.u.l. series resistance',
                'inductance_title': r'P.u.l. series inductance',
                'series_to_plot': self.scc_params,
                'comsol_series_to_plot': _generate_comsol_series(comsol_layout_template, 'vcoil_2')
            },
            'shunt_admittance_composition_core': {
                'suptitle': r'P.u.l. core self-admittance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'conductance_title': r'P.u.l. shunt conductance',
                'capacitance_title': r'P.u.l. shunt capacitance',
                'y_lim': {'conductance': (0, 20), 'capacitance': (0, 3)},
                'series_to_plot': self.scc_params,
            },
            'shunt_admittance_composition_sheath': {
                'suptitle': r'P.u.l. sheath self-admittance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'conductance_title': r'P.u.l. shunt conductance',
                'capacitance_title': r'P.u.l. shunt capacitance',
                'y_lim': {'conductance': (0, 20), 'capacitance': (0, 3)},
                'series_to_plot': self.scc_params,
            },
            'shunt_admittance_composition_core_sheath': {
                'suptitle': r'P.u.l. core-sheath mutual admittance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'conductance_title': r'P.u.l. shunt conductance',
                'capacitance_title': r'P.u.l. shunt capacitance',
                'y_lim': {'conductance': (0, 20), 'capacitance': (0, 3)},
                'series_to_plot': self.scc_params,
            },
            'potential_coefficients_composition_core': {
                'suptitle': r'P.u.l. core potential coefficient ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'norm_title': 'Absolute Value',
                'angle_title': 'Angle of Potential Coefficient',
                'series_to_plot': self.scc_params,
            },
            'potential_coefficients_composition_sheath': {
                'suptitle': r'P.u.l. sheath potential coefficient ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'norm_title': 'Absolute Value',
                'angle_title': 'Angle of Potential Coefficient',
                'series_to_plot': self.scc_params,
            },
             'potential_coefficients_composition_core_sheath': {
                'suptitle': r'P.u.l. core-sheath potential coefficient ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'norm_title': 'Absolute Value',
                'angle_title': 'Angle of Potential Coefficient',
                'series_to_plot': self.scc_params,
            },
            'earth_return_impedance': {
                'suptitle': r'P.u.l. earth-return impedance matrix of the single-core cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'resistance_title': r'P.u.l. series resistance',
                'inductance_title': r'P.u.l. series inductance',
                'series_to_plot': self.scc_series_params,
                'comsol_series_to_plot': _generate_comsol_series(comsol_layout_template, 'vcoil_1')
            },
            'earth_return_admittance': {
                'suptitle': r'P.u.l. earth-return admittance matrix of the single-core cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'conductance_title': r'P.u.l. shunt conductance',
                'capacitance_title': r'P.u.l. shunt capacitance',
                'series_to_plot': self.scc_shunt_params,
            },
            'earth_return_potential': {
                'suptitle': r'P.u.l. earth-return potential coefficient of the single-core cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
                'norm_title': 'Absolute Value',
                'angle_title': 'Angle of Potential Coefficient',
                'series_to_plot': self.scc_series_params,
            }, 
            'internal_impedance': {
                'suptitle': r'P.u.l. internal impedance matrix of the single-core cable',
                'resistance_title': r'P.u.l. series resistance',
                'inductance_title': r'P.u.l. series inductance',
            },
            'internal_admittance': {
                'suptitle': 'P.u.l. internal admittance matrix of the single-core cable',
                'conductance_title': r'P.u.l. shunt conductance',
                'capacitance_title': r'P.u.l. shunt capacitance',
            },
            'internal_potential': {
                'suptitle': 'P.u.l. internal potential coefficient of the single-core cable',
                'norm_title': 'Absolute Value',
                'angle_title': 'Angle of Potential Coefficient',
            }, 
        }

        self.nakagawa_series = [
            {'key': 'magalhaes_xue', 'label': 'Magalhães/Xue',  'color': 'black', 'linestyle': '-', 'linewidth': 2},
            {'key': 'sunde', 'label': 'Sunde', 'color': 'gray', 'linestyle': ':', 'linewidth': 2},
            {'key': 'pollaczek', 'label': 'Pollaczek', 'color': 'darkblue', 'linestyle': '-', 'linewidth': 2},
            {'key': 'ametani', 'label': 'Ametani', 'color': 'cyan', 'linestyle': '-', 'linewidth': 1},
            {'key': 'deconti', 'label': 'De Conti', 'color': 'red', 'linestyle': '--', 'linewidth': 1},
            {'key': 'saad', 'label': 'Saad', 'color': 'red', 'linestyle': '--', 'linewidth': 1},
        ]
        
        self.prysmian_plot_configs = {
            'core': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. parameters of core conductor, $z_{{cs}}$ [Ametani, 2015]',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'series_to_plot': [
                    {'key': 'zcs', 'label': r'$z_{cs} = z_{11} + z_{12} + z_{2i}$', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
                    {'key': 'z11', 'label': r'$z_{11}$: Internal impedance of core outer surface', 'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.0},
                    {'key': 'z12', 'label': r'$z_{12}$: Core outer insulator impedance', 'color': 'darkblue',  'linestyle': '--', 'linewidth': 1.0},
                    {'key': 'z2i', 'label': r'$z_{2i}$: Internal impedance of sheath inner surface', 'color': 'red',   'linestyle': '-.', 'linewidth': 1.0},
                ]
            },
            'sheath': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. parameters of sheath conductor, $z_{{s3}} [2]$',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'series_to_plot': [
                    {'key': 'zs3', 'label': r'$z_{s3} = z_{20} + z_{23}$', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
                    {'key': 'z20', 'label': r'$z_{20}$: Internal impedance of sheath outer surface', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
                    {'key': 'z2m', 'label': r'$z_{2m}$: Sheath mutual impedance', 'color': 'darkblue',  'linestyle': '-.', 'linewidth': 1.0},
                    {'key': 'z23', 'label': r'$z_{23}$: Sheath outer insulator impedance', 'color': 'darkgreen',  'linestyle': ':', 'linewidth': 1.0},
                ]
            },
            'core_sheath': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. Internal Impedance Matrix, $[z_{{i}}]$ [2]',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'series_to_plot': [
                    {'key': 'Zcc', 'label': r'$z_{cc} = z_{cs} + z_{s3} - 2z_{2m}$: Core self-impedance ', 'color': 'black', 'linestyle': '-', 'linewidth': 1.0},
                    {'key': 'Zcs', 'label': r'$z_{cs} = z_{20} + z_{23} - z_{2m}$: Mutual impedance between the core and sheath ', 'color': 'darkblue',  'linestyle': '--', 'linewidth': 1.0},
                    {'key': 'Zss', 'label': r'$z_{ss} = z_{20} + z_{23}$: Sheath self-impedance', 'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.0},
                ]
            },
            'internal_parameters': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. parameters of SCC [2]',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'series_to_plot': [
                    {'key': 'z11', 'label': r'$z_{11}$: Internal impedance of core outer surface', 'color': 'black', 'linestyle': '-', 'linewidth': 1.0},
                    {'key': 'z2m', 'label': r'$z_{2m}$: Sheath mutual impedance', 'color': 'darkblue', 'linestyle': '-', 'linewidth': 1.0},
                    {'key': 'z2i', 'label': r'$z_{2i}$: Internal impedance of sheath inner surface', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 1.0},
                    {'key': 'z20', 'label': r'$z_{20}$: Internal impedance of sheath outer surface', 'color': 'darkgray', 'linestyle': '-', 'linewidth': 1.0},
                ]
            },
            'ground_return_impedance': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. ground-return impedance for $\rho_1={self.rho_1:.0f} \;\Omega$ m and $\epsilon_{{r1}}={self.epsr_1:.0f}$',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'p': 0, 'q': 0,
                'series_to_plot': self.nakagawa_series
            },
        }

    def _internal_impedance_subplots(self, graph_key, impedance_data):
        """Generic plotting function for impedance-like data."""
        config = self.prysmian_plot_configs[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series_config in config['series_to_plot']:
            z = impedance_data[series_config['key']]
            style = {k: series_config[k] for k in ['label', 'color', 'linestyle', 'linewidth']}            
            ax1.plot(self.f, np.real(z) * 1e3, **style)
            ax2.plot(self.f, np.imag(z) / self.w * 1e6, **style)

        # Configure axes
        ax1.set_xscale('log')
        ax1.set_xlim(self.xlim)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'Resistance $(\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.legend(fontsize='small')
        
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'Inductance $(mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'impedance_parameters_{graph_key}')

    def _impedance_subplots(self, graph_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.xue_plot_configs[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        if graph_key in ['earth_return_impedance']:
            for series in config['series_to_plot']:
                Zg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['impedance_matrix']
                ax1.plot(self.f, np.real(Zg[:, 0, 0]) * 1e3, **series['type']['earth_return_term'])
                ax2.plot(self.f, np.imag(Zg[:, 0, 0]) / self.w * 1e6, **series['type']['earth_return_term'])

            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$R \, (\Omega/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['resistance_title'])
            
            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$L \, (mH/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['inductance_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if graph_key in ['internal_impedance']:
            Ri = self.pul_data['internal_matrices']['resistance_matrix']
            Li = self.pul_data['internal_matrices']['inductance_matrix']

            elements_to_plot = [
                {'index': (0, 0), 'label': 'Self-core', 'color': 'black', 'linestyle': '-'},
                {'index': (0, 1), 'label': 'Mutual core-sheath', 'color': 'red', 'linestyle': '--'},
                {'index': (1, 1), 'label': 'Self-Sheath', 'color': 'blue', 'linestyle': ':'},
            ]

            for value in elements_to_plot:
                plot_style = {'color': value['color'], 'linestyle': value['linestyle'], 'label': value['label']}
                ax1.plot(self.f, Ri[:, value['index'][0], value['index'][1]] * 1e3, **plot_style)
                ax2.plot(self.f, Li[:, value['index'][0], value['index'][1]] * 1e6, **plot_style)

            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$R \, (\Omega/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['resistance_title'])
            
            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$L \, (mH/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['inductance_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if graph_key in ['series_impedance_matrix']:            
            for series in config['series_to_plot']:
                Zs = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['series_impedance_matrix']
                ax1.plot(self.f, np.real(Zs[:, 0, 0]) * 1e3, **series['type']['self-core'])
                ax1.plot(self.f, np.real(Zs[:, 1, 1]) * 1e3, **series['type']['self-sheath'])
                ax1.plot(self.f, np.real(Zs[:, 0, 1]) * 1e3, **series['type']['mutual-core-sheath'])
                
                ax2.plot(self.f, np.imag(Zs[:, 0, 0]) / self.w * 1e6, **series['type']['self-core'])
                ax2.plot(self.f, np.imag(Zs[:, 1, 1]) / self.w * 1e6, **series['type']['self-sheath'])
                ax2.plot(self.f, np.imag(Zs[:, 0, 1]) / self.w * 1e6, **series['type']['mutual-core-sheath'])

            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$R \, (\Omega/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            
            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$L \, (mH/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if graph_key in ['series_impedance_composition_core', 'series_impedance_composition_sheath', 'series_impedance_composition_core_sheath']:
            if graph_key == 'series_impedance_composition_core':
                idx_i, idx_j = 0, 0
            elif graph_key == 'series_impedance_composition_sheath':
                idx_i, idx_j = 1, 1
            elif graph_key == 'series_impedance_composition_core_sheath':
                idx_i, idx_j = 0, 1
            else:
                raise ValueError("Invalid graph_key for series impedance.")
            
            Zi = self.pul_data['internal_matrices']['impedance_matrix']
            for series in config['series_to_plot']:
                Zg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['impedance_matrix']
                Zs = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['series_impedance_matrix']

                # Loop to build the block matrix for each frequency
                model = self.pul_data['scenarios'][series['key']]['mtl']
                M = model.num_conductors_per_scc
                Z0 = np.zeros_like(Zi, dtype=complex)
                for i in range(len(self.f)):
                    Z0[i, :, :] = np.kron(Zg[i, :, :], np.ones((M, M)))
                Zt = Zi + Z0

                ax1.plot(self.f, np.real(Zi[:, idx_i, idx_j]) * 1e3, **series['type']['internal_term'])
                ax1.plot(self.f, np.real(Zg[:, 0, 0]) * 1e3, **series['type']['earth_return_term'])
                ax1.plot(self.f, np.real(Zs[:, idx_i, idx_j]) * 1e3, **series['type']['series_term'])
                ax1.plot(self.f, np.real(Zt[:, idx_i, idx_j]) * 1e3, **series['type']['series_composition'])

                ax2.plot(self.f, np.imag(Zi[:, idx_i, idx_j]) / self.w * 1e6, **series['type']['internal_term'])
                ax2.plot(self.f, np.imag(Zg[:, 0, 0]) / self.w * 1e6, **series['type']['earth_return_term'])
                ax2.plot(self.f, np.imag(Zs[:, idx_i, idx_j]) / self.w * 1e6, **series['type']['series_term'])
                ax2.plot(self.f, np.imag(Zt[:, idx_i, idx_j]) / self.w * 1e6, **series['type']['series_composition'])

            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$R \, (\Omega/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['resistance_title'])
            
            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$L \, (mH/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['inductance_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')

    def _admittance_subplots(self, graph_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.xue_plot_configs[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        if graph_key == 'shunt_admittance_matrix':
            for series in config['series_to_plot']:
                Ysh = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['shunt_admittance_matrix']                
                ax1.plot(self.f, np.real(Ysh[:, 0, 0]) * 1e3, **series['type']['self-core'])
                ax1.plot(self.f, np.real(Ysh[:, 1, 1]) * 1e3, **series['type']['self-sheath'])
                ax1.plot(self.f, np.real(Ysh[:, 0, 1]) * 1e3, **series['type']['mutual-core-sheath'])
                
                ax2.plot(self.f, np.imag(Ysh[:, 0, 0]) / self.w * 1e9, **series['type']['self-core'])
                ax2.plot(self.f, np.imag(Ysh[:, 1, 1]) / self.w * 1e9, **series['type']['self-sheath'])
                ax2.plot(self.f, np.imag(Ysh[:, 0, 1]) / self.w * 1e9, **series['type']['mutual-core-sheath'])

            ax1.set_xscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$G \, (S/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            
            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$C \, (\mu F/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if graph_key == 'earth_return_admittance':
            for series in config['series_to_plot']:
                Yg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['admittance_matrix']
                ax1.plot(self.f, np.real(Yg[:, 0, 0]) * 1e3, **series['type']['shunt_term'])
                ax2.plot(self.f, np.imag(Yg[:, 0, 0]) / self.w * 1e9, **series['type']['shunt_term'])
            
            ax1.set_xscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$G \, (S/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['conductance_title'])
            
            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$C \, (\mu F/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['capacitance_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if graph_key == 'internal_admittance':
            Yi = self.pul_data['internal_matrices']['shunt_admittance_matrix']
            
            elements_to_plot = [
                {'index': (0, 0), 'label': 'Self-core', 'color': 'black', 'linestyle': '-'},
                {'index': (0, 1), 'label': 'Mutual core-sheath', 'color': 'red', 'linestyle': '--'},
                {'index': (1, 1), 'label': 'Self-Sheath', 'color': 'blue', 'linestyle': ':'},
            ]

            for value in elements_to_plot:
                cond = np.real(Yi[:, value['index'][0], value['index'][1]]) * 1e3
                cap = np.imag(Yi[:, value['index'][0], value['index'][1]]) / self.w * 1e9

                plot_style = {'color': value['color'], 'linestyle': value['linestyle'], 'label': value['label']}
                ax1.plot(self.f, cond, **plot_style)
                ax2.plot(self.f, cap, **plot_style)

            ax1.set_xscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$G \, (S/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['conductance_title'])
            
            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$C \, (\mu F/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['capacitance_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if graph_key in ['shunt_admittance_composition_core', 'shunt_admittance_composition_sheath', 'shunt_admittance_composition_core_sheath']:
            if graph_key == 'shunt_admittance_composition_core':
                idx_i, idx_j = 0, 0
            elif graph_key == 'shunt_admittance_composition_sheath':
                idx_i, idx_j = 1, 1
            elif graph_key == 'shunt_admittance_composition_core_sheath':
                idx_i, idx_j = 0, 1
            else:
                raise ValueError("Invalid graph_key for shunt admittance.")
            
            Yi = self.pul_data['internal_matrices']['shunt_admittance_matrix']            
            for series in config['series_to_plot']:
                Yg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['admittance_matrix']
                Ysh = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['shunt_admittance_matrix']
                
                Yt = np.zeros_like(Yi, dtype=complex)
                model = self.pul_data['scenarios'][series['key']]['mtl']
                M = model.num_conductors_per_scc
                for i in range(len(self.f)):
                    inv_Yi = np.linalg.inv(Yi[i, :, :])
                    inv_Yg = np.linalg.inv(Yg[i, :, :]) * np.ones((M, M))
                    Yt[i, :, :] = np.linalg.inv(inv_Yi + inv_Yg)

                ax1.plot(self.f, np.real(Yi[:, idx_i, idx_j]) * 1e3, **series['type']['internal_term'])
                ax1.plot(self.f, np.real(Yg[:, 0, 0]) * 1e3, **series['type']['earth_return_term'])
                ax1.plot(self.f, np.real(Ysh[:, idx_i, idx_j]) * 1e3, **series['type']['shunt_term'])
                ax1.plot(self.f, np.real(Yt[:, idx_i, idx_j]) * 1e3, **series['type']['series_composition'])

                ax2.plot(self.f, np.imag(Yi[:, idx_i, idx_j]) / self.w * 1e9, **series['type']['internal_term'])
                ax2.plot(self.f, np.imag(Yg[:, 0, 0]) / self.w * 1e9, **series['type']['earth_return_term'])
                ax2.plot(self.f, np.imag(Ysh[:, idx_i, idx_j]) / self.w * 1e9, **series['type']['shunt_term'])
                ax2.plot(self.f, np.imag(Yt[:, idx_i, idx_j]) / self.w * 1e9, **series['type']['series_composition'])

            ax1.set_xscale('log')
            ax1.set_xlim(self.xlim)
            ax1.set_ylim(config['y_lim']['conductance'])
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$G \, (S/km)$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['conductance_title'])
            
            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.set_ylim(config['y_lim']['capacitance'])
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'$C \, (\mu F/km)$')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['capacitance_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')

    def _potential_subplots(self, graph_key):
        config = self.xue_plot_configs[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)
        
        if graph_key == 'earth_return_potential':
            for series in config['series_to_plot']:                
                Pg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['potential_coefficient']
                ax1.plot(self.f, np.abs(Pg[:, 0, 0]) * 1e-9, **series['type']['earth_return_term'])
                ax2.plot(self.f, np.angle(Pg[:, 0, 0], deg=True), **series['type']['earth_return_term'])  

            ax1.set_xscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$|P| \times 10^9 \, (\Omega m s^{{-1}})$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['norm_title'])

            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'Angle of $P$ (Degrees)')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['angle_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        elif graph_key == 'internal_potential':
            # Pega a matriz Pi, que é constante com a frequência
            Pi = self.pul_data['internal_matrices']['potential_coefficient_matrix']

            # 1. Define os elementos da matriz a serem plotados, com seus rótulos e estilos
            elements_to_plot = [
                {'index': (0, 0), 'label': '$P_{00}$', 'color': 'black', 'linestyle': '-'},
                {'index': (0, 1), 'label': '$P_{01}$', 'color': 'red',   'linestyle': '--'},
                {'index': (1, 1), 'label': '$P_{11}$', 'color': 'blue',  'linestyle': ':'},
            ]

            # 2. Itera sobre a lista para plotar cada elemento
            for value in elements_to_plot:
                # Cria arrays repetindo os valores para corresponder ao eixo de frequência
                Pi_norm = np.full(len(self.f), np.abs(Pi[value['index']])) * 1e-9
                Pi_angle = np.full(len(self.f), np.angle(Pi[value['index']], deg=True))

                # Plota os dados com o rótulo e estilo definidos
                plot_style = {'color': value['color'], 'linestyle': value['linestyle'], 'label': value['label']}
                ax1.plot(self.f, Pi_norm, **plot_style)
                ax2.plot(self.f, Pi_angle, **plot_style)

            # Configurações dos eixos (o restante do código permanece o mesmo)
            ax1.set_xscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')  # A legenda agora mostrará os rótulos definidos
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$|P| \times 10^9 \, (\Omega m s^{{-1}})$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax1.set_title(config['norm_title'])

            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')  # A legenda agora mostrará os rótulos definidos
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'Angle of $P$ (Degrees)')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax2.set_title(config['angle_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        elif graph_key in ['potential_coefficients_composition_core', 'potential_coefficients_composition_sheath', 'potential_coefficients_composition_core_sheath']:
            if graph_key == 'potential_coefficients_composition_core':
                idx_i, idx_j = 0, 0
            elif graph_key == 'potential_coefficients_composition_sheath':
                idx_i, idx_j = 1, 1
            elif graph_key == 'potential_coeff_core_sheath':
                idx_i, idx_j = 0, 1
            else:
                raise ValueError("Invalid graph_key for potential coefficients.")
            
            Pi = self.pul_data['internal_matrices']['potential_coefficient_matrix']
            for series in config['series_to_plot']:
                Pg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['potential_coefficient']
                Psh = self.pul_data['scenarios'][series['key']]['quasi_tem_matrices']['potential_coefficient']
                Pt = Pi + Pg

                ax1.plot(self.f, np.full(len(self.f), np.abs(Pi[idx_i, idx_j]))  * 1e-9, **series['type']['internal_term'])
                ax1.plot(self.f, np.abs(Pg[:, 0, 0]) * 1e-9, **series['type']['earth_return_term'])
                ax1.plot(self.f, np.abs(Psh[:, idx_i, idx_j]) * 1e-9, **series['type']['potential_term'])
                ax1.plot(self.f, np.abs(Pt[:, idx_i, idx_j]) * 1e-9, **series['type']['series_composition'])

                ax2.plot(self.f, np.full(len(self.f), np.angle(Pi[idx_i, idx_j], deg=True)), **series['type']['internal_term'])
                ax2.plot(self.f, np.angle(Pg[:, 0, 0], deg=True), **series['type']['earth_return_term'])
                ax2.plot(self.f, np.angle(Psh[:, idx_i, idx_j], deg=True), **series['type']['potential_term'])
                ax2.plot(self.f, np.angle(Pt[:, idx_i, idx_j], deg=True), **series['type']['series_composition'])

            ax1.set_xscale('log')
            ax1.set_xlim(self.xlim)
            ax1.legend(fontsize='small')
            ax1.set_xlabel('Frequency (Hz)')
            ax1.set_ylabel(fr'$|P| \times 10^9 \, (\Omega m s^{{-1}})$')
            ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
            # ax1.set_title(config['norm_title'])

            ax2.set_xscale('log')
            ax2.set_xlim(self.xlim)
            ax2.legend(fontsize='small')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.set_ylabel(fr'Angle of $P$ (Degrees)')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            # ax2.set_title(config['angle_title'])
            plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')

    def ground_return_impedance(self, graph_key='ground_return_impedance'):
        config = self.prysmian_plot_configs[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['impedance_matrix']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': series['linewidth']}
            ax1.plot(self.f, np.real(zg[:, p, q]), **style)
            ax2.plot(self.f, np.imag(zg[:, p, q]) / self.w * 1e6, **style)

        ax1.set(xscale='log', xlim=self.xlim, xlabel='Frequency (Hz)', ylabel=r'$R_g \, (\Omega/m)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.legend(fontsize='small')
        
        ax2.set(xscale='log', xlim=self.xlim, xlabel='Frequency (Hz)', ylabel=r'$L_g \, (\mu H/m)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.legend(fontsize='small')        
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'ground_return_impedance')

    def series_impedance_earth_return(self):
        self._impedance_subplots(graph_key='earth_return_impedance')
    
    def internal_impedance_parameters(self, graph_key):
        data = self.pul_data['internal_parameters']

        if graph_key == 'core':
            z11 = data['zcs']['z11']
            z12 = data['zcs']['z12']
            z2i = data['zcs']['z2i']
            impedance_data = {'zcs': z11 + z12 + z2i, 'z11': z11, 'z12': z12, 'z2i': z2i}

        elif graph_key == 'sheath':
            z20 = data['zs3']['z20']
            z23 = data['zs3']['z23']
            z2m = data['z2m']
            impedance_data = {'zs3': z20 + z23, 'z20': z20, 'z23': z23, 'z2m': z2m}

        elif graph_key == 'core_sheath':
            zcs = data['zcs']['z11'] + data['zcs']['z12'] + data['zcs']['z2i']
            zs3 = data['zs3']['z20'] + data['zs3']['z23']
            z2m = data['z2m']
            impedance_data = {
                'Zcc': zcs + zs3 - 2 * z2m,
                'Zss': zs3,
                'Zcs': zs3 - z2m
            }

        elif graph_key == 'internal_parameters':
            impedance_data = {
                'z11': data['zcs']['z11'],
                'z2m': data['z2m'],
                'z2i': data['zcs']['z2i'],
                'z20': data['zs3']['z20']
            }
        
        self._internal_impedance_subplots(graph_key, impedance_data)
    
    def series_impedance_internal(self,):
        self._impedance_subplots(graph_key='internal_impedance')

    def series_impedance_matrix(self):
        self._impedance_subplots(graph_key='series_impedance_matrix')

    def shunt_admittance_matrix(self):
        self._admittance_subplots(graph_key='shunt_admittance_matrix')

    def series_impedance_composition(self, condutor='core'):
        self._impedance_subplots(graph_key='series_impedance_composition_{}'.format(condutor))

    def shunt_admittance_composition(self, condutor='core'):
        self._admittance_subplots(graph_key='shunt_admittance_composition_{}'.format(condutor))

    def potential_coefficients_composition(self, condutor='core'):
        self._potential_subplots(graph_key='potential_coefficients_composition_{}'.format(condutor))

    def shunt_admittance_earth_return(self):
        self._admittance_subplots(graph_key='earth_return_admittance')
    
    def shunt_admittance_internal(self):
        self._admittance_subplots(graph_key='internal_admittance')

    def potential_coefficients_earth_return(self):
        self._potential_subplots(graph_key='earth_return_potential')

    def potential_coefficients_internal(self):
        self._potential_subplots(graph_key='internal_potential')
