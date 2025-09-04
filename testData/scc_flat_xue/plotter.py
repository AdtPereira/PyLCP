import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

class ModelPlotter:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    This version is adapted for the vectorized data structure.
    """
    def __init__(self, pul_parameters):
        self.pul_data = pul_parameters
        self.f = pul_parameters['frequencies']
        self.w = 2 * np.pi * self.f

        # ... (the rest of the __init__ method remains unchanged) ...
        self.xue_series = [
            {'key': 'p100_deconti',      'label': 'De Conti Approx.', 'color': 'red', 'linestyle': ':'},
            {'key': 'p100_er20_deconti', 'label': '', 'color': 'red', 'linestyle': ':'},
            {'key': 'p500_deconti',      'label': '', 'color': 'red', 'linestyle': ':'},
            {'key': 'p100',      'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20', 'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--'},
            {'key': 'p500',      'label': r'$\rho_e=500 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-.'},
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
            'fig419': {
                'suptitle': 'Figure 4.19: P.u.l. Self-impedance of phase - a sheath with Magalhães/Xue formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, # Sheath of phase-a
                'q': 1, # Sheath of phase-a
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (0.5, 2.0)},
                'series_to_plot': self.xue_series
            },
            'fig421': {
                'suptitle': 'Figure 4.21: P.u.l. Mutual impedance between phase - a and phase - b sheaths with Magalhães/Xue formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, # Sheath of phase-a
                'q': 3, # Sheath of phase-b
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (0.0, 1.5)},
                'series_to_plot': self.xue_series
            },
            'fig423': {
                'suptitle': 'Figure 4.23: P.u.l. Self-admittance of phase - a sheath with Magalhães/Xue formulation [1]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 1, # Sheath of phase-a
                'q': 1, # Sheath of phase-a
                'x_lim': {'conductance': (1E3, 1E7), 'capacitance': (1E3, 1E7)},
                'y_lim': {'conductance': (0.0, 20), 'capacitance': (0.0, 3.0)},
                'series_to_plot': self.xue_series
            },
        }

    def _plot_impedance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            # --- MODIFICATION START ---
            # Instead of a list comprehension, we now use direct NumPy slicing.
            # self.pul_data[series['key']] is a dictionary containing the 3D matrices.
            # We get the 3D matrix and slice it: [all_frequencies, row_p, col_q]
            zs_matrix_3d = self.pul_data[series['key']]['series_impedance_matrix']
            zs = zs_matrix_3d[:, p, q]
            # --- MODIFICATION END ---
            
            rs = np.real(zs) * 1e3  # Resistance in Ohm/km
            ls = np.imag(zs) / self.w * 1e6  # Inductance in mH/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': 1.0}
            ax1.plot(self.f, rs, **style)
            ax2.plot(self.f, ls, **style)

        # ... (the rest of the method for plotting axes remains unchanged) ...
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
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$Ls_{{{p+1}{q+1}}} \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _plot_admittance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            # --- MODIFICATION START ---
            # Using direct NumPy slicing for admittance as well.
            ysh_matrix_3d = self.pul_data[series['key']]['shunt_admittance_matrix']
            ysh = ysh_matrix_3d[:, p, q]
            # --- MODIFICATION END ---

            g_pq = np.real(ysh) * 1e3  # Conductance in S/km
            c_pq = np.imag(ysh) / self.w * 1e9  # Capacitance in uF/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, g_pq, **style)
            ax2.plot(self.f, c_pq, **style)

        # ... (the rest of the method for plotting axes remains unchanged) ...
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
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$C_{{{p+1}{q+1}}} \, (\mu F/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['capacitance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def plot_fig419(self):
        """Plots the data corresponding to Figure 4.19 from the reference."""
        self._plot_impedance_subplots('fig419')

    def plot_fig421(self):
        """Plots the data corresponding to Figure 4.21 from the reference."""
        self._plot_impedance_subplots('fig421')

    def plot_fig423(self):
        """Plots the data corresponding to Figure 4.23 from the reference."""
        self._plot_admittance_subplots('fig423')