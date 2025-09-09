import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

class XueModels:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    This version is adapted for the vectorized data structure.
    """
    def __init__(self, pul_parameters):
        self.pul_data = pul_parameters
        self.f = pul_parameters['frequencies']
        self.w = 2 * np.pi * self.f

        self.overhead_line = [
            {'key': 'p100',         'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20',    'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--'},
            {'key': 'p2000',        'label': r'$\rho_e=2000 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-.'}
        ]

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
            'fig42': {
                'suptitle': 'Figure 4.2: P.u.l. series impedance with Nakagawa formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 0, 'q': 0,
                'x_lim': {'resistance': (1E3, 1E9), 'inductance': (1E3, 1E9)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (1, 2.5)},
                'series_to_plot': self.overhead_line
            },
            'fig43': {
                'suptitle': 'Figure 4.3: P.u.l. series impedance comparison [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 0, 'q': 0,
                'x_lim': {'resistance': (1E3, 1E9), 'inductance': (1E3, 1E9)},
                'y_lim': {'resistance': (1E0, 1E4), 'inductance': (1.4, 2.2)},
                'series_to_plot': self.nakagawa_carson_series
            },
            'fig45': {
                'suptitle': 'Figure 4.5: P.u.l. shunt admittance with Nakagawa formulation [1]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 0, 'q': 0,
                'x_lim': {'conductance': (1E3, 1E9), 'capacitance': (1E3, 1E9)},
                'y_lim': {'conductance': (-0.3, 0.1), 'capacitance': (6.6, 7.4)},
                'series_to_plot': self.overhead_line
            },
            'fig46': {
                'suptitle': 'Figure 4.6: P.u.l. shunt admittance comparison [1]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 0, 'q': 0,
                'x_lim': {'conductance': (1E3, 1E9), 'capacitance': (1E3, 1E9)},
                'y_lim': {'conductance': (-0.06, 0.02), 'capacitance': (7.22, 7.32)},
                'series_to_plot': self.nakagawa_carson_series
            },
            'fig47': {
                'suptitle': 'Figure 4.7: Propagation constant with Nakagawa formulation [1]',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1E3, 1E9), 'phase_velocity': (1E3, 1E9)},
                'y_lim': {'attenuation': (1E-3, 1E1), 'phase_velocity': (0.8, 1.1)},
                'series_to_plot': self.overhead_line
            },
            'fig48': {
                'suptitle': r'Figure 4.8: Propagation constant comparison for $\varepsilon_r = 1$ [1]',
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized Phase velocity',
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1E3, 1E9), 'phase_velocity': (1E3, 1E9)},
                'y_lim': {'attenuation': (1E-3, 1E2), 'phase_velocity': (0.8, 1.1)},
                'series_to_plot': self.attenuation_constant_series
            },
            'fig419': {
                'suptitle': 'Figure 4.19: P.u.l. Self-impedance of phase - a sheath with Magalhães/Xue formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, 'q': 1,
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (0.5, 2.0)},
                'series_to_plot': self.xue_series
            },
            'fig421': {
                'suptitle': 'Figure 4.21: P.u.l. Mutual impedance between phase - a and phase - b sheaths with Magalhães/Xue formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, 'q': 3,
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (0.0, 1.5)},
                'series_to_plot': self.xue_series
            },
            'fig423': {
                'suptitle': 'Figure 4.23: P.u.l. Self-admittance of phase - a sheath with Magalhães/Xue formulation [1]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 1, 'q': 1,
                'x_lim': {'conductance': (1E3, 1E7), 'capacitance': (1E3, 1E7)},
                'y_lim': {'conductance': (0.0, 20), 'capacitance': (0.0, 3.0)},
                'series_to_plot': self.xue_series
            },
        }

    def _scc_impedance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zs_3d = self.pul_data[series['key']]['series_impedance_matrix']
            zs = zs_3d[:, p, q]
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': 1.0}
            ax1.plot(self.f, zs.real * 1e3, **style)
            ax2.plot(self.f, zs.imag / self.w * 1e6, **style)

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

    def _scc_admittance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            ysh_3d = self.pul_data[series['key']]['shunt_admittance_matrix']
            ysh = ysh_3d[:, p, q]
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, ysh.real * 1e3, **style)
            ax2.plot(self.f, ysh.imag / self.w * 1e9, **style)

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

    def _overhead_impedance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zs_3d = self.pul_data[series['key']]['series_impedance_matrix']
            zs = zs_3d[:, p, q]            
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, zs.real * 1e3, **style)
            ax2.plot(self.f, zs.imag / self.w * 1e6, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(config['x_lim']['resistance'])
        ax1.set_ylim(config['y_lim']['resistance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_s \, (\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['inductance'])
        ax2.set_ylim(config['y_lim']['inductance'])
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
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
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
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
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

    def plot_fig419(self):
        """Plots the data corresponding to Figure 4.19 from the reference."""
        self._scc_impedance_subplots('fig419')

    def plot_fig421(self):
        """Plots the data corresponding to Figure 4.21 from the reference."""
        self._scc_impedance_subplots('fig421')

    def plot_fig423(self):
        """Plots the data corresponding to Figure 4.23 from the reference."""
        self._scc_admittance_subplots('fig423')