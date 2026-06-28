import os
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from utils.case_utils import save_figure

# Chaves retornadas por PerUnitParameters.pul_matrices()
_SERIES_Z  = 'series_impedance_matrix'
_SHUNT_Y   = 'shunt_admittance_matrix'
_EARTH_Z   = 'earth-return_impedance_matrix'
_GAMMA_V   = 'propagation_voltage_matrix'
_ZC        = 'characteristic_impedance_matrix'


class OverheadLineModels:
    """
    Camada de dados e plotagem para linhas aéreas de transmissão (OHTL).
    Análogo a SingleCoreCableModels (scc_models.py) para cabos subterrâneos.

    Espera pul_data com formato FLAT:
        pul_data[scenario_key][matrix_key]  →  ndarray (n_freq, N, N)
        pul_data['frequencies']              →  ndarray (n_freq,)
        pul_data['comsol']                   →  dict | None

    Cada serie em series_to_plot tem a forma:
        {'key': <str>, 'type': {'main': {matplotlib kwargs}}}
    """

    def __init__(self, file_path: str, pul_data: dict, autoSave: bool = True):
        self.script_path = Path(file_path)
        self.autoSave    = autoSave
        self.pul_data    = pul_data
        self.f           = pul_data['frequencies']
        self.w           = 2 * np.pi * self.f
        self.cmsl        = pul_data.get('comsol')
        self.figsize     = (12, 5)
        self.xlim        = (self.f[0], self.f[-1])

        self.results_dir = os.path.join('testData', self.script_path.stem, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

        # ------------------------------------------------------------------
        # Séries genéricas — formulações de impedância de retorno por terra
        # (usadas em ohtl_single_xue, ohtl_single_lima, ohtl_single_deConti)
        # ------------------------------------------------------------------

        self.xue_series = [
            {
                'key': 'p100',
                'type': {
                    'main': {'label': r'Nakagawa — $\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=1$',
                             'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                }
            },
            {
                'key': 'p100_er20',
                'type': {
                    'main': {'label': r'Nakagawa — $\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=20$',
                             'color': 'black', 'linestyle': '--', 'linewidth': 1.5},
                }
            },
            {
                'key': 'p2000',
                'type': {
                    'main': {'label': r'Nakagawa — $\rho_e=2000\;\Omega\mathrm{m},\;\epsilon_r=1$',
                             'color': 'black', 'linestyle': '-.', 'linewidth': 1.5},
                }
            },
            {
                'key': 'p100_carson',
                'type': {
                    'main': {'label': r'Carson — $\rho_e=100\;\Omega\mathrm{m}$',
                             'color': 'red', 'linestyle': '-', 'linewidth': 1.0},
                }
            },
            {
                'key': 'p2000_carson',
                'type': {
                    'main': {'label': r'Carson — $\rho_e=2000\;\Omega\mathrm{m}$',
                             'color': 'red', 'linestyle': '-.', 'linewidth': 1.0},
                }
            },
        ]

        self.lima_series = [
            {
                'key': 'quasi_tem',
                'type': {
                    'main': {'label': 'Quasi-TEM (Integral Eq.)',
                             'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                }
            },
            {
                'key': 'quasi_tem_log',
                'type': {
                    'main': {'label': 'Quasi-TEM (Approx. Log.)',
                             'color': 'blue', 'linestyle': '--', 'linewidth': 1.5},
                }
            },
            {
                'key': 'nakagawa',
                'type': {
                    'main': {'label': 'Wise/Nakagawa (1981)',
                             'color': 'red', 'linestyle': '-', 'linewidth': 1.5},
                }
            },
            {
                'key': 'sunde',
                'type': {
                    'main': {'label': 'Sunde (1968)',
                             'color': 'darkgreen', 'linestyle': '-', 'linewidth': 1.5},
                }
            },
            {
                'key': 'carson',
                'type': {
                    'main': {'label': 'Carson (1926)',
                             'color': 'black', 'linestyle': ':', 'linewidth': 1.5},
                }
            },
        ]

        self.nakagawa_resistivity = [
            {
                'key': 'nakagawa_1u',
                'type': {
                    'main': {'label': r'$\rho_g = 1\;\mu\Omega\mathrm{m}$',
                             'color': 'black', 'linestyle': '--', 'linewidth': 1.5},
                }
            },
            {
                'key': 'nakagawa',
                'type': {
                    'main': {'label': r'$\rho_g = 200\;\Omega\mathrm{m}$',
                             'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                }
            },
            {
                'key': 'nakagawa_5k',
                'type': {
                    'main': {'label': r'$\rho_g = 5000\;\Omega\mathrm{m}$',
                             'color': 'red', 'linestyle': ':', 'linewidth': 1.5},
                }
            },
        ]

        # ------------------------------------------------------------------
        # Configurações de gráficos
        # ------------------------------------------------------------------

        self.xue_plot_configs = {
            'series_impedance_matrix': {
                'suptitle': r'P.u.l. series impedance — single overhead line [Xue, 2012]',
                'resistance_title': r'P.u.l. series resistance',
                'inductance_title': r'P.u.l. series inductance',
                'matrix_key': _SERIES_Z,
                'p': 0, 'q': 0,
                'series_to_plot': self.xue_series,
            },
            'shunt_admittance_matrix': {
                'suptitle': r'P.u.l. shunt admittance — single overhead line [Xue, 2012]',
                'conductance_title': r'P.u.l. shunt conductance',
                'capacitance_title': r'P.u.l. shunt capacitance',
                'matrix_key': _SHUNT_Y,
                'p': 0, 'q': 0,
                'series_to_plot': self.xue_series,
            },
            'earth_return_impedance': {
                'suptitle': r'P.u.l. earth-return impedance — single overhead line [Xue, 2012]',
                'resistance_title': r'P.u.l. earth-return resistance',
                'inductance_title': r'P.u.l. earth-return inductance',
                'matrix_key': _EARTH_Z,
                'p': 0, 'q': 0,
                'series_to_plot': self.xue_series,
            },
        }

        self.lima_plot_configs = {
            'propagation_constant_forms': {
                'suptitle': (r'Propagation constant — single overhead line'
                             '\n'
                             r'$r_1=0.01\,\mathrm{m},\;h_1=10\,\mathrm{m},\;'
                             r'\sigma=6.496\times10^7\,\mathrm{S/m}$ [Lima, 2015]'),
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized phase velocity',
                'matrix_key': _GAMMA_V,
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1e2, 1e10), 'phase_velocity': (1e0, 1e10)},
                'y_lim': {'attenuation': (0, 4), 'phase_velocity': (0.5, 1.1)},
                'series_to_plot': self.lima_series,
            },
            'propagation_constant_nakagawa': {
                'suptitle': (r'Nakagawa (1981) propagation constant — single overhead line'
                             '\n'
                             r'$r_1=0.01\,\mathrm{m},\;h_1=10\,\mathrm{m},\;'
                             r'\sigma=6.496\times10^7\,\mathrm{S/m}$ [Lima, 2015]'),
                'attenuation_title': 'Attenuation constant',
                'phase_velocity_title': 'Normalized phase velocity',
                'matrix_key': _GAMMA_V,
                'p': 0, 'q': 0,
                'x_lim': {'attenuation': (1e2, 1e10), 'phase_velocity': (1e0, 1e10)},
                'y_lim': {'attenuation': (0, 1.6), 'phase_velocity': (0.5, 1.1)},
                'series_to_plot': self.nakagawa_resistivity,
            },
            'series_impedance_matrix': {
                'suptitle': r'P.u.l. series impedance — single overhead line [Lima, 2015]',
                'resistance_title': r'P.u.l. series resistance',
                'inductance_title': r'P.u.l. series inductance',
                'matrix_key': _SERIES_Z,
                'p': 0, 'q': 0,
                'series_to_plot': self.lima_series,
            },
        }

    # ------------------------------------------------------------------
    # Métodos de plotagem internos
    # ------------------------------------------------------------------

    def _impedance_subplots(self, graph_key: str, configs: dict):
        """Plota resistance (esq.) e inductance (dir.) a partir de matriz complexa."""
        cfg = configs[graph_key]
        p, q = cfg['p'], cfg['q']
        matrix_key = cfg['matrix_key']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)

        for series in cfg['series_to_plot']:
            key  = series['key']
            data = self.pul_data.get(key)
            if data is None or matrix_key not in data:
                continue
            Z   = data[matrix_key]
            sty = series['type']['main']
            ax1.plot(self.f, Z[:, p, q].real * 1e3, **sty)
            ax2.plot(self.f, Z[:, p, q].imag / self.w * 1e6, **sty)

        for ax, yl in [(ax1, cfg.get('resistance_title', r'$R_s\;(\Omega/\mathrm{km})$')),
                       (ax2, cfg.get('inductance_title', r'$L_s\;(\mathrm{mH/km})$'))]:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(self.xlim)
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel(yl)
            ax.legend(fontsize='small')
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def _admittance_subplots(self, graph_key: str, configs: dict):
        """Plota conductance (esq.) e capacitance (dir.) a partir de matriz complexa."""
        cfg = configs[graph_key]
        p, q = cfg['p'], cfg['q']
        matrix_key = cfg['matrix_key']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)

        for series in cfg['series_to_plot']:
            key  = series['key']
            data = self.pul_data.get(key)
            if data is None or matrix_key not in data:
                continue
            Y   = data[matrix_key]
            sty = series['type']['main']
            ax1.plot(self.f, Y[:, p, q].real * 1e3, **sty)
            ax2.plot(self.f, Y[:, p, q].imag / self.w * 1e9, **sty)

        for ax, yl in [(ax1, cfg.get('conductance_title', r'$G_s\;(\mathrm{mS/km})$')),
                       (ax2, cfg.get('capacitance_title', r'$C_s\;(\mathrm{nF/km})$'))]:
            ax.set_xscale('log')
            ax.set_xlim(self.xlim)
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel(yl)
            ax.legend(fontsize='small')
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def _propagation_subplots(self, graph_key: str, configs: dict):
        """Plota attenuation (esq.) e phase velocity normalizada (dir.)."""
        import scipy.constants as sc
        cfg = configs[graph_key]
        p, q = cfg['p'], cfg['q']
        matrix_key = cfg['matrix_key']
        x_lim = cfg.get('x_lim', {})
        y_lim = cfg.get('y_lim', {})

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)

        for series in cfg['series_to_plot']:
            key  = series['key']
            data = self.pul_data.get(key)
            if data is None or matrix_key not in data:
                continue
            gamma_v = data[matrix_key][:, p, q]
            sty     = series['type']['main']
            ax1.plot(self.f, gamma_v.real * 1e3, **sty)
            ax2.plot(self.f, self.w / gamma_v.imag / sc.c, **sty)

        ax1.set_xscale('log')
        ax1.set_xlim(x_lim.get('attenuation', self.xlim))
        if 'attenuation' in y_lim:
            ax1.set_ylim(y_lim['attenuation'])
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$\alpha_\nu\;(\mathrm{Np/km})$')
        ax1.set_title(cfg.get('attenuation_title', 'Attenuation constant'))
        ax1.legend(fontsize='small')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

        ax2.set_xscale('log')
        ax2.set_xlim(x_lim.get('phase_velocity', self.xlim))
        if 'phase_velocity' in y_lim:
            ax2.set_ylim(y_lim['phase_velocity'])
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$c_\nu/c_0$')
        ax2.set_title(cfg.get('phase_velocity_title', 'Normalized phase velocity'))
        ax2.legend(fontsize='small')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    # ------------------------------------------------------------------
    # Métodos públicos — Xue
    # ------------------------------------------------------------------

    def series_impedance_matrix(self):
        self._impedance_subplots('series_impedance_matrix', self.xue_plot_configs)

    def shunt_admittance_matrix(self):
        self._admittance_subplots('shunt_admittance_matrix', self.xue_plot_configs)

    def earth_return_impedance(self):
        self._impedance_subplots('earth_return_impedance', self.xue_plot_configs)

    # ------------------------------------------------------------------
    # Métodos públicos — Lima
    # ------------------------------------------------------------------

    def propagation_constant_forms(self):
        self._propagation_subplots('propagation_constant_forms', self.lima_plot_configs)

    def propagation_constant_nakagawa(self):
        self._propagation_subplots('propagation_constant_nakagawa', self.lima_plot_configs)

    def lima_series_impedance_matrix(self):
        self._impedance_subplots('series_impedance_matrix', self.lima_plot_configs)
