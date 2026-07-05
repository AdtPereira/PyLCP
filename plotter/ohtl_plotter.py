import numpy as np
import matplotlib.pyplot as plt
from utils.case_utils import save_figure
from .models_base import BasePlotter


class OHTLPlotter(BasePlotter):
    """
    Config-driven plotter for overhead transmission lines (OHTL).
    Analogous to SCCPlotter (scc_plotter.py) for underground cables.

    Expects pul_data in FLAT format:
        pul_data[scenario_key][matrix_key]  →  ndarray (n_freq, N, N)
        pul_data['frequencies']              →  ndarray (n_freq,)

    plot_config follows the same schema as BasePlotter, with additional keys:
        'series_to_plot': [{'key': <str>, 'type': {'main': {kwargs}}}]
        'matrix_key': <str>    — key inside pul_data[scenario_key]
        'p', 'q': <int>        — matrix indices
        'path': [<str>, ...]   — navigation inside pul_data[scenario_key] (optional,
                                  replaces matrix_key when the structure is deeper)
    """

    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)
        _freq = pul_data.get('frequencies')
        self._ohtl_freq = _freq
        self._ohtl_xlim = (_freq[0], _freq[-1]) if _freq is not None else self.xlim

    # ------------------------------------------------------------------
    # Internal plotting methods (analogous to _plot_scc_matrix)
    # ------------------------------------------------------------------

    def _plot_ohtl_matrix(self, graph_key: str):
        """
        Plots matricial parameters using the FLAT pul_data format.
        Supports 'matrix_key' (simple) or 'path' (arbitrary navigation).
        """
        cfg = self.plot_config[graph_key]
        p, q = cfg['p'], cfg['q']
        path = cfg.get('path') or [cfg['matrix_key']]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg, right_cfg = cfg['left_plot'], cfg['right_plot']

        freq = self._ohtl_freq
        w    = 2 * np.pi * freq

        for series in cfg.get('series_to_plot', []):
            scenario = self.pul_data.get(series['key'])
            if scenario is None:
                continue
            matrix = scenario
            try:
                for k in path:
                    matrix = matrix[k]
            except (KeyError, TypeError):
                print(f"Warning: path {path} not found for '{series['key']}'.")
                continue

            style = series['type']['main']
            ax1.plot(freq, self._calculate_plot_data(matrix[:, p, q], w, left_cfg),  **style)
            ax2.plot(freq, self._calculate_plot_data(matrix[:, p, q], w, right_cfg), **style)

        comsol_key = cfg.get('comsol_key')
        if comsol_key and self.cmsl is not None and comsol_key in self.cmsl:
            comsol_df = self.cmsl[comsol_key]
            comsol_col = cfg.get('comsol_column', 'coil_impedance')
            comsol_style = cfg.get('comsol_style', {
                'label': 'COMSOL (mf)', 'marker': 'o', 'facecolors': 'black', 's': 10, 'zorder': 2,
            })
            c_freq = comsol_df['freq'].to_numpy()
            c_w = 2 * np.pi * c_freq
            c_vals = comsol_df[comsol_col].to_numpy()
            ax1.scatter(c_freq, self._calculate_plot_data(c_vals, c_w, left_cfg),  **comsol_style)
            ax2.scatter(c_freq, self._calculate_plot_data(c_vals, c_w, right_cfg), **comsol_style)

        self._format_axis(ax1, self._ohtl_xlim, left_cfg)
        self._format_axis(ax2, self._ohtl_xlim, right_cfg)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def _plot_ohtl_propagation(self, graph_key: str):
        """
        Plots the propagation constant: attenuation (left) and phase velocity (right).
        """
        import scipy.constants as sc
        cfg  = self.plot_config[graph_key]
        p, q = cfg['p'], cfg['q']
        path = cfg.get('path') or [cfg.get('matrix_key', 'propagation_voltage_matrix')]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)

        freq = self._ohtl_freq
        w    = 2 * np.pi * freq
        x_lim = cfg.get('x_lim', {})
        y_lim = cfg.get('y_lim', {})
        y_scale = cfg.get('y_scale', {})

        for series in cfg.get('series_to_plot', []):
            scenario = self.pul_data.get(series['key'])
            if scenario is None:
                continue
            data = scenario
            try:
                for k in path:
                    data = data[k]
            except (KeyError, TypeError):
                continue

            gamma_v = data[:, p, q]
            style   = series['type']['main']
            ax1.plot(freq, gamma_v.real * 1e3, **style)
            ax2.plot(freq, w / gamma_v.imag / sc.c, **style)

        ax1.set_xscale('log')
        ax1.set_yscale(y_scale.get('attenuation', 'linear'))
        ax1.set_xlim(x_lim.get('attenuation', self._ohtl_xlim))
        if 'attenuation' in y_lim:
            ax1.set_ylim(y_lim['attenuation'])
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$\alpha_\nu\;(\mathrm{Np/km})$')
        ax1.set_title(cfg.get('attenuation_title', 'Attenuation constant'))
        ax1.legend(fontsize='small')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

        ax2.set_xscale('log')
        ax2.set_yscale(y_scale.get('phase_velocity', 'linear'))
        ax2.set_xlim(x_lim.get('phase_velocity', self._ohtl_xlim))
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
    # Public methods
    # ------------------------------------------------------------------

    def ohtl_series_impedance_matrix(self, graph_key: str = 'series_impedance_matrix'):
        self._plot_ohtl_matrix(graph_key)

    def ohtl_shunt_admittance_matrix(self, graph_key: str = 'shunt_admittance_matrix'):
        self._plot_ohtl_matrix(graph_key)

    def ohtl_earth_return_impedance(self, graph_key: str = 'earth_return_impedance'):
        self._plot_ohtl_matrix(graph_key)

    def ohtl_propagation_constant(self, graph_key: str = 'propagation_constant'):
        self._plot_ohtl_propagation(graph_key)

    def ohtl_characteristic_impedance(self, graph_key: str = 'characteristic_impedance'):
        self._plot_ohtl_matrix(graph_key)
