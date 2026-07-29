import numpy as np
import matplotlib.pyplot as plt
from utils.case_utils import save_figure
from .models_base import BasePlotter

class SCCPlotter(BasePlotter):
    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)
        _freq = pul_data.get('frequencies')
        self._scc_freq = _freq
        self._scc_xlim = (_freq[0], _freq[-1]) if _freq is not None else self.xlim

    def _plot_scc_matrix(self, graph_key):
        """Plots matricial parameters using the 'series_to_plot' + 'path' + 'p'/'q' config format."""
        cfg = self.plot_config[graph_key]
        p, q = cfg['p'], cfg['q']
        path = cfg['path']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg, right_cfg = cfg['left_plot'], cfg['right_plot']

        freq = self._scc_freq
        w = 2 * np.pi * freq
        for series in cfg.get('series_to_plot', []):
            scenario = self.pul_data['scenarios'].get(series['key'])
            if scenario is None:
                continue
            matrix = scenario
            for k in path:
                matrix = matrix[k]
            style = series['type']['main']
            ax1.plot(freq, self._calculate_plot_data(matrix[:, p, q], w, left_cfg), **style)
            ax2.plot(freq, self._calculate_plot_data(matrix[:, p, q], w, right_cfg), **style)

        if self.cmsl:
            comsol_freq = self.cmsl.get('frequencies')
            if comsol_freq is not None:
                comsol_w = 2 * np.pi * comsol_freq
                comsol_path = path[:-1] + [cfg.get('comsol_matrix_key', path[-1])]
                for series in cfg.get('comsol_series_to_plot', []):
                    key = series['key']
                    if key not in self.cmsl.get('scenarios', {}):
                        continue
                    try:
                        matrix = self.cmsl['scenarios'][key]
                        for k in comsol_path:
                            matrix = matrix[k]
                        style = {k: v for k, v in series.items() if k != 'key'}
                        ax1.scatter(comsol_freq, self._calculate_plot_data(matrix[:, p, q], comsol_w, left_cfg), **style)
                        ax2.scatter(comsol_freq, self._calculate_plot_data(matrix[:, p, q], comsol_w, right_cfg), **style)
                    except (KeyError, TypeError, IndexError):
                        print(f"Warning: COMSOL data not found for '{key}' with path {comsol_path}.")

        if self.mtlb:
            mtlb_freq = self.mtlb.get('frequencies')
            if mtlb_freq is not None:
                mtlb_w = 2 * np.pi * mtlb_freq
                mtlb_key = cfg.get('matlab_matrix_key', path[-1])
                for series in cfg.get('matlab_series_to_plot', []):
                    key = series['key']
                    if key not in self.mtlb.get('scenarios', {}):
                        continue
                    try:
                        matrix = self.mtlb['scenarios'][key][mtlb_key]
                        style = {k: v for k, v in series.items() if k != 'key'}
                        ax1.scatter(mtlb_freq, self._calculate_plot_data(matrix[:, p, q], mtlb_w, left_cfg), **style)
                        ax2.scatter(mtlb_freq, self._calculate_plot_data(matrix[:, p, q], mtlb_w, right_cfg), **style)
                    except (KeyError, TypeError, IndexError):
                        print(f"Warning: MATLAB data not found for '{key}' with key '{mtlb_key}'.")

        xlim = cfg.get('xlim', self._scc_xlim)
        self._format_axis(ax1, xlim, left_cfg)
        self._format_axis(ax2, xlim, right_cfg)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def _plot_scc_scalar(self, graph_key):
        """Plots scalar parameters (1D over frequency) using the 'series_to_plot' + 'path' config format."""
        cfg = self.plot_config[graph_key]
        path = cfg['path']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg, right_cfg = cfg['left_plot'], cfg['right_plot']

        freq = self._scc_freq
        w = 2 * np.pi * freq
        for series in cfg.get('series_to_plot', []):
            scenario = self.pul_data['scenarios'].get(series['key'])
            if scenario is None:
                continue
            data = scenario
            for k in path:
                data = data[k]
            style = series['type']['main']
            ax1.plot(freq, self._calculate_plot_data(data, w, left_cfg), **style)
            ax2.plot(freq, self._calculate_plot_data(data, w, right_cfg), **style)

        if self.cmsl:
            comsol_freq = self.cmsl.get('frequencies')
            if comsol_freq is not None:
                comsol_w = 2 * np.pi * comsol_freq
                comsol_key = cfg.get('comsol_matrix_key')
                for series in cfg.get('comsol_series_to_plot', []):
                    key = series['key']
                    if key not in self.cmsl.get('scenarios', {}):
                        continue
                    try:
                        data = self.cmsl['scenarios'][key]
                        for k in path[:-1]:
                            data = data[k]
                        data = data[comsol_key]
                        style = {k: v for k, v in series.items() if k != 'key'}
                        ax1.scatter(comsol_freq, self._calculate_plot_data(data, comsol_w, left_cfg), **style)
                        ax2.scatter(comsol_freq, self._calculate_plot_data(data, comsol_w, right_cfg), **style)
                    except (KeyError, TypeError, IndexError):
                        print(f"Warning: COMSOL data not found for '{key}' [{path[:-1]} + '{comsol_key}'].")

        self._format_axis(ax1, self._scc_xlim, left_cfg)
        self._format_axis(ax2, self._scc_xlim, right_cfg)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def _plot_scc_internal_matrix(self, graph_key):
        """Plots an analytical internal-only matrix (Zi or Yi) against the MATLAB reference matrix.

        Accepts either a single (p, q) pair via the 'p'/'q'/'internal_style' keys, or several
        overlaid on the same axes via a 'components' list of such dicts (each with its own 'p',
        'q', 'internal_style' and 'matlab_series_to_plot'). The internal matrix to plot is picked
        via 'internal_matrix_key' (defaults to 'impedance_matrix', i.e. Zi).
        """
        cfg = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg, right_cfg = cfg['left_plot'], cfg['right_plot']

        freq = self._scc_freq
        w = 2 * np.pi * freq

        components = cfg.get('components') or [{
            'p': cfg['p'], 'q': cfg['q'],
            'internal_style': cfg['internal_style'],
            'matlab_series_to_plot': cfg.get('matlab_series_to_plot', []),
        }]

        internal_matrices = self.pul_data.get('internal_matrices')
        internal_matrix_key = cfg.get('internal_matrix_key', 'impedance_matrix')
        mtlb_freq = self.mtlb.get('frequencies') if self.mtlb else None
        mtlb_w = 2 * np.pi * mtlb_freq if mtlb_freq is not None else None
        mtlb_key = cfg.get('matlab_matrix_key', 'internal_impedance_matrix')

        for comp in components:
            p, q = comp['p'], comp['q']
            if internal_matrices is not None:
                Mi = internal_matrices[internal_matrix_key]
                style = comp['internal_style']
                ax1.plot(freq, self._calculate_plot_data(Mi[:, p, q], w, left_cfg), **style)
                ax2.plot(freq, self._calculate_plot_data(Mi[:, p, q], w, right_cfg), **style)

            if self.mtlb and mtlb_freq is not None:
                for series in comp.get('matlab_series_to_plot', []):
                    key = series['key']
                    if key not in self.mtlb.get('scenarios', {}):
                        continue
                    try:
                        matrix = self.mtlb['scenarios'][key][mtlb_key]
                        style = {k: v for k, v in series.items() if k != 'key'}
                        ax1.scatter(mtlb_freq, self._calculate_plot_data(matrix[:, p, q], mtlb_w, left_cfg), **style)
                        ax2.scatter(mtlb_freq, self._calculate_plot_data(matrix[:, p, q], mtlb_w, right_cfg), **style)
                    except (KeyError, TypeError, IndexError):
                        print(f"Warning: MATLAB data not found for '{key}' with key '{mtlb_key}'.")

        xlim = cfg.get('xlim', self._scc_xlim)
        self._format_axis(ax1, xlim, left_cfg)
        self._format_axis(ax2, xlim, right_cfg)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def compare_complete_matrices(self, key_list):
        for key in key_list:
            self._plot_scc_matrix(key)

    def compare_internal_matrices(self, key_list):
        for key in key_list:
            self._plot_scc_internal_matrix(key)

    def scc_earth_propagation_constant(self):
        self._plot_scc_scalar('earth_propagation_constant')

class HDPEPlotter(BasePlotter):
    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)

    def hdpe_internal_impedance_matrix(self):
        for key in ['internal_impedance_matrix']:
            self._plot_matrix_parameters(key)

    def hdpe_internal_impedance_elements(self):
        for key in ['coaxial_cable', 'internal_impedance_elements']:
            self._plot_non_matrix_parameters(key)

class CoaxialCablePlotter(BasePlotter):
    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)

    def coaxial_cable_internal_impedance_matrix(self):
        for key in ['partial_internal_impedance']:
            self._plot_matrix_parameters(key)

    def coaxial_cable_impedance(self):
        self._plot_non_matrix_parameters('coaxial_cable_impedance')
