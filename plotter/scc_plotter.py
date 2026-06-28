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
        """Plota parâmetros matriciais usando o formato de config 'series_to_plot' + 'path' + 'p'/'q'."""
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
                        print(f"Aviso: Dados COMSOL não encontrados para '{key}' com caminho {comsol_path}.")

        self._format_axis(ax1, self._scc_xlim, left_cfg)
        self._format_axis(ax2, self._scc_xlim, right_cfg)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def _plot_scc_scalar(self, graph_key):
        """Plota parâmetros escalares (1D por frequência) usando o formato de config 'series_to_plot' + 'path'."""
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
                        print(f"Aviso: Dados COMSOL não encontrados para '{key}' [{path[:-1]} + '{comsol_key}'].")

        self._format_axis(ax1, self._scc_xlim, left_cfg)
        self._format_axis(ax2, self._scc_xlim, right_cfg)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def scc_series_impedance_matrix(self, graph_key_list):
        for key in graph_key_list:
            self._plot_scc_matrix(key)

    def scc_shunt_admittance_matrix(self, graph_key_list):
        for key in graph_key_list:
            self._plot_scc_matrix(key)

    def scc_earth_return_impedance_matrix(self):
        for key in ['earth_return_impedance_self', 'earth_return_impedance_mutual_ab', 'earth_return_impedance_mutual_ac']:
            self._plot_scc_matrix(key)

    def scc_earth_return_admittance_matrix(self):
        for key in ['earth_return_admittance_self', 'earth_return_admittance_mutual_ab', 'earth_return_admittance_mutual_ac']:
            self._plot_scc_matrix(key)

    def scc_earth_propagation_constant(self):
        self._plot_scc_scalar('earth_propagation_constant')


class HDPEPlotter(BasePlotter):
    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)

    def hdpe_internal_impedance_matrix(self):
        for key in ['internal_impedance_matrix']:
            self._plot_matricial_parameters(key)

    def hdpe_internal_impedance_elements(self):
        for key in ['coaxial_cable', 'internal_impedance_elements']:
            self._plot_non_matricial_parameters(key)

class CoaxialCablePlotter(BasePlotter):
    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)

    def coaxial_cable_internal_impedance_matrix(self):
        for key in ['partial_internal_impedance']:
            self._plot_matricial_parameters(key)

    def coaxial_cable_impedance(self):
        self._plot_non_matricial_parameters('coaxial_cable_impedance')
