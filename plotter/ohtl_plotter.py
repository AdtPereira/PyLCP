import numpy as np
import matplotlib.pyplot as plt
from utils.case_utils import save_figure
from .models_base import BasePlotter


class OHTLPlotter(BasePlotter):
    """
    Plotter config-driven para linhas aéreas de transmissão (OHTL).
    Análogo a SCCPlotter (scc_plotter.py) para cabos subterrâneos.

    Espera pul_data com formato FLAT:
        pul_data[scenario_key][matrix_key]  →  ndarray (n_freq, N, N)
        pul_data['frequencies']              →  ndarray (n_freq,)

    O plot_config segue o mesmo esquema de BasePlotter, com chaves adicionais:
        'series_to_plot': [{'key': <str>, 'type': {'main': {kwargs}}}]
        'matrix_key': <str>    — chave dentro de pul_data[scenario_key]
        'p', 'q': <int>        — índices matriciais
        'path': [<str>, ...]   — navegação dentro de pul_data[scenario_key] (opcional,
                                  substitui matrix_key quando a estrutura é mais profunda)
    """

    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        super().__init__(file_path, pul_data, plot_config, autoSave=autoSave)
        _freq = pul_data.get('frequencies')
        self._ohtl_freq = _freq
        self._ohtl_xlim = (_freq[0], _freq[-1]) if _freq is not None else self.xlim

    # ------------------------------------------------------------------
    # Métodos internos de plotagem (análogos a _plot_scc_matrix)
    # ------------------------------------------------------------------

    def _plot_ohtl_matrix(self, graph_key: str):
        """
        Plota parâmetros matriciais usando formato FLAT de pul_data.
        Suporta 'matrix_key' (simples) ou 'path' (navegação arbitrária).
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
                print(f"Aviso: caminho {path} não encontrado para '{series['key']}'.")
                continue

            style = series['type']['main']
            ax1.plot(freq, self._calculate_plot_data(matrix[:, p, q], w, left_cfg),  **style)
            ax2.plot(freq, self._calculate_plot_data(matrix[:, p, q], w, right_cfg), **style)

        self._format_axis(ax1, self._ohtl_xlim, left_cfg)
        self._format_axis(ax2, self._ohtl_xlim, right_cfg)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=graph_key)

    def _plot_ohtl_propagation(self, graph_key: str):
        """
        Plota constante de propagação: attenuation (esq.) e phase velocity (dir.).
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
        ax1.set_xlim(x_lim.get('attenuation', self._ohtl_xlim))
        if 'attenuation' in y_lim:
            ax1.set_ylim(y_lim['attenuation'])
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$\alpha_\nu\;(\mathrm{Np/km})$')
        ax1.set_title(cfg.get('attenuation_title', 'Attenuation constant'))
        ax1.legend(fontsize='small')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

        ax2.set_xscale('log')
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
    # Métodos públicos
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
