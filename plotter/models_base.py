import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from utils.case_utils import *

class BasePlotter:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    This version is adapted for the vectorized data structure.
    """
    def __init__(self, file_path: str, pul_data: dict, plot_config: dict, autoSave: bool = True):
        """
        A highly refactored class to handle plotting for the Xue model results.
        It uses a configuration-driven approach to generate complex subplot figures.
        This version is adapted for the vectorized data structure.
        """
        self.script_path = Path(file_path)
        self.pul_data = pul_data
        self.plot_config = plot_config
        self.autoSave = autoSave
        
        self.cmsl = self.pul_data.get('comsol', None)
        self.f = pul_data['frequencies']
        self.w = 2 * np.pi * self.f        

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.script_path.stem, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)
        
        self.xlim = tuple(pul_data['frequencies'][[0, -1]])
        self.figsize = (12, 5)

    def _plot_matricial_complex_quantity(self, graph_key):
        """
        Método genérico para plotar uma grandeza complexa em dois subplots 
        (ex: R/L ou G/C).
        """
        config = self.plot_config[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        # --- Configurações dos subplots (ex: 'left_plot', 'right_plot' em config)
        left_cfg = config['left_plot']
        right_cfg = config['right_plot']
        
        # --- Caminho para os dados (ex: 'data_path' em config)
        # ex: ['scenarios', '{key}', 'quasi_tem_matrices', 'series_impedance_matrix']
        data_path = config['data_path'] 
        
        for series in config['series_to_plot']:
            # Lógica para buscar os dados usando data_path e a 'key' da série
            # key = series['key']
            matrix = self.pul_data[data_path[0]][series['key']][data_path[2]][data_path[3]]
            
            # Plotar dados (Ex: Resistência)
            y1_data = self._calculate_plot_data(matrix[:, p, q], self.w, left_cfg['component'], left_cfg['scale'])
            ax1.plot(self.f, y1_data, **series['type']['main'])
            
            # Plotar dados (Ex: Indutância)
            y2_data = self._calculate_plot_data(matrix[:, p, q], self.w, right_cfg['component'], right_cfg['scale'])
            ax2.plot(self.f, y2_data, **series['type']['main'])

        # --- Bloco COMSOL Genérico ---
        if self.cmsl is not None and 'comsol_series_to_plot' in config:
            freq = self.cmsl['frequencies']
            w_cmsl = self.cmsl['angular_frequencies']
            
            for series in config['comsol_series_to_plot']:
                cmsl_key = series.get('base_key')
                cmsl_graph_type = config['comsol_matrix_key']
                
                if cmsl_key not in self.cmsl['scenarios']:
                    print(f"Aviso: Chave COMSOL '{cmsl_key}' não encontrada.")
                    continue

                if cmsl_graph_type in ['impedance_matrix', 'admittance_matrix']:
                    cmsl_matrix = self.cmsl['scenarios'][cmsl_key]['earth_return_parameters'][cmsl_graph_type]

                elif cmsl_graph_type in ['series_impedance_matrix', 'shunt_admittance_matrix']:
                    cmsl_matrix = self.cmsl['scenarios'][cmsl_key]['quasi_tem_matrices'][cmsl_graph_type]

                # Obter o estilo de plotagem
                style = {k: v for k, v in series.items() if k != 'base_key'}
                
                # Plotar dados do subplot esquerdo
                y1_data_cmsl = self._calculate_plot_data(cmsl_matrix[:, p, q], w_cmsl, left_cfg['component'], left_cfg['scale'])
                ax1.scatter(freq, y1_data_cmsl, **style)
                
                # Plotar dados do subplot direito
                y2_data_cmsl = self._calculate_plot_data(cmsl_matrix[:, p, q], w_cmsl, right_cfg['component'], right_cfg['scale'])
                ax2.scatter(freq, y2_data_cmsl, **style)

        # --- Formatação Genérica dos Eixos ---
        self._format_axis(ax1, self.f, self.xlim, left_cfg)
        self._format_axis(ax2, self.f, self.xlim, right_cfg)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')

    def _plot_scalar_complex_quantity(self, graph_key):
        """
        Método genérico para plotar uma grandeza complexa em dois subplots 
        (ex: R/L ou G/C).
        """
        config = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        # --- Configurações dos subplots (ex: 'left_plot', 'right_plot' em config)
        left_cfg = config['left_plot']
        right_cfg = config['right_plot']
        
        # --- Caminho para os dados (ex: 'data_path' em config)
        # ex: ['scenarios', '{key}', 'quasi_tem_matrices', 'series_impedance_matrix']
        data_path = config['data_path'] 
        
        for series in config['series_to_plot']:
            # Lógica para buscar os dados usando data_path e a 'key' da série
            item = self.pul_data[data_path[0]][series['key']][data_path[2]][data_path[3]]
            
            # Plotar dados do subplot esquerdo
            y1_data = self._calculate_plot_data(item, self.w, left_cfg['component'], left_cfg['scale'])
            ax1.plot(self.f, y1_data, **series['type']['main'])
            
            # Plotar dados do subplot direito
            y2_data = self._calculate_plot_data(item, self.w, right_cfg['component'], right_cfg['scale'])
            ax2.plot(self.f, y2_data, **series['type']['main'])

        # --- Bloco COMSOL Genérico ---
        if self.cmsl is not None and 'comsol_series_to_plot' in config:
            freq = self.cmsl['frequencies']
            w_cmsl = self.cmsl['angular_frequencies']
            
            for series in config['comsol_series_to_plot']:
                cmsl_key = series.get('base_key')
                cmsl_graph_type = config['comsol_matrix_key']
                
                if cmsl_key not in self.cmsl['scenarios']:
                    print(f"Aviso: Chave COMSOL '{cmsl_key}' não encontrada.")
                    continue

                if cmsl_graph_type in ['gamma_earth']:
                    cmsl_item = self.cmsl['scenarios'][cmsl_key]['earth_return_parameters'][cmsl_graph_type]

                # Obter o estilo de plotagem
                style = {k: v for k, v in series.items() if k != 'base_key'}
                
                # Plotar dados do subplot esquerdo
                y1_data_cmsl = self._calculate_plot_data(cmsl_item, w_cmsl, left_cfg['component'], left_cfg['scale'])
                ax1.scatter(freq, y1_data_cmsl, **style)
                
                # Plotar dados do subplot direito
                y2_data_cmsl = self._calculate_plot_data(cmsl_item, w_cmsl, right_cfg['component'], right_cfg['scale'])
                ax2.scatter(freq, y2_data_cmsl, **style)

        # --- Formatação Genérica dos Eixos ---
        self._format_axis(ax1, self.f, self.xlim, left_cfg)
        self._format_axis(ax2, self.f, self.xlim, right_cfg)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')
            
    def _calculate_plot_data(self, data, w, component_type, scale):
        if component_type == 'real':
            return np.real(data) * scale
        elif component_type == 'imag':
            return np.imag(data) * scale
        elif component_type == 'imag_div_w':
            return np.imag(data) / w * scale
        # ... outras transformações ...
        
    def _format_axis(self, ax, f, xlim, cfg):
        ax.set_xscale(cfg.get('xscale', 'log'))
        ax.set_yscale(cfg.get('yscale', 'linear'))
        ax.set_xlim(xlim)
        ax.legend(fontsize='small')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(cfg['label'])
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

