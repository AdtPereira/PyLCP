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
        
        self.mom_so = self.pul_data.get('mom_so', None)
        self.cmsl = self.pul_data.get('comsol', None)
        self.f = pul_data['frequencies']
        self.w = 2 * np.pi * self.f        

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.script_path.stem, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)
        
        self.xlim = tuple(pul_data['frequencies'][[0, -1]])
        self.figsize = (12, 5)

    def _plot_matricial_input(self, graph_key):
        """
        Método genérico para plotar uma grandeza complexa em dois subplots 
        (ex: R/L ou G/C).
        """
        cfg = self.plot_config[graph_key]
        p, q = cfg['p'], cfg['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg = cfg['left_plot']
        right_cfg = cfg['right_plot']
        
        # --- Bloco Analítico Genérico ---
        if self.pul_data['scenarios']:
            for series in cfg['series_to_plot']:
                matrix = self.pul_data['scenarios'][series['key']][cfg['path'][0]][cfg['path'][1]]
                
                # Plotar dados (Ex: Resistência)
                y1_data = self._calculate_plot_data(matrix[:, p, q], self.w, left_cfg['component'], left_cfg['scale'])
                ax1.plot(self.f, y1_data, **series['type']['main'])
                
                # Plotar dados (Ex: Indutância)
                y2_data = self._calculate_plot_data(matrix[:, p, q], self.w, right_cfg['component'], right_cfg['scale'])
                ax2.plot(self.f, y2_data, **series['type']['main'])

        # --- Bloco COMSOL Genérico ---
        if self.cmsl is not None and 'comsol_series_to_plot' in cfg:
            freq = self.cmsl['frequencies']
            w_cmsl = self.cmsl['angular_frequencies']
            
            for series in cfg['comsol_series_to_plot']:
                cmsl_key = series.get('key')
                cmsl_graph_type = cfg['comsol_matrix_key']
                
                if cmsl_key not in self.cmsl['scenarios']:
                    print(f"Aviso: Chave COMSOL '{cmsl_key}' não encontrada.")
                    continue

                if cmsl_graph_type in ['impedance_matrix', 'admittance_matrix']:
                    cmsl_matrix = self.cmsl['scenarios'][cmsl_key]['earth_return_parameters'][cmsl_graph_type]

                elif cmsl_graph_type in ['series_impedance_matrix', 'shunt_admittance_matrix']:
                    cmsl_matrix = self.cmsl['scenarios'][cmsl_key]['quasi_tem_matrices'][cmsl_graph_type]

                # Obter o estilo de plotagem
                style = {k: v for k, v in series.items() if k != 'key'}
                
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

    def _plot_matricial_upper_triangular(self, graph_key):
        """
        Método genérico para plotar uma grandeza complexa em dois subplots 
        (ex: R/L ou G/C).
        """
        cfg = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg = cfg['left_plot']
        right_cfg = cfg['right_plot']
        
        # --- Bloco Analítico Genérico ---
        if self.pul_data['scenarios']:
            for series in cfg['series_to_plot']:
                matrix = self.pul_data['scenarios'][series['key']][cfg['path'][0]][cfg['path'][1]]

                for value in series['type'].values():
                    p, q = value['p'], value['q']
                    
                    # Plotar dados do subplot esquerdo
                    y1_data = self._calculate_plot_data(matrix[:, p, q], self.w, left_cfg['component'], left_cfg['scale'])                
                    ax1.plot(self.f, y1_data, color=value['color'], linestyle=value['linestyle'], linewidth=value['linewidth'], label=value['label'])
                    
                    # Plotar dados do subplot direito
                    y2_data = self._calculate_plot_data(matrix[:, p, q], self.w, right_cfg['component'], right_cfg['scale'])                
                    ax2.plot(self.f, y2_data, color=value['color'], linestyle=value['linestyle'], linewidth=value['linewidth'], label=value['label'])
            
        # --- Bloco COMSOL Genérico ---
        if self.cmsl is not None and 'comsol_series_to_plot' in cfg:
            freq = self.cmsl['frequencies']
            w = self.cmsl['angular_frequencies']
            
            for series in cfg['comsol_series_to_plot']:
                key = series.get('key')
                graph_type = cfg['comsol_matrix_key']
                
                if key not in self.cmsl['scenarios']:
                    print(f"Aviso: Chave COMSOL '{key}' não encontrada.")
                    continue

                if graph_type in ['internal_impedance_matrix']:
                    cmsl_matrix = self.cmsl['scenarios'][key][graph_type]
                
                else:
                    continue

                for key, value in series['type'].items():
                    p, q = value['p'], value['q']
                    style = {k: v for k, v in value.items() if k not in ['p', 'q']}

                    y1_data = self._calculate_plot_data(cmsl_matrix[:, p, q], w, left_cfg['component'], left_cfg['scale'])
                    ax1.scatter(freq, y1_data, **style)

                    y2_data = self._calculate_plot_data(cmsl_matrix[:, p, q], w, right_cfg['component'], right_cfg['scale'])
                    ax2.scatter(freq, y2_data, **style)

        # --- Bloco MoM-SO Genérico ---
        if self.mom_so is not None and 'mom_so_series_to_plot' in cfg:
            freq = self.mom_so['frequencies']
            w = 2 * np.pi * freq
            
            for series in cfg['mom_so_series_to_plot']:
                key = series.get('key')
                graph_type = cfg['mom_so_matrix_key']

                if key not in self.mom_so['scenarios']:
                    print(f"Aviso: Chave MoM-SO '{key}' não encontrada.")
                    continue

                if graph_type in ['partial_internal_impedance']:
                    matrix = self.mom_so['scenarios'][key][graph_type]

                else:
                    continue

                for key, value in series['type'].items():
                    p, q = value['p'], value['q']
                    style = {k: v for k, v in value.items() if k not in ['p', 'q']}
                    
                    y1_data = self._calculate_plot_data(matrix[:, p, q], w, left_cfg['component'], left_cfg['scale'])
                    ax1.scatter(freq, y1_data, **style)

                    y2_data = self._calculate_plot_data(matrix[:, p, q], w, right_cfg['component'], right_cfg['scale'])
                    ax2.scatter(freq, y2_data, **style)

        # --- Formatação Genérica dos Eixos ---
        self._format_axis(ax1, self.f, self.xlim, left_cfg)
        self._format_axis(ax2, self.f, self.xlim, right_cfg)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')

    def _plot_non_matricial_list_parameter(self, graph_key):
        """
        Método genérico para plotar uma grandeza complexa em dois subplots 
        (ex: R/L ou G/C).
        """
        cfg = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg = cfg['left_plot']
        right_cfg = cfg['right_plot']
        
        # --- Bloco Analítico Genérico ---
        if self.pul_data['scenarios']:
            for series in cfg['series_to_plot']:
                matrix = self.pul_data['scenarios'][series['key']][cfg['path'][0]][cfg['path'][1]]

                for key, value in series['type'].items():
                    # Plotar dados do subplot esquerdo
                    y1_data = self._calculate_plot_data(matrix[key], self.w, left_cfg['component'], left_cfg['scale'])                
                    ax1.plot(self.f, y1_data, **value)
                    
                    # Plotar dados do subplot direito
                    y2_data = self._calculate_plot_data(matrix[key], self.w, right_cfg['component'], right_cfg['scale'])
                    ax2.plot(self.f, y2_data, **value)
            
        # --- Bloco COMSOL Genérico ---
        if self.cmsl is not None and 'comsol_series_to_plot' in cfg:
            freq = self.cmsl['frequencies']
            w = self.cmsl['angular_frequencies']
            
            for series in cfg['comsol_series_to_plot']:
                key = series.get('key')
                graph_type = cfg['comsol_matrix_key']
                
                if key not in self.cmsl['scenarios']:
                    print(f"Aviso: Chave COMSOL '{key}' não encontrada.")
                    continue

                if graph_type in ['impedance_matrix', 'admittance_matrix']:
                    matrix = self.cmsl['scenarios'][key]['earth_return_parameters'][graph_type]

                elif graph_type in ['series_impedance_matrix', 'shunt_admittance_matrix']:
                    matrix = self.cmsl['scenarios'][key]['quasi_tem_matrices'][graph_type]

                elif graph_type in ['coaxial_cable_impedance', 'internal_impedance_elements']:
                    matrix = self.cmsl['scenarios'][key][graph_type]

                else:
                    continue

                for key, value in series['type'].items():
                    y1_data_cmsl = self._calculate_plot_data(matrix[key], w, left_cfg['component'], left_cfg['scale'])
                    ax1.scatter(freq, y1_data_cmsl, **value)

                    y2_data_cmsl = self._calculate_plot_data(matrix[key], w, right_cfg['component'], right_cfg['scale'])
                    ax2.scatter(freq, y2_data_cmsl, **value)

        # --- Bloco MoM-SO Genérico ---
        if self.mom_so is not None and 'mom_so_series_to_plot' in cfg:
            freq = self.mom_so['frequencies']
            w = 2 * np.pi * freq
            
            for series in cfg['mom_so_series_to_plot']:
                key = series.get('key')
                graph_type = cfg['mom_so_matrix_key']

                if key not in self.mom_so['scenarios']:
                    print(f"Aviso: Chave MoM-SO '{key}' não encontrada.")
                    continue

                if graph_type in ['coaxial_cable_impedance']:
                    matrix = self.mom_so['scenarios'][key][graph_type]

                else:
                    continue

                for key, value in series['type'].items():
                    y1_data = self._calculate_plot_data(matrix[:, 0, 0], w, left_cfg['component'], left_cfg['scale'])
                    ax1.scatter(freq, y1_data, **value)

                    y2_data = self._calculate_plot_data(matrix[:, 0, 0], w, right_cfg['component'], right_cfg['scale'])
                    ax2.scatter(freq, y2_data, **value)
                    
        # --- Formatação Genérica dos Eixos ---
        self._format_axis(ax1, self.f, self.xlim, left_cfg)
        self._format_axis(ax2, self.f, self.xlim, right_cfg)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure_multiformat(fig, self.results_dir, base_filename=f'{graph_key}')

    def _plot_non_matricial_parameter(self, graph_key):
        """
        Método genérico para plotar uma grandeza complexa em dois subplots 
        (ex: R/L ou G/C).
        """
        cfg = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg = cfg['left_plot']
        right_cfg = cfg['right_plot']        
        
        if self.pul_data['scenarios']:
            for series in cfg['series_to_plot']:
                item = self.pul_data['scenarios'][series['key']][cfg['path'][0]][cfg['path'][1]]
                
                y1_data = self._calculate_plot_data(item, self.w, left_cfg['component'], left_cfg['scale'])
                ax1.plot(self.f, y1_data, **series['type']['main'])
                
                y2_data = self._calculate_plot_data(item, self.w, right_cfg['component'], right_cfg['scale'])
                ax2.plot(self.f, y2_data, **series['type']['main'])

        # --- Bloco COMSOL Genérico ---
        if self.cmsl is not None and 'comsol_series_to_plot' in cfg:
            freq = self.cmsl['frequencies']
            w_cmsl = self.cmsl['angular_frequencies']
            
            for series in cfg['comsol_series_to_plot']:
                cmsl_key = series.get('key')
                cmsl_graph_type = cfg['comsol_matrix_key']
                
                if cmsl_key not in self.cmsl['scenarios']:
                    print(f"Aviso: Chave COMSOL '{cmsl_key}' não encontrada.")
                    continue

                if cmsl_graph_type in ['gamma_earth']:
                    cmsl_item = self.cmsl['scenarios'][cmsl_key]['earth_return_parameters'][cmsl_graph_type]

                # Obter o estilo de plotagem
                style = {k: v for k, v in series.items() if k != 'key'}
                
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
        if 'y_lim' in cfg:
            ax.set_ylim(cfg['y_lim'])
        if cfg.get('legend', False):
           ax.legend(fontsize='small')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(cfg['label'])
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
