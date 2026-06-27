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
        self.analytical = self.pul_data.get('analytical', None)

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.script_path.stem, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)
        
        _freq_src = (pul_data.get('analytical')
                     or pul_data.get('mom_so')
                     or pul_data.get('comsol')
                     or pul_data)
        _freqs = _freq_src.get('frequencies', np.array([1.0, 1e6]))
        self.xlim = tuple(_freqs[[0, -1]])
        self.figsize = (12, 5)

    def plot_graph(self, graph_key_list):
        """
        Plota um ou mais gráficos com base na configuração.
        O método de plotagem real é lido do PLOT_CONFIG.
        """
        if not isinstance(graph_key_list, list):
            graph_key_list = [graph_key_list]
            
        for key in graph_key_list:
            if key not in self.plot_config:
                print(f"Aviso: Chave de gráfico '{key}' não encontrada no PLOT_CONFIG.")
                continue
                
            cfg = self.plot_config[key]
            
            # Busca o nome da função de plotagem na config
            plot_function_name = cfg.get('plot_function') 
            
            if not plot_function_name:
                print(f"Erro: Chave 'plot_function' não definida em PLOT_CONFIG['{key}'].")
                continue
                
            # Obtém o método real da classe e o chama
            plot_method = getattr(self, plot_function_name, None)
            
            if plot_method and callable(plot_method):
                plot_method(key) # Chama self._plot_matricial_upper_triangular(key)
            else:
                print(f"Erro: Método '{plot_function_name}' não encontrado em BasePlotter.")

    def _get_data_from_source(self, series_def):
        """
        Busca o conjunto de dados (matriz, vetor) e suas frequências
        correspondentes com base na fonte e no caminho.
        """
        base_obj, freq = None, None
        source_type = series_def['source']
        scenario_key = series_def['scenario_key']
        
        # 1. Obter o objeto base da fonte e as frequências
        if source_type == 'analytical':
            if self.analytical and self.analytical['scenarios'].get(scenario_key):
                base_obj = self.analytical['scenarios'][scenario_key]
                freq = self.analytical['frequencies']
        elif source_type == 'comsol':
            if self.cmsl and self.cmsl['scenarios'].get(scenario_key):
                base_obj = self.cmsl['scenarios'][scenario_key]
                freq = self.cmsl['frequencies']
        elif source_type == 'mom_so':
            if self.mom_so and self.mom_so['scenarios'].get(scenario_key):
                base_obj = self.mom_so['scenarios'][scenario_key]
                freq = self.mom_so['frequencies']
        
        if base_obj is None:
            print(f"Aviso: Chave '{scenario_key}' não encontrada para a fonte '{source_type}'.")
            return None, None

        # 2. Navegar no dicionário usando o data_path
        data = base_obj
        try:
            for key in series_def['data_path']:
                data = data[key]
            return data, freq
        except (KeyError, TypeError):
            print(f"Aviso: Caminho de dados {series_def['data_path']} não encontrado para {source_type}['{scenario_key}'].")
            return None, None

    def _plot_matricial_parameters(self, graph_key):
        """
        Método genérico (refatorado) para plotar dados matriciais 
        de múltiplas fontes.
        """
        cfg = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg = cfg['left_plot']
        right_cfg = cfg['right_plot']
        
        # --- Bloco Único de Plotagem ---
        for series_def in cfg.get('data_series', []):            
            # 1. Buscar os dados usando o helper
            matrix, freq = self._get_data_from_source(series_def)
            w = 2 * np.pi * freq
            
            if matrix is None:
                continue # Helper já emitiu o aviso

            plot_style_func_ax1 = ax1.plot if series_def.get('plot_style', 'line') == 'line' else ax1.scatter
            plot_style_func_ax2 = ax2.plot if series_def.get('plot_style', 'line') == 'line' else ax2.scatter

            # 2. Iterar sobre os elementos da matriz a plotar
            for key, value in series_def['series'].items():
                p, q = value['p'], value['q']
                
                # Copia o estilo, removendo chaves de controle
                style = {k: v for k, v in value.items() if k not in ['p', 'q']}

                # Plotar dados do subplot esquerdo
                y1_data = self._calculate_plot_data(matrix[:, p, q], w, left_cfg)
                plot_style_func_ax1(freq, y1_data, **style)
                
                # Plotar dados do subplot direito
                y2_data = self._calculate_plot_data(matrix[:, p, q], w, right_cfg)
                plot_style_func_ax2(freq, y2_data, **style)

        # --- Formatação Genérica dos Eixos ---
        self._format_axis(ax1, self.xlim, left_cfg)
        self._format_axis(ax2, self.xlim, right_cfg)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=f'{graph_key}')

    def _plot_non_matricial_parameters(self, graph_key):
        """
        Método genérico (refatorado) para plotar parâmetros 
        armazenados em listas, de múltiplas fontes.
        
        Ex: Plotar o 1º e 3º item de 'internal_impedance_elements'
        """
        cfg = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg = cfg['left_plot']
        right_cfg = cfg['right_plot']
        
        # --- Bloco Único de Plotagem ---
        for series_def in cfg.get('data_series', []):            
            # 1. Buscar os dados (deve ser uma lista/array de arrays)
            data_list, freq = self._get_data_from_source(series_def)
            
            if data_list is None:
                continue # Helper já emitiu o aviso

            w = 2 * np.pi * freq
            
            plot_style_func = ax1.plot if series_def.get('plot_style', 'line') == 'line' else ax1.scatter
            plot_style_func_ax2 = ax2.plot if series_def.get('plot_style', 'line') == 'line' else ax2.scatter

            # 2. Iterar sobre os elementos da lista a plotar
            # 'series' é um dict onde 'value' contém o 'idx' (índice da lista)
            # e o resto é o estilo de plotagem.
            for key, value in series_def['series'].items():
                
                # if 'idx' not in value:
                #     print(f"Aviso: 'idx' não encontrado em series['{key}'] para {graph_key}.")
                #     continue
                    
                # idx = value['idx']
                
                # Copia o estilo, removendo chaves de controle
                style = {k: v for k, v in value.items() if k != 'idx'}

                try:
                    # Acessa o item específico da lista/array.
                    # Se data_list for (n_freq, n_elements), pegamos [:, idx]
                    # Se data_list for uma lista de (n_freq,), pegamos [idx]
                    if isinstance(data_list, np.ndarray) and data_list.ndim == 2:
                        data_item = data_list[:, key]
                    else:
                        data_item = data_list[key] # Assume lista de arrays
                
                except (IndexError, TypeError, KeyError) as e:
                    print(f"Aviso: Não foi possível acessar o índice {key} para {series_def['source']}['{series_def['scenario_key']}']. Erro: {e}")
                    continue

                # Plotar dados do subplot esquerdo
                y1_data = self._calculate_plot_data(data_item, w, left_cfg)
                plot_style_func(freq, y1_data, **style)
                
                # Plotar dados do subplot direito
                y2_data = self._calculate_plot_data(data_item, w, right_cfg)
                plot_style_func_ax2(freq, y2_data, **style)

        # --- Formatação Genérica dos Eixos ---
        # Note a remoção de 'self.f' da chamada, para alinhar com o
        # método '_plot_matricial_upper_triangular' já revisado.
        self._format_axis(ax1, self.xlim, left_cfg)
        self._format_axis(ax2, self.xlim, right_cfg)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=f'{graph_key}')

    # def _plot_non_matricial_parameter(self, graph_key):
    #     """
    #     Método genérico para plotar uma grandeza complexa em dois subplots 
    #     (ex: R/L ou G/C).
    #     """
    #     cfg = self.plot_config[graph_key]
    #     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
    #     fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
    #     left_cfg = cfg['left_plot']
    #     right_cfg = cfg['right_plot']        
        
    #     if self.pul_data['scenarios']:
    #         for series in cfg['series_to_plot']:
    #             item = self.pul_data['scenarios'][series['key']][cfg['path'][0]][cfg['path'][1]]
                
    #             y1_data = self._calculate_plot_data(item, self.w, left_cfg['component'], left_cfg['scale'])
    #             ax1.plot(self.f, y1_data, **series['type']['main'])
                
    #             y2_data = self._calculate_plot_data(item, self.w, right_cfg['component'], right_cfg['scale'])
    #             ax2.plot(self.f, y2_data, **series['type']['main'])

    #     # --- Bloco COMSOL Genérico ---
    #     if self.cmsl is not None and 'comsol_series_to_plot' in cfg:
    #         freq = self.cmsl['frequencies']
    #         w_cmsl = self.cmsl['angular_frequencies']
            
    #         for series in cfg['comsol_series_to_plot']:
    #             cmsl_key = series.get('key')
    #             cmsl_graph_type = cfg['comsol_matrix_key']
                
    #             if cmsl_key not in self.cmsl['scenarios']:
    #                 print(f"Aviso: Chave COMSOL '{cmsl_key}' não encontrada.")
    #                 continue

    #             if cmsl_graph_type in ['gamma_earth']:
    #                 cmsl_item = self.cmsl['scenarios'][cmsl_key]['earth_return_parameters'][cmsl_graph_type]

    #             # Obter o estilo de plotagem
    #             style = {k: v for k, v in series.items() if k != 'key'}
                
    #             # Plotar dados do subplot esquerdo
    #             y1_data_cmsl = self._calculate_plot_data(cmsl_item, w_cmsl, left_cfg['component'], left_cfg['scale'])
    #             ax1.scatter(freq, y1_data_cmsl, **style)
                
    #             # Plotar dados do subplot direito
    #             y2_data_cmsl = self._calculate_plot_data(cmsl_item, w_cmsl, right_cfg['component'], right_cfg['scale'])
    #             ax2.scatter(freq, y2_data_cmsl, **style)

    #     # --- Formatação Genérica dos Eixos ---
    #     self._format_axis(ax1, self.f, self.xlim, left_cfg)
    #     self._format_axis(ax2, self.f, self.xlim, right_cfg)
        
    #     plt.tight_layout(rect=[0, 0, 1, 0.96])
    #     if self.autoSave:
    #         save_figure(fig, self.results_dir, base_filename=f'{graph_key}')
            
    def _calculate_plot_data(self, data, w, cfg):
        component_type = cfg.get('component', 'real')
        scale = cfg.get('scale', 1.0)
        if component_type == 'real':
            return np.real(data) * scale
        elif component_type == 'imag':
            return np.imag(data) * scale
        elif component_type == 'imag_div_w':
            return np.imag(data) / w * scale
        # ... outras transformações ...
        
    def _format_axis(self, ax, xlim, cfg):
        ax.set_xscale(cfg.get('xscale', 'log'))
        ax.set_yscale(cfg.get('yscale', 'linear'))
        ax.set_xlim(xlim)
        if 'y_lim' in cfg:
            ax.set_ylim(cfg['y_lim'])
        if cfg.get('legend', True):
           ax.legend(fontsize='small')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(cfg['label'])
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
