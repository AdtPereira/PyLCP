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
        Plots one or more graphs based on the configuration.
        The actual plotting method is read from PLOT_CONFIG.
        """
        if not isinstance(graph_key_list, list):
            graph_key_list = [graph_key_list]

        for key in graph_key_list:
            if key not in self.plot_config:
                print(f"Warning: graph key '{key}' not found in PLOT_CONFIG.")
                continue

            cfg = self.plot_config[key]

            # Look up the plotting function name in the config
            plot_function_name = cfg.get('plot_function')

            if not plot_function_name:
                print(f"Error: 'plot_function' key not defined in PLOT_CONFIG['{key}'].")
                continue

            # Get the actual method from the class and call it
            plot_method = getattr(self, plot_function_name, None)

            if plot_method and callable(plot_method):
                plot_method(key) # Calls self._plot_matricial_upper_triangular(key)
            else:
                print(f"Error: method '{plot_function_name}' not found in BasePlotter.")

    def _get_data_from_source(self, series_def):
        """
        Fetches the dataset (matrix, vector) and its corresponding
        frequencies based on the source and path.
        """
        base_obj, freq = None, None
        source_type = series_def['source']
        scenario_key = series_def['scenario_key']

        # 1. Get the base object from the source and the frequencies
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
            print(f"Warning: key '{scenario_key}' not found for source '{source_type}'.")
            return None, None

        # 2. Navigate the dictionary using data_path
        data = base_obj
        try:
            for key in series_def['data_path']:
                data = data[key]
            return data, freq
        except (KeyError, TypeError):
            print(f"Warning: data path {series_def['data_path']} not found for {source_type}['{scenario_key}'].")
            return None, None

    def _plot_matricial_parameters(self, graph_key):
        """
        Generic (refactored) method to plot matricial data from
        multiple sources.
        """
        cfg = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg = cfg['left_plot']
        right_cfg = cfg['right_plot']

        # --- Single Plotting Block ---
        any_plotted = False
        for series_def in cfg.get('data_series', []):
            # 1. Fetch the data using the helper
            matrix, freq = self._get_data_from_source(series_def)

            if matrix is None:
                continue # Helper already emitted the warning

            w = 2 * np.pi * freq

            plot_style_func_ax1 = ax1.plot if series_def.get('plot_style', 'line') == 'line' else ax1.scatter
            plot_style_func_ax2 = ax2.plot if series_def.get('plot_style', 'line') == 'line' else ax2.scatter

            # 2. Iterate over the matrix elements to plot
            for key, value in series_def['series'].items():
                p, q = value['p'], value['q']

                # Copy the style, removing control keys
                style = {k: v for k, v in value.items() if k not in ['p', 'q']}

                # Plot data on the left subplot
                y1_data = self._calculate_plot_data(matrix[:, p, q], w, left_cfg)
                plot_style_func_ax1(freq, y1_data, **style)

                # Plot data on the right subplot
                y2_data = self._calculate_plot_data(matrix[:, p, q], w, right_cfg)
                plot_style_func_ax2(freq, y2_data, **style)
                any_plotted = True

        if not any_plotted:
            plt.close(fig)
            print(f"Warning: no data available for '{graph_key}'. Graph skipped.")
            return

        # --- Generic Axis Formatting ---
        self._format_axis(ax1, self.xlim, left_cfg)
        self._format_axis(ax2, self.xlim, right_cfg)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=f'{graph_key}')

    def _plot_non_matricial_parameters(self, graph_key):
        """
        Generic (refactored) method to plot parameters stored in
        lists, from multiple sources.

        E.g.: Plot the 1st and 3rd items of 'internal_impedance_elements'
        """
        cfg = self.plot_config[graph_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(cfg['suptitle'], fontsize=12, y=0.98)
        left_cfg = cfg['left_plot']
        right_cfg = cfg['right_plot']

        # --- Single Plotting Block ---
        any_plotted = False
        for series_def in cfg.get('data_series', []):
            # 1. Fetch the data (should be a list/array of arrays)
            data_list, freq = self._get_data_from_source(series_def)

            if data_list is None:
                continue # Helper already emitted the warning

            w = 2 * np.pi * freq

            plot_style_func = ax1.plot if series_def.get('plot_style', 'line') == 'line' else ax1.scatter
            plot_style_func_ax2 = ax2.plot if series_def.get('plot_style', 'line') == 'line' else ax2.scatter

            for key, value in series_def['series'].items():
                style = {k: v for k, v in value.items() if k != 'idx'}

                try:
                    if isinstance(data_list, np.ndarray) and data_list.ndim == 2:
                        data_item = data_list[:, key]
                    else:
                        data_item = data_list[key]

                except (IndexError, TypeError, KeyError) as e:
                    print(f"Warning: could not access index {key} for {series_def['source']}['{series_def['scenario_key']}']. Error: {e}")
                    continue

                y1_data = self._calculate_plot_data(data_item, w, left_cfg)
                plot_style_func(freq, y1_data, **style)

                y2_data = self._calculate_plot_data(data_item, w, right_cfg)
                plot_style_func_ax2(freq, y2_data, **style)
                any_plotted = True

        if not any_plotted:
            plt.close(fig)
            print(f"Warning: no data available for '{graph_key}'. Graph skipped.")
            return

        # --- Generic Axis Formatting ---
        self._format_axis(ax1, self.xlim, left_cfg)
        self._format_axis(ax2, self.xlim, right_cfg)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=f'{graph_key}')

    def _calculate_plot_data(self, data, w, cfg):
        component_type = cfg.get('component', 'real')
        scale = cfg.get('scale', 1.0)
        if component_type == 'real':
            return np.real(data) * scale
        elif component_type == 'imag':
            return np.imag(data) * scale
        elif component_type == 'imag_div_w':
            return np.imag(data) / w * scale
        elif component_type == 'abs':
            return np.abs(data) * scale
        elif component_type == 'angle_deg':
            return np.angle(data, deg=True) * scale
        # ... other transformations ...
        
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
        if 'title' in cfg:
            ax.set_title(cfg['title'])
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
