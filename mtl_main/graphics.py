import os
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
from mtl_main.source import MulticonductorTransmissionLine
from utils.case_utils import *

class BaseMTLRepresentation:
    """ 
    Base class for creating graphical representations of an MTL model.
    It contains common functionalities and attributes shared by specialized
    representation classes.
    """
    def __init__(self, file_path: str, model: MulticonductorTransmissionLine, autoSave: bool = True, units: str = 'meter'):        
        
        self.script_path = Path(file_path)
        self.model = model
        self.autoSave = autoSave
        unit_info = UNITS_DATA.get(units, UNITS_DATA['meter'])
        active_conductors = list(self.model.mtl.values())
        
        self.num_sc_cables = len(set(cond['center_point'] for cond in active_conductors))
        self.unit_factor = unit_info['scale']
        self.label_unit = unit_info['label']
        self.figsize = (12, 5)
        self.color_map = {
            'core': 'gray',
            'XLPE': 'lightblue',
            'sheath': 'darkgreen',
            'HDPE': 'lightgreen',
            'core insulation': 'lightblue',
            'outer insulation': 'lightyellow',
            'enclosure': 'darkgray',
            'soil': 'tan',
            'return': 'brown',
            'default_conductor': 'lightgray',
            'default_insulation': 'white',
            'pipe': 'lightgray',
            'air': 'ghostwhite',
            'vacuum': 'ghostwhite',
        }

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.script_path.stem, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

    def _calculate_schematic_parameters(self) -> dict:
        """
        Calculates key parameters required for plotting the schematic.
        
        Returns:
            A dictionary containing parameters like max_radius, h_factor,
            title, and the primary core conductor data.
        """
        active_conductors = [c for c in self.model.mtl.values() if c.get('line_type') == 'active']
        core_conductor = active_conductors[0] if active_conductors else None
        
        # Calculate the average vertical position and find the deepest conductor
        unique_centers = list(set(tuple(c['center_point']) for c in active_conductors))
        y_avg, y_min = 0, 0
        deepest_conductor = core_conductor
        if unique_centers:
            y_coords = [center[1] for center in unique_centers]
            y_avg = sum(y_coords) / len(y_coords)
            y_min = min(y_coords)
            deepest_conductor = next((c for c in active_conductors if c['center_point'][1] == y_min), core_conductor)

        max_radius = 0
        if core_conductor:
            for conductor in self.model.mtl.values():
                radius = conductor['radius'][1]
                if 'insulation' in conductor and conductor['insulation'] is not None:
                    radius += conductor['insulation']['thickness']
                if 'enclosure' in conductor and conductor['enclosure'] is not None:
                    radius = max(radius, conductor['enclosure']['outer_radius'])
                if radius > max_radius:
                    max_radius = radius
        
        if self.model.mtl_type == 'overhead':
            title = 'Overhead Transmission Line'
            h_factor = 8
        elif self.model.mtl_type == 'scc':
            title = 'Buried Single-Core Cable'
            h_factor = -2.5
        elif self.model.mtl_type == 'hdpe':
            title = 'Single-Core Cable Buried in HDPE-Air Gap Enclosure'
            h_factor = -2
        elif self.model.mtl_type == 'shared-hdpe':
            title = 'Single-Core Cable Buried in HDPE-Air Gap Enclosure'
            h_factor = -2
        elif self.model.mtl_type == 'pipe':
            title = 'Pipe-Type Cable'
            h_factor = 1
        elif self.model.mtl_type == 'coaxial':
            title = 'Isolated Coaxial Cable'
            h_factor = 1
        else: 
            title = 'Multicondutor Transmission Line'
            h_factor = 1

        return {
            'core_conductor': core_conductor,
            'deepest_conductor': deepest_conductor,
            'max_radius': max_radius,
            'h_factor': h_factor,
            'title': title,
            'y_avg': y_avg,
            'y_min': y_min,
        }
    
    def _plot_conductor_graphic(self, ax, conductor_data, parameters, used_labels):
        """
        Plots a single cylindrical layer of a cable (as a Circle or Wedge).
        Handles the logic for avoiding duplicate legend entries.
        """
        conductor_label = conductor_data.get('conductor_name', 'Conductor')
        conductor_center = np.array(conductor_data['center_point']) if conductor_data['center_point'] else np.array([0, 0])
        insulation_data = conductor_data.get('insulation')
        insulation_center = np.array(insulation_data['center_point']) if insulation_data else np.copy(conductor_center)
        label_to_plot = conductor_label if conductor_label not in used_labels else None

        # Adjust y-position for ground-return systems while preserving relative heights
        if parameters['h_factor'] != 1:
            y_real_conductor = conductor_data['center_point'][1]
            y_avg = parameters['y_avg']
            schematic_y_avg = parameters['h_factor'] * parameters['max_radius']
            
            # Calculate the new schematic y-position for the conductor
            schematic_y_conductor = schematic_y_avg + (y_real_conductor - y_avg)
            
            # Handle insulation center, preserving any original offset
            y_offset_insulation = insulation_center[1] - conductor_center[1]
            
            conductor_center[1] = schematic_y_conductor
            insulation_center[1] = schematic_y_conductor + y_offset_insulation
        
        if label_to_plot:
            used_labels.add(label_to_plot)

        if 'insulation' in conductor_data and conductor_data['insulation'] is not None:
            insulation = conductor_data['insulation']
            insulation_type = insulation.get('type', 'Insulation')
            patch = Wedge(
                center=insulation_center * self.unit_factor,
                r=(conductor_data['radius'][1] + insulation['thickness']) * self.unit_factor,
                theta1=0, theta2=360,
                width=insulation['thickness'] * self.unit_factor,
                edgecolor='black',
                facecolor=self.color_map.get(insulation_type, self.color_map['default_insulation']),
                linestyle='solid',
                label=insulation_type if insulation_type not in used_labels else None,
                zorder=4)
            ax.add_patch(patch)
            
        # Wedge do material condutor
        patch = Wedge(
            center=conductor_center * self.unit_factor,
            r=conductor_data['radius'][1] * self.unit_factor,
            theta1=0, theta2=360,
            width=(conductor_data['radius'][1] - conductor_data['radius'][0]   ) * self.unit_factor,
            edgecolor='black',
            facecolor=self.color_map.get(conductor_label, self.color_map['default_conductor']),
            linestyle='solid',
            label=label_to_plot,
            zorder=5)
        ax.add_patch(patch)        

        # Preenche o interior oco do tipo 'pipe' or 'core'.
        if conductor_label.lower() in ['pipe', 'core'] and conductor_data['radius'][0] > 0:
            fill_patch = Circle(
                xy=conductor_center * self.unit_factor,
                radius=conductor_data['radius'][0] * self.unit_factor,
                fill=True,
                edgecolor='black',
                facecolor='ghostwhite',
                linestyle='solid',
                label='air' if 'air' not in used_labels else None,
                zorder=2)
            ax.add_patch(fill_patch)

class IsolatedMTLRepresentation(BaseMTLRepresentation):
    """
    Creates graphical representations for isolated MTL systems,
    such as single wires, coaxial cables, and pipe-type cables where
    the ground effect is not the primary focus of the schematic.
    """
    def __init__(self, file_path: str, model: MulticonductorTransmissionLine, autoSave=True, units='meter'):
        super().__init__(file_path, model, autoSave, units)

    def _finalize_plot(self, ax, title):
        """Applies final settings for an isolated system plot."""
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'x ({self.label_unit})')
        ax.set_ylabel(f'')
        ax.set_yticks([])
        if title:
            ax.set_title(title)
        ax.grid(True, linestyle='--', linewidth=0.5, zorder=0)

        ax.relim()
        ax.autoscale_view()
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        x_margin = (xlim[1] - xlim[0]) * 0.5
        y_margin = (ylim[1] - ylim[0]) * 0.2
        ax.set_xlim(xlim[0] - x_margin, xlim[1] + x_margin)
        ax.set_ylim(ylim[0] - y_margin, ylim[1] + y_margin)
        
        final_xlim, final_ylim = ax.get_xlim(), ax.get_ylim()
        text_x = final_xlim[0] + (final_xlim[1] - final_xlim[0]) * 0.05
        y_pos = final_ylim[0] + (final_ylim[1] - final_ylim[0]) * 0.5
        ax.text(text_x, y_pos, 'Air ($\\varepsilon_0$, $\\mu_0$)', 
                va='center', fontsize=10, style='italic', zorder=4)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            unique_labels = dict(zip(labels, handles))
            legend = ax.legend(unique_labels.values(), unique_labels.keys(), loc='best')
            legend.set_zorder(5)

    def system_schematic(self, base_filename='system_schematic') -> None:
        """
        Generates a schematic plot for isolated systems (e.g., pipe-type cables).
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        used_labels = set()
        parameters = self._calculate_schematic_parameters()

        for conductor_data in reversed(self.model.mtl.values()):                 
            self._plot_conductor_graphic(ax, conductor_data, parameters, used_labels)

        self._finalize_plot(ax, parameters['title'])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename)

class GroundReturnMTLRepresentation(BaseMTLRepresentation):
    """
    Creates graphical representations for ground-return MTL systems,
    such as buried cables or overhead lines, where the ground plane
    is an essential part of the schematic.
    """
    def __init__(self, file_path: str, model: MulticonductorTransmissionLine, autoSave: bool = True, units: str = 'meter'):
        super().__init__(file_path, model, autoSave, units)

    def system_schematic(self, base_filename='system_schematic') -> None:
        """
        Generates a schematic plot for ground-return systems.
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        used_labels = set()
        parameters = self._calculate_schematic_parameters()
        
        for conductor_data in reversed(self.model.mtl.values()):
            if 'enclosure' in conductor_data and conductor_data['enclosure'] is not None:
                self._plot_enclosure_graphic(ax, conductor_data, parameters, used_labels)
            self._plot_conductor_graphic(ax, conductor_data, parameters, used_labels)

        self._schematic_annotations(ax, parameters)        
        x_margin_scale = 0.2 if self.num_sc_cables > 2 else 1.0
        self._finalize_plot(ax, parameters['title'], x_margin_scale)

        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename)

    def _plot_enclosure_graphic(self, ax, conductor_data, parameters, used_labels):
        """
        Plots a single cylindrical layer of a cable (as a Circle or Wedge).
        Handles the logic for avoiding duplicate legend entries.
        """
        enclosure = conductor_data['enclosure']
        enclosure_center = np.array(enclosure['center_point']) if 'center_point' in enclosure else np.array([0, 0])
        enclosure_type = enclosure.get('type', 'HDPE')

        # --- INÍCIO DA CORREÇÃO ---
        # Ajusta a posição y para sistemas com retorno pelo solo (ground-return)
        # Esta lógica DEVE ser idêntica à de _plot_conductor_graphic
        if parameters['h_factor'] != 1:
            # Obter os parâmetros de escala
            y_real_enclosure = enclosure['center_point'][1]
            y_avg = parameters['y_avg']
            schematic_y_avg = parameters['h_factor'] * parameters['max_radius']
            
            # Calcular a nova posição y esquemática para o duto
            schematic_y_enclosure = schematic_y_avg + (y_real_enclosure - y_avg)
            
            # Aplicar a nova posição
            enclosure_center[1] = schematic_y_enclosure
        # --- FIM DA CORREÇÃO ---
        
        patch = Wedge(
            center=enclosure_center * self.unit_factor,
            r=enclosure['outer_radius'] * self.unit_factor,
            theta1=0, theta2=360,
            width=(enclosure['outer_radius'] - enclosure['inner_radius']) * self.unit_factor,
            edgecolor='black',
            facecolor=self.color_map.get(enclosure_type, self.color_map['default_conductor']),
            linestyle='solid',
            label=enclosure_type if enclosure_type not in used_labels else None,
            zorder=3)
        ax.add_patch(patch)

        fill_patch = Circle(
            xy=enclosure_center * self.unit_factor,
            radius=enclosure['inner_radius'] * self.unit_factor,
            fill=True,
            edgecolor='black',
            facecolor='ghostwhite',
            label='air' if 'air' not in used_labels else None,
            zorder=2)
        ax.add_patch(fill_patch)

    def _schematic_annotations(self, ax, params):
        """Draws annotations like the ground level and depth/height line."""
        deepest_cond = params['deepest_conductor']
        deepest_cond_x = deepest_cond['center_point'][0]
        
        # Calculate the schematic y-position of the deepest conductor
        y_real_deepest = params['y_min']
        y_avg = params['y_avg']
        schematic_y_avg = params['h_factor'] * params['max_radius']
        schematic_y_deepest = schematic_y_avg + (y_real_deepest - y_avg)
        
        schematic_y_scaled = schematic_y_deepest * self.unit_factor
        
        dim_x_arrow = (deepest_cond_x + params['max_radius'] * 3.5) * self.unit_factor
        ax.annotate('', xy=(dim_x_arrow, 0), xycoords='data', 
                    xytext=(dim_x_arrow, schematic_y_scaled), textcoords='data',
                    arrowprops=dict(arrowstyle='<->', color='black', lw=1, zorder=3))
        
        real_h = abs(y_real_deepest)
        label_text = f'h = {real_h:.2f} m'
        
        ax.text(dim_x_arrow, schematic_y_scaled / 2, label_text, ha='center', va='center', fontsize=9, 
                zorder=3, bbox=dict(boxstyle='square,pad=0.3', fc='white', ec='none', alpha=0.8))

        ax.axhline(y=0, color='darkgreen', linestyle=':', linewidth=1.5, zorder=2)

    def _finalize_plot(self, ax, title, x_margin_scale=2.0):
        """Applies final settings for a ground-return system plot."""
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'x ({self.label_unit})')
        ax.set_ylabel(f'')
        ax.set_yticks([])
        ax.set_title(title)
        ax.grid(True, linestyle='--', linewidth=0.5, zorder=0)

        ax.relim()
        ax.autoscale_view()
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        x_margin = (xlim[1] - xlim[0]) * x_margin_scale
        y_margin = (ylim[1] - ylim[0]) * 0.2
        ax.set_xlim(xlim[0] - x_margin, xlim[1] + x_margin)
        ax.set_ylim(ylim[0] - y_margin, ylim[1] + y_margin)
        
        final_xlim, final_ylim = ax.get_xlim(), ax.get_ylim()
        ax.fill_between(final_xlim, final_ylim[0], 0, color='saddlebrown', alpha=0.2, zorder=1)
        
        text_x = final_xlim[0] + (final_xlim[1] - final_xlim[0]) * 0.05
        y_offset = (final_ylim[1] - final_ylim[0]) * 0.03
        ax.text(text_x, y_offset, 'Air ($\\varepsilon_0$, $\\mu_0$)', 
                va='bottom', fontsize=10, style='italic', zorder=4)
        ax.text(text_x, -y_offset, 'Ground ($\\varepsilon_1$, $\\mu_1$, $\\sigma_1$)', 
                va='top', fontsize=10, style='italic', zorder=4)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            unique_labels = dict(zip(labels, handles))
            legend = ax.legend(unique_labels.values(), unique_labels.keys(), loc='best')
            legend.set_zorder(5)

