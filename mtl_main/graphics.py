import os
import numpy as np
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
    def __init__(self, model: MulticonductorTransmissionLine, case_name, units='meter'):        
        self.model = model
        self.case_name = case_name
        unit_info = UNITS_DATA.get(units, UNITS_DATA['meter'])
        active_conductors = list(self.model.mtl.values())
        self.num_sc_cables = len(set(cond['center_point'] for cond in active_conductors))
        self.unit_factor = unit_info['scale']
        self.label_unit = unit_info['label']
        self.color_map = {
            'core': 'gray',
            'XLPE': 'lightblue',
            'sheath': 'darkgreen',
            'HDPE': 'lightgreen',
            'core insulation': 'lightblue',
            'outer insulation': 'lightyellow',
            'insulation': 'white',
            'enclosure': 'darkgray',
            'soil': 'tan',
            'air_gap': 'white',
            'return': 'brown',
            'default': 'lightgray'
        }
        self.figsize = (12, 5)

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.case_name, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

    def _calculate_schematic_parameters(self) -> dict:
        """
        Calculates key parameters required for plotting the schematic.
        
        Returns:
            A dictionary containing parameters like max_radius, h_factor,
            title, and the primary core conductor data.
        """
        core_conductor = next((c for c in self.model.mtl.values() if c.get('line_type') == 'active'), None)
        
        max_radius = 0
        if core_conductor:
            for c in self.model.mtl.values():
                radius = c['radius'][1]
                if 'insulation' in c and c['insulation'] is not None:
                    radius += c['insulation']['thickness']
                if 'enclosure' in c and c['enclosure'] is not None:
                    radius = max(radius, c['enclosure']['outer_radius'])
                if radius > max_radius:
                    max_radius = radius
        
        if self.model.mtl_type == 'overhead':
            title = 'Overhead Transmission Line'
            h_factor = 8
        elif self.model.mtl_type == 'scc':
            title = 'Buried Single-Core Cable'
            h_factor = -2
        elif self.model.mtl_type == 'pipe':
            title = 'Pipe-Type Cable'
            h_factor = 1
        else: 
            title = 'Buried Single-Core Cable'
            h_factor = -2

        return {
            'core_conductor': core_conductor,
            'max_radius': max_radius,
            'h_factor': h_factor,
            'title': title
        }

    def _plot_cable_layer(self, ax, center, outer_radius_m, thickness_m, color, label, used_labels):
        """
        Plots a single cylindrical layer of a cable (as a Circle or Wedge).
        Handles the logic for avoiding duplicate legend entries.
        """
        outer_radius = outer_radius_m * self.unit_factor
        thickness = thickness_m * self.unit_factor
        
        label_to_plot = label if label not in used_labels else None
        if label_to_plot:
            used_labels.add(label_to_plot)
            
        if thickness < 1e-9 * self.unit_factor:
            patch = Circle(center, outer_radius, fill=True, edgecolor='black', 
                           facecolor=color, label=label_to_plot)
        else:
            patch = Wedge(center, outer_radius, 0, 360, width=thickness, 
                          edgecolor='black', facecolor=color, linestyle='solid',
                          label=label_to_plot)
        ax.add_patch(patch)
        
    def _finalize_plot(self, ax, title=""):
        """Applies final settings to a matplotlib axes object."""
        ax.relim()
        ax.autoscale_view()
        ax.margins(0.2)
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'x ({self.label_unit})')
        ax.set_ylabel(f'y ({self.label_unit})')
        ax.grid(True, linestyle='--', linewidth=0.5)
        if title:
            ax.set_title(title)
        
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            unique_labels = dict(zip(labels, handles))
            ax.legend(unique_labels.values(), unique_labels.keys(), loc='best')

# === Specialized Class for Isolated Systems ===

class IsolatedMTLRepresentation(BaseMTLRepresentation):
    """
    Creates graphical representations for isolated MTL systems,
    such as single wires, coaxial cables, and pipe-type cables where
    the ground effect is not the primary focus of the schematic.
    """
    def __init__(self, model: MulticonductorTransmissionLine, case_name, units):
        super().__init__(model, case_name, units)

    def wires(self, base_filename='system_schematic') -> None:
        """Plots the geometry of wires with insulating coating."""
        fig, ax = plt.subplots(figsize=self.figsize)

        for conductor in reversed(self.model.mtl.values()):
            conductor_radius = conductor['radius'][1] * self.unit_factor
            conductor_center = np.array(conductor['center_point']) * self.unit_factor

            if 'insulation' in conductor and conductor['insulation'] is not None:
                center = np.array(conductor['insulation']['center_point']) * self.unit_factor
                thickness = conductor['insulation']['thickness'] * self.unit_factor
                outer_radius = conductor_radius + thickness
                ax.add_patch(Wedge(center, outer_radius, 0, 360, width=thickness, edgecolor='black', facecolor='lightblue', linestyle='solid'))

            ax.add_patch(Circle(conductor_center, conductor_radius, fill=True, edgecolor='black', facecolor='darkgrey'))

        self._finalize_plot(ax, title='Multiconductor Transmission Line Cross-Section')
        save_figure_multiformat(fig, self.results_dir, base_filename)

    def coaxial_cable(self, base_filename='system_schematic') -> None:
        """Plots the geometry of a coaxial cable."""
        if 0 not in self.model.mtl or 1 not in self.model.mtl:
            print("Error: Coaxial data must contain keys for conductor 0 (sheath) and 1 (core).")
            return

        sheath = self.model.mtl[0]
        core = self.model.mtl[1]

        fig, ax = plt.subplots(figsize=self.figsize)
        center = np.array(core['center_point']) * self.unit_factor        
        core_outer_radius = core['radius'][1] * self.unit_factor
        sheath_inner_radius = sheath['radius'][0] * self.unit_factor
        sheath_outer_radius = sheath['radius'][1] * self.unit_factor
        
        sheath_thickness = sheath_outer_radius - sheath_inner_radius
        ax.add_patch(Wedge(center, sheath_outer_radius, 0, 360, width=sheath_thickness, 
                           edgecolor='black', facecolor='darkgrey', label='Sheath'))

        dielectric_thickness = sheath_inner_radius - core_outer_radius
        if dielectric_thickness > 0:
            ax.add_patch(Wedge(center, sheath_inner_radius, 0, 360, width=dielectric_thickness, 
                               edgecolor='black', facecolor='ivory', linestyle='--', label='Dielectric'))

        ax.add_patch(Circle(center, core_outer_radius, fill=True, edgecolor='black', 
                            facecolor='sandybrown', label='Core'))

        self._finalize_plot(ax, title='Coaxial Cable Cross-Section')
        save_figure_multiformat(fig, self.results_dir, base_filename)

    def system_schematic(self, base_filename='system_schematic') -> None:
        """
        Generates a schematic plot for isolated systems (e.g., pipe-type cables).
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        used_labels = set()
        params = self._calculate_schematic_parameters()

        if not params['core_conductor']:
            print("Error: Could not find an 'active' conductor in the model.")
            return

        for conductor_data in self.model.mtl.values():
            self._single_cable(ax, conductor_data, used_labels)

        self._finalize_plot(ax, params['title'])
        save_figure_multiformat(fig, self.results_dir, base_filename)

    def _single_cable(self, ax, conductor_data, used_labels):
        """Plots a single cable at its real coordinates."""
        center_point = np.array(conductor_data['center_point']) * self.unit_factor
        
        outer_radius_m = conductor_data['radius'][1]
        thickness_m = outer_radius_m - conductor_data['radius'][0]
        label = conductor_data.get('conductor_name', 'Conductor')
        color = self.color_map.get(label, self.color_map['default'])
        
        self._plot_cable_layer(ax, center_point, outer_radius_m, thickness_m, color, label, used_labels)
        
        if 'insulation' in conductor_data and conductor_data['insulation'] is not None:
            ins_data = conductor_data['insulation']
            ins_label = ins_data.get('type', 'Insulation')
            ins_color = self.color_map.get(ins_label, self.color_map['insulation'])
            self._plot_cable_layer(ax, center_point, outer_radius_m + ins_data['thickness'],
                                   ins_data['thickness'], ins_color, ins_label, used_labels)

    def _finalize_plot(self, ax, title):
        """Applies final settings for an isolated system plot."""
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'x ({self.label_unit})')
        ax.set_ylabel('')
        ax.set_yticks([])
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

# === Specialized Class for Ground-Return Systems ===

class GroundReturnMTLRepresentation(BaseMTLRepresentation):
    """
    Creates graphical representations for ground-return MTL systems,
    such as buried cables or overhead lines, where the ground plane
    is an essential part of the schematic.
    """
    def __init__(self, model: MulticonductorTransmissionLine, case_name, units):
        super().__init__(model, case_name, units)

    def system_schematic(self, base_filename='system_schematic') -> None:
        """
        Generates a schematic plot for ground-return systems.
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        used_labels = set()
        params = self._calculate_schematic_parameters()
        
        if not params['core_conductor']:
            print("Error: Could not find an 'active' conductor in the model.")
            return

        for conductor_data in self.model.mtl.values():
            self._scaled_single_cable(ax, conductor_data, params, used_labels)

        self._schematic_annotations(ax, params)
        
        x_margin_scale = 0.2 if self.num_sc_cables > 2 else 2.0
        self._finalize_plot(ax, params['title'], x_margin_scale)
        save_figure_multiformat(fig, self.results_dir, base_filename)

    def _scaled_single_cable(self, ax, conductor_data, params, used_labels):
        """
        Plots a single cable with schematic scaling for the y-position.
        """
        cable_center_point = np.array([
            conductor_data['center_point'][0], 
            params['h_factor'] * params['max_radius']
        ]) * self.unit_factor
        
        outer_radius_m = conductor_data['radius'][1]
        thickness_m = outer_radius_m - conductor_data['radius'][0]
        label = conductor_data.get('conductor_name', 'Conductor')
        color = self.color_map.get(label, self.color_map['default'])
        
        self._plot_cable_layer(ax, cable_center_point, outer_radius_m, thickness_m, color, label, used_labels)
        
        if 'insulation' in conductor_data and conductor_data['insulation'] is not None:
            ins_data = conductor_data['insulation']
            ins_label = ins_data.get('type', 'Insulation')
            ins_color = self.color_map.get(ins_label, self.color_map['insulation'])
            self._plot_cable_layer(ax, cable_center_point, outer_radius_m + ins_data['thickness'],
                                   ins_data['thickness'], ins_color, ins_label, used_labels)

        if 'enclosure' in conductor_data and conductor_data['enclosure'] is not None:
            enc_data = conductor_data['enclosure']
            enc_label = enc_data.get('type', 'HDPE')
            enc_color = self.color_map.get(enc_label, self.color_map['enclosure'])
            enc_thickness_m = enc_data['outer_radius'] - enc_data['inner_radius']
            
            real_cable_center = np.array(conductor_data['center_point'])
            real_enc_center = np.array(enc_data['center_point'])
            offset_vector = real_enc_center - real_cable_center
            enc_center_point_schematic = cable_center_point + (offset_vector * self.unit_factor)
            
            self._plot_cable_layer(ax, enc_center_point_schematic, enc_data['outer_radius'], 
                                   enc_thickness_m, enc_color, f'{enc_label}_enclosure', used_labels)

    def _schematic_annotations(self, ax, params):
        """Draws annotations like the ground level and depth/height line."""
        core_center_x = params['core_conductor']['center_point'][0]
        schematic_y = (params['h_factor'] * params['max_radius']) * self.unit_factor
        
        dim_x_arrow = (core_center_x + params['max_radius'] * 2.5) * self.unit_factor
        ax.annotate('', xy=(dim_x_arrow, 0), xycoords='data', 
                    xytext=(dim_x_arrow, schematic_y), textcoords='data',
                    arrowprops=dict(arrowstyle='<->', color='black', lw=1, zorder=3))
        
        real_h = abs(params['core_conductor']['center_point'][1])
        label_text = f'h = {real_h:.2f} m'
        
        ax.text(dim_x_arrow, schematic_y / 2, label_text, ha='center', va='center', fontsize=9, 
                zorder=3, bbox=dict(boxstyle='square,pad=0.3', fc='white', ec='none', alpha=0.8))

        ax.axhline(y=0, color='darkgreen', linestyle=':', linewidth=1.5, zorder=2)

    def _finalize_plot(self, ax, title, x_margin_scale=2.0):
        """Applies final settings for a ground-return system plot."""
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'x ({self.label_unit})')
        ax.set_ylabel('')
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

