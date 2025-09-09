import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge

from mtl_main.source import MulticonductorTransmissionLine
from utils.case_utils import UNITS_DATA

class MTLRepresentation:
    """ 
    This class creates a graphical representation of an MTL model.
    It operates on a MulticonductorTransmissionLine object.
    """
    def __init__(self, model: MulticonductorTransmissionLine, units='meter'):        
        self.model = model
        unit_info = UNITS_DATA.get(units, UNITS_DATA['meter'])
        active_conductors = [v for k, v in self.model.mtl.items()]
        self.num_sc_cables = len(set([cond['center_point'] for cond in active_conductors]))
        self.unit_factor = unit_info['scale']
        self.label_unit = unit_info['label']
        self.color_map = {
            'default': 'grey',
            'insulation': 'yellow',
            'core': 'darkgray',
            'sheath': 'darkgreen',
            'XLPE': 'lightblue',
            'HDPE': 'lightgreen',
            'pipe': 'lightgrey',
        }

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
            ax.legend()

    def isolated_wires(self):
        """
        This function plots the geometry of wires with insulating coating.
        """
        _, ax = plt.subplots(figsize=(8, 5))   

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

    def isolated_coaxial_cables(self):
        """
        This function plots the geometry of a coaxial cable,
        correctly interpreting a structure with separate core and sheath conductors.
        It elegantly defines the axis limits after plotting.
        """
        # Ensure the data structure has the expected conductors (0 and 1)
        if 0 not in self.model.mtl or 1 not in self.model.mtl:
            print("Error: The cable data must contain keys for conductor 0 (sheath) and 1 (core).")
            return

        # Assign core and sheath based on the provided structure
        # Conductor 0 is the sheath (return), Conductor 1 is the core (active)
        sheath = self.model.mtl[0]
        core = self.model.mtl[1]

        # Create the plot and axes
        _, ax = plt.subplots(figsize=(8, 5))
        
        # --- Define Radii and Center (applying scale factor) ---
        center = np.array(core['center_point']) * self.unit_factor        
        core_outer_radius = core['radius'][1] * self.unit_factor
        sheath_inner_radius = sheath['radius'][0] * self.unit_factor
        sheath_outer_radius = sheath['radius'][1] * self.unit_factor

        # --- Plotting from outside to inside ---

        # 1. Draw the outer sheath (conductor)
        sheath_thickness = sheath_outer_radius - sheath_inner_radius
        ax.add_patch(Wedge(
            center, sheath_outer_radius, 0, 360, 
            width=sheath_thickness, 
            edgecolor='black', 
            facecolor='darkgrey', 
            label='Sheath'
        ))

        # 2. Draw the dielectric insulator (the space between core and sheath)
        dielectric_thickness = sheath_inner_radius - core_outer_radius
        if dielectric_thickness > 0:
            ax.add_patch(Wedge(
                center, sheath_inner_radius, 0, 360, 
                width=dielectric_thickness, 
                edgecolor='black', 
                facecolor='ivory', 
                linestyle='--',
                label='Dielectric'
            ))

        # 3. Draw the inner core (conductor)
        ax.add_patch(Circle(
            center, core_outer_radius, 
            fill=True, 
            edgecolor='black', 
            facecolor='sandybrown', 
            label='Core'
        ))

        # --- Final plot adjustments ---
        self._finalize_plot(ax, title='Coaxial Cable Cross-Section')

    def isolated_systems(self):
        """
        Generates and displays a schematic plot for isolated systems 
        (e.g., single-core cables, pipe-type cables).
        
        This method acts as a coordinator, delegating tasks to specialized
        helper methods for plotting, annotation, and finalization.
        """
        # 1. Initialization
        _, ax = plt.subplots(figsize=(12, 6))
        used_labels = set()

        # 2. Calculate schematic parameters based on the model
        params = self._calculate_schematic_parameters()
        if not params['core_conductor']:
            print("Error: Could not find an 'active' conductor in the model.")
            return

        # 3. Plot each cable in the transmission line system
        for conductor_data in self.model.mtl.values():
            self._plot_single_cable(ax, conductor_data, params, used_labels)

        # 4. Add schematic annotations and finalize the plot
        self._finalize_isolated_plot(ax, params['title'])

    def ground_return_systems(self):
        """
        Generates and displays a schematic plot for ground-return systems 
        (e.g., buried cables or overhead lines).
        
        This method acts as a coordinator, delegating tasks to specialized
        helper methods for plotting, annotation, and finalization.
        """
        # 1. Initialization
        _, ax = plt.subplots(figsize=(12, 5))
        used_labels = set()

        # 2. Calculate schematic parameters based on the model
        params = self._calculate_schematic_parameters()
        if not params['core_conductor']:
            print("Error: Could not find an 'active' conductor in the model.")
            return

        # 3. Plot each cable in the transmission line system
        for conductor_data in self.model.mtl.values():
            self._plot_scaled_single_cable(ax, conductor_data, params, used_labels)

        # 4. Add schematic annotations and finalize the plot
        self._draw_schematic_annotations(ax, params)

        if self.num_sc_cables > 2:
            x_margin_scale = 0.2
        else:
            x_margin_scale = 2.0
        self._finalize_ground_return_plot(ax, params['title'], x_margin_scale)

    def _calculate_schematic_parameters(self):
        """
        Calculates key parameters required for plotting the schematic.
        
        Returns:
            dict: A dictionary containing parameters like max_radius, h_factor,
                  title, and the primary core conductor data.
        """
        core_conductor = next((c for c in self.model.mtl.values() if c.get('line_type') == 'active'), None)
        
        max_radius = 0
        if core_conductor:
            for c in self.model.mtl.values():
                radius = c['radius'][1]
                if 'insulation' in c and c['insulation'] is not None:
                    radius += c['insulation']['thickness']
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

    def _plot_single_cable(self, ax, conductor_data, params, used_labels):
        """
        Plots a single cable with all its constituent layers (conductor, insulation).
        
        Args:
            ax (matplotlib.axes.Axes): The axes to plot on.
            conductor_data (dict): The data dictionary for one conductor.
            params (dict): The pre-calculated schematic parameters.
            used_labels (set): A set of labels already used in the legend.
        """
        # Calculate the schematic center point for this specific cable
        center_point = np.array([
            conductor_data['center_point'][0], 
            conductor_data['center_point'][1]
        ]) * self.unit_factor
        
        # Plot the main conductor body
        outer_radius_m = conductor_data['radius'][1]
        thickness_m = outer_radius_m - conductor_data['radius'][0]
        label = conductor_data.get('conductor_name', 'Conductor')
        color = self.color_map.get(label, self.color_map['default'])
        
        self._plot_cable_layer(
            ax, center_point, outer_radius_m, thickness_m,
            color, label, used_labels
        )
        
        # Plot the insulation layer if it exists
        if 'insulation' in conductor_data and conductor_data['insulation'] is not None:
            ins_data = conductor_data['insulation']
            ins_label = ins_data.get('type', 'Insulation')
            ins_color = self.color_map.get(ins_label, self.color_map['insulation'])
            
            self._plot_cable_layer(
                ax, center_point, outer_radius_m + ins_data['thickness'],
                ins_data['thickness'], ins_color, ins_label, used_labels
            )

    def _plot_scaled_single_cable(self, ax, conductor_data, params, used_labels):
        """
        Plots a single cable with all its constituent layers (conductor, insulation).
        
        Args:
            ax (matplotlib.axes.Axes): The axes to plot on.
            conductor_data (dict): The data dictionary for one conductor.
            params (dict): The pre-calculated schematic parameters.
            used_labels (set): A set of labels already used in the legend.
        """
        # Calculate the schematic center point for this specific cable
        center_point = np.array([
            conductor_data['center_point'][0], 
            params['h_factor'] * params['max_radius']
        ]) * self.unit_factor
        
        # Plot the main conductor body
        outer_radius_m = conductor_data['radius'][1]
        thickness_m = outer_radius_m - conductor_data['radius'][0]
        label = conductor_data.get('conductor_name', 'Conductor')
        color = self.color_map.get(label, self.color_map['default'])
        
        self._plot_cable_layer(
            ax, center_point, outer_radius_m, thickness_m,
            color, label, used_labels
        )
        
        # Plot the insulation layer if it exists
        if 'insulation' in conductor_data and conductor_data['insulation'] is not None:
            ins_data = conductor_data['insulation']
            ins_label = ins_data.get('type', 'Insulation')
            ins_color = self.color_map.get(ins_label, self.color_map['insulation'])
            
            self._plot_cable_layer(
                ax, center_point, outer_radius_m + ins_data['thickness'],
                ins_data['thickness'], ins_color, ins_label, used_labels
            )

    def _plot_cable_layer(self, ax, center, outer_radius_m, thickness_m, color, label, used_labels):
        """
        Plots a single cylindrical layer of a cable (as a Circle or Wedge).
        Handles the logic for avoiding duplicate legend entries.
        """
        # Apply scaling
        outer_radius = outer_radius_m * self.unit_factor
        thickness = thickness_m * self.unit_factor
        
        # Determine if the label should be added to the legend
        label_to_plot = label if label not in used_labels else None
        if label_to_plot:
            used_labels.add(label_to_plot)
            
        # A thickness of 0 implies a solid core
        if thickness < 1e-9 * self.unit_factor: # Use a small tolerance for floating point
            patch = Circle(center, outer_radius, fill=True, edgecolor='black', 
                           facecolor=color, label=label_to_plot)
        else:
            patch = Wedge(center, outer_radius, 0, 360, width=thickness, 
                          edgecolor='black', facecolor=color, linestyle='solid',
                          label=label_to_plot)
        ax.add_patch(patch)

    def _draw_schematic_annotations(self, ax, params):
        """
        Draws annotations like the ground level, depth dimension line, and
        surrounding medium text.
        """
        core_center_x = params['core_conductor']['center_point'][0]
        schematic_y = (params['h_factor'] * params['max_radius']) * self.unit_factor
        
        # Dimension line for burial depth 'h'
        dim_x_arrow = (core_center_x + params['max_radius'] * 2.5) * self.unit_factor
        ax.annotate('', xy=(dim_x_arrow, 0), xycoords='data', 
                    xytext=(dim_x_arrow, schematic_y), textcoords='data',
                    arrowprops=dict(arrowstyle='<->', color='black', lw=1, zorder=3))
        
        real_h = abs(params['core_conductor']['center_point'][1])
        label_text = f'h = {real_h:.2f} m'
        
        # Adjust text position and zorder
        ax.text(dim_x_arrow, schematic_y / 2, label_text, 
                ha='center', va='center', fontsize=9, zorder=3,
                bbox=dict(boxstyle='square,pad=0.3', fc='white', ec='none', alpha=0.8))

        # Ground level line
        ax.axhline(y=0, color='darkgreen', linestyle=':', linewidth=1.5, zorder=2)

    def _finalize_ground_return_plot(self, ax, title, x_margin_scale = 2.0):
        """
        Applies final settings to the plot, including limits, labels,
        backgrounds, and legend formatting.
        """
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'x ({self.label_unit})')
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.set_title(title)
        ax.grid(True, linestyle='--', linewidth=0.5, zorder=0)

        # Set plot limits with generous margins
        ax.relim()
        ax.autoscale_view()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        x_margin = (xlim[1] - xlim[0]) * x_margin_scale
        y_margin = (ylim[1] - ylim[0]) * 0.2
        ax.set_xlim(xlim[0] - x_margin, xlim[1] + x_margin)
        ax.set_ylim(ylim[0] - y_margin, ylim[1] + y_margin)
        final_xlim = ax.get_xlim()
        final_ylim = ax.get_ylim()
        
        # Add background fills and text
        ax.fill_between(final_xlim, final_ylim[0], 0, color='saddlebrown', alpha=0.2, zorder=1)
        
        text_x = final_xlim[0] + (final_xlim[1] - final_xlim[0]) * 0.05
        y_offset = (final_ylim[1] - final_ylim[0]) * 0.03
        ax.text(text_x, y_offset, 'Air ($\\varepsilon_0$, $\\mu_0$)', va='bottom', fontsize=10, style='italic', zorder=4)
        ax.text(text_x, -y_offset, 'Ground ($\\varepsilon_1$, $\\mu_1$, $\\sigma_1$)', va='top', fontsize=10, style='italic', zorder=4)

        # Format legend: create it first, then set its zorder.
        handles, labels = ax.get_legend_handles_labels()
        if handles: 
            unique_labels = dict(zip(labels, handles))
            legend = ax.legend(unique_labels.values(), unique_labels.keys(), loc='best')
            legend.set_zorder(5)

    def _finalize_isolated_plot(self, ax, title):
        """
        Applies final settings to the plot, including limits, labels,
        backgrounds, and legend formatting.
        """
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'x ({self.label_unit})')
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.set_title(title)
        ax.grid(True, linestyle='--', linewidth=0.5, zorder=0)

        # Set plot limits with generous margins
        ax.relim()
        ax.autoscale_view()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        x_margin = (xlim[1] - xlim[0]) * 2.0
        y_margin = (ylim[1] - ylim[0]) * 0.2
        ax.set_xlim(xlim[0] - x_margin, xlim[1] + x_margin)
        ax.set_ylim(ylim[0] - y_margin, ylim[1] + y_margin)
        final_xlim = ax.get_xlim()
        final_ylim = ax.get_ylim()
        
        text_x = final_xlim[0] + (final_xlim[1] - final_xlim[0]) * 0.05
        y_offset = (final_ylim[1] - final_ylim[0]) * 0.35
        ax.text(text_x, y_offset, 'Air ($\\varepsilon_0$, $\\mu_0$)', va='bottom', fontsize=10, style='italic', zorder=4)

        # Format legend: create it first, then set its zorder.
        handles, labels = ax.get_legend_handles_labels()
        if handles: 
            unique_labels = dict(zip(labels, handles))
            legend = ax.legend(unique_labels.values(), unique_labels.keys(), loc='best')
            legend.set_zorder(5)