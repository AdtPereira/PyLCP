import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge

from mtl_main.source import MulticonductorTransmissionLine

UNITS_DATA = {
    'meter':      {'scale': 1,       'label': 'm'},
    'centimeter': {'scale': 100,     'label': 'cm'},
    'millimeter': {'scale': 1000,    'label': 'mm'},
    'mil':        {'scale': 39370.1, 'label': 'mil'},
}

class MTLRepresentation:
    """ 
    This class creates a graphical representation of an MTL model.
    It operates on a MulticonductorTransmissionLine object.
    """
    def __init__(self, model: MulticonductorTransmissionLine, units='meter'):        
        self.model = model
        unit_info = UNITS_DATA.get(units, UNITS_DATA['meter'])
        self.scale_factor = unit_info['scale']
        self.label_unit = unit_info['label']
        self.color_map = {
            'default': 'darkgrey',
            'primary_insulation': 'lightblue',
            'sheath': 'dimgray',
            'jacket': 'black'
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
            conductor_radius = conductor['radius'][1] * self.scale_factor
            conductor_center = np.array(conductor['center_point']) * self.scale_factor

            if 'insulation' in conductor and conductor['insulation'] is not None:
                center = np.array(conductor['insulation']['center_point']) * self.scale_factor
                thickness = conductor['insulation']['thickness'] * self.scale_factor
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
        center = np.array(core['center_point']) * self.scale_factor        
        core_outer_radius = core['radius'][1] * self.scale_factor
        sheath_inner_radius = sheath['radius'][0] * self.scale_factor
        sheath_outer_radius = sheath['radius'][1] * self.scale_factor

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

    def ground_return_systems(self):
        """
        This function plots a simplified schematic of a single-core buried cable,
        including any defined insulation layers around the conductors.
        The burial depth 'h' is graphically represented as h = 3*r_max, where r_max
        is the largest radius of the cable's outer layer, including insulation.
        """
        _, ax = plt.subplots(figsize=(8, 6))

        core_conductor = None
        for cond in self.model.mtl.values():
            if cond.get('line_type') == 'active':
                core_conductor = cond
                break
                
        if not core_conductor:
            print("Error: Could not find an 'active' conductor.")
            return

        max_radius = 0
        for conductor in self.model.mtl.values():
            if 'radius' in conductor:
                current_outer_radius = conductor['radius'][1]
                if 'insulation' in conductor and conductor['insulation'] is not None:
                    current_outer_radius += conductor['insulation']['thickness']
                if current_outer_radius > max_radius:
                    max_radius = current_outer_radius
        
        if self.model.mtl_type == 'overhead':
            mtl_title = 'Overhead Transmission Line'
            h_factor = 5
        elif self.model.mtl_type == 'scc':
            mtl_title = 'Buried Single-Core Cable'
            h_factor = -3

        center_x = core_conductor['center_point'][0] * self.scale_factor
        center_y_depth = h_factor * max_radius * self.scale_factor
        plot_center_point = np.array([center_x, center_y_depth])

        # --- Plot Conductor and Insulation Layers ---
        sorted_conductors = sorted(self.model.mtl.values(), key=lambda c: c.get('radius', [0,0])[1])
        for conductor in sorted_conductors:
            if 'radius' not in conductor:
                continue
            conductor_outer_radius = conductor['radius'][1] * self.scale_factor
            conductor_inner_radius = conductor['radius'][0] * self.scale_factor
            conductor_thickness = conductor_outer_radius - conductor_inner_radius
            label = conductor.get('conductor_name', 'Conductor').capitalize()
            color = self.color_map.get(label.lower(), self.color_map['default'])
            if conductor_thickness > 0:
                ax.add_patch(Wedge(plot_center_point, conductor_outer_radius, 0, 360, width=conductor_thickness,
                                edgecolor='black', facecolor=color, linestyle='solid', label=label))
            else:
                ax.add_patch(Circle(plot_center_point, conductor_outer_radius, fill=True, 
                                    edgecolor='black', facecolor=color, label=label))
            if 'insulation' in conductor and conductor['insulation'] is not None:
                insulation_data = conductor['insulation']
                insulation_thickness = insulation_data['thickness'] * self.scale_factor
                insulation_outer_radius = conductor_outer_radius + insulation_thickness
                insulation_label = insulation_data.get('name', 'Insulation').replace('_', ' ').capitalize()
                insulation_color = self.color_map.get(insulation_data.get('name'), 'cyan')
                ax.add_patch(Wedge(plot_center_point, insulation_outer_radius, 0, 360, width=insulation_thickness,
                    edgecolor='black', facecolor=insulation_color, linestyle='solid', label=insulation_label))

        # --- Finalize Plot and Add Schematic Ground ---
        ax.axhline(y=0, color='darkgreen', linestyle=':', linewidth=1.5, label='Ground Level')
        ax.relim()
        ax.autoscale_view()

        # 3. Manually calculate margins and set final, explicit limits
        tight_xlim = ax.get_xlim()
        x_margin = (tight_xlim[1] - tight_xlim[0]) * 0.45
        final_xmin = tight_xlim[0] - x_margin
        final_xmax = tight_xlim[1] + x_margin
        ax.set_xlim(final_xmin, final_xmax)

        tight_ylim = ax.get_ylim()
        y_margin = (tight_ylim[1] - tight_ylim[0]) * 0.2
        final_ymin = tight_ylim[0] - y_margin
        final_ymax = tight_ylim[1] + y_margin
        ax.set_ylim(final_ymin, final_ymax)

        # 4. Apply other plot settings
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'x ({self.label_unit})')
        ax.set_ylabel('')
        ax.grid(True, linestyle='--', linewidth=0.5, zorder=0)
        ax.set_title(mtl_title)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles[::-1], labels=labels[::-1], loc='upper right')
        ax.set_yticks([])

        # 5. Fill background and place text using the final limits
        ax.fill_between([final_xmin, final_xmax], final_ymin, 0, color='saddlebrown', alpha=0.2)
        plot_height = final_ymax - final_ymin
        text_x = final_xmin + (final_xmax - final_xmin) * 0.05
        y_offset = plot_height * 0.03
        ax.text(text_x, +y_offset, 'Air ($\\varepsilon_0$, $\\mu_0$)', verticalalignment='bottom', fontsize=10, style='italic')
        ax.text(text_x, -y_offset, 'Ground ($\\varepsilon_1$, $\\mu_1$, $\\sigma_1$)', verticalalignment='top', fontsize=10, style='italic')    
      