import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib import patches

from mtl_data.mtl import MulticonductorTransmissionLine

class MTLRepresentation(MulticonductorTransmissionLine):
    """ This class defines the coaxial cable with a sheath. """

    def bare_and_coated_wires(self):
        """
        Esta função plota a geometria de fios com revestimento isolante,
        utilizando uma abordagem elegante para definir os limites dos eixos.
        """
        _, ax = plt.subplots(figsize=(8, 5))

        # Itera e desenha todos os condutores e isolamentos
        for conductor in reversed(self.mtl.values()):
            conductor_radius = conductor['radius'][1]
            conductor_center = conductor['center_point']

            # Desenha a camada de isolamento se existir
            if 'sheath' in conductor and conductor['sheath'] is not None:
                center = conductor['sheath']['center_point']
                thickness = conductor['sheath']['thickness']
                outer_radius = conductor_radius + thickness
                ax.add_patch(Wedge(center, outer_radius, 0, 360, width=thickness, edgecolor='black', facecolor='lightblue', linestyle='solid'))

            # Desenha o condutor principal (núcleo)            
            ax.add_patch(patches.Circle(conductor_center, conductor_radius, fill=True, edgecolor='black', facecolor='darkgrey'))

        ax.relim()                      # Recalcula os limites dos dados para incluir todos os patches
        ax.autoscale_view()             # Ajusta a visualização para os novos limites
        ax.margins(0.2)                 # Aplica uma margem de 20% aos limites calculados
        ax.set_aspect('equal', 'box')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.grid(True, linestyle='--', linewidth=0.5)

    # def underground_system(self):
    #     """This function plots the geometry of the buried SCC."""
    #     _, ax = plt.subplots(figsize=(10, 5))
    #     all_points = []
    #     facecolor_list = ['black', 'lightblue', 'lightblue']

    #     # Plot the conductors and insulations
    #     for i, item in enumerate(reversed(self.mtl)):

    #         # Plot the insulation using Wedge
    #         ax.add_patch(Wedge(item['center_point'], item['radius'][1], 0, 360,
    #                            width=item['radius'][1] - item['radius'][0],
    #                            edgecolor='black',
    #                            facecolor=facecolor_list[i],
    #                            linestyle='solid'))

    #         # Add the center point to the list of points
    #         all_points.append((item['center_point'][0],
    #                            item['center_point'][1],
    #                            item['radius'][1]))

    #     # Set the aspect of the plot to be equal
    #     # ax.set_aspect('equal', 'box')
    #     ax.set_aspect('equal', adjustable='box')
    #     plt.xlabel('x')
    #     plt.ylabel('y')
    #     plt.title('Geometry of the underground system.')
    #     plt.grid(True, linestyle='--', linewidth=0.5)

    #     # Determine axis limits
    #     margin = 0.2  # 20 % margin for better visualization
    #     max_x = max(point[0] + point[2] for point in all_points)
    #     min_x = min(point[0] - point[2] for point in all_points)
    #     min_y = min(point[1] - point[2] for point in all_points)

    #     # Adjust the y-axis to show the ground level
    #     min_x_level = 10*(min_x - margin * abs(min_x))
    #     max_x_level = 10*(max_x + margin * abs(max_x))
    #     min_y_level = min_y - margin * abs(min_y)
    #     plt.xlim(min_x_level, max_x_level)
    #     plt.ylim(min_y_level, 0.2)

    #     # Draw the ground level reference line
    #     # plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)

    #     # Fill the ground and air regions and Add labels
    #     ax.fill_between([min_x_level, max_x_level], min_y_level, 0,
    #                     color='lightgrey', alpha=0.5)
    #     plt.text(min_x_level, 0.1,
    #              'air ($\\varepsilon_0$, $\\mu_0$)',
    #              verticalalignment='center', fontsize=10)
    #     plt.text(min_x_level, - 0.15,
    #              'ground ($\\varepsilon_0$, $\\mu_0$, $\\sigma_g$)',
    #              verticalalignment='center', fontsize=10)

    #     # plt.legend()
    #     plt.show()

    # def single_core_cable(self):
    #     """This function plots the geometry of the bifilar line or coaxial cable."""
    #     _, ax = plt.subplots(figsize=(10, 5))
    #     all_points = []

    #     # Plot the conductors and insulations
    #     for item in reversed(self.mtl):
    #         # Plot the insulation covering the conductors
    #         if item['insulation']['thickness'] is not None:
    #             # Calculate the insulation inner and outer radius
    #             ins_inner_radius = item['radius'][1]
    #             ins_outer_radius = item['radius'][1] + item['insulation']['thickness']

    #             # Add the insulation to the list of points
    #             all_points.append((item['center_point'][0], item['center_point'][1], ins_outer_radius))

    #             # Plot the insulation using Wedge
    #             wedge = Wedge(item['center_point'], ins_outer_radius, 0, 360,
    #                           width=ins_outer_radius - ins_inner_radius,
    #                           edgecolor='black', facecolor='lightblue', linestyle='solid')
    #             ax.add_patch(wedge)

    #         # Plot the outer conductor using Wedge
    #         outer_wedge = Wedge(item['center_point'], item['radius'][1], 0, 360,
    #                             width=item['radius'][1] - item['radius'][0],
    #                             edgecolor='black', facecolor='grey', linestyle='dotted')
    #         ax.add_patch(outer_wedge)

    #         # Plot the inner conductor as a circle
    #         inner_circle = patches.Circle(item['center_point'], item['radius'][0], fill=True,
    #                                       edgecolor='black', facecolor='white', linestyle='dotted')
    #         ax.add_patch(inner_circle)

    #         # Add the center point to the list of points
    #         all_points.append((item['center_point'][0], item['center_point'][1], item['radius'][1]))

    #     # Set the aspect of the plot to be equal
    #     # ax.set_aspect('equal', 'box')
    #     ax.set_aspect('equal', adjustable='box')
    #     plt.xlabel('x')
    #     plt.ylabel('y')
    #     plt.title('Geometry of the underground system.')
    #     plt.grid(True, linestyle='--', linewidth=0.5)

    #     # Determine axis limits
    #     margin = 0.2  # 20 % margin for better visualization
    #     max_x = max(point[0] + point[2] for point in all_points)
    #     min_x = min(point[0] - point[2] for point in all_points)
    #     min_y = min(point[1] - point[2] for point in all_points)

    #     # Adjust the y-axis to show the ground level
    #     min_x_level = min_x - margin * abs(min_x)
    #     max_x_level = max_x + margin * abs(max_x)
    #     min_y_level = min_y - margin * abs(min_y)
    #     plt.xlim(min_x_level, max_x_level)
    #     plt.ylim(min_y_level, 0.2)

    #     # Draw the ground level reference line
    #     # plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)

    #     # Fill the ground and air regions and Add labels
    #     ax.fill_between([min_x_level, max_x_level], min_y_level, 0, color='lightgrey', alpha=0.5)
    #     plt.text(min_x_level, 0.1, 'air ($\\varepsilon_0$, $\\mu_0$)', verticalalignment='center', fontsize=10)
    #     plt.text(min_x_level, - 0.15, 'ground ($\\varepsilon_0$, $\\mu_0$, $\\sigma_g$)', verticalalignment='center', fontsize=10)

    #     # plt.legend()
    #     plt.show()

    # def plot_trefoil(self):
    #     """This function extracts the conductors from the list."""
    #     _, ax = plt.subplots()
    #     all_points = []

    #     # Core subconductors list
    #     core_sublist = [conductor
    #                     for conductor in self.mtl
    #                     if conductor['conductor_name'] == 'core']

    #     # Screen subconductors list
    #     screen_sublist = [conductor
    #                       for conductor in self.mtl
    #                       if conductor['conductor_name'] == 'screen']

    #     # Screen subconductors list
    #     armor_sublist = [conductor
    #                      for conductor in self.mtl
    #                      if conductor['conductor_name'] == 'armor']

    #     # Jacket dimensions
    #     jack_thickness = [conductor['insulation']['thickness']
    #                       for conductor in self.original_mtl
    #                       if conductor['conductor_name'] == 'armor'][0]

    #     # Armor outer radius
    #     arm_out_radius = [conductor['radius'][1]
    #                       for conductor in self.original_mtl
    #                       if conductor['conductor_name'] == 'armor'][0]

    #     # Create the core and screen plots
    #     for i, core in enumerate((core_sublist)):
    #         # core conductor
    #         cx, cy = core['center_point']
    #         ax.add_patch(patches.Circle(core['center_point'], self.core_radius,
    #                                     fill=True, edgecolor='black', facecolor='darkgrey'))

    #         # Add the core number
    #         ax.text(cx, cy, f'Core {i+1}', ha='center',
    #                 va='center', fontsize=10)

    #         # Primary insulation
    #         prim_ins_radius = self.core_radius + self.prim_ins_thick
    #         ax.add_patch(patches.Circle((cx, cy), prim_ins_radius,
    #                                     fill=False, edgecolor='blue', linestyle='dashed'))

    #         # Screen wires
    #         for i, screen in enumerate((screen_sublist)):
    #             # core conductor
    #             ax.add_patch(patches.Circle(screen['center_point'], screen['radius'][1],
    #                                         fill=True, edgecolor='gray', facecolor='lightblue'))

    #             # Add the core number
    #             sc_x, sc_y = screen['center_point']
    #             ax.text(sc_x, sc_y, f'{i+1}',
    #                     ha='center', va='center', fontsize=6)

    #         # Secondary insulation
    #         ax.add_patch(patches.Circle((cx, cy), self.scc_out_radius,
    #                                     fill=False, edgecolor='black', linestyle='solid'))

    #     # Create the armor plots
    #     for i, armor in enumerate((armor_sublist)):
    #         cx, cy = armor['center_point']
    #         ax.add_patch(patches.Circle(armor['center_point'], armor['radius'][1],
    #                                     fill=True, edgecolor='gray', facecolor='lightgreen'))

    #         # Add the armor number
    #         ax.text(cx, cy, f'{i+1}', ha='center', va='center', fontsize=8)

    #     # Jacket and armor dimensions
    #     jack_inn_radius = self.trefoil_radius + self.scc_out_radius
    #     jack_out_radius = jack_inn_radius + jack_thickness
    #     ax.add_patch(patches.Circle((0, 0), jack_inn_radius,
    #                                 fill=False, edgecolor='black', linestyle='solid'))
    #     ax.add_patch(patches.Circle((0, 0), jack_out_radius,
    #                                 fill=False, edgecolor='blue', linestyle='dashed'))
    #     ax.add_patch(patches.Circle((0, 0), arm_out_radius,
    #                                 fill=False, edgecolor='black', linestyle='solid'))

    #     # Include the relevant points for the axis limits
    #     all_points.append((0, 0, arm_out_radius))

    #     # Determine axis limits
    #     margin = 0.05  # 5% margin for better visualization
    #     max_x = max(p[0] + p[2] for p in all_points)
    #     min_x = min(p[0] - p[2] for p in all_points)
    #     max_y = max(p[1] + p[2] for p in all_points)
    #     min_y = min(p[1] - p[2] for p in all_points)

    #     # Set the aspect of the plot to be equal
    #     ax.set_aspect('equal', 'box')
    #     plt.xlabel('x')
    #     plt.ylabel('y')
    #     plt.title('Trefoil Power Cable Geometry.')
    #     plt.grid(True)
    #     plt.xlim(min_x - margin * abs(min_x), max_x + margin * abs(max_x))
    #     plt.ylim(min_y - margin * abs(min_y), max_y + margin * abs(max_y))
    #     plt.show()

    # def plot_trefoil_from_original_data(self, conductors_list):
    #     """This function plots the geometry of the trefoil power cable."""
    #     _, ax = plt.subplots()

    #     def rotate_point(
    #             point, center_point, angle):
    #         """This function rotates a point around another point by a given angle."""
    #         angle_rad = np.radians(angle)

    #         x, y = point
    #         x0, y0 = center_point

    #         # Rotation x-value
    #         x_rot = x0 + (x - x0) * np.cos(angle_rad) - \
    #             (y - y0) * np.sin(angle_rad)

    #         # Rotation y-value
    #         y_rot = y0 + (x - x0) * np.sin(angle_rad) + \
    #             (y - y0) * np.cos(angle_rad)

    #         return (x_rot, y_rot)

    #     def plot_wires_subconductors(
    #             wires, layer_radius, wire_radius, base_center, fig, shift=0):
    #         """This function plots the geometry of the subconductors."""
    #         base_x, base_y = base_center

    #         for m in range(wires):
    #             # Calculate the angle for each subconductor
    #             angle_step = 2 * np.pi / wires

    #             # angle = m * angle_step + shift / 2 * angle_step
    #             angle = angle_step * (m + shift / 2)

    #             # Calculate the position of the subconductor
    #             cx = base_x + layer_radius * np.cos(angle)
    #             cy = base_y + layer_radius * np.sin(angle)

    #             # Plot the subconductor
    #             fig.add_patch(patches.Circle((cx, cy), wire_radius, fill=True,
    #                                          edgecolor='gray', facecolor='lightblue'))

    #             # Add the subconductor number
    #             fig.text(cx, cy, f'{m+1}', ha='center',
    #                      va='center', fontsize=6)

    #     mtl = conductors_list['data']
    #     all_points = []
    #     core = [c for c in mtl if c['conductor_name'] == 'core'][0]
    #     screen = [c for c in mtl if c['conductor_name'] == 'screen'][0]
    #     armor = [c for c in mtl if c['conductor_name'] == 'armor'][0]

    #     # primary insulation
    #     prim_insul_radius = core['radius'][1] + core['insulation']['thickness']

    #     # screen wires
    #     screen_radius = prim_insul_radius + \
    #         screen['subconductors']['wire_radius']

    #     # secondary insulation
    #     sec_insul_radius = prim_insul_radius + self.sec_ins_thick

    #     # Jacket dimensions
    #     jacket_inner_radius = self.trefoil_radius + self.scc_out_radius
    #     jacket_outer_radius = jacket_inner_radius + \
    #         armor['insulation']['thickness']

    #     # Define the core positions for the trefoil power cable
    #     core_positions = [(self.trefoil_radius, 0),
    #                       rotate_point((self.trefoil_radius, 0), (0, 0), 120),
    #                       rotate_point((self.trefoil_radius, 0), (0, 0), 240)]

    #     # Plot the core and screen subconductors
    #     for i, (cx, cy) in enumerate(core_positions):
    #         # core conductor
    #         ax.add_patch(patches.Circle((cx, cy), core['radius'][1],
    #                                     fill=True, edgecolor='black', facecolor='darkgrey'))

    #         ax.text(cx, cy, f'Core {i+1}', ha='center',
    #                 va='center', fontsize=10)

    #         # primary insulation
    #         ax.add_patch(patches.Circle((cx, cy), prim_insul_radius,
    #                                     fill=False, edgecolor='blue', linestyle='dashed'))

    #         # screen wires
    #         plot_wires_subconductors(screen['subconductors']['wires_per_layer'],
    #                                  screen_radius,
    #                                  screen['subconductors']['wire_radius'],
    #                                  (cx, cy), ax, shift=0)

    #         # secondary insulation
    #         ax.add_patch(patches.Circle((cx, cy), sec_insul_radius, fill=False,
    #                                     edgecolor='black', linestyle='solid'))

    #     # Jacket layer
    #     ax.add_patch(patches.Circle((0, 0), jacket_inner_radius,
    #                                 fill=False, edgecolor='black', linestyle='solid'))
    #     ax.add_patch(patches.Circle((0, 0), jacket_outer_radius,
    #                                 fill=False, edgecolor='blue', linestyle='dashed'))

    #     # Armoured layer
    #     ax.add_patch(patches.Circle((0, 0), armor['radius'][1],
    #                                 fill=False, edgecolor='black', linestyle='solid'))

    #     # Armor wires
    #     for k in range(armor['subconductors']['layers']):
    #         armor_radius = armor['radius'][1] - \
    #             (2*k + 1) * armor['subconductors']['wire_radius']

    #         plot_wires_subconductors(armor['subconductors']['wires_per_layer'],
    #                                  armor_radius,
    #                                  armor['subconductors']['wire_radius'],
    #                                  (0, 0), ax, shift=k)

    #     # Include the relevant points for the axis limits
    #     all_points.append((0, 0, armor['radius'][1]))
    #     margin = 0.05  # 5% margin for better visualization
    #     max_x = max(p[0] + p[2] for p in all_points)
    #     min_x = min(p[0] - p[2] for p in all_points)
    #     max_y = max(p[1] + p[2] for p in all_points)
    #     min_y = min(p[1] - p[2] for p in all_points)

    #     # Set the aspect of the plot to be equal
    #     ax.set_aspect('equal', 'box')
    #     plt.xlabel('x')
    #     plt.ylabel('y')
    #     plt.title('Trefoil Power Cable Geometry from Original Data')
    #     plt.grid(True)
    #     plt.xlim(min_x - margin * abs(min_x), max_x + margin * abs(max_x))
    #     plt.ylim(min_y - margin * abs(min_y), max_y + margin * abs(max_y))
    #     plt.show()

