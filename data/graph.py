"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

REFERENCES:
[1] 

"""
import copy
import matplotlib.pyplot as plt
from scipy.constants import mu_0, epsilon_0
from matplotlib.patches import Wedge
from matplotlib import patches

class GraphicRepresentation():
    """ This class plots the geometry of the system. """

    def __init__(self, mtl):
        """Initialize the MulticonductorTransmissionLine class.

        Key Points
        Deep Copy: When initializing the MulticonductorTransmissionLine class, 
        we use copy.deepcopy to ensure that self.conductors is an independent 
        copy of TREFOIL.

        Isolation of Modifications: Any modifications made to self.conductors 
        within the class will not affect the original TREFOIL object.
        """

        # Deep copy of the conductors list
        self.mtl = copy.deepcopy(mtl['data'])

    def underground_system(self):
        """This function plots the geometry of the buried SCC."""
        _, ax = plt.subplots(figsize=(10, 5))
        all_points = []
        facecolor_list = ['black', 'lightblue', 'lightblue']

        # Plot the conductors and insulations
        for i, item in enumerate(reversed(self.mtl)):

            # Plot the insulation using Wedge
            ax.add_patch(Wedge(item['center_point'], item['radius'][1], 0, 360,
                               width=item['radius'][1] - item['radius'][0],
                               edgecolor='black',
                               facecolor=facecolor_list[i],
                               linestyle='solid'))

            # Add the center point to the list of points
            all_points.append((item['center_point'][0],
                               item['center_point'][1],
                               item['radius'][1]))

        # Set the aspect of the plot to be equal
        # ax.set_aspect('equal', 'box')
        ax.set_aspect('equal', adjustable='box')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Geometry of the underground system.')
        plt.grid(True, linestyle='--', linewidth=0.5)

        # Determine axis limits
        margin = 0.2  # 20 % margin for better visualization
        max_x = max(point[0] + point[2] for point in all_points)
        min_x = min(point[0] - point[2] for point in all_points)
        min_y = min(point[1] - point[2] for point in all_points)

        # Adjust the y-axis to show the ground level
        min_x_level = 10*(min_x - margin * abs(min_x))
        max_x_level = 10*(max_x + margin * abs(max_x))
        min_y_level = min_y - margin * abs(min_y)
        plt.xlim(min_x_level, max_x_level)
        plt.ylim(min_y_level, 0.2)

        # Draw the ground level reference line
        # plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)

        # Fill the ground and air regions and Add labels
        ax.fill_between([min_x_level, max_x_level], min_y_level, 0,
                        color='lightgrey', alpha=0.5)
        plt.text(min_x_level, 0.1,
                 'air ($\\varepsilon_0$, $\\mu_0$)',
                 verticalalignment='center', fontsize=10)
        plt.text(min_x_level, - 0.15,
                 'ground ($\\varepsilon_0$, $\\mu_0$, $\\sigma_g$)',
                 verticalalignment='center', fontsize=10)

        # plt.legend()
        plt.show()

    def single_core_cable(self):
        """This function plots the geometry of the bifilar line or coaxial cable."""
        _, ax = plt.subplots(figsize=(10, 5))
        all_points = []

        # Plot the conductors and insulations
        for item in reversed(self.mtl):

            # Plot the insulation covering the conductors
            if item['insulation']['thickness'] is not None:
                # Calculate the insulation inner and outer radius
                ins_inner_radius = item['radius'][1]
                ins_outer_radius = item['radius'][1] + \
                    item['insulation']['thickness']

                # Add the insulation to the list of points
                all_points.append(
                    (item['center_point'][0], item['center_point'][1], ins_outer_radius))

                # Plot the insulation using Wedge
                wedge = Wedge(item['center_point'], ins_outer_radius, 0, 360,
                              width=ins_outer_radius - ins_inner_radius,
                              edgecolor='black', facecolor='lightblue', linestyle='solid')
                ax.add_patch(wedge)

            # Plot the outer conductor using Wedge
            outer_wedge = Wedge(item['center_point'], item['radius'][1], 0, 360,
                                width=item['radius'][1] - item['radius'][0],
                                edgecolor='black', facecolor='grey', linestyle='dotted')
            ax.add_patch(outer_wedge)

            # Plot the inner conductor as a circle
            inner_circle = patches.Circle(item['center_point'], item['radius'][0], fill=True,
                                          edgecolor='black', facecolor='white', linestyle='dotted')
            ax.add_patch(inner_circle)

            # Add the center point to the list of points
            all_points.append(
                (item['center_point'][0], item['center_point'][1], item['radius'][1]))

        # Set the aspect of the plot to be equal
        # ax.set_aspect('equal', 'box')
        ax.set_aspect('equal', adjustable='box')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Geometry of the underground system.')
        plt.grid(True, linestyle='--', linewidth=0.5)

        # Determine axis limits
        margin = 0.2  # 20 % margin for better visualization
        max_x = max(point[0] + point[2] for point in all_points)
        min_x = min(point[0] - point[2] for point in all_points)
        min_y = min(point[1] - point[2] for point in all_points)

        # Adjust the y-axis to show the ground level
        min_x_level = min_x - margin * abs(min_x)
        max_x_level = max_x + margin * abs(max_x)
        min_y_level = min_y - margin * abs(min_y)
        plt.xlim(min_x_level, max_x_level)
        plt.ylim(min_y_level, 0.2)

        # Draw the ground level reference line
        # plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)

        # Fill the ground and air regions and Add labels
        ax.fill_between([min_x_level, max_x_level], min_y_level,
                        0, color='lightgrey', alpha=0.5)
        plt.text(min_x_level, 0.1, 'air ($\\varepsilon_0$, $\\mu_0$)',
                 verticalalignment='center', fontsize=10)
        plt.text(min_x_level, - 0.15, 'ground ($\\varepsilon_0$, $\\mu_0$, $\\sigma_g$)',
                 verticalalignment='center', fontsize=10)

        # plt.legend()
        plt.show()

    def wires_and_cables(self, line_type):
        """This function plots the geometry of the bifilar line."""
        _, ax = plt.subplots()
        all_points = []
        # colors = ['lightgreen', 'lightblue']

        # Plot the conductors
        for _, item in enumerate((self.mtl)):
            all_points.append((item['center_point'][0],
                               item['center_point'][1],
                               item['radius'][1]))

            # Plot the inner and outer conductors
            ax.add_patch(patches.Circle(item['center_point'], item['radius'][1],
                                        fill=True, edgecolor='black',
                                        facecolor='lightblue', linestyle='dotted'))

            ax.add_patch(patches.Circle(item['center_point'], item['radius'][0],
                                        fill=True, edgecolor='black', facecolor='white',
                                        linestyle='dotted'))

        # Determine axis limits
        margin_x = 4  # 400% margin for better visualization
        margin_y = 0.01   # 1% margin for better visualization
        max_x = max(point[0] + point[2] for point in all_points)
        min_x = min(point[0] - point[2] for point in all_points)
        max_y = max(point[1] + point[2] for point in all_points)
        min_y = min(point[1] - point[2] for point in all_points)

        # Set the aspect of the plot to be equal
        ax.set_aspect('equal', 'box')
        plt.xlabel('x')
        plt.ylabel('y')
        if line_type == 'overhead':
            plt.title(
                'Geometry of the overhead transmission line.'
                '\nDimensions are in meters.')
        if line_type == 'bifilar':
            plt.title(
                'Geometry of the bifilar transmission line.'
                '\nDimensions are in meters.')
        elif line_type == 'scc':
            plt.title(
                'Geometry of the coaxial (scc) transmission line.'
                '\nDimensions are in meters.')
        elif line_type == 'enclosure-gib':
            plt.title(
                'Cross section of three-phase enclosure-type GIB.'
                '\nDimensions are in meters.')
        plt.grid(True)
        plt.xlim(min_x - margin_x * abs(min_x), max_x + margin_x * abs(max_x))
        plt.ylim(min_y - margin_y * abs(min_y), max_y + margin_y * abs(max_y))
        plt.show()
