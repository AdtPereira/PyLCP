"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

Below is a high-level overview of the script components:

Imports and Global Variables:

Required libraries and global variables are imported and defined.

BIFILAR_TL: A list containing properties of the two conductors.
Classes:

Geometry: Handles basic geometry calculations, such as distance matrices between conductor centers.
ParametersWithFrequency: Extends Geometry to include frequency-dependent parameters.
GreensMatrices: Uses the geometry to compute Green's matrices, which are essential for the method 
of moments.
MoMSuperficialOperator: Implements the method of moments, calculating matrices like U, Ys, G, and Z.
AnalyticalFormulation: Provides analytical formulations for high-frequency resistance, external 
inductance, and impedance.
Plotter: Handles plotting of series resistance and inductance against frequency.

Functions:

clear_screen: Clears the console screen.
main: The main function orchestrates the scattering calculations and plotting. It performs the 
following steps:
Clears the screen.
Initializes objects for the method of moments and analytical formulations.
Computes series resistance, external inductance, and impedance over a range of frequencies.
Plots the results using the Plotter class.
Detailed Class and Function Explanations
Geometry
__init__: Initializes the geometry of the system based on conductor properties.
distance_matrices: Calculates matrices for distances and angles between conductor centers.
ParametersWithFrequency
__init__: Extends the Geometry class to include frequency-dependent parameters such as 
conductivity, permeability, and permittivity.
ynp_operator: Calculates the surface admittance operator for a conductor.
GreensMatrices
__init__: Initializes Green's matrices using the geometry of the system.
dissertation and ieee_paper: Calculate Green's functions using different methods.
sub_matrices: Generates Green's sub-matrices.
MoMSuperficialOperator
__init__: Extends ParametersWithFrequency to initialize the method of moments parameters.
matrix_u: Constructs matrix U.
matrix_ys: Constructs matrix Ys.
matrix_g and matrix_g_ieee: Constructs matrix G using different methods.
matrix_z: Computes the impedance matrix Z.
AnalyticalFormulation
__init__: Initializes the analytical formulation based on the conductor properties and frequency.
pul_parameters: Calculates high-frequency resistance, external inductance, and impedance.
Plotter
__init__: Initializes the plotting class with frequency and impedance data.
series_resistance: Plots series resistance against frequency.
series_inductance: Plots series inductance against frequency.

REFERENCES:
[1] PATEL, Utkarsh R. A Surface Admittance Approach For Fast Calculation of the 
    Series Impedance of Cables Including Skin, Proximity, and Ground Return Effects.
    2014. University of Toronto, Graduate Department of The Edward S. Rogers Sr. 
    Department of Electrical & Computer Engineering. 

[2] U. R. Patel, B. Gustavsen and P. Triverio, "An Equivalent Surface Current Approach
    for the Computation of the Series Impedance of Power Cables with Inclusion of Skin
    and Proximity Effects," in IEEE Transactions on Power Delivery, vol. 28, no. 4, pp.
    2474-2482, Oct. 2013, doi: 10.1109/TPWRD.2013.2267098.

[3] U. R. Patel, B. Gustavsen and P. Triverio, "Application of the MoM-SO Method for 
    Accurate Impedance Calculation of Single-Core Cables Enclosed by a Conducting Pipe," 
    Proc. International Conference on Power Systems Transients (IPST 2013), Vancouver, 
    Canada July 18-20, 2013. https://www.ipstconf.org/Proc_IPST2013.php

"""

import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib import patches


class MulticonductorTransmissionLine():
    """ This class defines the coaxial cable with a sheath. """

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
        self.original_mtl = copy.deepcopy(mtl['data'])
        self.mtl = copy.deepcopy(mtl['data'])

        # Coaxial Cable Model with screen wires
        if mtl['type'] == 'scc_screen_wires':
            # Include the subconductors in the list
            self.scc_screen_wires()

        # Trefoil power cable Model
        elif mtl['type'] == 'trefoil':
            # Core radius
            self.core_radius = [c['radius'][1]
                                for c in self.original_mtl
                                if c['conductor_name'] == 'core'][0]

            # Primary insulation thickness
            self.prim_ins_thick = [c['insulation']['thickness']
                                   for c in self.original_mtl
                                   if c['conductor_name'] == 'core'][0]

            # Secondary insulation thickness
            self.sec_ins_thick = [c['insulation']['thickness']
                                  for c in self.original_mtl
                                  if c['conductor_name'] == 'screen'][0]

            # Outer radius of the core
            self.scc_out_radius = (
                self.core_radius + self.prim_ins_thick + self.sec_ins_thick)

            # radius of the circumscribed circle of the trefoil triangle
            self.trefoil_radius = 2 * self.scc_out_radius * np.sqrt(3) / 3

            # Include the core subconductors in the list
            self.trefoil_core_subconductors()

            # Include the screen subconductors in the list
            self.trefoil_screen_subconductors()

            # Include the armor subconductors in the list
            self.trefoil_armor_subconductors()

        # Underground [cable-hole] system Model
        elif mtl['type'] == 'cable_hole':
            # Subconductors list [list]
            self.subconductors = [
                c for c in self.mtl if c['line_type'] == 'active']

            # Hole list [list]
            self.holes = [c for c in self.mtl if c['line_type'] == 'hole']

            # Include the subconductors and insulation (hole) in the list
            # print('Underground [cable-hole] system model...')

    def wire_list_of_dict(
            self, conductor, wires, layer_radius, wire_radius, base_center, shift=0):
        """This function calculates the number of wires per layer in a coaxial cable."""
        base_x, base_y = base_center

        to_add = []
        # Generate the new subconductors dictionary
        for m in range(wires):
            # Calculate the angle for each subconductor
            angle_step = 2 * np.pi / wires

            # angle = m * angle_step + shift / 2 * angle_step
            angle = angle_step * (m + shift / 2)

            # Calculate the position of the subconductor
            cx = base_x + layer_radius * np.cos(angle)
            cy = base_y + layer_radius * np.sin(angle)

            # Generate the subconductor dictionary
            subconductor = {
                'line_id': conductor['line_id'],
                'conductor_name': conductor['conductor_name'],
                'subconductor_id': m+1,
                'line_type': conductor['line_type'],
                'line_return': conductor['line_return'],
                'center_point': (cx, cy),
                'radius': [0, wire_radius],
                'thickness': conductor['thickness'],
                'conductivity': conductor['conductivity'],
                'subconductors': {'status': 'no'},
                'relative_permeability': conductor['relative_permeability'],
                'relative_permittivity': conductor['relative_permittivity'],
                'relative_permittivity_out': conductor['relative_permittivity_out'],
                'surface_points': conductor['surface_points'],
            }

            # Substitute the subconductor to the multiconductor list
            to_add.append(subconductor)

        return to_add

    def scc_screen_wires(self):
        """
        Coletar Itens para Remover e Adicionar: Em vez de modificar a lista diretamente 
        dentro do loop, coletamos os índices dos itens a serem removidos em to_remove e
        os subcondutores a serem adicionados em to_add.

        Remover Itens: Após a iteração, removemos os itens da lista self.conductors
        usando pop em ordem reversa dos índices. Isso impede problemas de índice ao
        remover múltiplos itens.

        Adicionar Itens: Após remover os itens, adicionamos os novos subconductores
        à lista self.conductors.

        Atualização do Campo subconductors: Definimos subconductors['status'] como
        'no' para os subconductores recém-criados, conforme a estrutura original.

        *** Atenção: O código foi modificado para acomodar a nova estrutura de dados. ***
        for idx in sorted(to_remove, reverse=True):
            self.conductors.pop(idx)

            Nesta parte do código, queremos remover múltiplos itens da lista self.conductors.
            Esses itens foram identificados durante a iteração anterior e seus índices foram
            armazenados na lista to_remove.

            Integridade dos Índices: Remover itens do final para o início garante que os
            índices dos itens restantes não sejam alterados durante o processo de remoção.
            Simplicidade: Evita a necessidade de ajustar manualmente os índices ou criar
            cópias intermediárias da lista.
        """

        to_remove = []
        to_add = []

        # Check if the conductors list contains subconductors
        for idx, conductor in enumerate(self.mtl):
            if conductor['conductor_layers'] != 0:

                # Remove the conductor from the list
                to_remove.append(idx)

                # Add the conductor to the subconductors list
                wire_radius = (
                    conductor['radius'][1] - conductor['radius'][0]) / 2
                screen_radius = conductor['radius'][1] - wire_radius
                wires = int(
                    np.floor((2 * np.pi * screen_radius) / (2 * wire_radius)))
                # wires = 2

                # Generate the new subconductors dictionary
                subconductor_list = self.wire_list_of_dict(
                    conductor, wires, screen_radius, wire_radius, base_center=(0, 0))

                # Substitute the subconductor to the multiconductor list
                to_add.extend(subconductor_list)

        # Remove the conductors from the list
        for idx in sorted(to_remove, reverse=True):
            self.mtl.pop(idx)

        # Add the subconductors to the list
        self.mtl.extend(to_add)

    def trefoil_core_subconductors(self):
        """
        Coletar Itens para Remover e Adicionar: Em vez de modificar a lista diretamente 
        dentro do loop, coletamos os índices dos itens a serem removidos em to_remove e
        os subcondutores a serem adicionados em to_add.

        Remover Itens: Após a iteração, removemos os itens da lista self.conductors
        usando pop em ordem reversa dos índices. Isso impede problemas de índice ao
        remover múltiplos itens.

        Adicionar Itens: Após remover os itens, adicionamos os novos subconductores
        à lista self.conductors.

        Atualização do Campo subconductors: Definimos subconductors['status'] como
        'no' para os subconductores recém-criados, conforme a estrutura original.

        *** Atenção: O código foi modificado para acomodar a nova estrutura de dados. ***
        for idx in sorted(to_remove, reverse=True):
            self.conductors.pop(idx)

            Nesta parte do código, queremos remover múltiplos itens da lista self.conductors.
            Esses itens foram identificados durante a iteração anterior e seus índices foram
            armazenados na lista to_remove.

            Integridade dos Índices: Remover itens do final para o início garante que os
            índices dos itens restantes não sejam alterados durante o processo de remoção.
            Simplicidade: Evita a necessidade de ajustar manualmente os índices ou criar
            cópias intermediárias da lista.
        """

        # Include the core subconductors in the list
        to_remove = []
        to_add = []

        for idx, conductor in enumerate(self.mtl):
            if conductor['conductor_name'] == 'core':

                # Remove the conductor from the list
                to_remove.append(idx)

                # Number of core subconductors in a trefoil power cable
                num_subconductors = 3

                # Generate the core subconductors dictionary
                subconductor_list = self.wire_list_of_dict(
                    conductor, num_subconductors, self.trefoil_radius,
                    conductor['radius'][1], base_center=(0, 0))

                # Substitute the subconductor to the multiconductor list
                to_add.extend(subconductor_list)

        # Remove the conductors from the list
        for idx in sorted(to_remove, reverse=True):
            self.mtl.pop(idx)

        # Add the subconductors to the list
        self.mtl.extend(to_add)

    def trefoil_screen_subconductors(self):
        """
        Coletar Itens para Remover e Adicionar: Em vez de modificar a lista diretamente 
        dentro do loop, coletamos os índices dos itens a serem removidos em to_remove e
        os subcondutores a serem adicionados em to_add.

        Remover Itens: Após a iteração, removemos os itens da lista self.conductors
        usando pop em ordem reversa dos índices. Isso impede problemas de índice ao
        remover múltiplos itens.

        Adicionar Itens: Após remover os itens, adicionamos os novos subconductores
        à lista self.conductors.

        Atualização do Campo subconductors: Definimos subconductors['status'] como
        'no' para os subconductores recém-criados, conforme a estrutura original.

        *** Atenção: O código foi modificado para acomodar a nova estrutura de dados. ***
        for idx in sorted(to_remove, reverse=True):
            self.conductors.pop(idx)

            Nesta parte do código, queremos remover múltiplos itens da lista self.conductors.
            Esses itens foram identificados durante a iteração anterior e seus índices foram
            armazenados na lista to_remove.

            Integridade dos Índices: Remover itens do final para o início garante que os
            índices dos itens restantes não sejam alterados durante o processo de remoção.
            Simplicidade: Evita a necessidade de ajustar manualmente os índices ou criar
            cópias intermediárias da lista.
        """

        # Include the core subconductors in the list
        to_remove = []
        to_add = []

        # Core subconductors list
        core_base_center_list = [conductor['center_point']
                                 for conductor in self.mtl
                                 if conductor['conductor_name'] == 'core']

        for idx, conductor in enumerate(self.mtl):
            if conductor['conductor_name'] == 'screen':

                # Remove the conductor from the list
                to_remove.append(idx)

                # Number of core subconductors in a trefoil power cable
                num_wires = conductor['subconductors']['wires_per_layer']

                # screen radius
                wire_radius = conductor['subconductors']['wire_radius']
                screen_radius = self.core_radius + self.prim_ins_thick + wire_radius

                # Generate the screen subconductors dictionary
                for _, base_center in enumerate(core_base_center_list):
                    wire_list = self.wire_list_of_dict(
                        conductor, num_wires, screen_radius,
                        wire_radius, base_center)

                    # Substitute the subconductor to the multiconductor list
                    to_add.extend(wire_list)

        # Remove the conductors from the list
        for idx in sorted(to_remove, reverse=True):
            self.mtl.pop(idx)

        # Add the subconductors to the list
        self.mtl.extend(to_add)

    def trefoil_armor_subconductors(self):
        """
        Coletar Itens para Remover e Adicionar: Em vez de modificar a lista diretamente 
        dentro do loop, coletamos os índices dos itens a serem removidos em to_remove e
        os subcondutores a serem adicionados em to_add.

        Remover Itens: Após a iteração, removemos os itens da lista self.conductors
        usando pop em ordem reversa dos índices. Isso impede problemas de índice ao
        remover múltiplos itens.

        Adicionar Itens: Após remover os itens, adicionamos os novos subconductores
        à lista self.conductors.

        Atualização do Campo subconductors: Definimos subconductors['status'] como
        'no' para os subconductores recém-criados, conforme a estrutura original.

        *** Atenção: O código foi modificado para acomodar a nova estrutura de dados. ***
        for idx in sorted(to_remove, reverse=True):
            self.conductors.pop(idx)

            Nesta parte do código, queremos remover múltiplos itens da lista self.conductors.
            Esses itens foram identificados durante a iteração anterior e seus índices foram
            armazenados na lista to_remove.

            Integridade dos Índices: Remover itens do final para o início garante que os
            índices dos itens restantes não sejam alterados durante o processo de remoção.
            Simplicidade: Evita a necessidade de ajustar manualmente os índices ou criar
            cópias intermediárias da lista.
        """

        # Include the core subconductors in the list
        to_remove = []
        to_add = []

        for idx, conductor in enumerate(self.mtl):
            if conductor['conductor_name'] == 'armor':

                # Remove the conductor from the list
                to_remove.append(idx)

                # Number of wire subconductors and layers in a trefoil power cable
                num_wires = conductor['subconductors']['wires_per_layer']
                num_layers = conductor['subconductors']['layers']

                # armor radius
                wire_radius = conductor['subconductors']['wire_radius']

                # Cable radius
                cable_radius = conductor['radius'][1]

                # Generate the armor subconductors dictionary
                for k in range(num_layers):
                    # Calculate the radius of the armor layer
                    armor_radius = cable_radius - (2*k + 1) * wire_radius

                    # Generate the armor subconductors dictionary
                    subconductor_list = self.wire_list_of_dict(
                        conductor, num_wires, armor_radius,
                        wire_radius, base_center=(0, 0), shift=k)

                    # Substitute the subconductor to the multiconductor list
                    to_add.extend(subconductor_list)

        # Remove the conductors from the list
        for idx in sorted(to_remove, reverse=True):
            self.mtl.pop(idx)

        # Add the subconductors to the list
        self.mtl.extend(to_add)


class GraphicRepresentation(MulticonductorTransmissionLine):
    """ This class defines the coaxial cable with a sheath. """

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
            all_points.append((item['center_point'][0], item['center_point'][1],
                               item['radius'][1]))

            # Plot the inner and outer conductors
            ax.add_patch(patches.Circle(item['center_point'], item['radius'][1],
                                        fill=True, edgecolor='black',
                                        facecolor='lightblue', linestyle='dotted'))

            ax.add_patch(patches.Circle(item['center_point'], item['radius'][0],
                                        fill=True, edgecolor='black', facecolor='white',
                                        linestyle='dotted'))

        # Determine axis limits
        margin = 0.05  # 5% margin for better visualization
        max_x = max(point[0] + point[2] for point in all_points)
        min_x = min(point[0] - point[2] for point in all_points)
        max_y = max(point[1] + point[2] for point in all_points)
        min_y = min(point[1] - point[2] for point in all_points)

        # Set the aspect of the plot to be equal
        ax.set_aspect('equal', 'box')
        plt.xlabel('x')
        plt.ylabel('y')
        if line_type == 'bifilar':
            plt.title('Geometry of the bifilar transmission line.')
        elif line_type == 'scc':
            plt.title('Geometry of the coaxial (scc) transmission line.')
        elif line_type == 'enclosure-gib':
            plt.title('Cross section of three-phase enclosure-type GIB')
        plt.grid(True)
        plt.xlim(min_x - margin * abs(min_x), max_x + margin * abs(max_x))
        plt.ylim(min_y - margin * abs(min_y), max_y + margin * abs(max_y))
        plt.show()

    def plot_trefoil(self):
        """This function extracts the conductors from the list."""
        _, ax = plt.subplots()
        all_points = []

        # Core subconductors list
        core_sublist = [conductor
                        for conductor in self.mtl
                        if conductor['conductor_name'] == 'core']

        # Screen subconductors list
        screen_sublist = [conductor
                          for conductor in self.mtl
                          if conductor['conductor_name'] == 'screen']

        # Screen subconductors list
        armor_sublist = [conductor
                         for conductor in self.mtl
                         if conductor['conductor_name'] == 'armor']

        # Jacket dimensions
        jack_thickness = [conductor['insulation']['thickness']
                          for conductor in self.original_mtl
                          if conductor['conductor_name'] == 'armor'][0]

        # Armor outer radius
        arm_out_radius = [conductor['radius'][1]
                          for conductor in self.original_mtl
                          if conductor['conductor_name'] == 'armor'][0]

        # Create the core and screen plots
        for i, core in enumerate((core_sublist)):
            # core conductor
            cx, cy = core['center_point']
            ax.add_patch(patches.Circle(core['center_point'], self.core_radius,
                                        fill=True, edgecolor='black', facecolor='darkgrey'))

            # Add the core number
            ax.text(cx, cy, f'Core {i+1}', ha='center',
                    va='center', fontsize=10)

            # Primary insulation
            prim_ins_radius = self.core_radius + self.prim_ins_thick
            ax.add_patch(patches.Circle((cx, cy), prim_ins_radius,
                                        fill=False, edgecolor='blue', linestyle='dashed'))

            # Screen wires
            for i, screen in enumerate((screen_sublist)):
                # core conductor
                ax.add_patch(patches.Circle(screen['center_point'], screen['radius'][1],
                                            fill=True, edgecolor='gray', facecolor='lightblue'))

                # Add the core number
                sc_x, sc_y = screen['center_point']
                ax.text(sc_x, sc_y, f'{i+1}',
                        ha='center', va='center', fontsize=6)

            # Secondary insulation
            ax.add_patch(patches.Circle((cx, cy), self.scc_out_radius,
                                        fill=False, edgecolor='black', linestyle='solid'))

        # Create the armor plots
        for i, armor in enumerate((armor_sublist)):
            cx, cy = armor['center_point']
            ax.add_patch(patches.Circle(armor['center_point'], armor['radius'][1],
                                        fill=True, edgecolor='gray', facecolor='lightgreen'))

            # Add the armor number
            ax.text(cx, cy, f'{i+1}', ha='center', va='center', fontsize=8)

        # Jacket and armor dimensions
        jack_inn_radius = self.trefoil_radius + self.scc_out_radius
        jack_out_radius = jack_inn_radius + jack_thickness
        ax.add_patch(patches.Circle((0, 0), jack_inn_radius,
                                    fill=False, edgecolor='black', linestyle='solid'))
        ax.add_patch(patches.Circle((0, 0), jack_out_radius,
                                    fill=False, edgecolor='blue', linestyle='dashed'))
        ax.add_patch(patches.Circle((0, 0), arm_out_radius,
                                    fill=False, edgecolor='black', linestyle='solid'))

        # Include the relevant points for the axis limits
        all_points.append((0, 0, arm_out_radius))

        # Determine axis limits
        margin = 0.05  # 5% margin for better visualization
        max_x = max(p[0] + p[2] for p in all_points)
        min_x = min(p[0] - p[2] for p in all_points)
        max_y = max(p[1] + p[2] for p in all_points)
        min_y = min(p[1] - p[2] for p in all_points)

        # Set the aspect of the plot to be equal
        ax.set_aspect('equal', 'box')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Trefoil Power Cable Geometry.')
        plt.grid(True)
        plt.xlim(min_x - margin * abs(min_x), max_x + margin * abs(max_x))
        plt.ylim(min_y - margin * abs(min_y), max_y + margin * abs(max_y))
        plt.show()

    def plot_trefoil_from_original_data(self, conductors_list):
        """This function plots the geometry of the trefoil power cable."""
        _, ax = plt.subplots()

        def rotate_point(
                point, center_point, angle):
            """This function rotates a point around another point by a given angle."""
            angle_rad = np.radians(angle)

            x, y = point
            x0, y0 = center_point

            # Rotation x-value
            x_rot = x0 + (x - x0) * np.cos(angle_rad) - \
                (y - y0) * np.sin(angle_rad)

            # Rotation y-value
            y_rot = y0 + (x - x0) * np.sin(angle_rad) + \
                (y - y0) * np.cos(angle_rad)

            return (x_rot, y_rot)

        def plot_wires_subconductors(
                wires, layer_radius, wire_radius, base_center, fig, shift=0):
            """This function plots the geometry of the subconductors."""
            base_x, base_y = base_center

            for m in range(wires):
                # Calculate the angle for each subconductor
                angle_step = 2 * np.pi / wires

                # angle = m * angle_step + shift / 2 * angle_step
                angle = angle_step * (m + shift / 2)

                # Calculate the position of the subconductor
                cx = base_x + layer_radius * np.cos(angle)
                cy = base_y + layer_radius * np.sin(angle)

                # Plot the subconductor
                fig.add_patch(patches.Circle((cx, cy), wire_radius, fill=True,
                                             edgecolor='gray', facecolor='lightblue'))

                # Add the subconductor number
                fig.text(cx, cy, f'{m+1}', ha='center',
                         va='center', fontsize=6)

        mtl = conductors_list['data']
        all_points = []
        core = [c for c in mtl if c['conductor_name'] == 'core'][0]
        screen = [c for c in mtl if c['conductor_name'] == 'screen'][0]
        armor = [c for c in mtl if c['conductor_name'] == 'armor'][0]

        # primary insulation
        prim_insul_radius = core['radius'][1] + core['insulation']['thickness']

        # screen wires
        screen_radius = prim_insul_radius + \
            screen['subconductors']['wire_radius']

        # secondary insulation
        sec_insul_radius = prim_insul_radius + self.sec_ins_thick

        # Jacket dimensions
        jacket_inner_radius = self.trefoil_radius + self.scc_out_radius
        jacket_outer_radius = jacket_inner_radius + \
            armor['insulation']['thickness']

        # Define the core positions for the trefoil power cable
        core_positions = [(self.trefoil_radius, 0),
                          rotate_point((self.trefoil_radius, 0), (0, 0), 120),
                          rotate_point((self.trefoil_radius, 0), (0, 0), 240)]

        # Plot the core and screen subconductors
        for i, (cx, cy) in enumerate(core_positions):
            # core conductor
            ax.add_patch(patches.Circle((cx, cy), core['radius'][1],
                                        fill=True, edgecolor='black', facecolor='darkgrey'))

            ax.text(cx, cy, f'Core {i+1}', ha='center',
                    va='center', fontsize=10)

            # primary insulation
            ax.add_patch(patches.Circle((cx, cy), prim_insul_radius,
                                        fill=False, edgecolor='blue', linestyle='dashed'))

            # screen wires
            plot_wires_subconductors(screen['subconductors']['wires_per_layer'],
                                     screen_radius,
                                     screen['subconductors']['wire_radius'],
                                     (cx, cy), ax, shift=0)

            # secondary insulation
            ax.add_patch(patches.Circle((cx, cy), sec_insul_radius, fill=False,
                                        edgecolor='black', linestyle='solid'))

        # Jacket layer
        ax.add_patch(patches.Circle((0, 0), jacket_inner_radius,
                                    fill=False, edgecolor='black', linestyle='solid'))
        ax.add_patch(patches.Circle((0, 0), jacket_outer_radius,
                                    fill=False, edgecolor='blue', linestyle='dashed'))

        # Armoured layer
        ax.add_patch(patches.Circle((0, 0), armor['radius'][1],
                                    fill=False, edgecolor='black', linestyle='solid'))

        # Armor wires
        for k in range(armor['subconductors']['layers']):
            armor_radius = armor['radius'][1] - \
                (2*k + 1) * armor['subconductors']['wire_radius']

            plot_wires_subconductors(armor['subconductors']['wires_per_layer'],
                                     armor_radius,
                                     armor['subconductors']['wire_radius'],
                                     (0, 0), ax, shift=k)

        # Include the relevant points for the axis limits
        all_points.append((0, 0, armor['radius'][1]))
        margin = 0.05  # 5% margin for better visualization
        max_x = max(p[0] + p[2] for p in all_points)
        min_x = min(p[0] - p[2] for p in all_points)
        max_y = max(p[1] + p[2] for p in all_points)
        min_y = min(p[1] - p[2] for p in all_points)

        # Set the aspect of the plot to be equal
        ax.set_aspect('equal', 'box')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Trefoil Power Cable Geometry from Original Data')
        plt.grid(True)
        plt.xlim(min_x - margin * abs(min_x), max_x + margin * abs(max_x))
        plt.ylim(min_y - margin * abs(min_y), max_y + margin * abs(max_y))
        plt.show()
