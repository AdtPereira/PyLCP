""" This module defines the MulticonductorTransmissionLine class.
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
from scipy.constants import mu_0, epsilon_0


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

        # Conductor surfaces dictionary
        self.surfaces = []

        # Check if the conductor is hollow or solid
        for conductor in self.mtl:
            # Check if the conductor is hollow
            if conductor['radius'][0] != 0:
                for radius in conductor['radius']:
                    self.surfaces.append({
                        'center': conductor['center_point'],
                        'fourier_order': conductor['fourier_order'],
                        'radius': radius
                    })

            # Then, the conductor is solid
            else:
                self.surfaces.append({
                    'center': conductor['center_point'],
                    'fourier_order': conductor['fourier_order'],
                    'radius': conductor['radius'][1]
                })

        # Dimension N
        # Equation (2.36) [1]
        self.N = sum([2*Np['fourier_order']+1 for Np in self.surfaces])

        # Conductors Permeability [np.array]
        self.mu = np.array([mu_0 * c['relative_permeability'] for c in self.mtl]) 

        # Conductors Permittivity [np.array]
        self.epsilon = np.array([epsilon_0 * c['relative_permittivity'] for c in self.mtl]) 

        # Conductors conductivity [np.array]
        self.sigma = np.array([c['conductivity'] for c in self.mtl])

        # Free-Space Permittivity [np.array]
        self.epsilon_out = np.array([epsilon_0 * c['relative_permittivity_out'] for c in self.mtl])

        # Distance matrices [np.array]
        self.D_pq, self.dx_pq, self.dy_pq, self.theta_pq = self.distance_matrices()

        # # Coaxial Cable Model with screen wires
        # if mtl['type'] == 'scc_screen_wires':
        #     # Include the subconductors in the list
        #     self.scc_screen_wires()

        # # Trefoil power cable Model
        # elif mtl['type'] == 'trefoil':
        #     # Core radius
        #     self.core_radius = [c['radius'][1]
        #                         for c in self.original_mtl
        #                         if c['conductor_name'] == 'core'][0]

        #     # Primary insulation thickness
        #     self.prim_ins_thick = [c['insulation']['thickness']
        #                             for c in self.original_mtl
        #                             if c['conductor_name'] == 'core'][0]

        #     # Secondary insulation thickness
        #     self.sec_ins_thick = [c['insulation']['thickness']
        #                             for c in self.original_mtl
        #                             if c['conductor_name'] == 'screen'][0]

        #     # Outer radius of the core
        #     self.scc_out_radius = (
        #         self.core_radius + self.prim_ins_thick + self.sec_ins_thick)

        #     # radius of the circumscribed circle of the trefoil triangle
        #     self.trefoil_radius = 2 * self.scc_out_radius * np.sqrt(3) / 3

        #     # Include the core subconductors in the list
        #     self.trefoil_core_subconductors()

        #     # Include the screen subconductors in the list
        #     self.trefoil_screen_subconductors()

        #     # Include the armor subconductors in the list
        #     self.trefoil_armor_subconductors()

        # # Underground [cable-hole] system Model
        # elif mtl['type'] == 'cable_hole':
        #     # Subconductors list [list]
        #     self.subconductors = [
        #         c for c in self.mtl if c['line_type'] == 'active']

        #     # Hole list [list]
        #     self.holes = [c for c in self.mtl if c['line_type'] == 'hole']

        #     # Include the subconductors and insulation (hole) in the list
        #     # print('Underground [cable-hole] system model...')

    # Contour Position Vector [np.array]
    # Equation (2.2) [1]
    def contour_vector_position(self, surfaces_list, theta, p):
        """ The boundary cp can be traced by the position vector rp(ap, θp) """

        xp, yp = surfaces_list[p]['center']
        ap = surfaces_list[p]['radius']

        return np.array([xp + ap * np.cos(theta), yp + ap * np.sin(theta)])

    # Auxiliary Vector distance [np.array]
    # Equation B.29 [1]
    # Equation B.39 [1]
    def distance_vector_dqp(self, p_surface_list, p, q_surface_list, q):
        """
        This function calculates the vector distance between the 
        center points of the conductor surfaces.

        Returns:
        numpy.ndarray: The distance matrix.
        """

        # Distance between the center points of the conductors
        x_qp = q_surface_list[q]['center'][0] - p_surface_list[p]['center'][0]
        y_qp = q_surface_list[q]['center'][1] - p_surface_list[p]['center'][1]

        # Angle between the center points of the conductors
        theta_qp = np.arctan2(y_qp, x_qp)

        # Distance between the center points of the conductors
        p = np.array(p_surface_list[p]['center'])
        q = np.array(q_surface_list[q]['center'])
        d_qp = np.linalg.norm(q - p)

        return d_qp, x_qp, y_qp, theta_qp

    # Distance matrices [np.array]
    def distance_matrices(self):
        """ Calculates the distance matrix between the center points of the conductors. """

        # Principal Dimension: Number of conductor surfaces
        N = len(self.mtl)  
        d_pq = np.zeros((N, N))
        x_pq = np.zeros_like(d_pq)
        y_pq = np.zeros_like(d_pq)
        theta_pq = np.zeros_like(d_pq)

        for p, conductor_p in enumerate(self.mtl):
            cp = conductor_p['center_point']

            for q, conductor_q in enumerate(self.mtl):
                cq = conductor_q['center_point']

                # Distance between the center points of the conductors
                x_pq[p, q] = cp[0] - cq[0]
                y_pq[p, q] = cp[1] - cq[1]

                # Angle between the center points of the conductors
                theta_pq[p, q] = np.arctan2(y_pq[p, q], x_pq[p, q])

                # Distance between the center points of the conductors
                d_pq[p, q] = np.linalg.norm(np.array(cp) - np.array(cq))

        return d_pq, x_pq, y_pq, theta_pq

    # Matrices Distance dnm and Dnm [np.array]
    # Equation 2.48 and 2.49 [2]
    def conductor_distances(self, n, m):
        """ Calculates the distance matrix Dnm between the center points of the conductors. """

        cn = self.mtl[n]['center_point']
        cm = self.mtl[m]['center_point']
        rn = self.mtl[n]['radius'][1]

        # Horizontal separation
        # Self-elements
        if n == m:
            dn_dm = rn
        else:
            dn_dm = cn[0] - cm[0]

        # Vertical separation
        hn_hm = cn[1] + cm[1]

        # Matrix Distance dnm [np.array]
        d = np.sqrt(dn_dm ** 2 + (cn[1] - cm[1]) ** 2)

        # Matrix Distance dnm [np.array]
        D = np.sqrt(dn_dm ** 2 + hn_hm ** 2)

        return d, D, hn_hm, dn_dm

    # Create a diagonal matrix by concatenating a list of matrices along the diagonal
    def create_diagonal_matrix(self, matrix_list):
        """
        Create a diagonal matrix by concatenating a list of matrices along the diagonal.

        Parameters:
        matrix_list (list): A list of numpy arrays representing matrices.

        Returns:
        numpy.ndarray: A diagonal matrix created by concatenating the matrices along the diagonal.
        """
        # Check if all elements in the list are square matrices
        for i, matrix in enumerate(matrix_list):
            if matrix.shape[0] != matrix.shape[1]:
                raise ValueError(f"Matrix at index {i} is not square: shape {matrix.shape}")

        # Determine the size of the resulting matrix Y
        size_y = sum(mat.shape[0] for mat in matrix_list)

        # Create a zero matrix of appropriate size
        matrix_y = np.zeros((size_y, size_y), dtype=complex)

        # Variables to keep track of the position in the matrix Y
        row_start = 0
        col_start = 0

        for matrix in matrix_list:
            matrix_size = matrix.shape[0]

            # Insert the matrix into the diagonal block of Y
            matrix_y[row_start:row_start + matrix_size,
                     col_start:col_start + matrix_size] = matrix

            # Update the position variables
            row_start += matrix_size
            col_start += matrix_size

        return matrix_y

    # Create a block matrix from a list of sub-matrices
    def create_block_matrix(self, matrix_list, num_rows, num_cols):
        """
        Create a block matrix from a list of sub-matrices.

        Args:
            matrix_list (list): A list of numpy arrays representing the sub-matrices.
            numRows (int): The number of rows in the block matrix.
            numColumns (int): The number of columns in the block matrix.

        Returns:
            numpy.ndarray: The block matrix created from the sub-matrices.

        Raises:
            ValueError: If the number of sub-matrices does not match numRows * numColumns.

        """

        # Verify the total number of sub-matrices matches the product of numRows and numColumns
        n = len(matrix_list)
        if num_rows * num_cols != n:
            raise ValueError(
                "The number of sub-matrices must match numRows * numColumns")

        # Determine the row and column sizes for the final matrix
        row_sizes = [matrix_list[i * num_cols].shape[0]
                     for i in range(num_rows)]
        col_sizes = [matrix_list[i].shape[1] for i in range(num_cols)]

        # Create the final matrix filled with zeros
        matrix = np.zeros((sum(row_sizes), sum(col_sizes)),
                          dtype=matrix_list[0].dtype)

        # Fill the final matrix with sub-matrices
        current_row = 0
        for i in range(num_rows):
            current_col = 0
            for j in range(num_cols):
                sub_matrix = matrix_list[i * num_cols + j]
                rows, cols = sub_matrix.shape
                matrix[current_row:current_row + rows,
                       current_col:current_col + cols] = sub_matrix
                current_col += cols
            current_row += row_sizes[i]

        return matrix

    # def wire_list_of_dict(self, conductor, wires, layer_radius, wire_radius, base_center, shift=0):
    #     """This function calculates the number of wires per layer in a coaxial cable."""
    #     base_x, base_y = base_center

    #     to_add = []
    #     # Generate the new subconductors dictionary
    #     for m in range(wires):
    #         # Calculate the angle for each subconductor
    #         angle_step = 2 * np.pi / wires

    #         # angle = m * angle_step + shift / 2 * angle_step
    #         angle = angle_step * (m + shift / 2)

    #         # Calculate the position of the subconductor
    #         cx = base_x + layer_radius * np.cos(angle)
    #         cy = base_y + layer_radius * np.sin(angle)

    #         # Generate the subconductor dictionary
    #         subconductor = {
    #             'line_id': conductor['line_id'],
    #             'conductor_name': conductor['conductor_name'],
    #             'subconductor_id': m+1,
    #             'line_type': conductor['line_type'],
    #             'line_return': conductor['line_return'],
    #             'center_point': (cx, cy),
    #             'radius': [0, wire_radius],
    #             'thickness': conductor['thickness'],
    #             'conductivity': conductor['conductivity'],
    #             'subconductors': {'status': 'no'},
    #             'relative_permeability': conductor['relative_permeability'],
    #             'relative_permittivity': conductor['relative_permittivity'],
    #             'relative_permittivity_out': conductor['relative_permittivity_out'],
    #             'surface_points': conductor['surface_points'],
    #         }

    #         # Substitute the subconductor to the multiconductor list
    #         to_add.append(subconductor)

    #     return to_add

    # def scc_screen_wires(self):
    #     """
    #     Coletar Itens para Remover e Adicionar: Em vez de modificar a lista diretamente 
    #     dentro do loop, coletamos os índices dos itens a serem removidos em to_remove e
    #     os subcondutores a serem adicionados em to_add.

    #     Remover Itens: Após a iteração, removemos os itens da lista self.conductors
    #     usando pop em ordem reversa dos índices. Isso impede problemas de índice ao
    #     remover múltiplos itens.

    #     Adicionar Itens: Após remover os itens, adicionamos os novos subconductores
    #     à lista self.conductors.

    #     Atualização do Campo subconductors: Definimos subconductors['status'] como
    #     'no' para os subconductores recém-criados, conforme a estrutura original.

    #     *** Atenção: O código foi modificado para acomodar a nova estrutura de dados. ***
    #     for idx in sorted(to_remove, reverse=True):
    #         self.conductors.pop(idx)

    #         Nesta parte do código, queremos remover múltiplos itens da lista self.conductors.
    #         Esses itens foram identificados durante a iteração anterior e seus índices foram
    #         armazenados na lista to_remove.

    #         Integridade dos Índices: Remover itens do final para o início garante que os
    #         índices dos itens restantes não sejam alterados durante o processo de remoção.
    #         Simplicidade: Evita a necessidade de ajustar manualmente os índices ou criar
    #         cópias intermediárias da lista.
    #     """

    #     to_remove = []
    #     to_add = []

    #     # Check if the conductors list contains subconductors
    #     for idx, conductor in enumerate(self.mtl):
    #         if conductor['conductor_layers'] != 0:

    #             # Remove the conductor from the list
    #             to_remove.append(idx)

    #             # Add the conductor to the subconductors list
    #             wire_radius = (
    #                 conductor['radius'][1] - conductor['radius'][0]) / 2
    #             screen_radius = conductor['radius'][1] - wire_radius
    #             wires = int(
    #                 np.floor((2 * np.pi * screen_radius) / (2 * wire_radius)))
    #             # wires = 2

    #             # Generate the new subconductors dictionary
    #             subconductor_list = self.wire_list_of_dict(
    #                 conductor, wires, screen_radius, wire_radius, base_center=(0, 0))

    #             # Substitute the subconductor to the multiconductor list
    #             to_add.extend(subconductor_list)

    #     # Remove the conductors from the list
    #     for idx in sorted(to_remove, reverse=True):
    #         self.mtl.pop(idx)

    #     # Add the subconductors to the list
    #     self.mtl.extend(to_add)

    # def trefoil_core_subconductors(self):
    #     """
    #     Coletar Itens para Remover e Adicionar: Em vez de modificar a lista diretamente 
    #     dentro do loop, coletamos os índices dos itens a serem removidos em to_remove e
    #     os subcondutores a serem adicionados em to_add.

    #     Remover Itens: Após a iteração, removemos os itens da lista self.conductors
    #     usando pop em ordem reversa dos índices. Isso impede problemas de índice ao
    #     remover múltiplos itens.

    #     Adicionar Itens: Após remover os itens, adicionamos os novos subconductores
    #     à lista self.conductors.

    #     Atualização do Campo subconductors: Definimos subconductors['status'] como
    #     'no' para os subconductores recém-criados, conforme a estrutura original.

    #     *** Atenção: O código foi modificado para acomodar a nova estrutura de dados. ***
    #     for idx in sorted(to_remove, reverse=True):
    #         self.conductors.pop(idx)

    #         Nesta parte do código, queremos remover múltiplos itens da lista self.conductors.
    #         Esses itens foram identificados durante a iteração anterior e seus índices foram
    #         armazenados na lista to_remove.

    #         Integridade dos Índices: Remover itens do final para o início garante que os
    #         índices dos itens restantes não sejam alterados durante o processo de remoção.
    #         Simplicidade: Evita a necessidade de ajustar manualmente os índices ou criar
    #         cópias intermediárias da lista.
    #     """

    #     # Include the core subconductors in the list
    #     to_remove = []
    #     to_add = []

    #     for idx, conductor in enumerate(self.mtl):
    #         if conductor['conductor_name'] == 'core':

    #             # Remove the conductor from the list
    #             to_remove.append(idx)

    #             # Number of core subconductors in a trefoil power cable
    #             num_subconductors = 3

    #             # Generate the core subconductors dictionary
    #             subconductor_list = self.wire_list_of_dict(
    #                 conductor, num_subconductors, self.trefoil_radius,
    #                 conductor['radius'][1], base_center=(0, 0))

    #             # Substitute the subconductor to the multiconductor list
    #             to_add.extend(subconductor_list)

    #     # Remove the conductors from the list
    #     for idx in sorted(to_remove, reverse=True):
    #         self.mtl.pop(idx)

    #     # Add the subconductors to the list
    #     self.mtl.extend(to_add)

    # def trefoil_screen_subconductors(self):
    #     """
    #     Coletar Itens para Remover e Adicionar: Em vez de modificar a lista diretamente 
    #     dentro do loop, coletamos os índices dos itens a serem removidos em to_remove e
    #     os subcondutores a serem adicionados em to_add.

    #     Remover Itens: Após a iteração, removemos os itens da lista self.conductors
    #     usando pop em ordem reversa dos índices. Isso impede problemas de índice ao
    #     remover múltiplos itens.

    #     Adicionar Itens: Após remover os itens, adicionamos os novos subconductores
    #     à lista self.conductors.

    #     Atualização do Campo subconductors: Definimos subconductors['status'] como
    #     'no' para os subconductores recém-criados, conforme a estrutura original.

    #     *** Atenção: O código foi modificado para acomodar a nova estrutura de dados. ***
    #     for idx in sorted(to_remove, reverse=True):
    #         self.conductors.pop(idx)

    #         Nesta parte do código, queremos remover múltiplos itens da lista self.conductors.
    #         Esses itens foram identificados durante a iteração anterior e seus índices foram
    #         armazenados na lista to_remove.

    #         Integridade dos Índices: Remover itens do final para o início garante que os
    #         índices dos itens restantes não sejam alterados durante o processo de remoção.
    #         Simplicidade: Evita a necessidade de ajustar manualmente os índices ou criar
    #         cópias intermediárias da lista.
    #     """

    #     # Include the core subconductors in the list
    #     to_remove = []
    #     to_add = []

    #     # Core subconductors list
    #     core_base_center_list = [conductor['center_point']
    #                              for conductor in self.mtl
    #                              if conductor['conductor_name'] == 'core']

    #     for idx, conductor in enumerate(self.mtl):
    #         if conductor['conductor_name'] == 'screen':

    #             # Remove the conductor from the list
    #             to_remove.append(idx)

    #             # Number of core subconductors in a trefoil power cable
    #             num_wires = conductor['subconductors']['wires_per_layer']

    #             # screen radius
    #             wire_radius = conductor['subconductors']['wire_radius']
    #             screen_radius = self.core_radius + self.prim_ins_thick + wire_radius

    #             # Generate the screen subconductors dictionary
    #             for _, base_center in enumerate(core_base_center_list):
    #                 wire_list = self.wire_list_of_dict(
    #                     conductor, num_wires, screen_radius,
    #                     wire_radius, base_center)

    #                 # Substitute the subconductor to the multiconductor list
    #                 to_add.extend(wire_list)

    #     # Remove the conductors from the list
    #     for idx in sorted(to_remove, reverse=True):
    #         self.mtl.pop(idx)

    #     # Add the subconductors to the list
    #     self.mtl.extend(to_add)

    # def trefoil_armor_subconductors(self):
    #     """
    #     Coletar Itens para Remover e Adicionar: Em vez de modificar a lista diretamente 
    #     dentro do loop, coletamos os índices dos itens a serem removidos em to_remove e
    #     os subcondutores a serem adicionados em to_add.

    #     Remover Itens: Após a iteração, removemos os itens da lista self.conductors
    #     usando pop em ordem reversa dos índices. Isso impede problemas de índice ao
    #     remover múltiplos itens.

    #     Adicionar Itens: Após remover os itens, adicionamos os novos subconductores
    #     à lista self.conductors.

    #     Atualização do Campo subconductors: Definimos subconductors['status'] como
    #     'no' para os subconductores recém-criados, conforme a estrutura original.

    #     *** Atenção: O código foi modificado para acomodar a nova estrutura de dados. ***
    #     for idx in sorted(to_remove, reverse=True):
    #         self.conductors.pop(idx)

    #         Nesta parte do código, queremos remover múltiplos itens da lista self.conductors.
    #         Esses itens foram identificados durante a iteração anterior e seus índices foram
    #         armazenados na lista to_remove.

    #         Integridade dos Índices: Remover itens do final para o início garante que os
    #         índices dos itens restantes não sejam alterados durante o processo de remoção.
    #         Simplicidade: Evita a necessidade de ajustar manualmente os índices ou criar
    #         cópias intermediárias da lista.
    #     """

    #     # Include the core subconductors in the list
    #     to_remove = []
    #     to_add = []

    #     for idx, conductor in enumerate(self.mtl):
    #         if conductor['conductor_name'] == 'armor':

    #             # Remove the conductor from the list
    #             to_remove.append(idx)

    #             # Number of wire subconductors and layers in a trefoil power cable
    #             num_wires = conductor['subconductors']['wires_per_layer']
    #             num_layers = conductor['subconductors']['layers']

    #             # armor radius
    #             wire_radius = conductor['subconductors']['wire_radius']

    #             # Cable radius
    #             cable_radius = conductor['radius'][1]

    #             # Generate the armor subconductors dictionary
    #             for k in range(num_layers):
    #                 # Calculate the radius of the armor layer
    #                 armor_radius = cable_radius - (2*k + 1) * wire_radius

    #                 # Generate the armor subconductors dictionary
    #                 subconductor_list = self.wire_list_of_dict(
    #                     conductor, num_wires, armor_radius,
    #                     wire_radius, base_center=(0, 0), shift=k)

    #                 # Substitute the subconductor to the multiconductor list
    #                 to_add.extend(subconductor_list)

    #     # Remove the conductors from the list
    #     for idx in sorted(to_remove, reverse=True):
    #         self.mtl.pop(idx)

    #     # Add the subconductors to the list
    #     self.mtl.extend(to_add)

