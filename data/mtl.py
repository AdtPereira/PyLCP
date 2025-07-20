"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

REFERENCES:
[1] 

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

        # Conductors Permeability [np.array]
        self.mu = np.array([mu_0 * c['relative_permeability']for c in self.mtl]) 

        # Conductors Permittivity [np.array]
        self.epsilon = np.array([epsilon_0 * c['relative_permittivity'] for c in self.mtl]) 

        # Conductors conductivity [np.array]
        self.sigma = np.array([c['conductivity'] for c in self.mtl])

        # Free-Space Permittivity [np.array]
        self.epsilon_out = np.array([epsilon_0 * c['relative_permittivity_out'] for c in self.mtl])

        # Distance matrices [np.array]
        self.D_pq, self.dx_pq, self.dy_pq, self.theta_pq = self.distance_matrices()

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

        vector_distance = {'norm': d_qp,
                           'angle': theta_qp,
                           'u_x': x_qp,
                           'u_y': y_qp}

        return vector_distance

    # Distance matrices [np.array]
    def distance_matrices(self):
        """ Calculates the distance matrix between the center points of the conductors. """

        # Principal Dimension: Number of conductor surfaces
        N = len(self.mtl)  # pylint: disable=invalid-name
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
        D = np.sqrt(dn_dm ** 2 + hn_hm ** 2)  # pylint: disable=invalid-name

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
