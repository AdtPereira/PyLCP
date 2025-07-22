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
"""

import numpy as np
from scipy.constants import epsilon_0
from .multiconductor import MulticonductorTransmissionLine as MTL


class FreeSpace(MTL):
    """ This class contains the basic geometry of the system. """

    def __init__(self, mtl):
        super().__init__(mtl)

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


class UndergroundSystem(MTL):
    """ This class contains the basic geometry of the system. """

    def __init__(self, mtl):
        super().__init__(mtl)

        # Conductor surfaces dictionary
        self.conductor_surfaces = []

        # Hole surfaces dictionary
        self.hole_surfaces = []

        # Ground dictionary
        self.ground = []

        # Number of conductor and hole surfaces
        self.surfaces_type = []

        # Verify the items in the dictionary
        for item in self.mtl:
            # Check if the item is conductor
            if item['line_type'] == 'active':
                self.surfaces_type.append(1)
                # Check if the conductor is hollow
                if item['radius'][0] != 0:
                    for radius in item['radius']:
                        self.conductor_surfaces.append({
                            'center': item['center_point'],
                            'fourier_order': item['surface_points'],
                            'radius': radius
                        })

                # Then, the conductor is solid
                else:
                    self.conductor_surfaces.append({
                        'center': item['center_point'],
                        'fourier_order': item['surface_points'],
                        'radius': item['radius'][1]
                    })

            # Check if the item is a hole
            elif item['line_type'] == 'hole':
                self.surfaces_type.append(0)
                self.hole_surfaces.append({
                    'center': item['center_point'],
                    'fourier_order': item['surface_points'],
                    'radius': item['radius'][1]
                })

            # Check if the item is ground
            elif item['line_type'] == 'return':
                self.ground.append({
                    'conductivity': item['conductivity'],
                    'permittivity': item['relative_permittivity'] * epsilon_0
                })

        # Dimension N
        # Equation (2.36) - PAG. 32 [1]
        self.N = sum([2*Np['fourier_order']+1 for Np in self.conductor_surfaces])

        # Dimension N_hat
        # Equation (3.3) - PAG. 59 [1]
        self.Nhat = sum([2*Nh['fourier_order']+1 for Nh in self.hole_surfaces])


class AuxiliaryGeometry():
    """ This class contains the basic geometry of the system. """

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
    def distance_matrices(self, surfaces_list):
        """
        This function calculates the distance matrix between the center points of the conductors.

        Returns:
        numpy.ndarray: The distance matrix.
        """
        # Principal Dimension: Number of conductor surfaces
        num_s = len(surfaces_list)
        d_pq = np.zeros((num_s, num_s))
        x_pq = np.zeros_like(d_pq)
        y_pq = np.zeros_like(d_pq)
        theta_pq = np.zeros_like(d_pq)

        for idx_p, p_cond in enumerate(surfaces_list):
            for idx_q, q_cond in enumerate(surfaces_list):

                # Distance between the center points of the conductors
                x_pq[idx_p, idx_q] = p_cond['center_point'][0] - \
                    q_cond['center_point'][0]
                y_pq[idx_p, idx_q] = p_cond['center_point'][1] - \
                    q_cond['center_point'][1]

                # Angle between the center points of the conductors
                theta_pq[idx_p, idx_q] = np.arctan2(
                    y_pq[idx_p, idx_q], x_pq[idx_p, idx_q])

                # Matrix Distance [np.array]
                if idx_p != idx_q:
                    d_pq[idx_p, idx_q] = np.linalg.norm(
                        np.array(p_cond['center_point']) - np.array(q_cond['center_point']))

        return d_pq, x_pq, y_pq, theta_pq

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
