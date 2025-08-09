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
        self.mtl = {key: value for key, value in mtl.items() if isinstance(key, int)}

        # MTL Type
        self.mtl_type = mtl.get('type', 'unknown')

        # MTL Conductor Reference Index
        self.idx_ref = mtl.get('idx_ref_conductor', 0)

        key_conductors = sorted(self.mtl.keys())
        key_expected = list(range(len(key_conductors)))
        
        assert key_conductors == key_expected, "As tags dos condutores devem ser uma sequência de inteiros começando em 0 (ex: 0, 1, 2, ...)."
        assert len(self.mtl) > 1, "The MulticonductorTransmissionLine must have at least two conductors."
        assert self.mtl[0]['line_type'] == 'return', "The conductor index '0' must be the return path."
        assert self.idx_ref is not None, "The reference conductor index must be defined."
        assert self.idx_ref in self.mtl, f"The reference conductor index {self.idx_ref} must be in the MTL dictionary."

        # Conductor surfaces dictionary
        self.surfaces = []

        for key, conductor in self.mtl.items():
            # Check if the conductor is hollow
            if conductor['radius'][0] != 0:
                for radius in conductor['radius']:
                    self.surfaces.append({
                        'type': 'conductor',
                        'tag': key,
                        'radius': radius,
                        'center_point': conductor['center_point'],
                        'fourier_order': conductor['fourier_order'],
                    })

            # Then, the conductor is solid
            else:
                self.surfaces.append({
                    'type': 'conductor',
                    'tag': key,
                    'radius': conductor['radius'][1],
                    'center_point': conductor['center_point'],
                    'fourier_order': conductor['fourier_order'],
                })

            # Check if the conductor has a sheath
            if conductor['sheath'] is not None:
                # Add the sheath surface
                self.surfaces.append({
                    'type': 'sheath',
                    'tag': key,
                    'radius': conductor['radius'][1] + conductor['sheath']['thickness'],
                    'center_point': conductor['sheath']['center_point'],
                    'fourier_order': conductor['sheath']['fourier_order'],
                    'relative_permittivity': conductor['sheath']['relative_permittivity'],
                })

        # Dimension N
        # Equation (2.36) [1]
        self.NF_List = [2 * surface['fourier_order'] + 1 for surface in self.surfaces]
        self.N = sum(self.NF_List)

        # Número de coeficientes harmônicos de Fourier por condutor
        self.NF = self.NF_List[self.idx_ref]

        # Conductors Permeability [np.array]
        self.mu = np.array([mu_0 * cond['relative_permeability'] for cond in self.mtl.values()]) 

        # Conductors Permittivity [np.array]
        self.epsilon = np.array([epsilon_0 * cond['relative_permittivity'] for cond in self.mtl.values()]) 

        # Conductors conductivity [np.array]
        self.sigma = np.array([cond['conductivity'] for cond in self.mtl.values()])

        # Free-Space Permittivity [np.array]
        self.epsilon_out = np.array([epsilon_0 * cond['relative_permittivity_out'] for cond in self.mtl.values()])

        # Distance matrices [np.array]
        self.D_pq, self.dx_pq, self.dy_pq, self.theta_pq = self.conductors_center_distance_matrix()

    # Contour Position Vector [np.array]
    # Equation (2.2) [1]
    def contour_vector_position(self, surfaces_list, theta, p):
        """ The boundary cp can be traced by the position vector rp(ap, θp) """

        xp, yp = surfaces_list[p]['center']
        ap = surfaces_list[p]['radius']

        return np.array([xp + ap * np.cos(theta), yp + ap * np.sin(theta)])

    # Distance matrices [np.array]
    def conductors_center_distance_matrix(self):
        """ Calculates the distance matrix between the center points of the conductors. """

        # Principal Dimension: Number of conductor surfaces
        N = len(self.mtl)  
        d_pq = np.zeros((N, N))
        x_pq = np.zeros_like(d_pq)
        y_pq = np.zeros_like(d_pq)
        theta_pq = np.zeros_like(d_pq)

        for p, conductor_p in self.mtl.items():
            cp = conductor_p['center_point']

            for q, conductor_q in self.mtl.items():
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

