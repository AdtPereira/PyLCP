"""
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
    Canada July 18-20, 2013. https://www.ipstconf.org/Proc_IPST2013.php.

[4] PAUL, Clayton R. Analysis of multiconductor transmission lines. 2. ed. Hoboken,
    N.J.: John Wiley & Sons, Inc., c2008.

[5] XUE, Haoyan. General Formulation and Accurate Evaluation of Earth-Return Parameters
    for Overhead / Underground Cables. PhD thesis, Department of Electrical Engineering,
    École Polytechnique de Montréal, Université de Montréal, August 2018. 

"""

import numpy as np
import scipy.constants as sc
from scipy.special import jv, jvp, h1vp, h2vp
from scipy.linalg import lu_factor, lu_solve, inv
from mtl_main.source import MulticonductorTransmissionLine

class HomogeneousLosslessMedium():
    """
    This class contains the frequency-dependent parameters of the system.
    It is designed to perform vectorized calculations over an array of frequencies.
    """

    def __init__(self, model: MulticonductorTransmissionLine, frequencies: np.ndarray):
        self.model = model

        # Ensure frequencies is a numpy array
        self.frequencies = np.asarray(frequencies)
        
        # Angular frequency [np.ndarray of shape (n_freqs,)]
        self.w = 2 * np.pi * self.frequencies

        # To enable broadcasting, we reshape frequency-dependent arrays to (n_freqs, 1)
        # and conductor-dependent arrays to (1, n_conductors).
        w_col = self.w[:, np.newaxis]
        mu_row = self.model.mu[np.newaxis, :]
        epsilon_row = self.model.epsilon[np.newaxis, :]
        sigma_row = self.model.sigma[np.newaxis, :]
        epsilon_out_row = self.model.epsilon_out[np.newaxis, :]

        # Conductors wave-number [np.ndarray of shape (n_freqs, n_conductors)]
        self.k = np.sqrt(w_col * mu_row * (w_col * epsilon_row - 1j * sigma_row))

        # Free-space wave-number [np.ndarray of shape (n_freqs, n_conductors)]
        self.kout = w_col * np.sqrt(mu_row * epsilon_out_row)

    # Surface admittance operator [np.array]
    # Equation (2.20) [1]
    def ynp_for_solid(self, n, p):
        """
        This method calculates the surface admittance operator for a solid conductor, Yn(p).
        It is vectorized to compute for all frequencies at once.
        """
        outer_radius = np.array([c['radius'][1] for c in self.model.mtl.values()])

        # k and kout are now 2D arrays (freqs, conductors). We select the column for conductor p.
        k_ap = self.k[:, p] * outer_radius[p]
        k0_ap = self.kout[:, p] * outer_radius[p]
        mu = self.model.mu[p]
        n = np.abs(n)
        
        # scipy.special functions are ufuncs and work element-wise on numpy arrays.
        term_1 = k_ap * jvp(n, k_ap) / mu / jv(n, k_ap)
        term_2 = k0_ap * jvp(n, k0_ap) / sc.mu_0 / jv(n, k0_ap)

        # self.w is a 1D array, so the operation is vectorized.
        return 2 * np.pi / (1j * self.w) * (term_1 - term_2)

    # chi_n function [int]
    # Equation (11) [3]
    def chi_n(self, n, alfa, beta):
        """
        This function calculates the qui_n function.

        Parameters:
            n (int): The order of the Bessel function.
            alfa (float): The alfa parameter.
            beta (float): The beta parameter.

        Returns:
            float: The value of the qui_n function.
        """

        n = np.abs(n)
        term_1 = h1vp(n, beta) * h2vp(n, alfa, 0)
        term_2 = h2vp(n, beta) * h1vp(n, alfa, 0)
        return beta * (term_1 - term_2)

    # eme_n function [int]
    # Equation (12) [3]
    def eme_n(self, n, alfa, beta):
        """
        This function calculates the qui_n function.

        Parameters:
            n (int): The order of the Bessel function.
            alfa (float): The alfa parameter.
            beta (float): The beta parameter.

        Returns:
            float: The value of the qui_n function.
        """

        n = np.abs(n)
        term_1 = h1vp(n, alfa, 0) * h2vp(n, beta, 0)
        term_2 = h1vp(n, beta, 0) * h2vp(n, alfa, 0)
        return term_1 - term_2

    # Surface admittance operator for hollow conductors [np.array]
    # Equation (2.31) [1]
    # Equation (10) [3]
    def ynp_for_hollow(self, n, cp):
        inner_radius = np.array([c['radius'][0] for c in self.model.mtl.values()])
        outer_radius = np.array([c['radius'][1] for c in self.model.mtl.values()])
        ap = outer_radius[cp]
        bp = inner_radius[cp]

        # Arguments are now 1D arrays of shape (n_freqs,)
        kap = self.k[:, cp] * ap
        kbp = self.k[:, cp] * bp
        kout_ap = self.kout[:, cp] * ap
        kout_bp = self.kout[:, cp] * bp

        mu = self.model.mu[cp]

        # Matrix elements will be 1D arrays
        y11_n = self.chi_n(n, kap, kbp) / self.eme_n(n, kap, kbp) / mu - (
            self.chi_n(n, kout_ap, kout_bp) / self.eme_n(n, kout_ap, kout_bp) / sc.mu_0)
        y12_n = self.chi_n(n, kout_bp, kout_bp) / self.eme_n(n, kout_ap, kout_bp) / sc.mu_0 - (
            self.chi_n(n, kbp, kbp) / self.eme_n(n, kap, kbp) / mu)
        y21_n = self.chi_n(n, kout_ap, kout_ap) / self.eme_n(n, kout_ap, kout_bp) / sc.mu_0 - (
            self.chi_n(n, kap, kap) / self.eme_n(n, kap, kbp) / mu)
        y22_n = self.chi_n(n, kbp, kap) / self.eme_n(n, kap, kbp) / mu - (
            self.chi_n(n, kout_bp, kout_ap) / self.eme_n(n, kout_ap, kout_bp) / sc.mu_0)

        # Stacking results to get a (2, 2, n_freqs) matrix
        matrix = np.array([[y11_n, y12_n], [y21_n, y22_n]])
        
        # Reshape self.w for broadcasting
        return (2 * np.pi / (1j * self.w))[:, np.newaxis, np.newaxis] * np.transpose(matrix, (2, 0, 1))

    # Matrix U [np.array]
    # Equation (2.39) [1]
    def u_matrix(self):
        """ This function calculates the matrix U."""

        # Initialize U matrix with zeros
        u_mtx = np.zeros((self.model.N, len(self.model.mtl.values())))

        # Initialize row index
        row_index = 0

        for p, conductor in enumerate(self.model.mtl.values()):

            # Number of surface points for the p-conductor
            surf_points = conductor['fourier_order']

            if conductor['radius'][0] == 0:
                # Central row for the solid conductor
                central_row = row_index + surf_points

                # Set the values in the U matrix
                u_mtx[central_row, p] = 1

                # Skip rows corresponding to this solid conductor
                row_index += 2 * surf_points + 1

            elif conductor['radius'][0] > 0:
                # Central row for the hollow conductor
                mid_row_1 = row_index + surf_points
                mid_row_2 = row_index + (2 * surf_points + 1) + (surf_points)

                u_mtx[mid_row_1, p] = 1
                u_mtx[mid_row_2, p] = 1

                # Skip rows corresponding to this hollow conductor
                row_index += 2 * (2 * surf_points + 1)

        # print("Matrix U: \n", u_mtx)
        return u_mtx

    # Matrix Ys [np.array]
    # Equation (2.38) [1]
    def ys_matrix(self):
        """
        Builds the block-diagonal Ys matrix for all frequencies at once.
        The resulting matrix has the shape (n_frequencies, N, N).
        """
        n_freqs = len(self.frequencies)
        ys_3d = np.zeros((n_freqs, self.model.N, self.model.N), dtype=np.complex128)
        
        current_idx = 0
        for p, conductor in enumerate(self.model.mtl.values()):
            Np = conductor['fourier_order']
            
            if conductor['radius'][0] == 0:  # Solid conductor
                block_size = 1
                num_blocks = 2 * Np + 1
                for n_idx, n in enumerate(range(-Np, Np + 1)):
                    ynp_vec = self.ynp_for_solid(n, p) # shape (n_freqs,)
                    start = current_idx + n_idx
                    # Assign the vector to the diagonal for all frequencies
                    ys_3d[:, start, start] = ynp_vec
                current_idx += num_blocks
            
            elif conductor['radius'][0] > 0:  # Hollow conductor
                block_size = 2
                num_blocks = 2 * Np + 1
                for n_idx, n in enumerate(range(-Np, Np + 1)):
                    # ynp_mat has shape (n_freqs, 2, 2)
                    ynp_mat = self.ynp_for_hollow(n, p)
                    start = current_idx + n_idx * block_size
                    end = start + block_size
                    # Assign the 2x2 block for all frequencies
                    ys_3d[:, start:end, start:end] = ynp_mat
                current_idx += num_blocks * block_size
                
            else:
                raise ValueError("Unknown conductor type.")

        return ys_3d

    # Matrix Z [np.array]
    # Equation (2.61) [1]
    def z_partial(self, G):
        """
        Calculates the partial impedance matrix Z for all frequencies.
        
        NOTE: A loop is necessary here because scipy.linalg.lu_solve does not
        support solving a stack of matrices in a vectorized manner.
        
        Returns:
            np.ndarray: A 3D array of shape (n_frequencies, n_conductors, n_conductors).
        """
        n = len(self.frequencies)
        N = len(self.model.mtl.values())
        
        # Pre-calculate frequency-dependent and independent matrices
        I = np.eye(self.model.N)
        ys = self.ys_matrix()           # 3D: (ns, N, N)
        U = self.u_matrix()             # 2D: (N, n)
        jwu0 = 1j * self.w * sc.mu_0    # 1D: (n,)
        
        # Pre-allocate result array
        z_partial = np.zeros((n, N, N), dtype=np.complex128)
        
        # Loop over each frequency
        for i in range(n):
            Ys = ys[i, :, :]
            M = I - jwu0[i] * (Ys @ G)
            
            # Compute A = U.T @ (M**-1) @ (Ys @ U)
            A = U.T @ lu_solve(lu_factor(M), Ys @ U)

            # Compute Z = A**-1
            z_partial[i, :, :] = lu_solve(lu_factor(A), np.eye(A.shape[0]))
            
        return z_partial

    # Generalized Capacitance Matrix [np.array]
    def generalized_capacitance_matrix(self, green_matrix):
        """
        Computes the physical capacitance matrix (n x n) from the generalized
        capacitance matrix ((n+1) x (n+1)), following Eq. 5.21 of Clayton Paul.

        The implemented formula is:
        C_ij = c_ij - ( (sum of row i of c) * (sum of column j of c) ) / (total sum of c)

        Where 'c' is the generalized matrix and 'C' is the resulting physical matrix.
        It is assumed that the conductor at index 0 of the generalized matrix is the
        reference one and is being eliminated.

        Args:
            green_matrix (np.ndarray): The symmetric generalized capacitance
                                            matrix of order (n+1) x (n+1).

        Returns:
            np.ndarray: The physical capacitance matrix of order n x n.

        Raises:
            ValueError: If the input matrix is not square or if the sum of
                        its elements is zero.
        """

        u = self.u_matrix()
        # Solve the linear system Gx = U and calculate U^T * G^{-1} * U
        uT_gInv_u = u.T @ lu_solve(lu_factor(green_matrix), u)

        # Generalized Capacitance Matrix [1]
        return -1 * self.model.epsilon[0] * uT_gInv_u

    # Maxwellian Capacitance Matrix [np.array]
    def maxwellian_capacitance_matrix(self, generalized_capacitance_matrix: np.ndarray) -> np.ndarray:
        """
        Calculates the physical (Maxwellian) capacitance matrix from the
        generalized matrix using a vectorized approach based on Eq. 5.21 of Clayton Paul.

        The formula C_ij = c_ij - ( (sum of row i) * (sum of col j) ) / (total sum of c)
        is implemented using NumPy slicing and outer product for efficiency.
        """
        gc = generalized_capacitance_matrix

        # --- Input validation ---
        assert isinstance(gc, np.ndarray), "Input must be a NumPy array."
        total_sum = np.sum(gc)
        assert total_sum != 0, "The total sum of the generalized matrix cannot be zero."
        assert gc.ndim == 2 and gc.shape[0] == gc.shape[1], "Input must be a square 2D matrix."
        assert gc.shape[0] >= 2, "The generalized matrix must be at least 2x2."

        # --- Vectorized Calculation ---
        
        # Extract the submatrix c_ij (excluding the reference conductor at index 0)
        c_ij_submatrix = gc[1:, 1:]

        # Calculate sum of rows and columns, excluding the reference conductor
        row_sums = np.sum(gc, axis=1)[1:]
        col_sums = np.sum(gc, axis=0)[1:]

        # Calculate the outer product of the row and column sums to create the adjustment matrix
        adjustment_matrix = np.outer(row_sums, col_sums) / total_sum
        
        # Apply the formula in a single vectorized operation
        matrix_c = c_ij_submatrix - adjustment_matrix

        return matrix_c

class LosslessPostProcessing():
    """ This class contains the post-processing parameters for the system. """

    def __init__(self, model: MulticonductorTransmissionLine):
        self.model = model

        # List of every line_id present in the system.
        self.line_id = [conductor['line_id'] for conductor in self.model.mtl.values()]

        # List of dictionaries for the conductors of type 'active'.
        self.active_lines = [line for line in self.model.mtl.values() if line['line_type'] == 'active']

    # Incident Matrix Q [np.array]
    # Equation (A.2) [1]
    def q_incident_matrix(self):
        """
        On the i-th row of Q, '1's are present in the columns
        corresponding to the conductors which are part of the 
        i-th line, and all other columns are zero.
        Reference: PAG. 122 [1]
        """

        matrix_q = np.zeros((len(set(self.line_id)), len(self.model.mtl)))

        for x in range(len(self.line_id)):
            for y in range(len(self.model.mtl)):
                if x == self.line_id[y]:
                    matrix_q[x, y] = 1

        # print("Incident Matrix Q: \n", matrix_q)
        return matrix_q

    # Incident Matrix S [np.array]
    # Equation (A.10) [1]
    def s_incident_matrix(self):
        """
        This function calculates the incident matrix S.

        S is made up of ‘1’s, ‘0’s, and ‘-1’s. In the 
        i-th row, we have a “1” in the column corresponding 
        to the active line’s line number, and “-1” in the 
        column corresponding to the line number of its return 
        line. 

        # Reference: PAG. 124 [1]

        Returns:
        numpy.ndarray: The incident matrix S.
        """

        matrix_s_transpose = np.zeros((len(self.active_lines), len(set(self.line_id))))

        for i, line in enumerate(self.active_lines):
            active_line_id = line['line_id']
            return_line_id = line['line_return']

            if return_line_id != 'null':
                matrix_s_transpose[i, active_line_id] = 1
                matrix_s_transpose[i, return_line_id] = -1

        # print("Incident Matrix S.T: \n", matrix_s_transpose)
        return matrix_s_transpose.T

    # Matrix Z_line [np.array]
    def z_line_matrix(self, z_partial):
        """
        This function calculates the impedance matrix Z_line.

        The impedance matrix Z_line is the impedance matrix of the lines.

        Returns:
        numpy.ndarray: The impedance matrix Z_line.
        """

        q = self.q_incident_matrix()

        # Perform LU factorization of the Z matrix
        lu, piv = lu_factor(z_partial)

        # Solve the linear system Zx = b for Q^T
        # Calculate the matrix QZ^{-1}Q^T
        qz_inv_qt = np.dot(q, lu_solve((lu, piv), q.T))

        return inv(qz_inv_qt)

    # Matrix Z_full [np.array]
    # Equation (A.12) [1]
    def z_total(self, z_partial):
        """
        Calculates the full impedance matrix Zs for all frequencies.

        Parameters:
            z_partial_stack (np.ndarray): 3D array of partial impedances from HomogeneousLosslessMedium.
        
        Returns:
            np.ndarray: A 3D array of shape (n_frequencies, n_active_lines, n_active_lines).
        """
        n = z_partial.shape[0]
        n_act = len(self.active_lines)
        z_total = np.zeros((n, n_act, n_act), dtype=np.complex128)

        S = self.s_incident_matrix()
        Q = self.q_incident_matrix()

        # Loop over each frequency's z_partial matrix
        for i in range(n):
            Zp = z_partial[i, :, :]
            qz_inv_qt = Q @ Zp @ Q.T
            z_total[i, :, :] = S.T @ qz_inv_qt @ S

        return z_total

    # Rs matrix [np.array]
    def rs_matrix(self, z_total_stack):
        return np.real(z_total_stack)

    # Ls matrix [np.array]
    def ls_matrix(self, z_total_stack, frequencies):
        w = 2 * np.pi * np.asarray(frequencies)
        # Reshape w to (n_freqs, 1, 1) for broadcasting with (n_freqs, N, N) matrix
        return np.imag(z_total_stack) / w[:, np.newaxis, np.newaxis]
