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
from mtl_main.source import MulticonductorTransmissionLine as MTL


class HomogeneousLosslessMedium(MTL):
    """ This class contains the frequency-dependent parameters of the system. """

    def __init__(self, mtl, frequency):
        super().__init__(mtl)

        # Angular frequency [float]
        self.w = 2 * np.pi * frequency

        # Conductors wave-number [float]
        self.k = np.sqrt(self.w * self.mu * (self.w * self.epsilon - 1j * self.sigma))

        # Free-space wave-number [float]
        self.kout = self.w * np.sqrt(self.mu * self.epsilon_out)

    # Surface admittance operator [np.array]
    # Equation (2.20) [1]
    def ynp_for_solid(self, n, p):
        """
        This method calculates the surface admittance operator for a solid conductor, Yn(p).

        Parameters:
            n (int): The order of the Bessel function.
            p (int): The index of the conductor.

        Returns:
            list: The surface admittance operator Yn_p.
        """
        outer_radius = np.array([c['radius'][1] for c in self.mtl.values()])

        k_ap = self.k[p] * outer_radius[p]
        k0_ap = self.kout[p] * outer_radius[p]
        mu = self.mu[p]
        n = np.abs(n)
        term_1 = k_ap * jvp(n, k_ap) / mu / jv(n, k_ap)
        term_2 = k0_ap * jvp(n, k0_ap) / sc.mu_0 / jv(n, k0_ap)

        return 2 * np.pi / 1j / self.w * (term_1 - term_2)

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

    # mu_n function [int]
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
        """
        This method calculates the surface admittance operator 
        for a hollow conductor, Yn(p).

        Parameters:
            n (int): The order of the Bessel function.
            p (int): The index of the conductor.

        Returns:
            list: The surface admittance operator Yn_p.
        """

        # inner and outer radii [float]
        inner_radius = np.array([c['radius'][0] for c in self.mtl.values()])
        outer_radius = np.array([c['radius'][1] for c in self.mtl.values()])
        ap = outer_radius[cp]
        bp = inner_radius[cp]

        # Bessel arguments [float]
        kap = self.k[cp] * ap
        kbp = self.k[cp] * bp
        kout_ap = self.kout[cp] * ap
        kout_bp = self.kout[cp] * bp

        # Permeability of the medium [float]
        mu = self.mu[cp]

        # Matrix elements [float]
        #try:
        y11_n = self.chi_n(n, kap, kbp) / self.eme_n(n, kap, kbp) / mu - (
            self.chi_n(n, kout_ap, kout_bp) / self.eme_n(n, kout_ap, kout_bp) / sc.mu_0)

        y12_n = self.chi_n(n, kout_bp, kout_bp) / self.eme_n(n, kout_ap, kout_bp) / sc.mu_0 - (
            self.chi_n(n, kbp, kbp) / self.eme_n(n, kap, kbp) / mu)

        y21_n = self.chi_n(n, kout_ap, kout_ap) / self.eme_n(n, kout_ap, kout_bp) / sc.mu_0 - (
            self.chi_n(n, kap, kap) / self.eme_n(n, kap, kbp) / mu)

        y22_n = self.chi_n(n, kbp, kap) / self.eme_n(n, kap, kbp) / mu - (
            self.chi_n(n, kout_bp, kout_ap) / self.eme_n(n, kout_ap, kout_bp) / sc.mu_0)

        # Matrix ynp [np.array]
        matrix = np.array([[y11_n, y12_n], [y21_n, y22_n]])
        return 2 * np.pi / 1j / self.w * matrix

    # Matrix U [np.array]
    # Equation (2.39) [1]
    def u_matrix(self):
        """ This function calculates the matrix U."""

        # Initialize U matrix with zeros
        u_mtx = np.zeros((self.N, len(self.mtl.values())))

        # Initialize row index
        row_index = 0

        for p, conductor in enumerate(self.mtl.values()):

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
        Fill the Ys matrix based on p (conductor), n (order of filling), and conductor type.

        Parameters:
        p (int): Conductor
        n_max (int): Maximum order of filling
        conductor_type (str): Type of conductor ('solid' or 'hollow')

        Returns:
        np.ndarray: Ys matrix
        """

        blocks = []

        for p, conductor in enumerate(self.mtl.values()):
            Np = conductor['fourier_order'] # pylint: disable=invalid-name

            if conductor['radius'][0] == 0:
                for n in range(-Np, Np + 1):
                    ynp = np.array([[self.ynp_for_solid(n, p)]])
                    blocks.append(ynp)

            elif conductor['radius'][0] > 0:
                for n in range(-Np, Np + 1):
                    ynp = self.ynp_for_hollow(n, p)
                    blocks.append(ynp)

            else:
                raise ValueError(
                    "Unknown conductor type. Use 'solid' or 'hollow'.")

        return self.create_diagonal_matrix(blocks)

    # Matrix Z [np.array]
    # Equation (2.61) [1]
    def z_partial(self, green_matrix):
        """
        This function calculates the matrix Z.

        The matrix Z is the impedance matrix of the system.

        Returns:
        numpy.ndarray: The matrix Z.
        """

        jwu0 = 1j * self.w * sc.mu_0
        ys = self.ys_matrix()
        u = self.u_matrix()

        # Perform LU factorization of the matrix M
        matrix = np.eye(self.N) - jwu0 * (ys @ green_matrix)
        lu, piv = lu_factor(matrix)

        # Solve the linear system Mx = b for Ys
        solution = lu_solve((lu, piv), ys @ u)

        return u.T @ solution

    # Generalized Capacitance Matrix [np.array]
    def generalized_capacitance_matrix(self, green_matrix):
        """
        Calcula a matriz de capacitância física (n x n) a partir da matriz de
        capacitância generalizada ((n+1) x (n+1)), seguindo a Eq. 5.21 de Clayton Paul.

        A fórmula implementada é:
        C_ij = c_ij - ( (soma da linha i de c) * (soma da coluna j de c) ) / (soma total de c)

        Onde 'c' é a matriz generalizada e 'C' é a matriz física resultante.
        Assume-se que o condutor de índice 0 da matriz generalizada é o de referência
        e está sendo eliminado.

        Args:
            matriz_generalizada (np.ndarray): A matriz de capacitância generalizada
                                            simétrica de ordem (n+1) x (n+1).

        Returns:
            np.ndarray: A matriz de capacitância física de ordem n x n.
            
        Raises:
            ValueError: Se a matriz de entrada não for quadrada ou se a soma de
                        seus elementos for zero.
        """

        u = self.u_matrix()
        # Solve the linear system Gx = U and calculate U^T * G^{-1} * U
        uT_gInv_u = u.T @ lu_solve(lu_factor(green_matrix), u)

        # Generalized Capacitance Matrix [1]
        return - self.epsilon[0] * uT_gInv_u

    # Maxwellian Capacitance Matrix [np.array]
    def maxwellian_capacitance_matrix(self, generalized_capacitance_matrix):
        """
        Calcula a matriz de capacitância física (n x n) a partir da matriz de
        capacitância generalizada ((n+1) x (n+1)), seguindo a Eq. 5.21 de Clayton Paul.

        A fórmula implementada é:
        C_ij = c_ij - ( (soma da linha i de c) * (soma da coluna j de c) ) / (soma total de c)

        Onde 'c' é a matriz generalizada e 'C' é a matriz física resultante.
        Assume-se que o condutor de índice 0 da matriz generalizada é o de referência
        e está sendo eliminado.

        Args:
            matriz_generalizada (np.ndarray): A matriz de capacitância generalizada
                                            simétrica de ordem (n+1) x (n+1).

        Returns:
            np.ndarray: A matriz de capacitância física de ordem n x n.
            
        Raises:
            ValueError: Se a matriz de entrada não for quadrada ou se a soma de
                        seus elementos for zero.
        """

        gc = generalized_capacitance_matrix

        # --- Validação da entrada com assert ---
        assert isinstance(gc, np.ndarray), "A entrada deve ser um array NumPy."
        assert np.sum(gc) != 0, "A soma total dos elementos da matriz generalizada não pode ser zero."
        assert gc.ndim == 2, "A entrada deve ser uma matriz 2D (array de 2 dimensões)."
        assert gc.shape[0] == gc.shape[1], "A entrada deve ser uma matriz quadrada."
        assert gc.shape[0] >= 2, "A matriz generalizada deve ser de ordem mínima 2x2."

        # Ordem da matriz generalizada (N = n+1)
        N = gc.shape[0]

        # 2. Numerador: Soma de cada linha e de cada coluna
        # Para uma matriz simétrica, as somas das linhas e colunas são iguais.
        row_sum = np.sum(gc, axis=1)     # axis=1 soma ao longo das colunas
        column_sum = np.sum(gc, axis=0)  # axis=0 soma ao longo das linhas

        # assert np.equal(row_sum, column_sum).all(), "As somas das linhas e colunas devem ser iguais."

        # Inicializa a matriz de capacitância física n x n com zeros
        matrix_c = np.zeros((N - 1, N - 1), dtype=gc.dtype)

        # Itera sobre os índices da matriz física (de 1 a n na matriz original)
        # Condutor de índice 0 é o de referência e não é incluído na matriz física
        for i in range(1, N):
            for j in range(1, N):
                c_ij = gc[i, j]
                row_i_sum = row_sum[i]
                column_j_sum = column_sum[j]
                
                # Eq. 5.21 [2]
                matrix_c[i - 1, j - 1] = c_ij - (row_i_sum * column_j_sum) / np.sum(gc)

        return matrix_c


class LosslessPostProcessing(MTL):
    """ This class contains the post-processing parameters for the system. """

    def __init__(self, mtl):
        super().__init__(mtl)

        # # Line ID [list] [int]
        # self.line_id = [conductor['line_id'] for conductor in self.mtl]

        # # Active lines [list] [int]
        # self.active_lines = [line for line in self.mtl if line['line_type'] == 'active']

        # Lista de todos os line_id's presentes no sistema.
        self.line_id = [conductor['line_id'] for conductor in self.mtl.values()]

        # Lista de dicionários dos condutores que são do tipo 'active'.
        self.active_lines = [line for line in self.mtl.values() if line['line_type'] == 'active']

    # Incident Matrix Q [np.array]
    # Equation (A.2) [1]
    def q_incident_matrix(self):
        """
        On the i-th row of Q, '1's are present in the columns
        corresponding to the conductors which are part of the 
        i-th line, and all other columns are zero.
        Reference: PAG. 122 [1]
        """

        matrix_q = np.zeros((len(set(self.line_id)), len(self.mtl)))

        for x in range(len(self.line_id)):
            for y in range(len(self.mtl)):
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
        This function calculates the full impedance matrix Z.

        The full impedance matrix Z is the impedance matrix of the system.

        Returns:
        numpy.ndarray: The full impedance matrix Z.
        """

        s = self.s_incident_matrix()
        q = self.q_incident_matrix()

        # Perform LU factorization of the Z matrix
        lu, piv = lu_factor(z_partial)

        # Solve the linear system Zx = b for Q^T
        # Calculate the matrix QZ^{-1}Q^T
        qz_inv_qt = q @ lu_solve((lu, piv), q.T)

        return s.T @ qz_inv_qt @ s
    

    # Matriz Rs [np.array]
    def rs_matrix(self, z_total):
        """
        This function calculates the series resistance matrix Rs.
        The series resistance matrix Rs is the real part of the total series impedance matrix.

        Returns:
        numpy.ndarray: The series resistance matrix Rs.
        """
        return np.real(z_total)


    # Matriz Ls [np.array]
    def ls_matrix(self, z_total, frequency):
        """
        This function calculates the series inductance matrix Ls.
        The series inductance matrix Ls is the imaginary part of the total series impedance matrix.

        Returns:
        numpy.ndarray: The series inductance matrix Ls.
        """
        return np.imag(z_total) / (2 * np.pi * frequency)
