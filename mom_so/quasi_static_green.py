"""
This script contains the Quasi-Static Green's matrices of the Lossless Homogeneous Medium.  

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

import math
import numpy as np
from scipy.integrate import quad
from mtl_main.source import MulticonductorTransmissionLine

class QuasiStatic():
    """ This class contains the Quasi-Static Green's matrix G of the system. """

    def __init__(self, model: MulticonductorTransmissionLine):
        self.model = model

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
        x_qp = q_surface_list[q]['center_point'][0] - p_surface_list[p]['center_point'][0]
        y_qp = q_surface_list[q]['center_point'][1] - p_surface_list[p]['center_point'][1]

        # Angle between the center points of the conductors
        theta_qp = np.arctan2(y_qp, x_qp)

        # Distance between the center points of the conductors
        p = np.array(p_surface_list[p]['center_point'])
        q = np.array(q_surface_list[q]['center_point'])
        d_qp = np.linalg.norm(q - p)

        return d_qp, x_qp, y_qp, theta_qp    

    # Equation (41) [2]
    def function_fn_theta(self, theta, n, p, q):
        """ This function calculates the function f_n(theta) """

        def integrand(theta_prime, theta, n, p, q):
            """ This function calculates the integrand of the function f_n(theta) """

            fn = 0
            rp = self.model.contour_vector_position(self.model.surfaces, theta, p)
            rq = self.model.contour_vector_position(self.model.surfaces, theta_prime, q)
            R = np.linalg.norm(rp - rq)

            if R != 0:
                fn = np.log(R) * np.exp(1j * n * theta_prime)

            return fn

        real_fn, _ = quad(lambda theta_prime: integrand(
            theta_prime, theta, n, p, q).real, 0, 2 * np.pi)

        imag_fn, _ = quad(lambda theta_prime: integrand(
            theta_prime, theta, n, p, q).imag, 0, 2 * np.pi)

        return (real_fn + 1j * imag_fn) / (2 * np.pi)

    # Analytical Green's sub-matrix G_{n',n}^{p,q} [np.array]
    # Equation (40) [2]
    def gpq_numeric(self, n_prime, n, p, q):
        """ Calculate the Greens' matrix G_{n',n}^{p,q} """

        def integrand(theta, n_prime, n, p, q):
            """ Calculate the integrand of Greens' matrix G_{n',n}^{p,q} """
            fn_theta = self.function_fn_theta(theta, n, p, q)

            return fn_theta * np.exp(-1j * n_prime * theta)

        real_g, _ = quad(lambda theta: integrand(
            theta, n_prime, n, p, q).real, 0, 2 * np.pi)

        imag_g, _ = quad(lambda theta: integrand(
            theta, n_prime, n, p, q).imag, 0, 2 * np.pi)

        return (real_g + 1j * imag_g) / (2 * np.pi)**2

    # Analytical Green's sub-matrix G_{n',n}^{p,q} [np.array]
    # Equation (40) [2]
    def gpq_master_thesis(self, n_prime, n, p, q):
        """
        This function calculates the Green's function G_{n',0}^{p,q}.

        Parameters:
        n (int): The order of the Green's function.
        p (int): The index of the first conductor.
        q (int): The index of the second conductor.

        Returns:
        float: The Green's function G_{n,0}.

        Reference: PAG. 127 [1]
        """

        # Surfaces radii [float]
        ap = self.model.surfaces[p]['radius']
        aq = self.model.surfaces[q]['radius']

        # Auxiliary vector of distances and angles [np.array]
        # PAG. 125 [1]
        dqp, x_qp, y_qp, theta_qp = self.distance_vector_dqp(self.model.surfaces, p, self.model.surfaces, q)

        # Auxiliary variable alfa [int]
        # PAG. 126 [1]
        alfa = np.abs(aq / (ap - dqp))

        # Auxiliary variable alfa_hat [int]
        # PAG. 127 [1]
        alfa_hat = dqp / ap

        # Green's function G_{n',n}^{p,q} [np.array]
        # 1st case: alfa < 1
        # PAG. 126 [1]
        if alfa < 1:
            # 1st sub-case: alfa < 1 and n == 0
            # PAG. 126 [1]
            if n == 0:
                # n' == 0
                if n_prime == 0:
                    common_term_1 = 1 / 2 / np.pi
                    # Equation (B.8) [1]
                    if alfa_hat < 1:
                        g_npn_pq = common_term_1 * np.log(ap)

                    # Equation (B.9) [1]
                    else:
                        g_npn_pq = common_term_1 * np.log(dqp)

                # n' != 0
                else:
                    common_term_2 = - np.exp(-1j * n_prime * theta_qp) / \
                        4 / np.pi / np.abs(n_prime)
                    # Equation (B.8) [1]
                    if alfa_hat < 1:
                        g_npn_pq = (common_term_2 *
                                    (- np.abs(dqp) / ap) ** np.abs(n_prime))

                    # Equation (B.9) [1]
                    else:
                        g_npn_pq = (common_term_2 * (- ap / dqp)
                                    ** np.abs(n_prime))

            # 2nd sub-case: alfa < 1 and n > 0
            elif n > 0:
                # Equation (B.12) [1]
                if alfa_hat > 1:
                    if n_prime >= 1:
                        g_npn_pq = 0

                    elif n_prime <= 0:
                        frac_1 = - np.pi / ((2 * np.pi) ** 2 * n)
                        frac_2 = (aq ** n) / (ap ** n_prime)
                        frac_3 = (-1) ** n_prime
                        frac_4 = math.comb(n - n_prime - 1, -n_prime)
                        frac_5 = (x_qp - 1j * y_qp) ** (-n + n_prime)

                        g_npn_pq = (frac_1 * frac_2 *
                                    frac_3 * frac_4 * frac_5)

                # Equation (B.13) [1]
                elif alfa_hat == 0:
                    if n_prime != n:
                        g_npn_pq = 0

                    else:
                        frac_6 = - aq ** n / (4 * np.pi * n * ap ** n_prime)
                        g_npn_pq = frac_6

                elif alfa_hat < 1:
                    # Equation (B.14) [1]
                    if n_prime < 1:
                        frac_6 = - aq ** n / (4 * np.pi * n * ap ** n_prime)
                        frac_4 = math.comb(n - n_prime - 1, -n_prime)
                        frac_7 = math.comb(n - n_prime - 1, n - 1)
                        frac_8 = (x_qp - 1j * y_qp) ** (-n + n_prime)
                        frac_9 = (-x_qp + 1j * y_qp) ** (n_prime - n)

                        g_npn_pq = (frac_6 * (
                            (-1) ** (-n_prime) * frac_4 * frac_8 + (
                                (-1) ** (n - 1) * frac_7 * frac_9))
                        )

                    # Equation (B.15) [1]
                    if n_prime >= 1:
                        if n_prime >= n and n_prime >= 0:
                            frac_6 = - aq ** n / \
                                (4 * np.pi * n * ap ** n_prime)
                            frac_9 = (-x_qp + 1j * y_qp) ** (n_prime - n)
                            frac_10 = math.comb(n_prime - 1, n - 1)

                            g_npn_pq = frac_6 * frac_9 * frac_10

                        elif n_prime >= 1 and n_prime < n:
                            g_npn_pq = 0

        # 2nd case: alfa > 1
        # PAG. 128 [1]
        elif alfa > 1:
            # 1st sub-case: alfa > 1 and n == 0
            # PAG. 128 [1]
            if n == 0:
                # Equation (B.17) [1]
                if n_prime == 0:
                    g_npn_pq = 1 / 2 / np.pi * np.log(aq)

                else:
                    g_npn_pq = 0

            # 2nd sub-case: alfa > 1 and n > 0
            # PAG. 129 [1]
            if n > 0:
                # Equation (B.20) [1]
                if n >= n_prime and n_prime >= 0:
                    frac_11 = - ap ** n_prime / (4 * np.pi * n * aq ** n)
                    frac_12 = math.comb(n, n_prime)
                    frac_13 = (x_qp + 1j * y_qp) ** (n - n_prime)

                    g_npn_pq = frac_11 * frac_12 * frac_13

                else:
                    g_npn_pq = 0

        # 3rd case: alfa == 1
        # PAG. 2481 [2] IEEE Paper
        elif alfa == 1:
            # 1st sub-case: alfa == 1 and n == 0
            # Equation (52) [2]
            if n == 0:
                if n_prime == 0:
                    g_npn_pq = 1 / 2 / np.pi * np.log(ap)

                else:
                    g_npn_pq = 0

            # 2st sub-case: alfa == 1 and n != 0
            # Equation (53) [2]
            elif n != 0:
                if n_prime == n:
                    g_npn_pq = - 1 / 4 / np.pi / np.abs(n)

                else:
                    g_npn_pq = 0

        return g_npn_pq

    # Analytical Green's sub-matrix G_{n',n}^{p,q} [np.array]
    # Equation (40) [2]
    def gpq_ieee_paper(self, n_prime, n, p, q, ap, aq):
        """
        This function calculates the Green's function G_{n',0}^{p,q}.

        Parameters:
        n (int): The order of the Green's function.
        p (int): The index of the first conductor.
        q (int): The index of the second conductor.

        Returns:
        float: The Green's function G_{n,0}.

        Reference: PAG. 127 [1]
        """

        # Auxiliary vector of distances and angles [np.array]
        # PAG. 125 [1]
        dqp, xqp, yqp, theta_qp = self.distance_vector_dqp(self.model.surfaces, p, self.model.surfaces, q)

        # Green's function G_{n',0}^{p,q} [float]
        if p != q and n == 0:
            # Expression (45) [2]
            if n_prime == 0:
                g_npn_pq = 1 / 2 / np.pi * np.log(dqp)

            # Expression (46) [2]
            else:
                term_1 = - 1 / 4 / np.pi / np.abs(n_prime)
                term_2 = (ap / dqp) ** np.abs(n_prime)
                term_3 = (-(xqp - 1j * yqp) /
                          dqp) ** n_prime
                g_npn_pq = term_1 * term_2 * term_3

        # Expression (49) [2]
        elif p != q and n > 0 and n_prime >= 1:
            g_npn_pq = 0

        # Expression (50) [2]
        elif p != q and n > 0 and n_prime <= 1:
            term_4 = - np.pi * (aq ** n) / ((-ap) ** n_prime)
            term_5 = math.comb(n - n_prime - 1, -n_prime)
            term_6 = (xqp - 1j * yqp) ** (-n + n_prime) / (
                ((2 * np.pi) ** 2 * n))

            g_npn_pq = term_4 * term_5 * term_6

        # Expression (52) [2]
        elif p == q and n == 0:
            if n_prime == 0:
                g_npn_pq = 1 / 2 / np.pi * np.log(ap)
            else:
                g_npn_pq = 0

        # Expression (53) [2]
        elif p == q and n != 0:
            if n_prime == n:
                g_npn_pq = - 1 / 4 / np.pi / np.abs(n)
            else:
                g_npn_pq = 0

        return g_npn_pq

    # Green's Matrix G [np.array]
    # Equation (2.54) [1]
    def green_matrix(self, green_mode='Analytically'):
        """ Computes the full Green's matrix for the given set of conductors """

        # Define the full Green's matrix
        g_matrix = np.zeros((self.model.N, self.model.N), dtype=complex)

        # Create an offset array to keep track of the starting index for each block
        offsets = np.cumsum([0] + [2 * Np['fourier_order'] + 1 for Np in self.model.surfaces[:-1]])

        # Number of conductor surfaces [int]
        ns = len(self.model.surfaces)

        for p in range(ns):
            # Offset for the p-th conductor [int]
            off_p = offsets[p]

            # Number of surface points for the p-th surface [int]
            Np = self.model.surfaces[p]['fourier_order']

            for q in range(ns):

                # Offset for the q-th conductor [int]
                off_q = offsets[q]

                # Number of surface points for the q-th surface [int]
                Nq = self.model.surfaces[q]['fourier_order']

                for n_prime in range(-Np, Np + 1):
                    for n in range(-Nq, Nq + 1):

                        # Analytical Evaluation
                        if green_mode == 'Analytically':
                            if n < 0:
                                g_value = np.conjugate(
                                    self.gpq_master_thesis(-n_prime, -n, p, q))
                            else:
                                g_value = self.gpq_master_thesis(n_prime, n, p, q)

                        # Numeric Evaluation
                        else:
                            if n < 0:
                                g_value = np.conjugate(
                                    self.gpq_numeric(-n_prime, -n, p, q))
                            else:
                                g_value = self.gpq_numeric(n_prime, n, p, q)

                        # Assign the Green's value to the Full Green's matrix G
                        g_matrix[off_p + Np + n_prime, off_q + Nq + n] = g_value

        return g_matrix
