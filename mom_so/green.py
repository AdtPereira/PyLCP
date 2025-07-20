"""
This script contains the Green's matrices of the Multilayered Lossy system.  

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
from scipy.constants import mu_0
from scipy.integrate import quad, dblquad
from scipy.special import jv, hankel2
from scipy.linalg import lu_solve, lu_factor
from .geometry import FreeSpace, AuxiliaryGeometry
from .patel import MultilayeredLossyMedium


class QuasiStatic(FreeSpace, AuxiliaryGeometry):
    """ This class contains the Quasi-Static Green's matrix G of the system. """

    # Equation (41) [2]
    def function_fn_theta(self, theta, n, p, q):
        """ This function calculates the function f_n(theta) """

        def integrand(theta_prime, theta, n, p, q):
            """ This function calculates the integrand of the function f_n(theta) """

            fn = 0
            rp = self.contour_vector_position(self.surfaces, theta, p)
            rq = self.contour_vector_position(self.surfaces, theta_prime, q)
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
        ap = self.surfaces[p]['radius']
        aq = self.surfaces[q]['radius']

        # Auxiliary vector of distances and angles [np.array]
        # PAG. 125 [1]
        dqp_dict = AuxiliaryGeometry().distance_vector_dqp(
            self.surfaces, p, self.surfaces, q)

        dqp = dqp_dict['norm']
        theta_qp = dqp_dict['angle']
        x_pq = dqp_dict['u_x']
        y_pq = dqp_dict['u_y']

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
                        frac_5 = (x_pq - 1j * y_pq) ** (-n + n_prime)

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
                        frac_8 = (x_pq - 1j * y_pq) ** (-n + n_prime)
                        frac_9 = (-x_pq + 1j * y_pq) ** (n_prime - n)

                        g_npn_pq = (frac_6 * (
                            (-1) ** (-n_prime) * frac_4 * frac_8 + (
                                (-1) ** (n - 1) * frac_7 * frac_9))
                        )

                    # Equation (B.15) [1]
                    if n_prime >= 1:
                        if n_prime >= n and n_prime >= 0:
                            frac_6 = - aq ** n / \
                                (4 * np.pi * n * ap ** n_prime)
                            frac_9 = (-x_pq + 1j * y_pq) ** (n_prime - n)
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
                    frac_13 = (x_pq + 1j * y_pq) ** (n - n_prime)

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
        dqp_dict = AuxiliaryGeometry().distance_vector_dqp(
            self.surfaces, p, self.surfaces, q)

        dqp = dqp_dict['norm']
        xqp = dqp_dict['u_x']
        yqp = dqp_dict['u_y']

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
    def g_tanaka(self, mode):
        """ Computes the full Green's matrix for the given set of conductors """

        # Define the full Green's matrix
        g_matrix = np.zeros((self.N, self.N), dtype=complex)

        # Create an offset array to keep track of the starting index for each block
        offsets = np.cumsum(
            [0] + [2 * Np['fourier_order'] + 1 for Np in self.surfaces[:-1]])

        # Number of conductor surfaces [int]
        ns = len(self.surfaces)

        for p in range(ns):
            # Offset for the p-th conductor [int]
            off_p = offsets[p]

            # Number of surface points for the p-th surface [int]
            Np = self.surfaces[p]['fourier_order']

            for q in range(ns):

                # Offset for the q-th conductor [int]
                off_q = offsets[q]

                # Number of surface points for the q-th surface [int]
                Nq = self.surfaces[q]['fourier_order']

                for n_prime in range(-Np, Np + 1):
                    for n in range(-Nq, Nq + 1):

                        # Analytical Evaluation
                        if mode == 'Analytically':
                            if n < 0:
                                g_value = np.conjugate(
                                    self.gpq_master_thesis(-n_prime, -n, p, q))
                            else:
                                g_value = self.gpq_master_thesis(
                                    n_prime, n, p, q)

                        # Numeric Evaluation
                        else:
                            if n < 0:
                                g_value = np.conjugate(
                                    self.gpq_numeric(-n_prime, -n, p, q))
                            else:
                                g_value = self.gpq_numeric(n_prime, n, p, q)

                        # Assign the Green's value to the Full Green's matrix G
                        g_matrix[off_p + Np + n_prime,
                                 off_q + Nq + n] = g_value

        return g_matrix


class FullWaveAnalytically(MultilayeredLossyMedium):
    """ This class contains the Analytically Green's matrices of the Multilayered Lossy system """

    # def __init__(self, conductors_list, frequency):
    #     super().__init__(conductors_list, frequency)

    # Green's function Gchat^{p,q}_{n',n} [np.array]
    # Equation (B.23) [1]
    def ghatc_pq(self, n_prime, n, p, q):
        """ Computes the full Green's matrix for the given set of conductors """

        # Auxiliary vector of distances and angles [np.array]
        dqp_dict = AuxiliaryGeometry().distance_vector_dqp(
            self.conductor_surfaces, p, self.conductor_surfaces, q)

        ap = self.conductor_surfaces[p]['radius']
        aq = self.conductor_surfaces[q]['radius']
        dqp = dqp_dict['norm']
        theta_qp = dqp_dict['angle']
        khat = self.khat[0]

        # 1st case: |a_q| < |a_p - d_qp|
        if aq < abs(ap - dqp):

            # Equation B.30 [1]
            if abs(dqp) > ap:
                bessel_1st = jv(n, khat * aq)
                hankel = hankel2(n - n_prime, khat * dqp)
                bessel_2nd = jv(-n_prime, khat * ap)
                exp_1st = np.exp(1j * n * (np.pi + theta_qp))
                exp_2nd = np.exp(-1j * n_prime * theta_qp)
                gchat_pq = 1 / 4 / 1j * bessel_1st * hankel * bessel_2nd * exp_1st * exp_2nd

            # Equation B.31 [1]
            elif abs(dqp) < ap:
                bessel_1st = jv(n, khat * aq)
                hankel = hankel2(n_prime, khat * ap)
                bessel_2nd = jv(n_prime - n, khat * dqp)
                exp_1st = np.exp(-1j * (n_prime - n) * theta_qp)
                gchat_pq = 1 / 4 / 1j * bessel_1st * hankel * bessel_2nd * exp_1st

        # 2nd case: |a_q| > a_p + d_qp
        elif aq > (ap + dqp):

            # Equation B.34 [1]
            if abs(dqp) > ap:
                hankel = hankel2(-n, khat * aq)
                bessel_1st = jv(-n + n_prime, khat * dqp)
                bessel_2nd = jv(n_prime, khat * ap)
                exp_1st = np.exp(-1j * n_prime * theta_qp)
                exp_2nd = np.exp(1j * n * (np.pi + theta_qp))
                gchat_pq = 1 / 4 / 1j * hankel * bessel_1st * bessel_2nd * exp_1st * exp_2nd

            # Equation B.35 [1]
            elif abs(dqp) < ap:
                hankel = hankel2(-n, khat * aq)
                bessel_1st = jv(-n_prime, khat * ap)
                bessel_2nd = jv(n - n_prime, khat * dqp)
                exp_1st = np.exp(1j * (n - n_prime) * theta_qp)
                gchat_pq = 1 / 4 / 1j * hankel * bessel_1st * bessel_2nd * exp_1st

        # 3rd case: |a_q| = |a_p - d_qp|
        # Green's Quasi-static formulation
        else:
            # 1st sub-case: alfa == 1 and n == 0
            # Equation (52) [2]
            if n == 0:
                if n_prime == 0:
                    gchat_pq = 1 / 2 / np.pi * np.log(ap)

                else:
                    gchat_pq = 0

            # 2st sub-case: alfa == 1 and n != 0
            # Equation (53) [2]
            elif n != 0:
                if n_prime == n:
                    gchat_pq = - 1 / 4 / np.pi / np.abs(n)

                else:
                    gchat_pq = 0

        return gchat_pq

    # Green's Matrix Gchat [np.array]
    # Equation (B.23) [1]
    def ghatc(self):
        """ Computes the full Green's matrix for the given set of conductors """

        # Define the full Green's matrix
        gchat = np.zeros((self.N, self.N), dtype=complex)

        # Create an offset array to keep track of the starting index for each block
        offsets = np.cumsum(
            [0] + [2 * Np['fourier_order'] + 1 for Np in self.conductor_surfaces[:-1]])

        # Number of conductor surfaces [int]
        ns = len(self.conductor_surfaces)

        for p in range(ns):
            # Offset for the p-th conductor [int]
            off_p = offsets[p]

            # Number of surface points for the p-th surface [int]
            Np = self.conductor_surfaces[p]['fourier_order']

            for q in range(ns):
                # Offset for the q-th conductor [int]
                off_q = offsets[q]

                # Number of surface points for the q-th surface [int]
                Nq = self.conductor_surfaces[q]['fourier_order']

                for n_prime in range(-Np, Np + 1):
                    for n in range(-Nq, Nq + 1):

                        if n < 0:
                            g_value = np.conjugate(
                                self.ghatc_pq(-n_prime, -n, p, q))

                        else:
                            g_value = self.ghatc_pq(n_prime, n, p, q)

                        # Assign the Green's value to the Full Green's matrix G
                        gchat[off_p + Np + n_prime, off_q + Nq + n] = g_value

        return gchat

    # Green's function G0hat^{q}_{n',n} [np.array]
    # Green's function G0til^{q}_{n',n} [np.array]
    # Equation (B.38) [1]
    # Equation (B.42) [1]
    def g0_q(self, n_prime, n, q):
        """ Computes the full Green's matrix for the given set of conductors """

        # Derivative of the Hankel function of the second kind
        def deriv_hankel2(n, z):
            """ Calculates the derivative of the Hankel function of the second kind. """
            hn = hankel2(n, z)
            hn_minus_1 = hankel2(n-1, z)
            return hn_minus_1 - (n / z) * hn

        # Auxiliary vector of distances [np.array]
        # Equation (B.39) - PAG. 135 [1]
        dqp_dict = AuxiliaryGeometry().distance_vector_dqp(
            self.hole_surfaces, 0, self.conductor_surfaces, q)

        a_hat = self.hole_surfaces[0]['radius']
        aq = self.conductor_surfaces[q]['radius']
        dqp = dqp_dict['norm']
        theta_qp = dqp_dict['angle']
        khat = self.khat[0]

        # Bessel and Hankel functions [np.array]
        bessel_1st = jv(n, khat * aq)
        hankel_1st = hankel2(n_prime, khat * a_hat)
        hankel_2nd = deriv_hankel2(n_prime, khat * a_hat)
        bessel_2nd = jv(n_prime - n, khat * dqp)

        # Exponencial term [np.array]
        exp_term = np.exp(-1j * (n_prime - n) * theta_qp)

        # Equation (B.38) [1]
        g0hat_q = 1 / 4 / 1j * bessel_1st * hankel_1st * bessel_2nd * exp_term

        # Equation (B.42) [1]
        g0til_q = khat / 4 / 1j * bessel_1st * hankel_2nd * bessel_2nd * exp_term

        return g0hat_q, g0til_q

    # Green's Matrix [Gchat] [np.array]
    # Equation (B.23) [1]
    def ghat_c(self):
        """ Computes the full Green's matrix for the given set of conductors """

        # Define the full Green's matrix
        gchat = np.zeros((self.N, self.N), dtype=complex)

        # Create an offset array to keep track of the starting index for each block
        offsets = np.cumsum(
            [0] + [2 * Np['fourier_order'] + 1 for Np in self.conductor_surfaces[:-1]])

        # Number of conductor surfaces [int]
        ns = len(self.conductor_surfaces)

        for p in range(ns):
            # Offset for the p-th conductor [int]
            off_p = offsets[p]

            # Number of surface points for the p-th surface [int]
            Np = self.conductor_surfaces[p]['fourier_order']

            for q in range(ns):
                # Offset for the q-th conductor [int]
                off_q = offsets[q]

                # Number of surface points for the q-th surface [int]
                Nq = self.conductor_surfaces[q]['fourier_order']

                for n_prime in range(-Np, Np + 1):
                    for n in range(-Nq, Nq + 1):

                        if n < 0:
                            g_value = np.conjugate(
                                self.ghatc_pq(-n_prime, -n, p, q))

                        else:
                            g_value = self.ghatc_pq(n_prime, n, p, q)

                        # Assign the Green's value to the Full Green's matrix G
                        gchat[off_p + Np + n_prime, off_q + Nq + n] = g_value

        return gchat

    # Green's Matrix [G0hat] [np.array]
    # Green's Matrix [G0til] [np.array]
    # Equation (3.16) - PAGE 62 [1]
    def g0_matrices(self):
        """ Computes the full Green's matrix for the given set of conductors """

        # Define the full Green's matrix
        g0hat_list = []
        g0til_list = []

        # Number of surface points for the p-th surface hole [int]
        Nh = self.hole_surfaces[0]['fourier_order']

        # Number of conductor surfaces [int]
        for p, item in enumerate(self.conductor_surfaces):

            # Number of surface points for the p-th surface conductor [int]
            Np = item['fourier_order']

            # Green's matrix G^(p, q) [np.array]
            g0hat_p = np.zeros((2*Nh+1, 2*Np+1), dtype=complex)
            g0til_p = np.zeros_like(g0hat_p)

            for n_prime in range(-Nh, Nh + 1):
                for n in range(-Np, Np + 1):

                    # Assign the Green's value to the Full Green's matrices G0
                    g0hat_p[n_prime + Nh, n +
                            Np] = self.g0_q(n_prime, n, p)[0]
                    g0til_p[n_prime + Nh, n +
                            Np] = self.g0_q(n_prime, n, p)[1]

            # Append Green's matrices G0 [np.array]
            g0hat_list.append(g0hat_p)
            g0til_list.append(g0til_p)

        return np.block([g0hat_list]), np.block([g0til_list])

    # Transformation matrix [T]
    # Equation (3.24) [1]
    def t_matrix(self):
        """ Computes the Transformation matrix T.

        Transformation matrix T which maps the equivalent currents Js(p)(θp)
        on conductor boundaries to the equivalent current Jbs(θˆ) on the hole boundary.

        """

        ahat = self.hole_surfaces[0]['radius']
        g0hat = self.g0_matrices()[0]
        g0til = self.g0_matrices()[1]
        d2 = self.d_matrices()[1]

        return 2 * np.pi * ahat * (g0til - d2 @ g0hat)

    # PSI matrix [Ψ] [np.array]
    # Equation (3.41) [1]
    def psi(self):
        """ Ψ is the Green's matrix of a lossy layered medium with holes inside it. """

        # Green's matrices [np.array]
        mu_gc = self.mu_hat[0] * self.ghat_c()
        mu_g0 = self.mu_hat[0] * self.g0_matrices()[0]
        gg = np.zeros((self.Nhat, self.Nhat), dtype=complex)

        # Matrices Products
        h_d1 = self.hhat() @ self.d_matrices()[0]
        gg_ys = np.eye(self.Nhat) + mu_0 * (gg @ self.yhats())

        # LU Decomposition of the matrix gg_ys
        lu, piv = lu_factor(gg_ys)
        gg_ys_inv = lu_solve((lu, piv), np.eye(gg_ys.shape[0]))
        bracket = mu_0 * gg_ys_inv @ gg @ self.t_matrix() - mu_g0

        return 1 / mu_0 * (h_d1 @ bracket + mu_gc)


class FullWaveNumerically(MultilayeredLossyMedium, AuxiliaryGeometry):
    """ This class contains the Numerically Green's matrices of the Multilayered Lossy system """

    # Equation (41) [2]
    def func_gc_theta(self, theta, n_prime, n, p, q):
        """ This function calculates the function f_n(theta) """

        def integrand(theta_prime, theta, n_prime, n, p, q):
            """ This function calculates the integrand of the function f_n(theta) """

            gc = 0
            rp = self.contour_vector_position(
                self.conductor_surfaces, theta, p)
            rq = self.contour_vector_position(
                self.conductor_surfaces, theta_prime, q)
            R = np.linalg.norm(rp - rq)
            khat_R = self.khat[0] * R
            exponent = np.exp(1j * (n * theta_prime - n_prime * theta))

            if R != 0:
                gc = hankel2(0, khat_R) * exponent

            return gc

        real_gc, _ = quad(lambda theta_prime: integrand(
            theta_prime, theta, n_prime, n, p, q).real, 0, 2 * np.pi)

        imag_gc, _ = quad(lambda theta_prime: integrand(
            theta_prime, theta, n_prime, n, p, q).imag, 0, 2 * np.pi)

        return real_gc + 1j * imag_gc

    # Equation (41) [2]
    def func_g0hat_theta(self, theta, n_prime, n, q):
        """ This function calculates the function f_n(theta) """

        def integrand(theta_prime, theta, n_prime, n, q):
            """ This function calculates the integrand of the function f_n(theta) """

            g0_matrix = 0
            r_hat = self.contour_vector_position(
                self.hole_surfaces, theta, p=0)
            rq = self.contour_vector_position(
                self.conductor_surfaces, theta_prime, q)
            R = np.linalg.norm(rq - r_hat)
            khat_R = self.khat[0] * R
            exponent = np.exp(1j * (n * theta_prime - n_prime * theta))

            if R != 0:
                g0_matrix = hankel2(0, khat_R) * exponent

            return g0_matrix

        real_g0, _ = quad(lambda theta_prime: integrand(
            theta_prime, theta, n_prime, n, q).real, 0, 2 * np.pi)

        imag_g0, _ = quad(lambda theta_prime: integrand(
            theta_prime, theta, n_prime, n, q).imag, 0, 2 * np.pi)

        return real_g0 + 1j * imag_g0

    # Equation (40) [2]
    def ghatc_pq(self, n_prime, n, p, q):
        """ Calculate the Greens' matrix Gc_hat_{n',n}^{p,q} """

        def integrand(theta, n_prime, n, p, q):
            """ Calculate the integrand of Greens' matrix G_{n',n}^{p,q} """
            return self.func_gc_theta(theta, n_prime, n, p, q)

        real_gchat, _ = quad(lambda theta: integrand(
            theta, n_prime, n, p, q).real, 0, 2 * np.pi)

        imag_gchat, _ = quad(lambda theta: integrand(
            theta, n_prime, n, p, q).imag, 0, 2 * np.pi)

        return (real_gchat + 1j * imag_gchat) / (2 * np.pi)**2 / (4 * 1j)

    # Green's Matrix Gchat [np.array]
    # Equation (B.23) [1]
    def ghatc(self):
        """ Computes the full Green's matrix for the given set of conductors """

        # Define the full Green's matrix
        gchat = np.zeros((self.N, self.N), dtype=complex)

        # Create an offset array to keep track of the starting index for each block
        offsets = np.cumsum(
            [0] + [2 * Np['fourier_order'] + 1 for Np in self.conductor_surfaces[:-1]])

        # Number of conductor surfaces [int]
        ns = len(self.conductor_surfaces)

        for p in range(ns):
            # Offset for the p-th conductor [int]
            off_p = offsets[p]

            # Number of surface points for the p-th surface [int]
            Np = self.conductor_surfaces[p]['fourier_order']

            for q in range(ns):
                # Offset for the q-th conductor [int]
                off_q = offsets[q]

                # Number of surface points for the q-th surface [int]
                Nq = self.conductor_surfaces[q]['fourier_order']

                for n_prime in range(-Np, Np + 1):
                    for n in range(-Nq, Nq + 1):

                        if n < 0:
                            g_value = np.conjugate(
                                self.ghatc_pq(-n_prime, -n, p, q))

                        else:
                            g_value = self.ghatc_pq(n_prime, n, p, q)

                        # Assign the Green's value to the Full Green's matrix G
                        gchat[off_p + Np + n_prime, off_q + Nq + n] = g_value

        return gchat

    # Equation (B.40) [1]
    def g0til_q(self, n_prime, n, q):
        """ This function calculates the function f_n(theta) """

        # Radius of the conductor [float]
        aq = self.conductor_surfaces[q]['radius']

        # Radius of the hole [float]
        ahat = self.hole_surfaces[0]['radius']

        # Auxiliary variables distances and angles [int]
        def distance(a_q, theta_q, rho_hat, theta_hat):
            return np.sqrt(a_q**2 + rho_hat**2 - 2 * a_q * rho_hat * np.cos(theta_q - theta_hat))

        # Integrand function
        def integrand(theta_q, theta_hat, rho_hat, a_q, n, n_prime):
            R = distance(a_q, theta_q, rho_hat, theta_hat)
            khat_R = self.khat[0] * R
            exp_term = np.exp(1j * (n * theta_q - n_prime * theta_hat))
            hankel_term = hankel2(0, khat_R)
            return hankel_term * exp_term

        # Função principal para calcular as partes real e imaginária de G_0
        def calculate_g0(a_q, n, n_prime, rho_hat):
            # Parte real da integral
            real_integral, _ = dblquad(
                lambda theta_hat, theta_q: integrand(
                    theta_q, theta_hat, rho_hat, a_q, n, n_prime).real,
                0, 2 * np.pi, lambda _: 0, lambda _: 2 * np.pi
            )

            # Parte imaginária da integral
            imag_integral, _ = dblquad(
                lambda theta_hat, theta_q: integrand(
                    theta_q, theta_hat, rho_hat, a_q, n, n_prime).imag,
                0, 2 * np.pi, lambda _: 0, lambda _: 2 * np.pi
            )

            return real_integral, imag_integral

        # Função para derivar numericamente G_0 com relação a rho_hat
        def derivative_g0(a_q, n, n_prime, rho_hat, h=1e-5):
            # Aproximação por diferenças finitas centradas
            g0_plus_real, g0_plus_imag = calculate_g0(
                a_q, n, n_prime, rho_hat + h)
            g0_minus_real, g0_minus_imag = calculate_g0(
                a_q, n, n_prime, rho_hat - h)

            deriv_real = (g0_plus_real - g0_minus_real) / (2 * h)
            deriv_imag = (g0_plus_imag - g0_minus_imag) / (2 * h)

            return deriv_real, deriv_imag

        g0_real, g0_imag = derivative_g0(aq, n, n_prime, ahat)

        return 1j / 16 / np.pi**2 * g0_real + 1j * g0_imag

    # Equation (B.37) [1]
    def g0hat_q(self, n_prime, n, q):
        """ Calculate the Greens' matrix G0_hat_{n',n}^{q} """

        def integrand(theta, n_prime, n, q):
            """ Calculate the integrand of Greens' matrix G_{n',n}^{p,q} """
            return self.func_g0hat_theta(theta, n_prime, n, q)

        re_g0hat, _ = quad(lambda theta: integrand(
            theta, n_prime, n, q).real, 0, 2 * np.pi)

        imag_g0hat, _ = quad(lambda theta: integrand(
            theta, n_prime, n, q).imag, 0, 2 * np.pi)

        return 1j * (re_g0hat + 1j * imag_g0hat) / (2 * np.pi)**2 / 4
