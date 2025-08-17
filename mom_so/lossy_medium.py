# import numpy as np
# from scipy.constants import mu_0, epsilon_0
# from scipy.special import jv, jvp, h1vp, h2vp
# from scipy.linalg import lu_factor, lu_solve, inv
# from mtl_data.mtl import MulticonductorTransmissionLine as MTL

# class MultilayeredLossyMedium(UndergroundSystem):
#     """ This class contains the frequency-dependent parameters of the system. """

#     def __init__(self, conductors_list, frequency):
#         super().__init__(conductors_list)

#         # Angular frequency [float]
#         self.w = 2 * np.pi * frequency

#         ### CONDUCTORS PARAMETERS ###

#         # Permeability [np.array]
#         self.mu = np.array([mu_0 * c['relative_permeability']
#                             for c in self.subconductors])

#         # Permittivity [np.array]
#         epsilon = np.array([epsilon_0 * c['relative_permittivity']
#                             for c in self.subconductors])

#         # Conductivity [np.array]
#         sigma = np.array([c['conductivity'] for c in self.subconductors])

#         # Wavenumber [float]
#         self.k = np.sqrt(self.w * self.mu * (self.w * epsilon - 1j * sigma))

#         ### HOLES PARAMETERS ###

#         # Medium permeability [np.array]
#         self.mu_hat = np.array([mu_0 * c['relative_permeability']
#                                 for c in self.holes])

#         # Permittivity [np.array]
#         epsilon_hat = np.array([epsilon_0 * c['relative_permittivity']
#                                 for c in self.holes])

#         # Wavenumber [float]
#         self.khat = self.w * np.sqrt(self.mu_hat * epsilon_hat)

#         ### GROUND PARAMETERS ###

#         # Ground conductivity [float]
#         sigma_g = np.array([c['conductivity'] for c in self.ground])

#         # Ground permittivity [float]
#         epsilon_g = np.array([c['permittivity'] for c in self.ground])

#         # Ground wavenumber [float]
#         self.kg = np.sqrt(self.w * mu_0 * (self.w * epsilon_g - 1j * sigma_g))

#     # Conductor Surface admittance operator [np.array]
#     # Equation (2.20) [1]
#     def ys_entries(self, n, p):
#         """ Calculates the surface admittance operator for a solid conductor, Yn(p). """

#         # Conductor Surface radius [float]
#         ap = self.conductor_surfaces[p]['radius']

#         # Conductor permeability [float]
#         mu = self.mu[p]

#         # Conductor wave number \times Conductor radius [float]
#         kap = self.k[p] * ap

#         # Hole wave number \times Conductor radius [float]
#         khat_ap = self.khat[0] * ap

#         # Hole permeability [float]
#         muhat = self.mu_hat[0]

#         # Bessel and Bessel derivative functions [np.array]
#         n = np.abs(n)
#         jn_kap = jv(n, kap)
#         jnp_kap = jvp(n, kap)
#         jn_khatap = jv(n, khat_ap)
#         jnp_khatap = jvp(n, khat_ap)

#         frac_1 = kap * jnp_kap / jn_kap / mu
#         frac_2 = khat_ap * jnp_khatap / jn_khatap / muhat

#         return 2 * np.pi / 1j / self.w * (frac_1 - frac_2)

#     # Ys Matrix [np.array]
#     # Equation (3.1) [1]
#     def ys(self):
#         """ Computes the full Green's matrix for the given set of conductors.
         
#         The surface admittance operator Ys was derived in Chapter 2, and its
#         diagonal entries are given by (2.20) with k0 replaced by 
#         k_hat = ω√(µ_hat . ε_hat), which is the wavenumber inside the hole. 
        
#         """

#         blocks = []

#         for p, conductor in enumerate(self.conductor_surfaces):

#             # Number of surface points for the p-th surface [int]
#             Np = conductor['fourier_order'] # pylint: disable=invalid-name

#             for n in range(-Np, Np + 1):
#                 ynp = np.array([[self.ys_entries(n, p)]])
#                 blocks.append(ynp)

#         return self.create_diagonal_matrix(blocks)

#     # Green's function G0hat^{q}_{n',n} [np.array]
#     # Equation (B.47) [1]
#     def hhat_p(self, n_prime, n, p):
#         """ Computes the full Green's matrix for the given set of conductors """

#         # Auxiliary vector of distances [np.array]
#         # Equation (B.45) - PAG. 137 [1]
#         dqp, xqp, yqp, theta_qp = self.distance_vector_dqp(self.conductor_surfaces, p, self.hole_surfaces, 0)

#         ap = self.conductor_surfaces[p]['radius']
#         khat = self.khat[0]

#         # Bessel and Hankel functions [np.array]
#         bessel_1st = jv(n_prime - n, khat * dqp)
#         bessel_2nd = jv(n_prime, khat * ap)

#         # Exponential term [np.array]
#         exp_term = np.exp(1j * (n - n_prime) * theta_qp)

#         # Equation (B.47) [1]
#         hhat_p = bessel_1st * bessel_2nd * exp_term

#         if n < 0:
#             hhat_p = hhat_p * (-1)**n

#         return hhat_p

#     # H_hat Matrix [np.array]
#     # Equation (3.38) [1]
#     def hhat(self):
#         """ Computes the full Green's matrix for the given set of conductors """

#         # Define the full Green's matrix
#         hhat_t_list = []

#         # Number of surface points for the p-th surface hole [int]
#         Nh = self.hole_surfaces[0]['fourier_order'] # pylint: disable=invalid-name

#         # Number of conductor surfaces [int]
#         for p, conductor in enumerate(self.conductor_surfaces):

#             # Number of surface points for the p-th surface [int]
#             Np = conductor['fourier_order'] # pylint: disable=invalid-name

#             # Green's matrix H_hat [np.array]
#             hhat_p = np.zeros((2*Nh+1, 2*Np+1), dtype=complex)

#             for n in range(-Nh, Nh + 1):
#                 for n_prime in range(-Np, Np + 1):

#                     # Assign the Green's value to the Full Green's matrix G
#                     hhat_p[n + Nh, n_prime + Np] = self.hhat_p(n_prime, n, p)

#             # Append Green's matrices H_hat [np.array]
#             hhat_t_list.append(hhat_p)

#         return np.block([hhat_t_list]).T

#     # D1 Matrix [np.array]
#     # Equation (3.15) [1]
#     def d_matrices(self):
#         """ Computes the full Green's matrix for the given set of conductors """

#         # Green's matrix H_hat [np.array]
#         d1 = np.zeros((self.Nhat, self.Nhat), dtype=complex)
#         d2 = np.zeros_like(d1)

#         # Number of surface points for the p-th surface hole [int]
#         Nh = self.hole_surfaces[0]['fourier_order'] # pylint: disable=invalid-name

#         # Wave number \times Hole radius [float]
#         ka_hat = self.khat[0] * self.hole_surfaces[0]['radius']

#         for n in range(-Nh, Nh + 1):

#             # Bessel and Bessel derivative functions [np.array]
#             bessel = jv(np.abs(n), ka_hat)
#             bessel_prime = jvp(np.abs(n), ka_hat)

#             # Assign the Green's value to the Full Green's matrix G
#             d1[n+Nh, n+Nh] = 1 / bessel

#             # Assign the Green's value to the Full Green's matrix G
#             d2[n+Nh, n+Nh] = self.khat[0] * bessel_prime / bessel

#         return d1, d2

#     # Hole Surface admittance operator [np.array]
#     # Equation (2.20) [1]
#     def yhats_entries(self, n):
#         """ Calculates the surface admittance operator for a solid conductor, Yn(p). """

#         # Hole Surface radius [float]
#         a_hat = self.hole_surfaces[0]['radius']

#         # Hole permeability [float]
#         mu_hat = self.mu_hat[0]

#         # Hole wave number \times Hole radius [float]
#         k_ahat = self.khat[0] * a_hat

#         # Ground wave number \times Hole radius [float]
#         kg_ahat = self.kg[0] * a_hat

#         # Bessel and Bessel derivative functions [np.array]
#         n = np.abs(n)
#         jn_kga = jv(n, kg_ahat)
#         jnp_kga = jvp(n, kg_ahat)
#         jn_ka = jv(n, k_ahat)
#         jnp_ka = jvp(n, k_ahat)

#         frac_1 = self.kg[0] * jnp_kga / jn_kga / mu_0
#         frac_2 = self.khat[0] * jnp_ka / jn_ka / mu_hat

#         return 2 * np.pi * a_hat * (frac_1 - frac_2)

#     # Yhat_s Matrix [np.array]
#     # Equation (3.23) [1]
#     def yhats(self):
#         """ Computes the full Green's matrix for the given set of conductors """

#         yhat_s = np.zeros((self.Nhat, self.Nhat), dtype=complex)

#         # Number of surface points for the p-th surface hole [int]
#         Nh = self.hole_surfaces[0]['fourier_order'] # pylint: disable=invalid-name

#         for n in range(-Nh, Nh + 1):

#             # Assign the Green's value to the Full Green's matrix G
#             yhat_s[n+Nh, n+Nh] = self.yhats_entries(n)

#         return yhat_s

#     # Matrix Z [np.array]
#     # Equation (3.43 | 3.44) [1]
#     def z_partial(self, psi):
#         """
#         This function calculates the matrix Z.

#         The matrix Z is the impedance matrix of the system.

#         Returns:
#         numpy.ndarray: The matrix Z.
#         """

#         jw_u0 = 1j * self.w * mu_0
#         ys = self.ys()
#         u = np.array([[1]])

#         # Perform LU factorization of the matrix M
#         matrix = np.eye(self.N) - jw_u0 * np.dot(ys, psi)
#         lu, piv = lu_factor(matrix)

#         # Solve the linear system Mx = b for Ys
#         solution = lu_solve((lu, piv), ys @ u)

#         return u.T @ solution
