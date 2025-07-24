import numpy as np
from scipy.special import jv, jvp
from scipy.constants import mu_0
from scipy.linalg import lu_factor, lu_solve

from mom_so.mtl import MulticonductorTransmissionLine


class WiresHomogeneousMedia(MulticonductorTransmissionLine):
    """
    This class contains the analytical formulation of the system. 

    References:
    [1] Clayton R. Paul, "Introduction to Electromagnetic Compatibility", 2nd Edition, Wiley, 2007.
        4.2.2 Per-Unit-Length Inductance and Capacitance for Wire-Type Lines
        5.2.1 Wide-Separation Approximations for Wires in Homogeneous Media 
    
    """

    def __init__(self, mtl):
        super().__init__(mtl)

        # Kelvin Functions
        self.kelvin_exp = np.exp(1j * 3 * np.pi / 4)

    def ber(self, xi):
        """ Kelvin ber(xi) function """
        return np.real(jv(0, xi * self.kelvin_exp))

    def bei(self, xi):
        """ Kelvin bei(xi) function """
        return np.imag(jv(0, xi * self.kelvin_exp))

    def ber_prime(self, xi):
        """ Kelvin ber'(xi) derivative function """
        return np.real(self.kelvin_exp * jvp(0, xi * self.kelvin_exp, 1))

    def bei_prime(self, xi):
        """ Kelvin bei'(xi) derivative function """
        return np.imag(self.kelvin_exp * jvp(0, xi * self.kelvin_exp, 1))

    def bifilar_pul_series_impedance(self, f):
        """
        This function calculates the series resistance of the system using 
        the high frequency approximation.

        Returns:
        tuple: A tuple containing the high frequency resistance, external inductance, 
        and matrix impedance.
        """
        # Skin Depth [np.array]
        w = 2 * np.pi * f
        delta = np.sqrt(1 / (w / 2 * mu_0 * self.sigma))

        # Surface Resistance [np.array]
        Rs = 1 / (self.sigma * delta) 

        # Outer Radii of the conductors [np.array]
        ap = np.array([cp['radius'][1] for cp in self.mtl])

        # Matrix Distance [np.array]
        D = self.D_pq

        # High Frequency Resistance and External Inductance [np.array]
        n = len(self.mtl) - 1
        Rhf = np.zeros((n, n))
        Lext = np.zeros_like(Rhf)
        Zi = np.zeros_like(Rhf, dtype=complex)

        # Constant Term and Bessel argument
        Xi = np.sqrt(2) * ap / delta
        constant_term = 1 / (np.sqrt(2) * np.pi * ap * self.sigma * delta)

        # 1st solution: High Frequency Approximation
        # These formulas account for proximity effect
        # only at high frequencies.These formulas
        # account for proximity effect only at high
        # frequencies.
        p = 0 

        # Common fraction term
        s_2rw = D[p][p+1] / 2 / ap[p]

        # Surface resistance
        Rs_pia = Rs[p] / np.pi / ap[p]

        # High Frequency Resistance (Ω/m)
        # Equation (2.64) [1]
        Rhf[p] = Rs_pia * s_2rw / np.sqrt(s_2rw ** 2 - 1)

        # Exact External Inductance (H/m)
        # Equation (2.65) [1]
        # Equation (4.41) [2]
        Lext[p] = mu_0 / np.pi * np.arccosh(s_2rw)

        # 2nd solution: Internal Impedance Matrix, z_int (Ω/m)
        # These formulas captures skin effect, but proximity
        # effect is neglected
        # Equation (2.67) [1]
        ber_bei = self.ber(Xi[p]) + 1j * self.bei(Xi[p])
        beip_berp = self.bei_prime(Xi[p]) - 1j * self.ber_prime(Xi[p])

        # Internal Impedance Matrix, Zi (Ω/m)
        Zi[p] = constant_term[p] * ber_bei / beip_berp

        # Matrix Impedance, Zs (Ω/m)
        # Equation (2.68) [1]
        Zs = 2 * Zi + 1j * w * Lext

        return Zs, Rhf

    def bifilar_pul_inductance_capacitance(self):
        """
        Calcula a capacitância por unidade de comprimento para a linha bifilar.

        Este método calcula tanto a fórmula exata (4.39) quanto a aproximada (4.21)
        do livro "Introduction to Electromagnetic Compatibility" de Clayton Paul.

        Fórmula Exata (4.39):
        c = (2 * pi * epsilon) / arccosh((s^2 - r_w1^2 - r_w2^2) / (2 * r_w1 * r_w2))

        Fórmula Aproximada (4.21):
        c = (2 * pi * epsilon) / ln(s^2 / (r_w1 * r_w2))

        Retorna:
            dict: Um dicionário contendo a capacitância 'exact' e 'approximate'
                  em Farads por metro (F/m).
        """
        # Assume que o meio externo é homogêneo (ex: ar).
        # self.epsilon_out é um array; para uma linha bifilar em um único meio,
        # o primeiro elemento é suficiente.
        epsilon = self.epsilon_out[0]

        # s: distância entre os centros dos condutores.
        # Da classe pai, D é a matriz de distância d_pq.
        s = self.D_pq[0, 1]

        # r_w1, r_w2: raios dos dois condutores.
        radii = np.array([conductor['radius'][1] for conductor in self.mtl])
        r_w1, r_w2 = radii
        assert len(radii) == 2, "A linha bifilar deve ter exatamente dois condutores."
        
        # --- Cálculo da Capacitância Aproximada (Eq. 4.21 [1]) ---
        pi2_epsilon = 2 * np.pi * epsilon
        den_approx = np.log((s**2) / (r_w1 * r_w2))
        capacitance_approx = pi2_epsilon / den_approx

        # --- Cálculo da Capacitância Exata (Eq. 4.39 [1]) ---
        arg_arccosh = (s**2 - r_w1**2 - r_w2**2) / (2 * r_w1 * r_w2) #
        den_exact = np.arccosh(arg_arccosh) #
        capacitance_exact = pi2_epsilon / den_exact #

        # --- Cálculo da Indutância  Aproximada (Eq. -) ---
        inductance_approx = mu_0 * epsilon / capacitance_approx

        # --- Cálculo da Indutância Exata (Eq. 4.41 [1]) ---
        inductance_exact = mu_0 * epsilon / capacitance_exact

        # Retorna um dicionário com ambos os resultados
        return {
            'capacitante': {
                'exact': capacitance_exact,
                'approximate': capacitance_approx
            },
            'indutância': {
                'exact': inductance_exact,
                'approximate': inductance_approx
            }
        }

    def n_wires_inductance_matrix(self):
        """
        This function calculates the series resistance of the system using 
        the high frequency approximation.

        Returns:
        tuple: A tuple containing the high frequency resistance, external inductance, 
        and matrix impedance.
        """
        # Conductors number
        n = len(self.mtl) - 1
        # High Frequency Resistance and External Inductance [np.array]
        Lext = np.zeros((n, n))        
       
        rw = np.array([cp['radius'][1] for cp in self.mtl])  # Outer Radii of the conductors [np.array]
        D = self.D_pq       # Matrix Distance [np.array]
        mu_2pi = self.mu / (2 * np.pi)
        rw0 = rw[-1]  # Reference conductor radius (outer radius of the last conductor)

        for i in range(n):
            di0 = D[i, n] # Distance to conductor i to the conductor n (reference)
            for j in range(n):
                dj0 = D[j, n] # Distance to conductor j to the conductor n (reference)
                dij = D[i, j]  # Distance between conductor i and j

                if i == j:
                    # Self Inductance (H/m)
                    Lext[i, j] = mu_2pi[i] * np.log(di0 ** 2 / (rw0 * rw[i]))
                else:
                    # Mutual Inductance (H/m)
                    Lext[i, j] = mu_2pi[i] * np.log(di0 * dj0 / (rw0 * dij))

        return Lext

    def n_wires_capacitance_matrix(self, L):
        """
        Calcula a matriz de capacitância por unidade de comprimento (C) a partir da
        matriz de indutância (L) para um meio homogêneo, seguindo a Eq. 5.24a.

        A fórmula implementada é: C = μ * ε * L⁻¹

        A inversa de L (L⁻¹) é calculada de forma numericamente estável
        resolvendo o sistema L @ X = I, onde I é a matriz identidade.

        Args:
            inductance_matrix (np.ndarray): A matriz de indutância L (n x n).
            permeability (float): A permeabilidade magnética do meio (μ).
            permittivity (float): A permissividade elétrica do meio (ε).

        Returns:
            np.ndarray: A matriz de capacitância C (n x n) resultante.
            
        Raises:
            ValueError: Se a matriz de indutância não for quadrada.
        """
        assert L.shape[0] == L.shape[1], "A matriz de indutância deve ser quadrada."
        
        # Identity matrix for the number of conductors    
        I = np.eye(L.shape[0])
        
        # Assuming a homogeneous medium, use the first value 
        mu = self.mu[0]  
        epsilon = self.epsilon_out[0]

        return mu * epsilon * lu_solve(lu_factor(L), I)

