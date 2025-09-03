"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

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

[4] PAUL, Clayton R. Analysis of multiconductor transmission lines. 2. ed. Hoboken,
    N.J.: John Wiley & Sons, Inc., c2008.

"""

import numpy as np
from scipy.special import iv, kv, jv, jvp
from scipy.linalg import lu_factor, lu_solve
from mtl_main.source import MulticonductorTransmissionLine

# 4.2.2 Per-Unit-Length Inductance and Capacitance for Wire-Type Lines [4]
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
        self.kelvin_exp = np.exp(1j * 3 * np.pi / 4)
        self.s = self.D_pq[0, 1]

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
        Calcula a impedância série do sistema bifilar usando a
        aproximação de alta frequência e funções de Bessel para a impedância interna.
        Assume um sistema com um condutor ativo (tag=1) e um de retorno (tag=0).

        Args:
            f (float): Frequência em Hertz.

        Returns:
            tuple: Uma tupla contendo a impedância série total (complexa, Zs) e
                a resistência de alta frequência (Rhf).
        """
        # 1. Obter dados dos condutores diretamente pelo seu tag
        cond_active = self.mtl[1]
        cond_ref = self.mtl[0]

        # 2. Calcular distância dinamicamente, eliminando self.D_pq
        center_active = np.array(cond_active['center_point'])
        center_ref = np.array(cond_ref['center_point'])
        D10 = np.linalg.norm(center_active - center_ref)

        # 3. Extrair parâmetros e calcular valores intermediários
        w = 2 * np.pi * f
        ap = cond_active['radius'][1]
        sigma = cond_active['conductivity']

        # Profundidade pelicular e Resistência superficial
        delta = np.sqrt(2 / (w * self.mu[0] * sigma))
        Rs = 1 / (sigma * delta)

        # 4. Calcular Resistência de Alta Frequência (Rhf) e Indutância Externa (Lext)
        # Termo comum para as equações
        s_2rw = D10 / (2 * ap)

        # Resistência de Alta Frequência (Ω/m) - Equação (2.64)
        Rs_pia = Rs / (np.pi * ap)
        Rhf_loop = Rs_pia * s_2rw / np.sqrt(s_2rw**2 - 1)

        # Indutância Externa (H/m) - Equação (2.65)
        Lext_loop = self.mu[0] / np.pi * np.arccosh(s_2rw)

        # 5. Calcular Impedância Interna (Zi) com Funções de Bessel
        # Equação (2.67)
        Xi = np.sqrt(2) * ap / delta
        constant_term = 1 / (np.sqrt(2) * np.pi * ap * sigma * delta)
        
        # As funções self.ber, self.bei, etc. devem estar definidas na classe
        ber_bei = self.ber(Xi) + 1j * self.bei(Xi)
        beip_berp = self.bei_prime(Xi) - 1j * self.ber_prime(Xi)
        Zi = constant_term * ber_bei / beip_berp

        # 6. Calcular a Impedância Série Total do Laço (Zs)
        # Equação (2.68)
        # O fator 2 em Zi assume que os condutores ativo e de retorno são idênticos.
        Zs_loop = 2 * Zi + 1j * w * Lext_loop

        return Zs_loop, Rhf_loop

    def bifilar_pul_inductance_and_capacitance(self):
        """
        Calcula a capacitância e indutância por unidade de comprimento para a
        linha bifilar, usando as fórmulas exata e aproximada.

        Returns:
            dict: Um dicionário contendo a capacitância e indutância
                ('exact', 'approximate').
        """
        # 1. Validação da estrutura de dados para o caso bifilar
        assert len(self.mtl) == 2, "Este método é específico para sistema bifilar (2 condutores)."

        # 2. Acesso direto e explícito aos dados dos condutores
        c0, c1 = self.mtl[0], self.mtl[1]
        rw0, rw1 = c0['radius'][1], c1['radius'][1]
        
        # 3. Cálculo dinâmico da distância 's'
        s = np.linalg.norm(np.array(c0['center_point']) - np.array(c1['center_point']))

        # 4. Validação do meio externo e definição de epsilon
        # A fórmula assume um meio dielétrico externo único e homogêneo.
        eps_out_0, eps_out_1 = self.epsilon_out
        pi2e = 2 * np.pi * eps_out_0
        assert np.isclose(eps_out_0, eps_out_1), "O meio externo deve ser homogêneo (permissividade externa igual para ambos os condutores)."
        
        # 5. Cálculo da Capacitância
        den_approx = np.log((s**2) / (rw0 * rw1))
        capacitance_approx = pi2e / den_approx
        arg_arccosh = (s**2 - rw0**2 - rw1**2) / (2 * rw0 * rw1)
        den_exact = np.arccosh(arg_arccosh)
        capacitance_exact = pi2e / den_exact

        # 6. Cálculo da Indutância (válido para meio não magnético)
        inductance_approx = self.mu[0] * eps_out_0 / capacitance_approx
        inductance_exact = self.mu[0] * eps_out_0 / capacitance_exact

        # 7. Retorno dos resultados em um dicionário estruturado
        return {
            'capacitance': {
                'exact': capacitance_exact,
                'approximate': capacitance_approx,
            },
            'inductance': {
                'exact': inductance_exact,
                'approximate': inductance_approx,
            }
        }

    def n_wires_inductance_matrix(self):
        """
        Calcula a matriz de indutância externa, retornando uma matriz NumPy.
        Aproveita a garantia de que as tags são sequenciais (0, 1, 2, ...).

        Returns:
            np.ndarray: A matriz de indutância (N-1 x N-1).
        """
        # O tamanho da matriz de indutância é para os condutores ativos.
        n_plus_1 = len(self.mtl)
        Lext = np.zeros((n_plus_1-1, n_plus_1-1), dtype=float)

        # Dados do condutor de referência.
        ref_center = np.array(self.mtl[0]['center_point'])
        rw0 = self.mtl[0]['radius'][1]
        mu_2pi = self.mu[0] / (2 * np.pi)

        # As tags ativas vão de 1 a n_active.
        for i in range(1, n_plus_1):
            center_i = np.array(self.mtl[i]['center_point'])
            rw_i = self.mtl[i]['radius'][1]
            di0 = np.linalg.norm(center_i - ref_center)

            for j in range(1, n_plus_1):
                if i == j:  # Autoindutância
                    Lext[i - 1, j - 1] = mu_2pi * np.log(di0 ** 2 / (rw0 * rw_i))
                
                else:  # Indutância Mútua
                    center_j = np.array(self.mtl[j]['center_point'])
                    dj0 = np.linalg.norm(center_j - ref_center)
                    dij = np.linalg.norm(center_i - center_j)
                    Lext[i - 1, j - 1] = mu_2pi * np.log(di0 * dj0 / (rw0 * dij))

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


# 4.2.2 Per-Unit-Length Inductance and Capacitance for Wire-Type Lines [4]
class CoaxialCable(MulticonductorTransmissionLine):
    """ This class contains the analytical formulation of the system. """

    def __init__(self, mtl):
        """
        Initialize the AnalyticalFormulation class.

        Parameters:
        conductor (list): List of dictionaries containing the properties of the conductors.
        frequency (float): The frequency of the system.
        """
        super().__init__(mtl)

        # Coaxial Cable radii
        for key, conductor in self.mtl.items():
            if isinstance(key, int):  # Ensures the key is an integer
                if conductor['conductor_name'] == 'core':
                    self.a = conductor['radius'][1]
                elif conductor['conductor_name'] == 'sheath':
                    self.b, self.c = conductor['radius']

    # Equation 2.70 [1] and 4.51 [4]
    def external_inductance(self):
        """ Calculate the external inductance for a lossless coaxial cable, L'. """
        return self.mu[0] / (2 * np.pi) * np.log(self.b / self.a)

    # Equation 2.71 [1]
    def internal_impedance(self, frequency):
        """ Calculate the internal impedance of the inner conductor Za (omega). """
        # Angular frequency, rad/s [float]
        jw = 1j * 2 * np.pi * frequency

        # Propagation Constant [np.array]
        gamma = np.sqrt(jw * self.mu[0] * self.sigma)

        # Propagation Constant of the inner conductor
        gama_a = gamma[0] * self.a

        # Intrinsic Impedance of the inner conductor
        eta = np.sqrt(jw * self.mu[0] / self.sigma)[0]

        # Internal Impedance of the inner conductor
        za = eta / (2 * np.pi * self.a) * iv(0, gama_a) / iv(1, gama_a)

        return za

    # Equation 2.72 [1]
    def external_impedance(self, frequency):
        """ Calculate the external impedance of the inner conductor Zb (omega). """
        # Angular frequency, rad/s [float]
        jw = 1j * 2 * np.pi * frequency

        # Propagation Constant [np.array]
        gamma = np.sqrt(jw * self.mu[0] * self.sigma)

        # Propagation Constant of the inner conductor
        gama_b = gamma[0] * self.b
        gama_c = gamma[0] * self.c

        # Intrinsic Impedance of the inner conductor
        eta = np.sqrt(jw * self.mu[0] / self.sigma)[0]

        numerator = iv(0, gama_b) * kv(1, gama_c) + (kv(0, gama_b) * iv(1, gama_c))
        denominator = iv(1, gama_c) * kv(1, gama_b) - (iv(1, gama_b) * kv(1, gama_c))

        # External Impedance of the inner conductor
        zb = eta / (2 * np.pi * self.b) * numerator / denominator

        return zb

    # Equation (2.69) [1]
    def pul_parameters(self, frequency):
        """
        This function calculates the series resistance of the system using 
        the high frequency approximation.

        Returns:
            tuple: A tuple containing the high frequency resistance, external inductance, 
            and matrix impedance.
        """
        # Angular frequency, rad/s [float]
        jw = 1j * 2 * np.pi * frequency

        # Matrix Impedance, z (Ω/m)
        l_ext = self.external_inductance()
        za = self.internal_impedance(frequency)
        zb = self.external_impedance(frequency)
        zs = jw * l_ext + za + zb

        return zs

# Obsolete. Substitute with scc.InternalPerUnitParameters class
class Ametani(MulticonductorTransmissionLine):
    """ This class contains the analytical formulation of the system. """

    def __init__(self, mtl):
        """
        Initialize the AnalyticalFormulation class.

        Parameters:
        conductor (list): List of dictionaries containing the properties of the conductors.
        frequency (float): The frequency of the system.
        """
        super().__init__(mtl)

        # Call the function to configure the parameters
        for key, conductor in self.mtl.items():
            if isinstance(key, int):  # Ensures the key is an integer
                if conductor['conductor_name'] == 'core':
                    self.a, self.b = conductor['radius']
                elif conductor['conductor_name'] == 'sheath':
                    self.b_prime, self.c = conductor['radius']

    def parameter_m(self, frequency, mu, sigma):
        """ Calculate the parameter m for the two-layered conductor. """
        # Angular frequency, rad/s [float]
        jw = 1j * 2 * np.pi * frequency
        return np.sqrt(jw * mu * sigma)

    def impedance_two_layered_conductor(self, frequency):
        """ Calculate the impedance of a two-layered conductor. """
        # Angular frequency, rad/s [float]
        jw = 1j * 2 * np.pi * frequency

        # Intermediate variables
        m1 = self.parameter_m(frequency, self.mu[0], self.sigma[0])
        m2 = self.parameter_m(frequency, self.mu[1], self.sigma[1])
        x1 = m1 * self.a
        x2 = m1 * self.b
        x3 = m2 * self.b_prime
        x4 = m2 * self.c

        # Intermediate variables A, B, E, F, R
        aa = kv(1, x1) * iv(0, x2) + iv(1, x1) * kv(0, x2)
        bb = kv(1, x1) * iv(1, x2) - iv(1, x1) * kv(1, x2)
        ee = iv(0, x3) * kv(1, x4) + kv(0, x3) * iv(1, x4)
        ff = kv(1, x3) * iv(1, x4) - iv(1, x3) * kv(1, x4)
        rr = kv(1, x3) * iv(0, x4) + iv(1, x3) * kv(0, x4)

        # Calculating z10, z2i, z2m, z20, z12
        rho1 = 1 / self.sigma[0]
        rho2 = 1 / self.sigma[1]

        # Solid conductor case
        if x1 == 0:
            z10 = (m1 * rho1 / (2 * np.pi * self.b)) * iv(0, x2) / iv(1, x2)
        else:
            z10 = (m1 * rho1 / (2 * np.pi * self.b)) * aa / bb

        z2i = (m2 * rho2 / (2 * np.pi * self.b_prime)) * ee / ff
        z2m = rho2 / (2 * np.pi * self.b_prime * self.c * ff)
        z20 = (m2 * rho2 / (2 * np.pi * self.c)) * rr / ff
        z12 = jw * (self.mu[0] / 2 / np.pi) * np.log(self.b_prime / self.b)

        # Calculating Z11, Z12, Z22
        zz22 = z20
        zz12 = zz22 - z2m
        zz11 = z10 + z12 + z2i - z2m + zz12

        return zz11, zz12, zz22
