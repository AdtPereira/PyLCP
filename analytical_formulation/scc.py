"""
REFERENCES:
[1] XUE, Haoyan. General Formulation and Accurate Evaluation of Earth-Return Parameters
    for Overhead / Underground Cables. PhD thesis, Department of Electrical Engineering,
    École Polytechnique de Montréal, Université de Montréal, August 2018. 

[2] A. Ametani, T. Ohno and N. Nagaoka, Cable System Transients: Theory, Modeling and 
    Simulation, Wiley-IEEE Press, 2015.
"""

import numpy as np
import scipy.special as ss
import scipy.constants as sc
from scipy.integrate import quad
from mtl_main.source import MulticonductorTransmissionLine
from scipy import linalg

def sommerfeld_quasi_tem_approx(hnm, dnm, ke2, ka2, s_form='s1', type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

    # soil refractive index
    if s_form == 's1':
        n = 1
    elif s_form == 's2':
        n = np.sqrt(ke2 / ka2)

    # Define the real and imaginary parts of the integrand
    def _int_quad_real(x, hnm, dnm, ke2, ka2, n):
        sqrt_x2_ke2 = np.sqrt(x**2 - ke2)            
        sqrt_x2_ka2 = np.sqrt(x**2 - ka2)
        num = np.exp(hnm * sqrt_x2_ke2) * np.cos(dnm * x)
        den = sqrt_x2_ka2 + sqrt_x2_ke2
        return (num/den).real

    def _int_quad_imag(x, hnm, dnm, ke2, ka2, n):
        sqrt_x2_ke2 = np.sqrt(x**2 - ke2)
        sqrt_x2_ka2 = np.sqrt(x**2 - ka2)
        num = np.exp(hnm * sqrt_x2_ke2) * np.cos(dnm * x)
        den = sqrt_x2_ka2 + sqrt_x2_ke2
        return (num/den).imag

    def _int_transformed_real(t, hnm, dnm, ke2, ka2, n):
        x = np.tan(t)
        sqrt_x2_ke2 = np.sqrt(x**2 - ke2)
        sqrt_x2_ka2 = np.sqrt(x**2 - ka2)
        num = np.exp(hnm * sqrt_x2_ke2) * np.cos(dnm * x)
        den = sqrt_x2_ka2 + sqrt_x2_ke2
        return (num/den).real * (1 / np.cos(t)**2)

    def _int_transformed_imag(t, hnm, dnm, ke2, ka2, n):
        x = np.tan(t)
        sqrt_x2_ke2 = np.sqrt(x**2 - ke2)
        sqrt_x2_ka2 = np.sqrt(x**2 - ka2)
        num = np.exp(hnm * sqrt_x2_ke2) * np.cos(dnm * x)
        den = sqrt_x2_ka2 + sqrt_x2_ke2
        return (num/den).imag * (1 / np.cos(t)**2)

    def _int_gauss_legendre(func, a, b, n, *args):
        """
        Integrates the function `func` over the interval [a, b] using the Gauss-Legendre method.
        - func: function to be integrated.
        - a, b: integration limits.
        - n: number of quadrature points.
        - args: additional arguments to be passed to `func`.
        """
        [x, w] = np.polynomial.legendre.leggauss(n)
        t = 0.5 * (x + 1) * (b - a) + a
        return np.sum(w * func(t, *args)) * 0.5 * (b - a)

    # Calculate the real and imaginary parts of the integral using Gauss-Legendre
    if type_form == 'gauss_legendre':
        gauss_legendre_real = _int_gauss_legendre(_int_transformed_real, 0, np.pi/2, pts, hnm, dnm, ke2, ka2, n)
        gauss_legendre_imag = _int_gauss_legendre(_int_transformed_imag, 0, np.pi/2, pts, hnm, dnm, ke2, ka2, n)
        Sn = gauss_legendre_real + 1j * gauss_legendre_imag

    # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
    else:
        quad_real, _ = quad(_int_quad_real, 0, np.inf, args=(hnm, dnm, ke2, ka2, n))
        quad_imag, _ = quad(_int_quad_imag, 0, np.inf, args=(hnm, dnm, ke2, ka2, n))
        Sn = quad_real + 1j * quad_imag

    return Sn

def sommerfeld_ametani_approx(hnm, dnm, ke2, type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

    # Define the real and imaginary parts of the integrand
    def _int_quad_real(x, hnm, dnm, ke2):
        num = np.exp(hnm * x) * np.cos(dnm * x)
        den = x + np.sqrt(x**2 - ke2)  
        return (num/den).real

    def _int_quad_imag(x, hnm, dnm, ke2):
        num = np.exp(hnm * x) * np.cos(dnm * x)
        den = x + np.sqrt(x**2 - ke2)  
        return (num/den).imag

    def _int_transformed_real(t, hnm, dnm, ke2):
        x = np.tan(t)
        num = np.exp(hnm * x) * np.cos(dnm * x)
        den = x + np.sqrt(x**2 - ke2)  
        return (num/den).real * (1 / np.cos(t)**2)

    def _int_transformed_imag(t, hnm, dnm, ke2):
        x = np.tan(t)
        num = np.exp(hnm * x) * np.cos(dnm * x)
        den = x + np.sqrt(x**2 - ke2)  
        return (num/den).imag * (1 / np.cos(t)**2)

    def _int_gauss_legendre(func, a, b, n, *args):
        """
        Integrates the function `func` over the interval [a, b] using the Gauss-Legendre method.
        - func: function to be integrated.
        - a, b: integration limits.
        - n: number of quadrature points.
        - args: additional arguments to be passed to `func`.
        """
        [x, w] = np.polynomial.legendre.leggauss(n)
        t = 0.5 * (x + 1) * (b - a) + a
        return np.sum(w * func(t, *args)) * 0.5 * (b - a)

    # Calculate the real and imaginary parts of the integral using Gauss-Legendre
    if type_form == 'gauss_legendre':
        gauss_legendre_real = _int_gauss_legendre(_int_transformed_real, 0, np.pi/2, pts, hnm, dnm, ke2)
        gauss_legendre_imag = _int_gauss_legendre(_int_transformed_imag, 0, np.pi/2, pts, hnm, dnm, ke2)
        Sn = gauss_legendre_real + 1j * gauss_legendre_imag

    # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
    else:
        quad_real, _ = quad(_int_quad_real, 0, np.inf, args=(hnm, dnm, ke2))
        quad_imag, _ = quad(_int_quad_imag, 0, np.inf, args=(hnm, dnm, ke2))
        Sn = quad_real + 1j * quad_imag

    return Sn

def t_sommerfeld(hnm, dnm, ke2, ka2, type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

    # Define the real and imaginary parts of the integrand
    def _int_quad_real(x, hnm, dnm, k_e2, k_a2):
        # soil refractive index
        n = np.sqrt(k_e2 / k_a2)

        exp_1 = np.exp(-hnm * x)
        exp_2 = np.exp(-0.5 * hnm * x)
        u2 = np.sqrt(x**2 + k_a2 - k_e2)
        num = u2 * (exp_2 - exp_1) * np.cos(dnm * x)
        den = (n * x)**2 + x * u2

        return (num/den).real

    def _int_quad_imag(x, hnm, dnm, k_e2, k_a2):
        # soil refractive index
        n = np.sqrt(k_e2 / k_a2)

        exp_1 = np.exp(-hnm * x)
        exp_2 = np.exp(-0.5 * hnm * x)
        u2 = np.sqrt(x**2 + k_a2 - k_e2)
        num = u2 * (exp_2 - exp_1) * np.cos(dnm * x)
        den = (n * x)**2 + x * u2

        return (num/den).imag

    def _int_transformed_real(t, hnm, dnm, k_e2, k_a2):
        x = np.tan(t)

        # soil refractive index
        n = np.sqrt(k_e2 / k_a2)

        exp_1 = np.exp(-hnm * x)
        exp_2 = np.exp(-0.5 * hnm * x)
        u2 = np.sqrt(x**2 + k_a2 - k_e2)
        num = u2 * (exp_2 - exp_1) * np.cos(dnm * x)
        den = (n * x)**2 + x * u2

        return (num/den).real * (1 / np.cos(t)**2)

    def _int_transformed_imag(t, hnm, dnm, k_e2, k_a2):
        x = np.tan(t)

        # soil refractive index
        n = np.sqrt(k_e2 / k_a2)

        exp_1 = np.exp(-hnm * x)
        exp_2 = np.exp(-0.5 * hnm * x)
        u2 = np.sqrt(x**2 + k_a2 - k_e2)
        num = u2 * (exp_2 - exp_1) * np.cos(dnm * x)
        den = (n * x)**2 + x * u2

        return (num/den).imag * (1 / np.cos(t)**2)

    def _int_gauss_legendre(func, a, b, n, *args):
        """
        Integrates the function `func` over the interval [a, b] using the Gauss-Legendre method.
        - func: function to be integrated.
        - a, b: integration limits.
        - n: number of quadrature points.
        - args: additional arguments to be passed to `func`.
        """
        [x, w] = np.polynomial.legendre.leggauss(n)
        t = 0.5 * (x + 1) * (b - a) + a
        return np.sum(w * func(t, *args)) * 0.5 * (b - a)

    # Calculate the real and imaginary parts of the integral using Gauss-Legendre
    if type_form == 'gauss_legendre':
        gauss_legendre_real = _int_gauss_legendre(_int_transformed_real, 0, np.pi/2, pts, hnm, dnm, ke2, ka2)
        gauss_legendre_imag = _int_gauss_legendre(_int_transformed_imag, 0, np.pi/2, pts, hnm, dnm, ke2, ka2)
        T = gauss_legendre_real + 1j * gauss_legendre_imag

    # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
    else:
        quad_real, _ = quad(_int_quad_real, 0, np.inf, args=(hnm, dnm, ke2, ka2))
        quad_imag, _ = quad(_int_quad_imag, 0, np.inf, args=(hnm, dnm, ke2, ka2))
        T = quad_real + 1j * quad_imag

    return T

class PerUnitParameters:    
    """ This class calculates PUL parameters using an MTL geometry model. """

    def __init__(self, model: MulticonductorTransmissionLine, f: float):
        # MTL Geometry Model
        self.model = model

        # Create a sorted list of CORE conductors to ensure consistent ordering
        self.core_conductors = {key: value for key, value in model.mtl.items() if value.get('conductor_name') in ['core']}.items()
        
        # Soil Relative Permittivity
        self.e1 = model.mtl_ref[0]['relative_permittivity'] * sc.epsilon_0

        # Soil conductivity (S/m)
        self.sigma_1 = model.mtl_ref[0]['conductivity']

        # Soil Relative Permeability
        self.mu1 = model.mtl_ref[0]['relative_permeability'] * sc.mu_0

        # External Conductance (S/m)
        # self.ge = model.mtl_ref[0]['external_conductance']

        # Soil resistivity (ohm.m)
        self.rho_1 = 1 / self.sigma_1

        # Angular frequency (rad/s)
        self.jw = 1j * 2 * np.pi * f
        self.jw_mu0_2pi = self.jw * sc.mu_0 / (2 * np.pi)
        self.jw_2pi_e0 = self.jw * 2 * np.pi * sc.epsilon_0

        # Air wave number - Equation (2.15) [1]
        self.k_air2 = - self.jw * sc.mu_0 * self.jw * sc.epsilon_0

        # Earth wave number - Equation (2.15) [1]
        self.k_earth2 = - self.jw * self.mu1 * (self.sigma_1 + self.jw * self.e1)

    # 2.2.1 Impedance of Single-core Coaxial Cable (SC Cable) [2]
    def internal_parameters_by_bessel(self):
        """
        Calculates the internal impedance matrix [zi] for a single-core cable (SCC)
        with a core, sheath, and optional armor and jacket.

        The method uses formulas for tubular conductors, which involve modified
        Bessel functions, to account for skin and proximity effects within the
        conductors. Special cases for solid conductors are also handled.

        The formulas are based on the component impedances for each layer:
        - z11: internal impedance of core outer surface
        - z12: core outer insulator impedance
        - z2i: internal impedance of sheath inner surface
        - z20: internal impedance of sheath outer surface
        - z2m: sheath mutual impedance
        - z23: sheath outer insulator impedance
        - z3i: internal impedance of armor inner surface
        - z30: internal impedance of armor outer surface
        - z34: armor outer insulator impedance

        Reference:  A. Ametani, T. Ohno and N. Nagaoka, Cable System Transients: Theory, Modeling and
                    Simulation, Wiley-IEEE Press, 2015.
        """
        s = self.jw
        scc = self.model.scc
        two_pi = 2 * np.pi

        # --- Initialize all impedance components ---
        z11, z12 = (None,) * 2
        z2i, z20, z23, z2m = (None,) * 4
        z3i, z30, z34, z3m = (None,) * 4

        # --- z11: internal impedance of core outer surface ---
        if 'core_outer_radius' in scc:
            rho1, mu1 = scc['core_resistivity'], scc['core_permeability']
            r1, r2 = scc['core_inner_radius'], scc['core_outer_radius']
            m_core = np.sqrt(s * mu1 / rho1)
            x1, x2 = m_core * r1, m_core * r2

            # Case 1: Solid core (r1 = 0)
            if np.isclose(r1, 0):
                with np.errstate(divide='ignore', invalid='ignore'):
                    z11 = (m_core * rho1 / (two_pi * r2)) * (ss.iv(0, x2) / ss.iv(1, x2)) if not np.isclose(x2, 0) else np.complex(0, np.inf)
            
            # Case 2: Tubular core (r1 > 0)
            else:
                # D1 = I1(x2) * K1(x1) - I1(x1) * K1(x2)
                D1 = ss.iv(1, x2) * ss.kv(1, x1) - ss.iv(1, x1) * ss.kv(1, x2)
                
                # N1 = I0(x2) * K1(x1) + K0(x2) * I1(x1)
                N1 = ss.iv(0, x2) * ss.kv(1, x1) + ss.kv(0, x2) * ss.iv(1, x1)

                with np.errstate(divide='ignore', invalid='ignore'):
                    z11 = (s * mu1 / two_pi) * (1 / (x2 * D1)) * N1 if not (np.isclose(x2, 0) or np.isclose(D1, 0)) else np.complex(0, np.inf)

        # --- z12: Core outer insulator impedance ---
        if 'core_insulation_outer_radius' in scc:
            mui1 = scc['core_insulation_permeability']
            r3 = scc['core_insulation_outer_radius']
            z12 = (s * mui1 / two_pi) * np.log(r3 / r2) if not np.isclose(r3, r2) else 0
        
        # --- Sheath Impedance (z20, z2i, z2m) ---
        if 'sheath_outer_radius' in scc:
            rho2, mu2 = scc['sheath_resistivity'], scc['sheath_permeability']
            r3, r4 = scc['sheath_inner_radius'], scc['sheath_outer_radius']
            m_sheath = np.sqrt(s * mu2 / rho2)
            x3, x4 = m_sheath * r3, m_sheath * r4

            # --- Common terms for sheath calculations ---
            # D2 = I1(x4) * K1(x3) - I1(x3) * K1(x4)
            D2 = ss.iv(1, x4) * ss.kv(1, x3) - ss.iv(1, x3) * ss.kv(1, x4)

            # --- z2m: sheath mutual impedance ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z2m = rho2 / (two_pi * r3 * r4 * D2) if not (np.isclose(r3, 0) or np.isclose(r4, 0) or np.isclose(D2, 0)) else np.complex(0, np.inf)
            
            # --- z2i: Internal impedance of sheath inner surface ---
            # N2i = I0(x3) * K1(x4) + K0(x3) * I1(x4)
            N2i = ss.iv(0, x3) * ss.kv(1, x4) + ss.kv(0, x3) * ss.iv(1, x4)
            with np.errstate(divide='ignore', invalid='ignore'):
                z2i = (s * mu2 / two_pi) * (1 / (x3 * D2)) * N2i if not (np.isclose(x3, 0) or np.isclose(D2, 0)) else np.complex(0, np.inf)

            # --- z20: Internal impedance of sheath outer surface ---
            # N20 = I0(x4) * K1(x3) + K0(x4) * I1(x3)
            N20 = ss.iv(0, x4) * ss.kv(1, x3) + ss.kv(0, x4) * ss.iv(1, x3)
            with np.errstate(divide='ignore', invalid='ignore'):
                z20 = (s * mu2 / two_pi) * (1 / (x4 * D2)) * N20 if not (np.isclose(x4, 0) or np.isclose(D2, 0)) else np.complex(0, np.inf)

        # --- z23: Sheath outer insulator impedance ---
        if 'sheath_insulation_outer_radius' in scc:
            mui2 = scc['sheath_insulation_permeability']
            r5 = scc['sheath_insulation_outer_radius']
            z23 = (s * mui2 / two_pi) * np.log(r5 / r4) if not np.isclose(r5, r4) else 0

        # --- Armor Impedances (z3i, z30) ---
        if 'armor_outer_radius' in scc:
            rho3, mu3 = scc['armor_resistivity'], scc['armor_permeability']
            r5, r6 = scc['armor_inner_radius'], scc['armor_outer_radius']
            m_armor = np.sqrt(s * mu3 / rho3)
            x5, x6 = m_armor * r5, m_armor * r6

            # --- z3m: armor mutual impedance ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z3m = rho3 / (two_pi * r5 * r6 * D3) if not (np.isclose(r5, 0) or np.isclose(r6, 0) or np.isclose(D3, 0)) else np.complex(0, np.inf)

            # --- Common terms for armor calculations ---
            # D3 = I1(x6) * K1(x5) - I1(x5) * K1(x6)
            # N3 = I0(x5) * K1(x6) + K0(x5) * I1(x6)
            D3 = ss.iv(1, x6) * ss.kv(1, x5) - ss.iv(1, x5) * ss.kv(1, x6)
            N3 = ss.iv(0, x5) * ss.kv(1, x6) + ss.kv(0, x5) * ss.iv(1, x6)

            # --- z3i: Internal impedance of armor inner surface ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z3i = (s * mu3 / two_pi) * (1 / (x5 * D3)) * N3 if not (np.isclose(x5, 0) or np.isclose(D3, 0)) else np.complex(0, np.inf)

            # --- z30: internal impedance of armor outer surface ---
            # N30 = I0(x6) * K1(x5) + K0(x6) * I1(x5)
            N30 = ss.iv(0, x6) * ss.kv(1, x5) + ss.kv(0, x6) * ss.iv(1, x5)
            with np.errstate(divide='ignore', invalid='ignore'):
                z30 = (s * mu3 / two_pi) * (1 / (x6 * D3)) * N30 if not (np.isclose(x6, 0) or np.isclose(D3, 0)) else np.complex(0, np.inf)

        # --- z34: armor outer insulator impedance ---
        if 'armor_insulation_outer_radius' in scc:
            mui3 = scc['armor_insulation_permeability']
            r7 = scc['armor_insulation_outer_radius']
            z34 = (s * mui3 / two_pi) * np.log(r7/r6) if not np.isclose(r7, r6) else 0

        return {
            'zcs': {'z11': z11, 'z12': z12, 'z2i': z2i},    # Eq. (2.10a) [2]
            'zsa': {'z20': z20, 'z23': z23, 'z3i': z3i},    # Eq. (2.10b) [2]
            'za4': {'z30': z30, 'z34': z34},                # Eq. (2.10c) [2]
            'zs3': {'z20': z20, 'z23': z23},                # Eq. (2.12a) [2]
            'z2m': z2m,                                     # sheath mutual impedance
            'z3m': z3m                                      # armor mutual impedance
        }
    
    # 3.3.2.2 Earth-return impedance and admittance formulas based on quasi-TEM assumption [1]
    def ground_return_parameters(self, form='magalhaes_xue'):
        """ This method calculates the impedance matrix of the earth return path. """
        N = len(self.core_conductors)
        S1c = np.zeros((N, N), dtype=complex)
        S2c = np.zeros((N, N), dtype=complex)
        Tc = np.zeros((N, N), dtype=complex)

        # Distance between the conductors (m)
        d_matrix = self.model.d_matrix_ground_return
        D_matrix = self.model.D_matrix_ground_return
        hnm = self.model.vertical_separation_matrix
        dnm = self.model.horizontal_separation_matrix

        k_air2, k_earth2 = self.k_air2, self.k_earth2
        
        if form in ['sunde', 'deconti_sunde']:
            k_air2 = 0
        
        elif form in ['pollaczek', 'ametani', 'saad', 'wedepohl']:
            k_air2 = 0
            k_earth2 = - self.jw * self.mu1 * self.sigma_1        
        
        K0_jke_dnm = ss.kv(0,  1j * np.sqrt(k_earth2) * d_matrix)
        K0_jke_Dnm = ss.kv(0,  1j * np.sqrt(k_earth2) * D_matrix)

        # Wedepohl e Wilcox Approximation Expression
        if form in ['wedepohl']:
            yg = 1j * np.sqrt(k_earth2)
            ln_term = np.log(0.5 * np.euler_gamma * yg * d_matrix)
            S1c = - ln_term + 0.5 + (2/3) * yg * hnm
            zg = self.jw_mu0_2pi * S1c

        # De Conti Approximation Expressions
        elif form in ['deconti', 'deconti_sunde', 'saad']:
            y0 = 1j * np.sqrt(k_air2)
            yg = 1j * np.sqrt(k_earth2)
            term_1 = (yg - y0) / (yg + y0)
            exp_term = np.exp(hnm * yg)
            term_2 = 2 / (4 + (yg**2 * dnm**2))
            zg = self.jw_mu0_2pi * (K0_jke_dnm + term_1 * exp_term * term_2)

        # Integral Expressions
        elif form in ['magalhaes_xue', 'sunde', 'pollaczek', 'ametani']:
            for n in range(N):
                for m in range(N):
                    hnm = self.model.vertical_separation_matrix[n, m]
                    dnm = self.model.horizontal_separation_matrix[n, m]                    
                    
                    # Ametani Integral Equation
                    if form == 'ametani':
                        S1c[n, m] = 2 * sommerfeld_ametani_approx(hnm, dnm, ke2=k_earth2)                    
                    
                    # Quasi-TEM Integral Equation
                    else:
                        S1c[n, m] = 2 * sommerfeld_quasi_tem_approx(hnm, dnm, ke2=k_earth2, ka2=k_air2)

            # Earth-return impedance based on quasi-TEM assumption [1]
            zg = self.jw_mu0_2pi * (K0_jke_dnm - K0_jke_Dnm + S1c)
        
        # ysh = self.jw_2pi_e0 * np.linalg.inv(M - Tc)
        # else:
        #     zg = self.jw_mu0_2pi * (S1c)
        #     ysh = self.jw_2pi_e0 * np.linalg.inv(M + S2c)

        # # Series Impedance Matrix
        # zs = zi + ze + zg 

        # # Currents and Voltages Propagation Constants Matrix
        # gamma_i = linalg.sqrtm(ysh @ zs)
        # gamma_v = linalg.sqrtm(zs @ ysh)

        # # Inverse of Y using LU decomposition
        # # Solve the system Y * Y_inv = I to find Y_inv
        # y_inverse = linalg.lu_solve(linalg.lu_factor(ysh), np.identity(ysh.shape[0]))

        # # Inverse of Z using LU decomposition
        # # Solve the system Z * Z_inv = I to find Z_inv
        # z_inverse = linalg.lu_solve(linalg.lu_factor(zs), np.identity(zs.shape[0]))

        # # Characteristic impedance Zc = Y_inv * sqrt(Y*Z)
        # zc = y_inverse @ gamma_i 

        # # Characteristic admittance Yc = Z_inv * sqrt(Z*Y)
        # yc = z_inverse @ gamma_v

        # return {
        #     'zi': zi, 'ze': ze, 'zg': zg, 'zs': zs, 'ysh': ysh,
        #     'gamma_i': gamma_i, 'gamma_v': gamma_v, 'zc': zc, 'yc': yc
        # }

        return {
            'zg': zg,
        }