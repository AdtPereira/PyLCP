"""
REFERENCES:
[1] XUE, Haoyan. General Formulation and Accurate Evaluation of Earth-Return Parameters
    for Overhead / Underground Cables. PhD thesis, Department of Electrical Engineering,
    École Polytechnique de Montréal, Université de Montréal, August 2018. 
"""

import numpy as np
import scipy.special as ss
import scipy.constants as sc
from scipy.integrate import quad
from mtl_main.source import MulticonductorTransmissionLine
from scipy import linalg

def sommerfeld_quasi_tem_approx(hnm, dnm, ke2, ka2, s_form='s1', type_form='quad', pts=150):
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

        # Soil Relative Permittivity
        self.er_1 = model.mtl_ref[0]['relative_permittivity']

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

        # Constant vacuum terms
        self.jw_mu0_2pi = self.jw * sc.mu_0 / (2 * np.pi)
        self.jw_2pi_e0 = self.jw * 2 * np.pi * sc.epsilon_0

        # Air wave number - Equation (2.15) [1]
        self.k_air2 = - self.jw * sc.mu_0 * self.jw * sc.epsilon_0

        # Earth wave number - Equation (2.15) [1]
        self.k_earth2 = - self.jw * self.mu1 * (self.sigma_1 + self.jw * self.er_1 * sc.epsilon_0)

    # 3.3.2.2 Earth-return impedance and admittance formulas based on quasi-TEM assumption [1]
    def pul_extended_theory(self, form='magalhaes_xue'):
        """ This method calculates the impedance matrix of the earth return path. """
        N = len(self.model.surfaces)
        S1c = np.zeros((N, N), dtype=complex)
        S2c = np.zeros((N, N), dtype=complex)
        Tc = np.zeros((N, N), dtype=complex)

        k_air2, k_earth2 = self.k_air2, self.k_earth2
        if form in ['sunde', 'deconti_sunde']:
            k_air2 = 0
        if form in ['pollaczek', 'ametani', 'saad', 'wedepohl']:
            k_air2 = 0
            k_earth2 = - self.jw * self.mu1 * self.sigma_1

        for conductor_n in self.model.surfaces:
            for conductor_m in self.model.surfaces:
                # Get conductor tags
                n, m = conductor_n['tag']-1, conductor_m['tag']-1

                # Distance between the conductors (m)
                d_matrix = self.model.d_matrix_ground_return
                D_matrix = self.model.D_matrix_ground_return
                hnm = self.model.vertical_separation_matrix[n, m]
                dnm = self.model.horizontal_separation_matrix[n, m]

                # Wedepohl e Wilcox Approximation Expression
                if form == 'wedepohl':
                    yg = 1j * np.sqrt(k_earth2)
                    ln_term = np.log(0.5 * np.euler_gamma * yg * d_matrix)
                    S1c[n, m] = - ln_term + 0.5 + (2/3) * yg * hnm

                # De Conti Approximation Expression
                elif form in ['deconti', 'deconti_sunde', 'saad']:
                    y0 = 1j * np.sqrt(k_air2)
                    yg = 1j * np.sqrt(k_earth2)
                    term_1 = (yg - y0) / (yg + y0)
                    exp_term = np.exp(hnm * yg)
                    term_2 = 2 / (4 + (yg**2 * dnm**2))
                    S1c[n, m] = term_1 * exp_term * term_2

                # Ametani Integral Equation
                elif form == 'ametani':
                    S1c[n, m] = 2 * sommerfeld_ametani_approx(hnm, dnm, ke2=k_earth2)                    
                
                # Quasi-TEM Integral Equation
                else:
                    S1c[n, m] = 2 * sommerfeld_quasi_tem_approx(hnm, dnm, ke2=k_earth2, ka2=k_air2)
                    
        # Approx. Internal Impedance Matrix
        # zi = self.internal_impedance_elements_solid_wires(ExactlyForms=False)['Zi_approx']

        # Earth return impedance and admittance (ohm/m)
        K0_jke_dnm = ss.kv(0,  1j * np.sqrt(k_earth2) * d_matrix)
        K0_jke_Dnm = ss.kv(0,  1j * np.sqrt(k_earth2) * D_matrix)

        if form in ['wedepohl']:
            zg = self.jw_mu0_2pi * S1c

        elif form in ['deconti', 'deconti_sunde', 'saad']:
            zg = self.jw_mu0_2pi * (K0_jke_dnm + S1c)

        else:
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