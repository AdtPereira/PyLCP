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

def sn_sommerfeld(hnm, dnm, ke2, ka2, s_form='s1', type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

    # soil refractive index
    if s_form == 's1':
        n = 1
    elif s_form == 's2':
        n = np.sqrt(ke2 / ka2)

    # Define the real and imaginary parts of the integrand
    def int_quad_real(x, hn_hm, dn_dm, ke2, ka2, n):
        num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
        den = np.sqrt(x**2 + ka2 - ke2) + (n**2 * x)
        return (num/den).real

    def int_quad_imag(x, hn_hm, dn_dm, ke2, ka2, n):
        num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
        den = np.sqrt(x**2 + ka2 - ke2) + (n**2 * x)
        return (num/den).imag

    def int_transformed_real(t, hn_hm, dn_dm, k_e2, k_a2, n):
        x = np.tan(t)
        num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
        den = np.sqrt(x**2 + k_a2 - k_e2) + (n**2 * x)
        return (num/den).real * (1 / np.cos(t)**2)

    def int_transformed_imag(t, hn_hm, dn_dm, k_e2, k_a2, n):
        x = np.tan(t)
        num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
        den = np.sqrt(x**2 + k_a2 - k_e2) + (n**2 * x)
        return (num/den).imag * (1 / np.cos(t)**2)

    def int_gauss_legendre(func, a, b, n, *args):
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
        gauss_legendre_real = int_gauss_legendre(
            int_transformed_real, 0, np.pi/2, pts, hnm, dnm, ke2, ka2, n)
        gauss_legendre_imag = int_gauss_legendre(
            int_transformed_imag, 0, np.pi/2, pts, hnm, dnm, ke2, ka2, n)
        S_n = gauss_legendre_real + 1j * gauss_legendre_imag

    else:
        # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
        quad_real, _ = quad(int_quad_real, 0, np.inf,
                            args=(hnm, dnm, ke2, ka2, n))
        quad_imag, _ = quad(int_quad_imag, 0, np.inf,
                            args=(hnm, dnm, ke2, ka2, n))
        S_n = quad_real + 1j * quad_imag

    return S_n

def t_sommerfeld(hnm, dnm, ke2, ka2, type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

    # Define the real and imaginary parts of the integrand
    def _int_quad_real(x, hn_hm, dn_dm, k_e2, k_a2):
        # soil refractive index
        n = np.sqrt(k_e2 / k_a2)

        exp_1 = np.exp(-hn_hm * x)
        exp_2 = np.exp(-0.5 * hn_hm * x)
        u2 = np.sqrt(x**2 + k_a2 - k_e2)
        num = u2 * (exp_2 - exp_1) * np.cos(dn_dm * x)
        den = (n * x)**2 + x * u2

        return (num/den).real

    def _int_quad_imag(x, hn_hm, dn_dm, k_e2, k_a2):
        # soil refractive index
        n = np.sqrt(k_e2 / k_a2)

        exp_1 = np.exp(-hn_hm * x)
        exp_2 = np.exp(-0.5 * hn_hm * x)
        u2 = np.sqrt(x**2 + k_a2 - k_e2)
        num = u2 * (exp_2 - exp_1) * np.cos(dn_dm * x)
        den = (n * x)**2 + x * u2

        return (num/den).imag

    def _int_transformed_real(t, hn_hm, dn_dm, k_e2, k_a2):
        x = np.tan(t)

        # soil refractive index
        n = np.sqrt(k_e2 / k_a2)

        exp_1 = np.exp(-hn_hm * x)
        exp_2 = np.exp(-0.5 * hn_hm * x)
        u2 = np.sqrt(x**2 + k_a2 - k_e2)
        num = u2 * (exp_2 - exp_1) * np.cos(dn_dm * x)
        den = (n * x)**2 + x * u2

        return (num/den).real * (1 / np.cos(t)**2)

    def _int_transformed_imag(t, hn_hm, dn_dm, k_e2, k_a2):
        x = np.tan(t)

        # soil refractive index
        n = np.sqrt(k_e2 / k_a2)

        exp_1 = np.exp(-hn_hm * x)
        exp_2 = np.exp(-0.5 * hn_hm * x)
        u2 = np.sqrt(x**2 + k_a2 - k_e2)
        num = u2 * (exp_2 - exp_1) * np.cos(dn_dm * x)
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

    else:
        # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
        quad_real, _ = quad(_int_quad_real, 0, np.inf, args=(hnm, dnm, ke2, ka2))
        quad_imag, _ = quad(_int_quad_imag, 0, np.inf, args=(hnm, dnm, ke2, ka2))
        T = quad_real + 1j * quad_imag

    return T

def sn_sommerfeld_low_frequencies(hnm, dnm, ke2, n2=1, type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""  

    def int_quad_real(x, hn_hm, dn_dm, k_e2, n2):
        num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
        den = np.sqrt(x**2 - k_e2) + n2 * x
        return (num/den).real

    def int_quad_imag(x, hn_hm, dn_dm, k_e2, n2):
        num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
        den = np.sqrt(x**2 - k_e2) + n2 * x
        return (num/den).imag

    def int_transformed_real(t, hn_hm, dn_dm, k_e2, n2):
        x = np.tan(t)
        num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
        den = np.sqrt(x**2 - k_e2) + n2 * x
        return (num/den).real * (1 / np.cos(t)**2)

    def int_transformed_imag(t, hn_hm, dn_dm, k_e2, n2):
        x = np.tan(t)
        num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
        den = np.sqrt(x**2 - k_e2) + n2 * x
        return (num/den).imag * (1 / np.cos(t)**2)

    def int_gauss_legendre(func, a, b, n, *args):
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
        gauss_legendre_real = int_gauss_legendre(int_transformed_real, 0, np.pi/2, pts, hnm, dnm, ke2, n2)
        gauss_legendre_imag = int_gauss_legendre(int_transformed_imag, 0, np.pi/2, pts, hnm, dnm, ke2, n2)
        Sn = gauss_legendre_real + 1j * gauss_legendre_imag

    # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
    else:
        quad_real, _ = quad(int_quad_real, 0, np.inf, args=(hnm, dnm, ke2, n2))
        quad_imag, _ = quad(int_quad_imag, 0, np.inf, args=(hnm, dnm, ke2, n2))
        Sn = quad_real + 1j * quad_imag

    return Sn

def verify_kelvin_functions(q = 1.5):
    """
    This script verifies the output of scipy.special.kelvin() by comparing it
    against the fundamental mathematical definitions that relate Kelvin functions
    to the modified Bessel functions I_0 and I_1 with a complex argument.
    
    This version correctly unpacks the complex tuple returned by ss.kelvin()
    as per the official SciPy documentation:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kelvin.html
    """
    
    print(f"--- Verification for q = {q} ---")
    
    # --- 1. Call ss.kelvin() and correctly extract real values ---
    print("\n[Reference] Unpacking values from ss.kelvin(q) according to SciPy docs:")
    
    # ss.kelvin() returns a tuple of 4 complex numbers: (Be, Ke, Bep, Kep)
    Be_complex, Ke_complex, Bep_complex, Kep_complex = ss.kelvin(q)
    
    # The real-valued functions are the real/imaginary parts of the complex results.
    # ber and bei come from the first element (Be).
    ber_ref = Be_complex.real
    bei_ref = Be_complex.imag
    # ber' and bei' come from the third element (Bep).
    ber_p_ref = Bep_complex.real
    bei_p_ref = Bep_complex.imag
    
    print(f"ber(q)  = Re[kelvin(q)[0]] = {ber_ref:>10.6f}")
    print(f"bei(q)  = Im[kelvin(q)[0]] = {bei_ref:>10.6f}")
    print(f"ber'(q) = Re[kelvin(q)[2]] = {ber_p_ref:>10.6f}")
    print(f"bei'(q) = Im[kelvin(q)[2]] = {bei_p_ref:>10.6f}")
    
    # --- 2. Verification using the Bessel function I_0 ---
    # According to the definition: I_0(q*sqrt(j)) = ber(q) + j*bei(q)
    complex_arg = q * (1j**0.5)
    i0_complex = ss.iv(0, complex_arg)
    
    ber_from_bessel = i0_complex.real
    bei_from_bessel = i0_complex.imag
    
    print("\n[Test 1] Values derived from I_0(q * sqrt(j)):")
    print(f"Re[I_0] = {ber_from_bessel:>10.6f} -> Matches ber(q)? {np.isclose(ber_ref, ber_from_bessel)}")
    print(f"Im[I_0] = {bei_from_bessel:>10.6f} -> Matches bei(q)? {np.isclose(bei_ref, bei_from_bessel)}")

    # --- 3. Verification of the derivatives using Bessel function I_1 ---
    # According to the definition: sqrt(j)*I_1(q*sqrt(j)) = ber'(q) + j*bei'(q)
    i1_complex_term = (1j**0.5) * ss.iv(1, complex_arg)

    ber_p_from_bessel = i1_complex_term.real
    bei_p_from_bessel = i1_complex_term.imag

    print("\n[Test 2] Values derived from sqrt(j) * I_1(q * sqrt(j)):")
    print(f"Re[...] = {ber_p_from_bessel:>10.6f} -> Matches ber'(q)? {np.isclose(ber_p_ref, ber_p_from_bessel)}")
    print(f"Im[...] = {bei_p_from_bessel:>10.6f} -> Matches bei'(q)? {np.isclose(bei_p_ref, bei_p_from_bessel)}")
    print("\nConclusion: With correct unpacking, the outputs from ss.kelvin() are consistent.")
    
class PerUnitParameters:
    """ This class calculates PUL parameters using an MTL geometry model. """

    def __init__(self, model: MulticonductorTransmissionLine, f: float):
        # MTL Geometry Model
        self.model = model

        # Skin Depth (m)
        self.skin_depth = 1 / np.sqrt(self.model.mu * np.pi * f * self.model.sigma)

        # Soil Relative Permittivity
        self.er_1 = model.mtl_ref[0]['relative_permittivity']

        # Soil conductivity (S/m)
        self.sigma_1 = model.mtl_ref[0]['conductivity']

        # Soil Relative Permeability
        self.mur_1 = model.mtl_ref[0]['relative_permeability']

        # External Conductance (S/m)
        self.ge = model.mtl_ref[0]['external_conductance']

        # Angular frequency (rad/s)
        self.jw = 1j * 2 * np.pi * f

        # Constant vacuum terms
        self.jw_mu0_2pi = self.jw * sc.mu_0 / (2 * np.pi)
        self.jw_2pi_e0 = self.jw * 2 * np.pi * sc.epsilon_0

        # Air wave number - Equation (2.15) [1] (rad/m)
        self.k_air2 = - self.jw * sc.mu_0 * self.jw * sc.epsilon_0

        # Earth wave number - Equation (2.15) [1] (rad/m)
        self.k_earth2 = - self.jw * self.mur_1 * sc.mu_0 * (self.sigma_1 + self.jw * self.er_1 * sc.epsilon_0)

    def internal_impedance_elements_solid_wires(self):
        """ This method calculates the internal impedance of solid wires. """
        N = len(self.model.surfaces)
        Zi_approx = np.zeros((N, N), dtype=complex)
        Zi_bessel = np.zeros_like(Zi_approx)
        Zi_kelvin = np.zeros_like(Zi_approx)
        Zi_nahman = np.zeros_like(Zi_approx)
        Ri_cc = np.zeros_like(Zi_approx)
        Zi_hf = np.zeros_like(Zi_approx)
        Li_cc = np.zeros_like(Zi_approx)

        for conductor in self.model.surfaces:
            p = conductor['tag'] - 1
            ro = conductor['radius']
            q = np.sqrt(2) * ro / self.skin_depth[p]
            jw_mu = self.jw * self.model.mu[p]
            sigma = self.model.sigma[p]

            # Approximation Closed-Form Expression
            ri_cc = 1.0/(sigma * np.pi * ro**2)
            li_cc = self.model.mu[p] / (8 * np.pi)
            zi_hf = 1.0/(2 * np.pi * ro) * np.sqrt(jw_mu / sigma)
            
            Li_cc[p, p] = li_cc
            Ri_cc[p, p] = ri_cc
            Zi_hf[p, p] = zi_hf
            Zi_approx[p, p] = np.sqrt(ri_cc ** 2 + zi_hf ** 2)
            Zi_nahman[p, p] = ri_cc + zi_hf

            # Exact Expression with Modified Bessel Functions
            bessel_arg = np.sqrt(jw_mu * sigma) * ro
            
            # Prevenir erro em DC (f=0), onde o argumento é 0
            if np.abs(bessel_arg) < 1e-9:
                Zi_bessel[p, p] = ri_cc
            else:
                Zi_bessel[p, p] = zi_hf * ss.iv(0, bessel_arg) / ss.iv(1, bessel_arg)

            # Caso DC (frequência zero)
            if q < 1e-6:
                Ri_val = ri_cc
                wLi_val = 0
            else:
                # ss.kelvin(q) retorna uma tupla de números complexos
                Be, Ke, Bep, Kep = ss.kelvin(q)
                scaling_factor = ri_cc * (q / 2) / (Bep.imag**2 + Bep.real**2)
                
                # Fórmula da Resistência (ca) - Parte Real
                numerador_R = Be.real * Bep.imag - Be.imag * Bep.real
                Ri_val = scaling_factor * numerador_R
                
                # Fórmula da Reatância (ca) - Parte Imaginária
                numerador_wL = Be.real * Bep.real + Be.imag * Bep.imag
                wLi_val = scaling_factor * numerador_wL

            # Reconstrói a impedância complexa
            Zi_kelvin[p, p] = Ri_val + 1j * wLi_val

        return {
            'Zi_bessel': Zi_bessel,
            'Zi_kelvin': Zi_kelvin,
            'Zi_approx': Zi_approx,
            'Zi_nahman': Zi_nahman,
            'Zi_hf': Zi_hf,
            'Ri_cc': Ri_cc,
            'Li_cc': Li_cc,
        }

    def internal_impedance_elements_tubular_wires(self):
        """ This method calculates the internal impedance of tubular wires. """
        N = len(self.model.surfaces)
        Zi_approx = np.zeros((N, N), dtype=complex)
        Zi_bessel_int = np.zeros_like(Zi_approx)
        Zi_bessel_ext = np.zeros_like(Zi_approx)
        Zi_kelvin = np.zeros_like(Zi_approx)
        Zi_nahman = np.zeros_like(Zi_approx)
        Ri_cc = np.zeros_like(Zi_approx)
        Zi_hf = np.zeros_like(Zi_approx)
        Li_cc = np.zeros_like(Zi_approx)

        for key, conductor in self.model.mtl.items():
            p = key - 1
            ri, ro = conductor['radius']
            q = np.sqrt(2) * ro / self.skin_depth[p]
            mu = self.model.mu[p]
            jw_mu = self.jw * mu
            mu_2pi = mu / (2 * np.pi)
            sigma = self.model.sigma[p]
            sqrt_jw = np.sqrt(jw_mu / sigma)

            # Approximation Closed-Form Expression
            ri_cc = 1.0/(sigma * np.pi * (ro**2 - ri**2))
            ro2, ri2 = ro**2, ri**2

            # Evita divisão por zero se o condutor não for tubular
            if ro > ri:
                ro2ri2 = ro2 - ri2
                term1 = (ro2 - 3 * ri2) / (4 * ro2ri2)
                term2 = (ri2**2 / ro2ri2**2) * np.log(ro / ri)
                li_cc = mu_2pi * (term1 + term2)
            else:
                li_cc = 0

            
            # --- IMPLEMENTAÇÃO DAS FÓRMULAS CORRIGIDAS ---
            zi_hf_ext = 1.0 / (2 * np.pi * ro) * sqrt_jw if ro > 0 else np.inf
            zi_hf_int = 1.0 / (2 * np.pi * ri) * sqrt_jw if ri > 0 else np.inf
            kappa = np.sqrt(jw_mu * sigma)
            
            # Caso de baixa frequência (DC)
            if np.abs(kappa * ro) < 1e-9:
                Zi_bessel_ext[p, p] = ri_cc
                Zi_bessel_int[p, p] = ri_cc
            else:
                kro, kri = kappa * ro, kappa * ri

                den =     ss.kv(1, kri) * ss.iv(1, kro) - ss.iv(1, kri) * ss.kv(1, kro)
                num_ext = ss.iv(0, kro) * ss.kv(1, kri) + ss.iv(1, kri) * ss.kv(0, kro)
                num_int = ss.iv(0, kri) * ss.kv(1, kro) + ss.iv(1, kro) * ss.kv(0, kri)

                if abs(den) < 1e-15:
                    Zi_bessel_ext[p, p] = zi_hf_ext
                    Zi_bessel_int[p, p] = zi_hf_int if ri > 0 else np.inf
                else:
                    # Z'_i com retorno externo
                    Zi_bessel_ext[p, p] = zi_hf_ext * (num_ext / den)

                    # Z'_i com retorno interno
                    if ri > 0:
                        Zi_bessel_int[p, p] = zi_hf_int * (num_int / den)
                    else:
                        Zi_bessel_int[p, p] = np.inf
            
            Li_cc[p, p] = li_cc
            Ri_cc[p, p] = ri_cc
            Zi_hf[p, p] = zi_hf_ext
            Zi_approx[p, p] = np.sqrt(ri_cc ** 2 + zi_hf_ext ** 2)
            Zi_nahman[p, p] = ri_cc + zi_hf_ext

            # Caso DC (frequência zero)
            if q < 1e-6:
                Ri_val = ri_cc
                wLi_val = 0
            else:
                # ss.kelvin(q) retorna uma tupla de números complexos
                Be, Ke, Bep, Kep = ss.kelvin(q)
                scaling_factor = ri_cc * (q / 2) / (Bep.imag**2 + Bep.real**2)
                
                # Fórmula da Resistência (ca) - Parte Real
                numerador_R = Be.real * Bep.imag - Be.imag * Bep.real
                Ri_val = scaling_factor * numerador_R
                
                # Fórmula da Reatância (ca) - Parte Imaginária
                numerador_wL = Be.real * Bep.real + Be.imag * Bep.imag
                wLi_val = scaling_factor * numerador_wL

            # Reconstrói a impedância complexa
            Zi_kelvin[p, p] = Ri_val + 1j * wLi_val

        return {
            'Zi_bessel': {'internal_return': Zi_bessel_int, 'external_return': Zi_bessel_ext},
            'Zi_kelvin': Zi_kelvin,
            'Zi_approx': Zi_approx,
            'Zi_nahman': Zi_nahman,
            'Zi_hf': Zi_hf,
            'Ri_cc': Ri_cc,
            'Li_cc': Li_cc,
        }

    def external_impedance_term(self):
        """ This method calculates the impedance matrix of the earth return path. """
        d_matrix = self.model.d_matrix_ground_return
        D_matrix = self.model.D_matrix_ground_return
        return (np.log(D_matrix) - np.log(d_matrix))

    def external_admittance(self, type_form='potentials'):
        """
        This method calculates the external admittance matrix of the system using a 
        vectorized approach for the Maxwell Potential Coefficients matrix (Pe).
        """

        N = len(self.model.surfaces)
        Ge = np.zeros((N, N))

        d_matrix = self.model.d_matrix_ground_return
        D_matrix = self.model.D_matrix_ground_return
        Pe = 1 / (2 * np.pi * sc.epsilon_0) * (np.log(D_matrix) - np.log(d_matrix))

        # External capacitance matrix
        if type_form == 'potentials':
            Ce = np.linalg.inv(Pe)

        elif type_form == 'indirect':
            Le = self.external_impedance_term() / self.jw
            Ce = np.linalg.inv(Le) * sc.mu_0 * sc.epsilon_0

        return Ge + self.jw * Ce

    def simplified_soil_admittance(self):
        """ This method calculates the external admittance matrix of the system. """
        Zg = self.earth_return_impedance(type_form='carson')
        return -self.k_earth2 * np.linalg.inv(Zg)

    def pul_extended_theory(self, form='nakagawa'):
        """ This method calculates the impedance matrix of the earth return path. """
        N = len(self.model.surfaces)
        S1 = np.zeros((N, N), dtype=complex)
        S2 = np.zeros((N, N), dtype=complex)
        T = np.zeros((N, N), dtype=complex)

        # Complex depth (m)
        if form == 'sunde_log':
            p_dot = 1 / np.sqrt(-self.k_earth2)

        if form == 'carson' or form == 'deri':
            # Soil wave number (rad/m)
            k_e2 = - self.jw * self.mur_1 * sc.mu_0 * self.sigma_1

            # A. Deri Complex depth (m)
            p_dot = 1 / np.sqrt(-k_e2)

        for conductor_n in self.model.surfaces:
            for conductor_m in self.model.surfaces:
                # Get conductor tags
                n, m = conductor_n['tag']-1, conductor_m['tag']-1

                # Distance between the conductors (m)
                hn_hm = self.model.vertical_separation_matrix[n, m]
                dn_dm = self.model.horizontal_separation_matrix[n, m]

                # Quasi-TEM Integral Equation
                if form == 'quasi_tem':
                    S1[n, m] = 2 * sn_sommerfeld(hn_hm, dn_dm, self.k_earth2, self.k_air2)
                    S2[n, m] = 2 * sn_sommerfeld(hn_hm, dn_dm, self.k_earth2, self.k_air2, s_form='s2')
                    T[n, m] = 2 * t_sommerfeld(hn_hm, dn_dm, self.k_earth2, self.k_air2)

                # Quasi-TEM Logarithmic Approximation
                elif form == 'quasi_tem_log':
                    eta = np.sqrt(self.k_air2 - self.k_earth2)
                    eta_sqrt = eta * np.sqrt(hn_hm**2 + dn_dm**2)
                    n2 = self.k_earth2 / self.k_air2

                    log_s1 = 1 + 2 / eta_sqrt
                    log_s2 = 1 + (n2 + 1) / eta_sqrt
                    log_t = 1 + 2*(n2 + 1) / eta_sqrt

                    S1[n, m] = np.log(log_s1)
                    S2[n, m] = 2 / (n2 + 1) * np.log(log_s2)
                    T[n, m] = 2 * np.log(2) + (2 * n2 / (n2 + 1)) * np.log(log_s2/log_t)

                # Nakagawa Integral Equation
                elif form == 'nakagawa':
                    # soil refractive index
                    n2_naka = self.k_earth2 / self.k_air2

                    # Earth wave number (rad/m)
                    ke2_naka = - self.jw * self.mur_1 * sc.mu_0 * \
                        (self.sigma_1 + self.jw * (self.er_1 - 1) * sc.epsilon_0)

                    S1[n, m] = 2 * sn_sommerfeld_low_frequencies(hn_hm, dn_dm, ke2_naka)
                    S2[n, m] = 2 * sn_sommerfeld_low_frequencies(hn_hm, dn_dm, ke2_naka, n2=n2_naka)

                # Sunde Integral Equation
                elif form == 'sunde':
                    S1[n, m] = 2 * sn_sommerfeld_low_frequencies(hn_hm, dn_dm, self.k_earth2)

                # Carson Integral Equation
                elif form == 'carson':
                    S1[n, m] = 2 * sn_sommerfeld_low_frequencies(hn_hm, dn_dm, k_e2)

                # Log. Approximation Closed-Form Expression
                else:
                    M = np.zeros((N, N), dtype=complex)
                    if n == m:
                        hn = conductor_n['center_point'][1]
                        S1[n, m] = np.log((hn + p_dot) / hn)

                    else:
                        num = np.sqrt((hn_hm + 2*p_dot)**2 + dn_dm**2)
                        den = np.sqrt(hn_hm**2 + dn_dm**2)
                        S1[n, m] = np.log(num / den)

        # Approx. Internal Impedance Matrix
        zi = self.internal_impedance_elements_solid_wires()['Zi_approx']

        # External Impedance term
        M = self.external_impedance_term()
        ze = self.jw_mu0_2pi * M

        # Earth return impedance and admittance (ohm/m)
        if form == 'quasi_tem' or form == 'quasi_tem_log':
            zg = self.jw_mu0_2pi * (S1 - (T + S2))
            ysh = self.jw_2pi_e0 * np.linalg.inv(M - T)
        else:
            zg = self.jw_mu0_2pi * (S1)
            ysh = self.jw_2pi_e0 * np.linalg.inv(M + S2)

        # Series Impedance Matrix
        zs = zi + ze + zg 

        # Currents and Voltages Propagation Constants Matrix
        gamma_i = linalg.sqrtm(ysh @ zs)
        gamma_v = linalg.sqrtm(zs @ ysh)

        # Inverse of Y using LU decomposition
        # Solve the system Y * Y_inv = I to find Y_inv
        y_inverse = linalg.lu_solve(linalg.lu_factor(ysh), np.identity(ysh.shape[0]))

        # Inverse of Z using LU decomposition
        # Solve the system Z * Z_inv = I to find Z_inv
        z_inverse = linalg.lu_solve(linalg.lu_factor(zs), np.identity(zs.shape[0]))

        # Characteristic impedance Zc = Y_inv * sqrt(Y*Z)
        zc = y_inverse @ gamma_i 

        # Characteristic admittance Yc = Z_inv * sqrt(Z*Y)
        yc = z_inverse @ gamma_v

        return {
            'zi': zi, 'ze': ze, 'zg': zg, 'zs': zs, 'ysh': ysh,
            'gamma_i': gamma_i, 'gamma_v': gamma_v, 'zc': zc, 'yc': yc
        }
