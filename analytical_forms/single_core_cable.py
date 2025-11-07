"""
REFERENCES:
[1] XUE, Haoyan. General Formulation and Accurate Evaluation of Earth-Return Parameters
    for Overhead / Underground Cables. PhD thesis, Department of Electrical Engineering,
    École Polytechnique de Montréal, Université de Montréal, August 2018. 

[2] A. Ametani, T. Ohno and N. Nagaoka, Cable System Transients: Theory, Modeling and 
    Simulation, Wiley-IEEE Press, 2015.

[3] Y. Yin, "Calculation of frequency-dependent parameters of underground power cables with 
        finite element method," PhD, Electrical and Computer Engineering, The University of 
        British Columbia, Ph.D. thesis, 1990

[4] I. Lafaia, N. Alatawneh, J. Mahseredjian, A. Ametani, M. T. Correia de Barros, I. Koçar and 
    A. Naud, "Modeling of an Underground Cable Installed in a Poly-Ethylene Tube for Transient
    Simulations," in IEEJ Power and Energy Society Conference, Nagoya, 2015. 
"""

import numpy as np
import copy
import scipy.special as ss
import scipy.constants as sc
from scipy.linalg import lu_factor, lu_solve
from mtl_main.source import MulticonductorTransmissionLine
from utils.case_utils import *

def sommerfeld_quasi_tem_approx_impedance(hnm, dnm, ke2, ka2, type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

    # --- Integrands are now broadcast-compatible ---
    def _int_transformed_real(t, hnm, dnm, ke2, ka2):
        # Reshape t to a column vector (150, 1) to broadcast with frequency vectors (40,)
        x = np.tan(t[:, np.newaxis])
        sqrt_x2_ke2 = np.sqrt(x**2 - ke2)
        sqrt_x2_ka2 = np.sqrt(x**2 - ka2)
        num = np.exp(hnm * sqrt_x2_ke2) * np.cos(dnm * x)
        den = sqrt_x2_ka2 + sqrt_x2_ke2
        return (num/den).real * (1 / np.cos(t[:, np.newaxis])**2)

    def _int_transformed_imag(t, hnm, dnm, ke2, ka2):
        # Reshape t to a column vector (150, 1)
        x = np.tan(t[:, np.newaxis])
        sqrt_x2_ke2 = np.sqrt(x**2 - ke2)
        sqrt_x2_ka2 = np.sqrt(x**2 - ka2)
        num = np.exp(hnm * sqrt_x2_ke2) * np.cos(dnm * x)
        den = sqrt_x2_ka2 + sqrt_x2_ke2
        return (num/den).imag * (1 / np.cos(t[:, np.newaxis])**2)

    # --- Integration helper is now vectorized ---
    def _int_gauss_legendre(func, a, b, n, *args):
        [x_gl, w] = np.polynomial.legendre.leggauss(n)
        t = 0.5 * (x_gl + 1) * (b - a) + a
        # Reshape weights w to a column vector (150, 1)
        # func(t, *args) returns a (150, 40) array
        # Sum along the integration axis (axis=0)
        return np.sum(w[:, np.newaxis] * func(t, *args), axis=0) * 0.5 * (b - a)

    # Calculate the real and imaginary parts of the integral using Gauss-Legendre
    if type_form == 'gauss_legendre':
        gauss_legendre_real = _int_gauss_legendre(_int_transformed_real, 0, np.pi/2, pts, hnm, dnm, ke2, ka2)
        gauss_legendre_imag = _int_gauss_legendre(_int_transformed_imag, 0, np.pi/2, pts, hnm, dnm, ke2, ka2)
        Sn = gauss_legendre_real + 1j * gauss_legendre_imag
    else:
        # Note: Scipy quad is not vectorized and would require a loop.
        raise NotImplementedError("Scipy Quad is not supported in vectorized mode.")

    return Sn

def sommerfeld_ametani_approx(hnm, dnm, ke2, type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

    def _int_transformed_real(t, hnm, dnm, ke2):
        x = np.tan(t[:, np.newaxis])
        num = np.exp(hnm * x) * np.cos(dnm * x)
        den = x + np.sqrt(x**2 - ke2)
        return (num/den).real * (1 / np.cos(t[:, np.newaxis])**2)

    def _int_transformed_imag(t, hnm, dnm, ke2):
        x = np.tan(t[:, np.newaxis])
        num = np.exp(hnm * x) * np.cos(dnm * x)
        den = x + np.sqrt(x**2 - ke2)
        return (num/den).imag * (1 / np.cos(t[:, np.newaxis])**2)

    def _int_gauss_legendre(func, a, b, n, *args):
        [x_gl, w] = np.polynomial.legendre.leggauss(n)
        t = 0.5 * (x_gl + 1) * (b - a) + a
        return np.sum(w[:, np.newaxis] * func(t, *args), axis=0) * 0.5 * (b - a)

    if type_form == 'gauss_legendre':
        gauss_legendre_real = _int_gauss_legendre(_int_transformed_real, 0, np.pi/2, pts, hnm, dnm, ke2)
        gauss_legendre_imag = _int_gauss_legendre(_int_transformed_imag, 0, np.pi/2, pts, hnm, dnm, ke2)
        Sn = gauss_legendre_real + 1j * gauss_legendre_imag
    else:
        raise NotImplementedError("Scipy Quad is not supported in vectorized mode.")

    return Sn

def sommerfeld_quasi_tem_approx_admittance(hnm, dnm, ke2, ka2, type_form='gauss_legendre', pts=150):
    """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

    def _int_transformed_real(t, hnm, dnm, ke2, ka2):
        x = np.tan(t[:, np.newaxis])
        sqrt_ke2 = np.sqrt(x**2 - ke2)
        sqrt_ka2 = np.sqrt(x**2 - ka2)
        num = sqrt_ka2 * np.exp(hnm * sqrt_ke2) * np.cos(dnm * x)
        den = sqrt_ke2 * (sqrt_ka2 + (ka2 / ke2) * sqrt_ke2)
        return (num/den).real * (1 / np.cos(t[:, np.newaxis])**2)

    def _int_transformed_imag(t, hnm, dnm, ke2, ka2):
        x = np.tan(t[:, np.newaxis])
        sqrt_ke2 = np.sqrt(x**2 - ke2)
        sqrt_ka2 = np.sqrt(x**2 - ka2)
        num = sqrt_ka2 * np.exp(hnm * sqrt_ke2) * np.cos(dnm * x)
        den = sqrt_ke2 * (sqrt_ka2 + (ka2 / ke2) * sqrt_ke2)
        return (num/den).imag * (1 / np.cos(t[:, np.newaxis])**2)

    def _int_gauss_legendre(func, a, b, n, *args):
        [x_gl, w] = np.polynomial.legendre.leggauss(n)
        t = 0.5 * (x_gl + 1) * (b - a) + a
        return np.sum(w[:, np.newaxis] * func(t, *args), axis=0) * 0.5 * (b - a)

    if type_form == 'gauss_legendre':
        gauss_legendre_real = _int_gauss_legendre(_int_transformed_real, 0, np.pi/2, pts, hnm, dnm, ke2, ka2)
        gauss_legendre_imag = _int_gauss_legendre(_int_transformed_imag, 0, np.pi/2, pts, hnm, dnm, ke2, ka2)
        Sn = gauss_legendre_real + 1j * gauss_legendre_imag
    else:
        raise NotImplementedError("Scipy Quad is not supported in vectorized mode.")

    return Sn

def coth(x):
    """
    Calculates the hyperbolic cotangent of x element-wise.

    The function is defined as: coth(x) = 1 / tanh(x)

    Args:
        x (float or numpy.ndarray): The input value or array.

    Returns:
        float or numpy.ndarray: The hyperbolic cotangent of the input.
    """
    return 1 / np.tanh(x)

def cosech(x):
    """
    Calculates the hyperbolic cosecant of x element-wise.

    The function is defined as: cosech(x) = 1 / sinh(x)

    Args:
        x (float or numpy.ndarray): The input value or array.

    Returns:
        float or numpy.ndarray: The hyperbolic cosecant of the input.
    """
    return 1 / np.sinh(x)

def capacitance_matrix_from_energy_method(energy_vector: np.ndarray, v0: float = 1.0) -> np.ndarray:
    """
    Calculates the capacitance matrix from the energy stored in the electric field.
    Application: Single-core cable with core and sheath conductors.

    Args:
        energy_vector (numpy.ndarray): A 1D array of energy values for each conductor (in Joules).
        v0 (float): The reference voltage (in Volts).

    Returns:
        numpy.ndarray: The capacitance matrix (C) in Farads.
    """
    # Ensure the energy vector has the correct number of elements
    if len(energy_vector) != 3:
        raise ValueError("The energy_vector must contain exactly three elements: [W11, W22, W12].")

    # The pre-calculated inverse of matrix A from Equation (4) [2]
    A_inv = np.array([
        [ 1.0,  0.0,  0.0],
        [ 0.0,  1.0,  0.0],
        [-0.5, -0.5,  0.5]
    ])

    # Define the right-hand side vector B from Equation (4) [2]
    B = (4 / v0**2) * np.asarray(energy_vector)

    # Calculate capacitance components directly using matrix-vector multiplication
    # C_components = A_inv * B
    capacitance_vector = A_inv @ B

    # Extract the individual capacitance values
    C11 = capacitance_vector[0]  # Self-capacitance of the core
    C22 = capacitance_vector[1]  # Self-capacitance of the sheath
    C12 = capacitance_vector[2]  # Magnitude of the mutual capacitance

    # Assemble the final 2x2 symmetric capacitance matrix.
    # The off-diagonal mutual capacitance terms are negative.
    return np.array([[C11, C12], [C12, C22]])   

class EquivalentRadiiSystems:
    def __init__(self, model: MulticonductorTransmissionLine):
        """
        Initializes the vectorized calculator.

        Args:
            model (MulticonductorTransmissionLine): The MTL geometry model.
        """
        self.model = model

        # Identify each enclosure by its name
        scc = model.scc if hasattr(model, 'scc') else None
        core_enclosure, sheath_enclosure, armor_enclosure = None, None, None
        for key, enclosure_data in scc['hdpe'].items():
            if key == 'core':
                core_enclosure = enclosure_data
            elif key == 'sheath':
                sheath_enclosure = enclosure_data
            elif key == 'armor':
                armor_enclosure = enclosure_data

        if 'core_outer_radius' in scc:
            print("Core enclosure detected.")
            self.rho1, self.mu1 = scc['core_resistivity'], scc['core_permeability']
            self.r1, self.r2 = scc['core_inner_radius'], scc['core_outer_radius']

        if 'core_insulation_outer_radius' in scc:
            print("Core insulation detected.")
            self.ei1, self.mui1 = scc['core_insulation_permittivity'], scc['core_insulation_permeability']
            self.r3 = scc['core_insulation_outer_radius']

        if 'sheath_outer_radius' in scc:
            print("Sheath enclosure detected.")
            self.rho2, self.mu2 = scc['sheath_resistivity'], scc['sheath_permeability']
            self.r3, self.r4 = scc['sheath_inner_radius'], scc['sheath_outer_radius']

        if 'sheath_insulation_outer_radius' in scc:
            print("Sheath insulation detected.")
            self.ei2, self.mui2 = scc['sheath_insulation_permittivity'], scc['sheath_insulation_permeability']
            self.r5 = scc['sheath_insulation_outer_radius']

        if sheath_enclosure is not None:
            print("Sheath air-gap enclosure detected.")
            self.r6, self.r7 = sheath_enclosure['inner_radius'], sheath_enclosure['outer_radius']
            self.ei3, self.mui3 = sheath_enclosure['permittivity'], sc.mu_0

        if 'armor_outer_radius' in scc:
            print("Armor enclosure detected.")
            self.rho3, self.mu3 = scc['armor_resistivity'], scc['armor_permeability']
            self.r5, self.r6 = scc['armor_inner_radius'], scc['armor_outer_radius']

        if 'armor_insulation_outer_radius' in scc:
            print("Armor insulation detected.")
            self.ei3, self.mui3 = scc['armor_insulation_permittivity'], scc['armor_insulation_permeability']
            self.r7 = scc['armor_insulation_outer_radius']
    
    def concentric_insulators(self):
        # A fórmula para a capacitância de um cilindro coaxial é C = 2*pi*epsilon / ln(r_externo / r_interno)

        # C_S: Capacitância da isolação da bainha do cabo
        C_S = (2 * np.pi * self.ei2) / np.log(self.r5 / self.r4) if self.r5 > self.r4 else np.inf

        # C_a: Capacitância do entreferro de ar
        C_a = (2 * np.pi * sc.epsilon_0) / np.log(self.r6 / self.r5) if self.r6 > self.r5 else np.inf

        # C_H: Capacitância do duto de PEAD
        C_H = (2 * np.pi * self.ei3) / np.log(self.r7 / self.r6) if self.r7 > self.r6 else np.inf

        # --- Etapa 3: Calcular a capacitância total C0 para a conexão em série ---
        # C0 = (1/C_S + 1/C_a + 1/C_H)^-1 Eq. 9 [4]
        inv_C_S = 1 / C_S if C_S != np.inf else 0
        inv_C_a = 1 / C_a if C_a != np.inf else 0
        inv_C_H = 1 / C_H if C_H != np.inf else 0
        C0 = 1 / (inv_C_S + inv_C_a + inv_C_H)

        # --- Etapa 4: Calcular a permissividade relativa equivalente eps_a ---
        eps_a_relative = C0 * np.log(self.r7 / self.r4) / (2 * np.pi * sc.epsilon_0)

        print("\n======== Equivalent Radii Systems (ERS) Results ========")
        print("r1:", self.r1, "r2:", self.r2, "r3:", self.r3, "r4:", self.r4, "r5:", self.r5, "r6:", self.r6, "r7:", self.r7)
        print("eri1:", self.ei1/sc.epsilon_0, "eri2:", self.ei2/sc.epsilon_0, "eri3:", self.ei3/sc.epsilon_0)
        print(f"C_S: {C_S:.4e}, C_a: {C_a:.4e}, C_H: {C_H:.4e}")
        print(f"Total Capacitance C0: {C0:.4e}")
        print(f"Equivalent Relative Permittivity eps_a: {eps_a_relative:.4e}")
        print("===================================================================\n")

        return {'total_capacitance': C0, 'equivalent_relative_permittivity': eps_a_relative}

    def equiv_rel_permittivity_epsr_area_weighted(self, ShowInfo: bool = True) -> float:
        """
        Calculates the equivalent relative permittivity for "Case 2" from Lafaia, 2015 [4].

        This method uses an area-weighted average of the permittivities of the
        insulating layers that are replaced in this simplified model. The model
        assumes a single, homogeneous insulator extending from the metallic sheath (r4)
        to the outer radius of the HDPE tube (r7).

        Returns:
            float: The calculated equivalent relative permittivity (eps_a).
        """
        # Relative permittivities of the materials
        eri2 = self.ei2 / sc.epsilon_0  # Cable's outer insulation
        eri3 = self.ei3 / sc.epsilon_0  # HDPE tube material

        # --- Cross-sectional area of each layer ---
        area_out_insul = self.r5**2 - self.r4**2
        area_air = self.r6**2 - self.r5**2
        area_hdpe = self.r7**2 - self.r6**2
        total_area = self.r7**2 - self.r4**2

        # --- Area-weighted average permittivity ---
        weighted_sum = (area_out_insul * eri2) + (area_air) + (area_hdpe * eri3)        
        eps_a = weighted_sum / total_area

        # --- Print a summary of the calculation for clarity ---
        if ShowInfo:
            print("\n======== Area-Weighted Average Permittivity Results ========")
            print(f"Replacing layers from r4={self.r4*1000:.2f} mm to r7={self.r7*1000:.2f} mm")
            print("--------------------------------------------------------------------")
            print(f"Layer 1 (Insulation): Area/pi =  {area_out_insul*1e6:.2f} mm^2, eps_r = {eri2}")
            print(f"Layer 2 (Air Gap):    Area/pi = {area_air*1e6:.2f} mm^2, eps_r = {1}")
            print(f"Layer 3 (HDPE Tube):  Area/pi = {area_hdpe*1e6:.2f} mm^2, eps_r = {eri3}")
            print("--------------------------------------------------------------------")
            print(f"Total Area/pi: {total_area*1e6:.2f} mm^2")
            print(f"Calculated Equivalent Relative Permittivity (eps_a): {eps_a:.4f}")
            print("====================================================================\n")

        return {
            'sheath_outer_radius': self.r4,
            'sheath_insulation_outer_radius': self.r5,
            'sheath_enclosure_inner_radius': self.r6,
            'sheath_enclosure_outer_radius': self.r7,
            'equivalent_relative_permittivity': eps_a
        }
    
    def equivalent_parameters_from_gmd(self, ShowInfo: bool = True):
        """
        Calculates equivalent parameters for a cable in an HDPE tube using the
        Geometric Mean Distance (GMD) method as described by Lafaia, 2015 [4].

        This method implements the logic from Equations (5) through (11) in the reference paper.
        It transforms the eccentric geometry into an equivalent concentric one, adjusts
        permittivity to preserve capacitance, and calculates a final equivalent
        permittivity for the entire outer insulation system.

        Args:
            D1 (float): Inner diameter of the HDPE tube [m].
            D2 (float): Outer diameter of the HDPE tube [m].

        Returns:
            dict: A dictionary containing the key calculated parameters, including
                  the equivalent radii, modified permittivity, total capacitance,
                  and the final equivalent relative permittivity (eps_a).
        """
        # --- Step 1: Get necessary parameters from the model ---
        two_pi_eo = 2 * np.pi * sc.epsilon_0
        r4, r5 = self.r4, self.r5
        D1, D2 = 2 * self.r6, 2 * self.r7
        
        # --- Equivalent outer radii r6 (the air gap) and r7 (HDPE) using GMD (Eq. 6 & 7) ---
        r6_gmd = np.exp(np.log(D1 / 2) + (2 * r5 / D1)**2 / 2 - 0.5)
        r7_gmd = np.exp(np.log(D2 / 2) + (2 * r5 / D2)**2 / 2 - 0.5)

        # --- Modify HDPE permittivity to preserve capacitance (Eq. 8) ---
        eps_3_prime = self.ei3 / sc.epsilon_0 * (np.log(r7_gmd / r6_gmd) / np.log(D2 / D1))

        # --- Total capacitance as a series of concentric layers (Eq. 9) ---
        # C_S: Capacitance of the cable's sheath insulator
        cs = 2 * np.pi * self.ei2 / np.log(r5 / r4)

        # C_a: Capacitance of the equivalent concentric air gap
        ca = two_pi_eo / np.log(r6_gmd / r5)

        # C_H: Capacitance of the equivalent concentric HDPE tube with modified permittivity
        ch = two_pi_eo * eps_3_prime / np.log(r7_gmd / r6_gmd)

        # Total series capacitance C0
        c0 = 1 / (1 / cs + 1 / ca + 1 / ch)

        # --- Final equivalent permittivity eps_a (Eq. 11) ---
        # This represents a single insulator from r4 to r5 with the total capacitance C0
        eps_a_31 = c0 * np.log(r5 / r4) / two_pi_eo
        eps_a_32 = c0 * np.log(r6_gmd / r4) / two_pi_eo
        eps_a_33 = c0 * np.log(r7_gmd / r4) / two_pi_eo

        # --- Print and return results for analysis ---
        if ShowInfo:
            print("\n======== GMD-Based Equivalent Systems Results (Lafaia, 2015) ========")
            print(f"Original HDPE Dims [mm]: D1={D1*1000:.2f}, D2={D2*1000:.2f}")
            print(f"Cable Outer Radius r5 [mm]: {r5*1000:.2f}")
            print("-------------------------------------------------------------------")
            print(f"GMD Air Gap Outer Radius r6 [mm]: {r6_gmd*1000:.2f}")
            print(f"GMD HDPE Tube Outer Radius r7 [mm]: {r7_gmd*1000:.2f}")
            print(f"Modified HDPE Relative Permittivity eps_3': {eps_3_prime:.4f}")
            print("-------------------------------------------------------------------")
            print(f"C_S: {cs:.4e} F/m, C_a: {ca:.4e} F/m, C_H: {ch:.4e} F/m")
            print(f"Total Series Capacitance C0: {c0:.4e} F/m")
            print(f"Final Equivalent Relative Permittivity, eps_a:")
            print(f"Case 3.1: {eps_a_31:.4f}")
            print(f"Case 3.2: {eps_a_32:.4f}")
            print(f"Case 3.3: {eps_a_33:.4f}")
            print("===================================================================\n")

        return {
            'air_gap_outer_radius': r6_gmd,
            'hdpe_outer_radius': r7_gmd,
            'modified_hdpe_rel_permittivity': eps_3_prime,
            'equivalent_outer_capacitance': c0,
            'equivalent_relative_permittivity': {
                'case 3.1': eps_a_31,
                'case 3.2': eps_a_32,
                'case 3.3': eps_a_33,
            }
        }

class InternalPerUnitParameters:
    """
    This class calculates vector-frequency PUL parameters using an MTL geometry model.
    The calculations are vectorized over the frequency axis for efficiency.
    """

    def __init__(self, model: MulticonductorTransmissionLine, frequencies: np.ndarray):
        """
        Initializes the vectorized calculator.

        Args:
            model (MulticonductorTransmissionLine): The MTL geometry model.
            f (np.ndarray): A NumPy array of frequencies to be calculated.
        """
        self.model = model
        self.f = np.asarray(frequencies)  # Ensure f is a NumPy array
        self.jw = 1j * 2 * np.pi * self.f

    def parameters_by_bessel(self):
        """
        Calculates the internal impedance matrix components for a single-core cable (SCC)
        over a vector of frequencies.

        The method uses formulas for tubular conductors, which involve modified
        Bessel functions, to account for skin and proximity effects within the
        conductors. The use of NumPy and SciPy's vectorized functions allows
        for efficient calculation across all frequencies simultaneously.

        The returned impedance components (z11, z2m, etc.) are NumPy arrays,
        where each element corresponds to a frequency in the input vector `f`.

        Reference: A. Ametani, T. Ohno and N. Nagaoka, Cable System Transients: Theory, Modeling and
                    Simulation, Wiley-IEEE Press, 2015.
        """
        s = self.jw
        scc = self.model.scc
        two_pi = 2 * np.pi

        # --- Initialize all impedance components ---
        z11, z12 = (None,) * 2
        z2i, z20, z23, z2m = (None,) * 4
        z3i, z30, z34, z3m = (None,) * 4
        pcj, psj, paj = (None,) * 3

        if 'core_outer_radius' in scc:
            rho1, mu1 = scc['core_resistivity'], scc['core_permeability']
            r1, r2 = scc['core_inner_radius'], scc['core_outer_radius']
            m_core = np.sqrt(s * mu1 / rho1)
            x1, x2 = m_core * r1, m_core * r2

            # --- z11: internal impedance of core outer surface ---
            if np.isclose(r1, 0):
                with np.errstate(divide='ignore', invalid='ignore'):
                    z11 = (m_core * rho1 / (two_pi * r2)) * (ss.iv(0, x2) / ss.iv(1, x2))
                    # Use the built-in complex type, not np.complex
                    invalid_mask = np.isclose(x2, 0) | np.isnan(x2)
                    if z11.ndim > 0:
                        z11[invalid_mask] = complex(np.inf, np.inf) # Use complex()
                    elif invalid_mask:
                        z11 = complex(np.inf, np.inf) # Use complex()
            
            else:
                with np.errstate(divide='ignore', invalid='ignore'):
                    D1 = ss.iv(1, x2) * ss.kv(1, x1) - ss.iv(1, x1) * ss.kv(1, x2)
                    N1 = ss.iv(0, x2) * ss.kv(1, x1) + ss.kv(0, x2) * ss.iv(1, x1)
                    z11 = (s * mu1 / two_pi) * (1 / (x2 * D1)) * N1
                    
                    # Create a mask for invalid conditions (nan or near-zero denominator)
                    invalid_mask = np.isnan(D1) | np.isclose(D1, 0)
                    if z11.ndim > 0:
                        z11[invalid_mask] = complex(np.inf, np.inf) # Use complex()
                    elif invalid_mask:
                        z11 = complex(np.inf, np.inf) # Use complex()

        if 'core_insulation_outer_radius' in scc:
            mui1 = scc['core_insulation_permeability']
            ei1 = scc['core_insulation_permittivity']
            r3 = scc['core_insulation_outer_radius']
            z12 = (s * mui1 / two_pi) * np.log(r3 / r2) if not np.isclose(r3, r2) else 0
            pcj = (1 / (two_pi * ei1)) * np.log(r3 / r2) if not np.isclose(r3, r2) else 0

        if 'sheath_outer_radius' in scc:
            rho2, mu2 = scc['sheath_resistivity'], scc['sheath_permeability']
            r3, r4 = scc['sheath_inner_radius'], scc['sheath_outer_radius']
            m_sheath = np.sqrt(s * mu2 / rho2)
            x3, x4 = m_sheath * r3, m_sheath * r4
            D2 = ss.iv(1, x4) * ss.kv(1, x3) - ss.iv(1, x3) * ss.kv(1, x4)
            with np.errstate(divide='ignore', invalid='ignore'):
                z2m = rho2 / (two_pi * r3 * r4 * D2)
            N2i = ss.iv(0, x3) * ss.kv(1, x4) + ss.kv(0, x3) * ss.iv(1, x4)
            with np.errstate(divide='ignore', invalid='ignore'):
                z2i = (s * mu2 / two_pi) * (1 / (x3 * D2)) * N2i
            N20 = ss.iv(0, x4) * ss.kv(1, x3) + ss.kv(0, x4) * ss.iv(1, x3)
            with np.errstate(divide='ignore', invalid='ignore'):
                z20 = (s * mu2 / two_pi) * (1 / (x4 * D2)) * N20

        if 'sheath_insulation_outer_radius' in scc:
            mui2 = scc['sheath_insulation_permeability']
            ei2 = scc['sheath_insulation_permittivity']
            r5 = scc['sheath_insulation_outer_radius']
            z23 = (s * mui2 / two_pi) * np.log(r5 / r4) if not np.isclose(r5, r4) else 0
            psj = (1 / (two_pi * ei2)) * np.log(r5 / r4) if not np.isclose(r5, r4) else 0

        if 'armor_outer_radius' in scc:
            rho3, mu3 = scc['armor_resistivity'], scc['armor_permeability']
            r5, r6 = scc['armor_inner_radius'], scc['armor_outer_radius']
            m_armor = np.sqrt(s * mu3 / rho3)
            x5, x6 = m_armor * r5, m_armor * r6
            D3 = ss.iv(1, x6) * ss.kv(1, x5) - ss.iv(1, x5) * ss.kv(1, x6)
            with np.errstate(divide='ignore', invalid='ignore'):
                z3m = rho3 / (two_pi * r5 * r6 * D3)
            N3 = ss.iv(0, x5) * ss.kv(1, x6) + ss.kv(0, x5) * ss.iv(1, x6)
            with np.errstate(divide='ignore', invalid='ignore'):
                z3i = (s * mu3 / two_pi) * (1 / (x5 * D3)) * N3
            N30 = ss.iv(0, x6) * ss.kv(1, x5) + ss.kv(0, x6) * ss.iv(1, x5)
            with np.errstate(divide='ignore', invalid='ignore'):
                z30 = (s * mu3 / two_pi) * (1 / (x6 * D3)) * N30

        if 'armor_insulation_outer_radius' in scc:
            mui3 = scc['armor_insulation_permeability']
            ei3 = scc['armor_insulation_permittivity']
            r7 = scc['armor_insulation_outer_radius']
            z34 = (s * mui3 / two_pi) * np.log(r7/r6) if not np.isclose(r7, r6) else 0
            paj = (1 / (two_pi * ei3)) * np.log(r7 / r6) if not np.isclose(r7, r6) else 0

        return {
            'zcs': {'z11': z11, 'z12': z12, 'z2i': z2i, 'Zcs': z11 + z12 + z2i},
            'zsa': {'z20': z20, 'z23': z23, 'z3i': z3i},
            'za4': {'z30': z30, 'z34': z34},
            'zs3': {'z20': z20, 'z23': z23},
            'z2m': z2m,
            'z3m': z3m,
            'potentials': {'pcj': pcj, 'psj': psj, 'paj': paj}
        }

    def parameters_approximation(self):
        """
        Calculates the internal impedance matrix components for a single-core cable (SCC)
        over a vector of frequencies.

        The method uses formulas for tubular conductors, which involve modified
        Bessel functions, to account for skin and proximity effects within the
        conductors. The use of NumPy and SciPy's vectorized functions allows
        for efficient calculation across all frequencies simultaneously.

        The returned impedance components (z11, z2m, etc.) are NumPy arrays,
        where each element corresponds to a frequency in the input vector `f`.

        Reference: A. Ametani, T. Ohno and N. Nagaoka, Cable System Transients: Theory, Modeling and
                    Simulation, Wiley-IEEE Press, 2015.
        """
        s = self.jw
        scc = self.model.scc
        two_pi = 2 * np.pi

        # --- Initialize all impedance components ---
        z11, z12 = (None,) * 2
        z2i, z20, z23, z2m = (None,) * 4
        z3i, z30, z34, z3m = (None,) * 4
        pcj, psj, paj = (None,) * 3

        if 'core_outer_radius' in scc:
            rho1, mu1 = scc['core_resistivity'], scc['core_permeability']
            r1, r2 = scc['core_inner_radius'], scc['core_outer_radius']
            m_core = np.sqrt(s * mu1 / rho1)
            x1, x2 = m_core * r1, m_core * r2
            pm = rho1 * m_core

            # --- z11: internal impedance of solid core outer surface ---
            if np.isclose(r1, 0):
                z11 = pm / (two_pi * r2) * coth(0.777 * x2) + 0.356 * rho1 / (np.pi * r2**2)
                
            # --- z11: internal impedance of tubular core outer surface ---
            else:
                with np.errstate(divide='ignore', invalid='ignore'):
                    # Fórmula análoga à de z20 (impedância externa da bainha)
                    z11 = pm / (two_pi * r2) * coth(m_core * (r2 - r1)) + rho1 / (two_pi * r2 * (r1 + r2))

        if 'core_insulation_outer_radius' in scc:
            mui1 = scc['core_insulation_permeability']
            ei1 = scc['core_insulation_permittivity']
            r3 = scc['core_insulation_outer_radius']

            # --- z12: core outer insulator impedance ---
            z12 = (s * mui1 / two_pi) * np.log(r3 / r2) if not np.isclose(r3, r2) else 0

            # --- pcj: core outer insulator potential coefficient ---
            pcj = (1 / (two_pi * ei1)) * np.log(r3 / r2) if not np.isclose(r3, r2) else 0

        if 'sheath_outer_radius' in scc:
            rho2, mu2 = scc['sheath_resistivity'], scc['sheath_permeability']
            r3, r4 = scc['sheath_inner_radius'], scc['sheath_outer_radius']
            m_sheath = np.sqrt(s * mu2 / rho2)
            pm = rho2 * m_sheath
                        
            # --- z2m: sheath mutual impedance ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z2m = pm / (np.pi * (r3 + r4)) * cosech(m_sheath * (r4 - r3))
            
            # --- z2i: internal impedance of sheath inner surface ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z2i = pm / (two_pi * r3) * coth(m_sheath * (r4 - r3)) - rho2 / (two_pi * r3 * (r3 + r4))
            
            # --- z20: internal impedance of sheath outer surface ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z20 = pm / (two_pi * r4) * coth(m_sheath * (r4 - r3)) + rho2 / (two_pi * r4 * (r3 + r4))

        if 'sheath_insulation_outer_radius' in scc:
            mui2 = scc['sheath_insulation_permeability']
            ei2 = scc['sheath_insulation_permittivity']
            r5 = scc['sheath_insulation_outer_radius']
            
            # --- z23: sheath outer insulator impedance ---
            z23 = (s * mui2 / two_pi) * np.log(r5 / r4) if not np.isclose(r5, r4) else 0

            # --- psj: sheath outer insulator potential coefficient ---
            psj = (1 / (two_pi * ei2)) * np.log(r5 / r4) if not np.isclose(r5, r4) else 0

        if 'armor_outer_radius' in scc:
            rho3, mu3 = scc['armor_resistivity'], scc['armor_permeability']
            r5, r6 = scc['armor_inner_radius'], scc['armor_outer_radius']
            m_armor = np.sqrt(s * mu3 / rho3)
            pm = rho3 * m_armor
            
            # --- z3m: armor mutual impedance ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z3m = pm / (np.pi * (r5 + r6)) * cosech(m_armor * (r6 - r5))
            
            # --- z3i: internal impedance of armor inner surface ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z3i = pm / (two_pi * r5) * coth(m_armor * (r6 - r5)) - rho3 / (two_pi * r5 * (r5 + r6))

            # --- z30: internal impedance of armor outer surface ---
            with np.errstate(divide='ignore', invalid='ignore'):
                z30 = pm / (two_pi * r6) * coth(m_armor * (r6 - r5)) + rho3 / (two_pi * r6 * (r5 + r6))

        if 'armor_insulation_outer_radius' in scc:
            mui3 = scc['armor_insulation_permeability']
            ei3 = scc['armor_insulation_permittivity']
            r7 = scc['armor_insulation_outer_radius']
            
            # --- z34: armor outer insulator impedance ---
            z34 = (s * mui3 / two_pi) * np.log(r7/r6) if not np.isclose(r7, r6) else 0

            # --- paj: armor outer insulator potential coefficient ---
            paj = (1 / (two_pi * ei3)) * np.log(r7 / r6) if not np.isclose(r7, r6) else 0

        return {
            'zcs': {'z11': z11, 'z12': z12, 'z2i': z2i, 'Zcs': z11 + z12 + z2i},
            'zsa': {'z20': z20, 'z23': z23, 'z3i': z3i},
            'za4': {'z30': z30, 'z34': z34},
            'zs3': {'z20': z20, 'z23': z23},
            'z2m': z2m,
            'z3m': z3m,
            'potentials': {'pcj': pcj, 'psj': psj, 'paj': paj}
        }

    def parameters_hybrid(self, transition_frequency=1e5):
        """
        Calculates internal parameters using a hybrid approach based on a single
        transition frequency.

        For frequencies below the transition_frequency, the results from the
        Bessel function formulation (`parameters_by_bessel`) are used.
        For frequencies at or above the transition_frequency, the results from
        the high-frequency approximation (`parameters_approximation`) are used.

        This provides a direct way to combine the accuracy of the Bessel model at
        low frequencies with the numerical stability of the approximation at high
        frequencies.

        Args:
            transition_frequency (float): The frequency (in Hz) at which to switch
                                          from the Bessel model to the approximation model.

        Returns:
            dict: A dictionary containing the calculated hybrid parameters, with the
                  same structure as the other methods.
        """
        # --- Step 1: Get results from both existing methods ---
        params_bessel = self.parameters_by_bessel()
        params_approx = self.parameters_approximation()

        # --- Step 2: Create a boolean mask based on the transition frequency ---
        # The mask is True for frequencies where the approximation should be used.
        use_approx_mask = self.f >= transition_frequency

        # --- Step 3: Create the hybrid results dictionary ---
        # Start with a copy of the Bessel results, then overwrite where necessary.
        params_hybrid = copy.deepcopy(params_bessel)

        # Iterate through the dictionary keys to apply the hybrid logic.
        # This approach is general and works for any cable configuration.
        for key, value in params_bessel.items():
            if isinstance(value, dict):
                # Handles nested dictionaries like 'zcs', 'zsa', etc.
                for sub_key, sub_value in value.items():
                    # Check if the item is a frequency-dependent numpy array
                    if isinstance(sub_value, np.ndarray) and sub_value.shape == self.f.shape:
                        params_hybrid[key][sub_key] = np.where(
                            use_approx_mask,
                            params_approx[key][sub_key], # Value if True (use approx)
                            sub_value                     # Value if False (use bessel)
                        )
            elif isinstance(value, np.ndarray) and value.shape == self.f.shape:
                # Handles top-level items like 'z2m', 'z3m'
                params_hybrid[key] = np.where(
                    use_approx_mask,
                    params_approx[key], # Value if True (use approx)
                    value               # Value if False (use bessel)
                )
        
        return params_hybrid
    
    def matrices(self, internal_form='hybrid'):
        """
        Assembles the full internal impedance [Zi] and shunt admittance [Ye] matrices
        for all specified frequencies.

        Returns:
            dict: A dictionary containing the calculated matrices.
                  'impedance_matrix' (Zi) and 'shunt_admittance_matrix' (Ye) are
                  3D NumPy arrays with shape (num_frequencies, num_total_conductors, num_total_conductors).
        """
        N, M = self.model.num_sc_cables, self.model.num_conductors_per_scc
        num_total_conductors = N * M
        
        if internal_form == 'bessel':
            zij = self.parameters_by_bessel()
        elif internal_form == 'approximation':
            zij = self.parameters_approximation()
        elif internal_form == 'hybrid':
            zij = self.parameters_hybrid(transition_frequency=1e5)
        else:
            raise ValueError("Invalid internal_form. Choose 'bessel' or 'approximation' or 'hybrid'.")

        # 1. SCC with core, core_insulation, sheath, sheath_insulation, armor, armor_insulation
        if 'armor_insulation_outer_radius' in self.model.scc:
            zcs = zij['zcs']['z11'] + zij['zcs']['z12'] + zij['zcs']['z2i']
            zsa = zij['zsa']['z20'] + zij['zsa']['z23'] + zij['zsa']['z3i']
            za4 = zij['za4']['z30'] + zij['za4']['z34']
            
            # core self-impedance
            Zcc_j = zcs + zsa + za4 - 2 * zij['z2m'] - 2 * zij['z3m']
            
            # sheath self-impedance
            Zss_j = zsa + za4 - 2 * zij['z3m']
            
            # armor self-impedance
            Zaa_j = za4
            
            # mutual impedance between the core and sheath
            Zcs_j = zsa + za4 - zij['z2m'] - 2 * zij['z3m']
            
            # mutual impedance between the core and armor
            Zca_j = za4 - zij['z3m']
            
            # mutual impedance between the sheath and armor
            Zsa_j = Zca_j
            
            # impedance matrix of the jth phase of SCC cable. Eq. (2.8) [2]
            Zij_values = np.array([[Zcc_j, Zcs_j, Zca_j],
                                   [Zcs_j, Zss_j, Zsa_j],
                                   [Zca_j, Zsa_j, Zaa_j]])
            
            # cable internal potential coefficient matrix. Eq. (2.19) [2]
            pcj, psj, paj = zij['potentials']['pcj'], zij['potentials']['psj'], zij['potentials']['paj']
            Pij = np.array([[pcj + psj + paj, psj + paj, paj],
                            [      psj + paj, psj + paj, paj],
                            [            paj,       paj, paj]])
            
        # 2. SCC with core, core_insulation, sheath, sheath_insulation, armor
        elif 'armor_outer_radius' in self.model.scc:
            pass
        
        # 3. SCC with core, core_insulation, sheath, sheath_insulation
        elif 'sheath_insulation_outer_radius' in self.model.scc:
            zcs = zij['zcs']['z11'] + zij['zcs']['z12'] + zij['zcs']['z2i']
            zs3 = zij['zs3']['z20'] + zij['zs3']['z23']
            
            Zcc_j = zcs + zs3 - 2 * zij['z2m']  # core self-impedance
            Zss_j = zs3                         # sheath self-impedance
            Zcs_j = zs3 - zij['z2m']            # mutual impedance between the core and sheath
            
            # impedance matrix of the j-th phase of SCC cable. Eq. (2.11) [2] 
            Zij_values = np.array([[Zcc_j, Zcs_j], [Zcs_j, Zss_j]])
            
            # cable internal potential coefficient matrix. Eq. (2.19) [2]
            pcj, psj = zij['potentials']['pcj'], zij['potentials']['psj']
            Pij = np.array([[pcj + psj, psj], [psj, psj]])
        
        # 4. SCC with core, core_insulation, sheath
        elif 'sheath_outer_radius' in self.model.scc:
            zcs = zij['zcs']['z11'] + zij['zcs']['z12'] + zij['zcs']['z2i']
            z20 = zij['zs3']['z20'] 
            
            Zcc_j = zcs + z20 - 2 * zij['z2m']   # core self-impedance
            Zss_j = z20                          # internal impedance of sheath outer surface
            Zcs_j = z20 - zij['z2m']             # mutual impedance between the core and sheath

            # impedance matrix of the j-th phase of SCC cable. Eq. (2.11) [2]
            Zij_values = np.array([[Zcc_j, Zcs_j],
                                   [Zcs_j, Zss_j]])
            
            # cable internal potential coefficient matrix. Eq. (2.19) [2]
            Pij = np.array([[zij['potentials']['pcj']]])

        # 5. SCC with core, core_insulation
        elif 'core_insulation_outer_radius' in self.model.scc:
            # core self-impedance
            Zcc_j = zij['zcs']['z11'] + zij['zcs']['z12']

            # impedance matrix of the j-th phase of SCC cable. Eq. (2.13) [2]
            Zij_values = np.array([[Zcc_j]])

            # cable internal potential coefficient matrix. Eq. (2.19) [2]
            Pij = np.array([[zij['potentials']['pcj']]])

        # 6. SCC with core
        elif 'core_outer_radius' in self.model.scc:
            # impedance matrix of the j-th phase of SCC cable.
            Zij_values = np.array([[zij['zcs']['z11']]])

        else:
            raise ValueError("Invalid SCC configuration. Please check the conductor layers.")

        # Transpose to get shape (num_freq, M, M)
        Zij = np.transpose(Zij_values, (2, 0, 1))

        # --- Assemble the full internal impedance matrix [Zi] ---
        # A loop is required here because np.kron does not operate on stacks of matrices.
        Zi = np.zeros((len(self.f), num_total_conductors, num_total_conductors), dtype=complex)
        Ri = np.zeros((len(self.f), num_total_conductors, num_total_conductors), dtype=float)
        Li = np.zeros((len(self.f), num_total_conductors, num_total_conductors), dtype=float)
        for i in range(len(self.f)):
            Zi[i, :, :] = np.kron(np.identity(N), Zij[i, :, :])
            Ri[i, :, :] = Zi[i, :, :].real
            Li[i, :, :] = Zi[i, :, :].imag / (2 * np.pi * self.f[i])

        # --- Assemble the full internal potential coefficient matrix [Pi] ---
        # This matrix is frequency-independent, so it's calculated only once.
        Pi = np.kron(np.identity(N), Pij)
        
        # --- Shunt Admittance Matrix, Ye = jw * Pi^-1 ---
        # inv_Pi = lu_solve(lu_factor(Pi), np.identity(Pi.shape[0]))
        inv_Pi = np.linalg.inv(Pi)
        
        # The result is multiplied by the jw vector using broadcasting.
        # jw[:, np.newaxis, np.newaxis] reshapes the 1D jw vector to (num_freq, 1, 1)
        # to correctly multiply with the 2D inv_Pi matrix.
        Ye = self.jw[:, np.newaxis, np.newaxis] * inv_Pi

        return {
            'impedance_matrix': Zi,             # 3D Array: (freq, cond, cond)
            'resistance_matrix': Ri,            # 3D Array: (freq, cond, cond)
            'inductance_matrix': Li,            # 3D Array: (freq, cond, cond)
            'shunt_admittance_matrix': Ye,      # 3D Array: (freq, cond, cond)
            'potential_coefficient_matrix': Pi, # 2D Array (freq-independent)
            'capacitance_matrix': inv_Pi        # 2D Array (freq-independent)
        }

class PerUnitParameters:
    """
    This class calculates vector-frequency PUL parameters using an MTL geometry model,
    including earth-return effects.
    """

    def __init__(self, model: MulticonductorTransmissionLine, frequencies: np.ndarray):
        """
        Initializes the vectorized calculator.

        Args:
            model (MulticonductorTransmissionLine): The MTL geometry model.
            f (np.ndarray): A NumPy array of frequencies to be calculated.
        """
        self.model = model
        self.f = np.asarray(frequencies)

        # Soil Relative Permittivity
        self.e1 = model.mtl_ref[0]['relative_permittivity'] * sc.epsilon_0
        
        # Soil conductivity (S/m)
        self.sigma_1 = model.mtl_ref[0]['conductivity']
        
        # Soil Relative Permeability
        self.mu1 = model.mtl_ref[0]['relative_permeability'] * sc.mu_0
        
        # Soil resistivity (ohm.m)
        self.rho_1 = 1 / self.sigma_1

        # Angular frequency (rad/s) is now a vector
        self.jw = 1j * 2 * np.pi * self.f
        self.jw_mu0_2pi = self.jw * sc.mu_0 / (2 * np.pi)
        self.jw_2pi_e0 = self.jw * 2 * np.pi * sc.epsilon_0
        self.jw_2pi_sg = self.jw / (2 * np.pi * (self.sigma_1 + self.jw * self.e1))

        # Wave numbers are now vectors
        self.k_air2 = -self.jw * sc.mu_0 * self.jw * sc.epsilon_0
        self.k_earth2 = -self.jw * self.mu1 * (self.sigma_1 + self.jw * self.e1)
        self.gamma_earth = np.sqrt(self.jw * self.mu1 * (self.sigma_1 + self.jw * self.e1))

    def earth_return_parameters(self, zg_form='magalhaes_xue', yg_form='magalhaes_xue'):
        """ Calculates Earth-return parameters over a vector of frequencies. """
        N = self.model.num_sc_cables
        d_matrix = self.model.d_matrix_ground_return
        D_matrix = self.model.D_matrix_ground_return
        hnm = self.model.images_vertical_distance_matrix
        dnm = self.model.horizontal_separation_matrix

        # Adjust wave numbers based on the formulation
        k_air2, k_earth2 = self.k_air2, self.k_earth2
        if zg_form in ['sunde', 'deconti_sunde']:
            k_air2 = np.zeros_like(self.f, dtype=complex)
        
        elif zg_form in ['pollaczek', 'ametani', 'saad', 'wedepohl']:
            k_air2 = np.zeros_like(self.f, dtype=complex)
            k_earth2 = -self.jw * self.mu1 * self.sigma_1

        # k_earth2 is (num_freq,), d_matrix is (N, N).
        # Result of sqrt() * d_matrix is (num_freq, N, N)
        # ss.kv operates element-wise, returning a (num_freq, N, N) array.
        K0_jke_dnm = ss.kv(0, 1j * np.sqrt(k_earth2)[:, np.newaxis, np.newaxis] * d_matrix)
        K0_jke_Dnm = ss.kv(0, 1j * np.sqrt(k_earth2)[:, np.newaxis, np.newaxis] * D_matrix)

        S1c = np.zeros((len(self.f), N, N), dtype=complex)
        S2c = np.zeros((len(self.f), N, N), dtype=complex)
        Pg = np.zeros((len(self.f), N, N), dtype=complex)
        Yg = np.zeros_like(Pg, dtype=complex)

        # Wedepohl-Wilcox Approximation
        if zg_form in ['wedepohl']:
            # yg is a 1D vector of shape (num_freq,)
            yg = 1j * np.sqrt(k_earth2)
            
            # Reshape 1D frequency vectors to (num_freq, 1, 1) for broadcasting
            # with 2D geometry matrices (N, N) to get a (num_freq, N, N) result.
            yg_3d = yg[:, np.newaxis, np.newaxis]
            jw_mu0_2pi = self.jw_mu0_2pi[:, np.newaxis, np.newaxis]
            
            ln_term = np.log(0.5 * np.euler_gamma * yg_3d * d_matrix)
            S1c = -ln_term + 0.5 + (2/3) * yg_3d * hnm
            Zg = jw_mu0_2pi * S1c

        # De Conti Closed-Form Approximation
        elif zg_form in ['deconti', 'deconti_sunde', 'saad']:
            # y0 and yg are 1D vectors of shape (num_freq,)
            y0 = 1j * np.sqrt(k_air2)
            yg = 1j * np.sqrt(k_earth2)

            # These operations are on 1D vectors, results are 1D vectors
            zg_term_1 = (yg - y0) / (yg + y0)
            yg_term_1 = (yg**2 - y0**2) / (yg**2 + y0**2)
            
            # Reshape for broadcasting with 2D geometry matrices
            yg_3d = yg[:, np.newaxis, np.newaxis]
            zg_term_1_3d = zg_term_1[:, np.newaxis, np.newaxis]
            yg_term_1_3d = yg_term_1[:, np.newaxis, np.newaxis]
            jw_mu0_2pi = self.jw_mu0_2pi[:, np.newaxis, np.newaxis]
            jw_2pi_sg = self.jw_2pi_sg[:, np.newaxis, np.newaxis]
            
            # Earth-return impedance based on quasi-TEM assumption [1]
            # All arrays are now (num_freq, N, N), allowing for element-wise operations
            term_2 = 2 / (4 + (yg_3d**2 * dnm**2))
            Zg = jw_mu0_2pi * (K0_jke_dnm + zg_term_1_3d * np.exp(hnm * yg_3d) * term_2)

        # Integral expressions based on quasi-TEM assumption [1] 
        elif zg_form in ['magalhaes_xue', 'sunde', 'pollaczek', 'ametani']:
            for n in range(N):
                for m in range(N):
                    if zg_form == 'ametani':
                        S1c[:, n, m] = 2 * sommerfeld_ametani_approx(hnm[n, m], dnm[n, m], ke2=k_earth2)
                    
                    else:
                        S1c[:, n, m] = 2 * sommerfeld_quasi_tem_approx_impedance(hnm[n, m], dnm[n, m], ke2=k_earth2, ka2=k_air2)
                        S2c[:, n, m] = 2 * sommerfeld_quasi_tem_approx_admittance(hnm[n, m], dnm[n, m], ke2=k_earth2, ka2=k_air2)
            
            # Earth-return impedance
            Zg = self.jw_mu0_2pi[:, np.newaxis, np.newaxis] * (K0_jke_dnm - K0_jke_Dnm + S1c)        

        # Earth-return admittance based on quasi-TEM assumption [1]
        if yg_form in ['deconti']:
            Pg = jw_2pi_sg * (K0_jke_dnm + yg_term_1_3d * K0_jke_Dnm)
            for i in range(len(self.f)):
                Yg[i, :, :] = self.jw[i] * np.linalg.inv(Pg[i, :, :])

        elif yg_form in ['magalhaes_xue']:
            Pg = self.jw_2pi_sg[:, np.newaxis, np.newaxis] * (K0_jke_dnm - K0_jke_Dnm + S2c)
            for i in range(len(self.f)):
                Yg[i, :, :] = self.jw[i] * np.linalg.inv(Pg[i, :, :])

        # Earth-Return Admittance based on Vance (1978) formulation
        elif yg_form in ['vance']:
            for i in range(len(self.f)):
                Yg[i, :, :] = (self.gamma_earth[i])**2 * np.linalg.inv(Zg[i, :, :])
                Pg[i, :, :] = self.jw[i] * np.linalg.inv(Yg[i, :, :])

        return {
            'impedance_matrix': Zg,
            'potential_coefficient': Pg,
            'admittance_matrix': Yg,
            'gamma_earth': self.gamma_earth,
        }

    def quasi_tem_approx_matrices(self, internal_matrices, earth_return_params):
        """
        Assembles the final PUL matrices for a vector of frequencies.
        """
        M = self.model.num_conductors_per_scc                   # M = 2; N = 3

        z0_jk = earth_return_params['impedance_matrix']         # Shape (num_freq, N, N)
        pg_jk = earth_return_params['potential_coefficient']    # Shape (num_freq, N, N)
        yg_jk = earth_return_params['admittance_matrix']        # Shape (num_freq, N, N)

        Zi = internal_matrices['impedance_matrix']              # Shape (num_freq, N*M, N*M)
        Yi = internal_matrices['shunt_admittance_matrix']       # Shape (num_freq, N*M, N*M)
        Pi = internal_matrices['potential_coefficient_matrix']  # Shape (N*M, N*M)

        # Loop to build the block matrix for each frequency
        Zg = np.zeros_like(Zi, dtype=complex)
        Pg = np.zeros_like(Zi, dtype=complex)
        Yg = np.zeros_like(Zi, dtype=complex)
        for i in range(len(self.f)):
            Zg[i, :, :] = np.kron(z0_jk[i, :, :], np.ones((M, M)))
            Pg[i, :, :] = np.kron(pg_jk[i, :, :], np.ones((M, M)))

        # Series impedance is a simple element-wise addition
        Zs = Zi + Zg

        # Shunt Admittance Matrix calculation
        # Pi is 2D, Pe is 3D. Use broadcasting to add them.
        Psh = Pi[np.newaxis, :, :] + Pg
        
        # The linear solve must be looped over the frequency axis
        Ysh = np.zeros_like(Psh, dtype=complex)
        for i in range(len(self.f)):
            Ysh[i, :, :] = self.jw[i] * np.linalg.inv(Psh[i, :, :])

        return {
            'earth_return_impedance_matrix': Zg,
            'earth_return_potential_coefficient': Pg,
            'earth_return_admittance_matrix': Yg,
            'potential_coefficient': Psh,
            'series_impedance_matrix': Zs,
            'shunt_admittance_matrix': Ysh,
        }