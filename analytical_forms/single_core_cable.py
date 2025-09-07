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
from scipy.linalg import lu_factor, lu_solve
from mtl_main.source import MulticonductorTransmissionLine

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

class InternalPerUnitParameters:
    """
    This class calculates vector-frequency PUL parameters using an MTL geometry model.
    The calculations are vectorized over the frequency axis for efficiency.
    """

    def __init__(self, model: MulticonductorTransmissionLine, f: np.ndarray):
        """
        Initializes the vectorized calculator.

        Args:
            model (MulticonductorTransmissionLine): The MTL geometry model.
            f (np.ndarray): A NumPy array of frequencies to be calculated.
        """
        self.model = model
        self.num_sc_cables, self.num_conductors_per_scc = model._count_scc_and_conductors()
        
        # Store the frequency array
        self.f = np.asarray(f)
        self.num_freq = len(self.f)

        # Angular frequency (rad/s) is now a vector
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
            'zcs': {'z11': z11, 'z12': z12, 'z2i': z2i},
            'zsa': {'z20': z20, 'z23': z23, 'z3i': z3i},
            'za4': {'z30': z30, 'z34': z34},
            'zs3': {'z20': z20, 'z23': z23},
            'z2m': z2m,
            'z3m': z3m,
            'potentials': {'pcj': pcj, 'psj': psj, 'paj': paj}
        }

    def internal_matrices(self):
        """
        Assembles the full internal impedance [Zi] and shunt admittance [Ye] matrices
        for all specified frequencies.

        Returns:
            dict: A dictionary containing the calculated matrices.
                  'impedance_matrix' (Zi) and 'shunt_admittance_matrix' (Ye) are
                  3D NumPy arrays with shape (num_frequencies, num_total_conductors, num_total_conductors).
        """
        N, M = self.num_sc_cables, self.num_conductors_per_scc
        num_total_conductors = N * M
        zij = self.parameters_by_bessel()
        
        # This part remains the same, but the variables are now 1D NumPy arrays
        # representing the value for each frequency.
        # --- (SCC configuration logic from original code) ---
        if 'armor_insulation_outer_radius' in self.model.scc:
            zcs = zij['zcs']['z11'] + zij['zcs']['z12'] + zij['zcs']['z2i']
            zsa = zij['zsa']['z20'] + zij['zsa']['z23'] + zij['zsa']['z3i']
            za4 = zij['za4']['z30'] + zij['za4']['z34']
            Zcc_j = zcs + zsa + za4 - 2 * zij['z2m'] - 2 * zij['z3m']
            Zss_j = zsa + za4 - 2 * zij['z3m']
            Zaa_j = za4
            Zcs_j = zsa + za4 - zij['z2m'] - 2 * zij['z3m']
            Zca_j = za4 - zij['z3m']
            Zsa_j = Zca_j
            Zij_values = np.array([[Zcc_j, Zcs_j, Zca_j], [Zcs_j, Zss_j, Zsa_j], [Zca_j, Zsa_j, Zaa_j]])
            pcj, psj, paj = zij['potentials']['pcj'], zij['potentials']['psj'], zij['potentials']['paj']
            Pij = np.array([[pcj + psj + paj, psj + paj, paj], [psj + paj, psj + paj, paj], [paj, paj, paj]])
        elif 'sheath_insulation_outer_radius' in self.model.scc:
            zcs = zij['zcs']['z11'] + zij['zcs']['z12'] + zij['zcs']['z2i']
            zs3 = zij['zs3']['z20'] + zij['zs3']['z23']
            Zcc_j = zcs + zs3 - 2 * zij['z2m']
            Zss_j = zs3
            Zcs_j = zs3 - zij['z2m']
            Zij_values = np.array([[Zcc_j, Zcs_j], [Zcs_j, Zss_j]])
            pcj, psj = zij['potentials']['pcj'], zij['potentials']['psj']
            Pij = np.array([[pcj + psj, psj], [psj, psj]])
        # ... (add other elif conditions as in the original code) ...
        else:
            raise ValueError("Invalid SCC configuration. Please check the conductor layers.")

        # Transpose to get shape (num_freq, M, M)
        Zij = np.transpose(Zij_values, (2, 0, 1))

        # --- Assemble the full internal impedance matrix [Zi] ---
        # A loop is required here because np.kron does not operate on stacks of matrices.
        Zi = np.zeros((self.num_freq, num_total_conductors, num_total_conductors), dtype=complex)
        identity_N = np.identity(N)
        for i in range(self.num_freq):
            Zi[i, :, :] = np.kron(identity_N, Zij[i, :, :])

        # --- Assemble the full internal potential coefficient matrix [Pi] ---
        # This matrix is frequency-independent, so it's calculated only once.
        Pi = np.kron(np.identity(N), Pij)
        
        # --- Shunt Admittance Matrix, Ye = jw * Pi^-1 ---
        # The expensive inversion is done only once.
        inv_Pi = lu_solve(lu_factor(Pi), np.identity(Pi.shape[0]))
        
        # The result is multiplied by the jw vector using broadcasting.
        # jw[:, np.newaxis, np.newaxis] reshapes the 1D jw vector to (num_freq, 1, 1)
        # to correctly multiply with the 2D inv_Pi matrix.
        Ye = self.jw[:, np.newaxis, np.newaxis] * inv_Pi

        return {
            'impedance_matrix': Zi, # 3D Array: (freq, cond, cond)
            'shunt_admittance_matrix': Ye, # 3D Array: (freq, cond, cond)
            'potential_coefficient_matrix': Pi, # 2D Array (freq-independent)
        }

class PerUnitParameters:
    """
    This class calculates vector-frequency PUL parameters using an MTL geometry model,
    including earth-return effects.
    """

    def __init__(self, model: MulticonductorTransmissionLine, f: np.ndarray):
        """
        Initializes the vectorized calculator.

        Args:
            model (MulticonductorTransmissionLine): The MTL geometry model.
            f (np.ndarray): A NumPy array of frequencies to be calculated.
        """
        # MTL Geometry Model
        self.model = model
        self.num_sc_cables, self.num_conductors_per_scc = model._count_scc_and_conductors()
        
        self.f = np.asarray(f)
        self.num_freq = len(self.f)

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

    def ground_return_parameters(self, zg_form='magalhaes_xue', yg_form='magalhaes_xue'):
        """
        Calculates Earth-return parameters over a vector of frequencies.
        """
        N = self.num_sc_cables
        d_matrix = self.model.d_matrix_ground_return
        D_matrix = self.model.D_matrix_ground_return
        hnm = self.model.vertical_separation_matrix
        dnm = self.model.horizontal_separation_matrix

        k_air2, k_earth2 = self.k_air2, self.k_earth2

        # Adjust wave numbers based on the formulation
        if zg_form in ['sunde', 'deconti_sunde']:
            k_air2 = np.zeros_like(self.f, dtype=complex)
        elif zg_form in ['pollaczek', 'ametani', 'saad', 'wedepohl']:
            k_air2 = np.zeros_like(self.f, dtype=complex)
            k_earth2 = -self.jw * self.mu1 * self.sigma_1

        # The arguments to ss.kv will broadcast correctly.
        # k_earth2 is (num_freq,), d_matrix is (N, N).
        # Result of sqrt() * d_matrix is (num_freq, N, N)
        # ss.kv operates element-wise, returning a (num_freq, N, N) array.
        arg_d = 1j * np.sqrt(k_earth2)[:, np.newaxis, np.newaxis] * d_matrix
        arg_D = 1j * np.sqrt(k_earth2)[:, np.newaxis, np.newaxis] * D_matrix
        K0_jke_dnm = ss.kv(0, arg_d)
        K0_jke_Dnm = ss.kv(0, arg_D)

        S1c = np.zeros((self.num_freq, N, N), dtype=complex)
        S2c = np.zeros((self.num_freq, N, N), dtype=complex)
        pg = np.zeros((self.num_freq, N, N), dtype=complex)

        # Wedepohl e Wilcox Approximation Expression
        if zg_form in ['wedepohl']:
            # yg is a 1D vector of shape (num_freq,)
            yg = 1j * np.sqrt(k_earth2)
            
            # Reshape 1D frequency vectors to (num_freq, 1, 1) for broadcasting
            # with 2D geometry matrices (N, N) to get a (num_freq, N, N) result.
            yg_3d = yg[:, np.newaxis, np.newaxis]
            jw_mu0_2pi_3d = self.jw_mu0_2pi[:, np.newaxis, np.newaxis]
            
            ln_term = np.log(0.5 * np.euler_gamma * yg_3d * d_matrix)
            S1c = -ln_term + 0.5 + (2/3) * yg_3d * hnm
            zg = jw_mu0_2pi_3d * S1c

        # De Conti Approximation Expressions
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
            jw_mu0_2pi_3d = self.jw_mu0_2pi[:, np.newaxis, np.newaxis]
            jw_2pi_sg_3d = self.jw_2pi_sg[:, np.newaxis, np.newaxis]

            exp_term = np.exp(hnm * yg_3d)
            term_2 = 2 / (4 + (yg_3d**2 * dnm**2))
            
            # Earth-return impedance based on quasi-TEM assumption [1]
            # All arrays are now (num_freq, N, N), allowing for element-wise operations
            zg = jw_mu0_2pi_3d * (K0_jke_dnm + zg_term_1_3d * exp_term * term_2)

            # Earth-return admittance based on quasi-TEM assumption [1]
            if yg_form in ['deconti']:
                pg = jw_2pi_sg_3d * (K0_jke_dnm + yg_term_1_3d * K0_jke_Dnm)
                
                # The linear solve part CANNOT be vectorized and requires a loop
                yg_matrix = np.zeros_like(pg, dtype=complex)
                identity_N = np.identity(N)
                for i in range(self.num_freq):
                    lu, piv = lu_factor(pg[i, :, :])
                    yg_matrix[i, :, :] = self.jw[i] * lu_solve((lu, piv), identity_N)
        
        # Integral expressions are also vectorized thanks to the custom Gauss-Legendre function
        elif zg_form in ['magalhaes_xue', 'sunde', 'pollaczek', 'ametani']:
            for n in range(N):
                for m in range(N):
                    if zg_form == 'ametani':
                        S1c[:, n, m] = 2 * sommerfeld_ametani_approx(hnm[n, m], dnm[n, m], ke2=k_earth2)
                    else:
                        S1c[:, n, m] = 2 * sommerfeld_quasi_tem_approx_impedance(hnm[n, m], dnm[n, m], ke2=k_earth2, ka2=k_air2)
                        if yg_form in ['magalhaes_xue']:
                            S2c[:, n, m] = 2 * sommerfeld_quasi_tem_approx_admittance(hnm[n, m], dnm[n, m], ke2=k_earth2, ka2=k_air2)
            
            # Use broadcasting for element-wise multiplication with the jw vectors
            zg = self.jw_mu0_2pi[:, np.newaxis, np.newaxis] * (K0_jke_dnm - K0_jke_Dnm + S1c)
            pg = self.jw_2pi_sg[:, np.newaxis, np.newaxis] * (K0_jke_dnm - K0_jke_Dnm + S2c)

        return {
            'earth-return_impedance_matrix': zg,
            'earth-return_potential_coefficient': pg,
            'k_earth2': k_earth2,
        }

    def quasi_tem_approximation(self, pul_internal, zg_form='magalhaes_xue', yg_form='magalhaes_xue'):
        """
        Assembles the final PUL matrices for a vector of frequencies.
        """
        N, M = self.num_sc_cables, self.num_conductors_per_scc
        num_total_conductors = N * M

        earth_return = self.ground_return_parameters(zg_form, yg_form)
        z0_jk = earth_return['earth-return_impedance_matrix']  # Shape (num_freq, N, N)
        pg_jk = earth_return['earth-return_potential_coefficient'] # Shape (num_freq, N, N)

        Zi = pul_internal['impedance_matrix'] # Shape (num_freq, N*M, N*M)
        Ye = pul_internal['shunt_admittance_matrix'] # Shape (num_freq, N*M, N*M)
        Pi = pul_internal['potential_coefficient_matrix'] # Shape (N*M, N*M)

        # Loop to build the block matrix for each frequency
        ones_MM = np.ones((M, M))
        Z0 = np.zeros_like(Zi, dtype=complex)
        Pe = np.zeros_like(Zi, dtype=complex)
        for i in range(self.num_freq):
            Z0[i, :, :] = np.kron(z0_jk[i, :, :], ones_MM)
            Pe[i, :, :] = np.kron(pg_jk[i, :, :], ones_MM)

        # Series impedance is a simple element-wise addition
        Zs = Zi + Z0

        # Shunt Admittance Matrix calculation
        # Pi is 2D, Pe is 3D. Use broadcasting to add them.
        P = Pi[np.newaxis, :, :] + Pe
        
        # The linear solve must be looped over the frequency axis
        Ysh = np.zeros_like(P, dtype=complex)
        identity_matrix = np.identity(num_total_conductors)
        for i in range(self.num_freq):
            lu, piv = lu_factor(P[i, :, :])
            Ysh[i, :, :] = self.jw[i] * lu_solve((lu, piv), identity_matrix)

        return {
            'internal_impedance_matrix': Zi,
            'earth-return_impedance_matrix': Z0,
            'series_impedance_matrix': Zs,
            'shunt_admittance_matrix': Ysh,
        }