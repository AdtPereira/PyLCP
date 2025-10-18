# analytical_forms/overhead_lines_vector.py

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
from scipy import linalg
from mtl_main.source import MulticonductorTransmissionLine
from numpy.lib import scimath

def sommerfeld(hnm, dnm, ke2, ka2, s_form='s1', pts=300):
    """ Vectorized calculation of the Sommerfeld integral (Form S1/S2). """
    if s_form == 's1':
        n = 1
    elif s_form == 's2':
        n = np.sqrt(ke2 / ka2)

    def _int_transformed(t, hnm, dnm, ke2, ka2, n):
        x = np.tan(t[:, np.newaxis])
        num = np.exp(-hnm * x) * np.cos(dnm * x)
        den = np.sqrt(x**2 + ka2 - ke2) + (n**2 * x)
        return (num / den) * (1 / np.cos(t[:, np.newaxis])**2)

    def _int_gauss_legendre(func, a, b, n_pts, *args):
        [x_gl, w] = np.polynomial.legendre.leggauss(n_pts)
        t = 0.5 * (x_gl + 1) * (b - a) + a
        result_matrix = func(t, *args)
        # Sum along the integration axis (axis=0)
        return np.sum(w[:, np.newaxis] * result_matrix, axis=0) * 0.5 * (b - a)

    return _int_gauss_legendre(_int_transformed, 0, np.pi/2, pts, hnm, dnm, ke2, ka2, n)

def t_sommerfeld(hnm, dnm, ke2, ka2, pts=300):
    """ Vectorized calculation of the T-type Sommerfeld integral. """
    
    def _int_transformed(t, hnm, dnm, ke2, ka2):
        x = np.tan(t[:, np.newaxis])
        # Refractive index 'n' is now a vector
        n = np.sqrt(ke2 / ka2)
        
        # All calculations are broadcasted
        exp_1 = np.exp(-hnm * x)
        exp_2 = np.exp(-0.5 * hnm * x)
        u2 = np.sqrt(x**2 + ka2 - ke2)
        num = u2 * (exp_2 - exp_1) * np.cos(dnm * x)
        den = (n * x)**2 + x * u2

        return (num / den) * (1 / np.cos(t[:, np.newaxis])**2)

    def _int_gauss_legendre(func, a, b, n_pts, *args):
        [x_gl, w] = np.polynomial.legendre.leggauss(n_pts)
        t = 0.5 * (x_gl + 1) * (b - a) + a
        result_matrix = func(t, *args)
        # Sum along the integration axis (axis=0)
        return np.sum(w[:, np.newaxis] * result_matrix, axis=0) * 0.5 * (b - a)

    return _int_gauss_legendre(_int_transformed, 0, np.pi/2, pts, hnm, dnm, ke2, ka2)

def sommerfeld_quasi_tem_approx(hnm, dnm, ke2, n2=1, pts=300):
    """ Vectorized calculation of the quasi-TEM Sommerfeld approximation. """

    def _int_transformed(t, hnm, dnm, ke2, n2):
        x = np.tan(t[:, np.newaxis])
        num = np.exp(-hnm * x) * np.cos(dnm * x)
        den = np.sqrt(x**2 - ke2) + n2 * x
        return (num / den) * (1 / np.cos(t[:, np.newaxis])**2)
        
    def _int_gauss_legendre(func, a, b, n_pts, *args):
        [x_gl, w] = np.polynomial.legendre.leggauss(n_pts)
        t = 0.5 * (x_gl + 1) * (b - a) + a
        result_matrix = func(t, *args)
        # Sum along the integration axis (axis=0)
        return np.sum(w[:, np.newaxis] * result_matrix, axis=0) * 0.5 * (b - a)

    return _int_gauss_legendre(_int_transformed, 0, np.pi/2, pts, hnm, dnm, ke2, n2)

def sommerfeld_adaptive(hnm, dnm, ke2, ka2, s_form='s1', limit=100, epsrel=1e-6,):
    """
    Adaptive calculation of the Sommerfeld integral (Form S1/S2)
    using scipy.integrate.quad for improved numerical stability.
    """
    if s_form == 's1':
        n = 1.0
    elif s_form == 's2':
        n = scimath.sqrt(ke2 / ka2)

    def _integrand(x, hnm, dnm, ke2, ka2, n):
        """ The integrand function in terms of the variable x. """
        numerator = np.exp(-hnm * x) * np.cos(dnm * x)
        denominator = scimath.sqrt(x**2 + ka2 - ke2) + (n**2 * x)
        if denominator == 0:
            return 0.0
        return numerator / denominator

    # Integrate with complex_func=True to ensure proper handling of complex values
    result, _ = quad(_integrand, 0, np.inf, args=(hnm, dnm, ke2, ka2, n),
                     limit=limit, epsrel=epsrel, complex_func=True)
    return result

def t_sommerfeld_adaptive(hnm, dnm, ke2, ka2, limit=100, epsrel=1e-6):
    """
    Adaptive calculation of the T-type Sommerfeld integral
    using scipy.integrate.quad for improved numerical stability.
    """
    n = scimath.sqrt(ke2 / ka2)

    def _integrand(x, hnm, dnm, ke2, ka2, n):
        """ The integrand function in terms of the variable x. """
        u2 = scimath.sqrt(x**2 + ka2 - ke2)
        exp_term = np.exp(-0.5 * hnm * x) - np.exp(-hnm * x)
        numerator = u2 * exp_term * np.cos(dnm * x)
        denominator = (n * x)**2 + x * u2
        if denominator == 0:
            return 0.0
        return numerator / denominator

    # Integrate with complex_func=True to ensure proper handling of complex values
    result, _ = quad(_integrand, 0, np.inf, args=(hnm, dnm, ke2, ka2, n),
                     limit=limit, epsrel=epsrel, complex_func=True)
    return result

def sommerfeld_quasi_tem_approx_adaptive(hnm, dnm, ke2, n2=1.0, limit=100, epsrel=1e-6):
    """
    Adaptive calculation of the quasi-TEM Sommerfeld approximation
    using scipy.integrate.quad for improved numerical stability.
    """
    def _integrand(x, hnm, dnm, ke2, n2):
        """ The integrand function in terms of the variable x. """
        numerator = np.exp(-hnm * x) * np.cos(dnm * x)
        denominator = scimath.sqrt(x**2 - ke2) + n2 * x
        if denominator == 0:
            return 0.0
        return numerator / denominator
        
    # Integrate with complex_func=True to ensure proper handling of complex values
    result, _ = quad(_integrand, 0, np.inf, args=(hnm, dnm, ke2, n2),
                     limit=limit, epsrel=epsrel, complex_func=True)
    return result

# =============================================================================
# CLASSE PARA PARÂMETROS INTERNOS
# =============================================================================
class InternalPerUnitParameters:
    """
    Calculates the vectorized internal PUL parameters for overhead lines.
    This version is modular and calculates multiple impedance formulations.
    """
    def __init__(self, model: MulticonductorTransmissionLine, f: np.ndarray):
        self.model = model
        self.f = np.asarray(f)
        self.num_freq = len(self.f)
        self.num_conductors = len(model.surfaces)
        self.jw = 1j * 2 * np.pi * self.f

        radii = [conductor['radius'] for conductor in model.mtl.values()]
        self.ri = np.array([r[0] if isinstance(r, (list, np.ndarray)) else 0 for r in radii])
        self.ro = np.array([r[1] if isinstance(r, (list, np.ndarray)) else r for r in radii])
        self.is_tubular_mask = self.ri > 0

    def _build_diag_matrix(self, diag_values):
        """Helper to build a 3D diagonal matrix from diagonal values."""
        N = self.num_conductors
        if diag_values.ndim == 2: # Frequency dependent
            matrix = np.zeros((self.num_freq, N, N), dtype=complex)
            matrix[:, np.arange(N), np.arange(N)] = diag_values
        else: # Frequency independent
            matrix = np.diag(diag_values)
        return matrix

    def dc_parameters(self):
        """Calculates frequency-independent DC resistance and inductance."""
        # DC Resistance
        ri_cc_solid = 1.0 / (self.model.sigma * np.pi * self.ro**2)
        
        # Use np.where para evitar divisão por zero se ro == ri em um tubo defeituoso
        ri_cc_tubular_denom = self.model.sigma * np.pi * (self.ro**2 - self.ri**2)
        ri_cc_tubular = np.divide(1.0, ri_cc_tubular_denom, where=ri_cc_tubular_denom!=0)        
        ri_cc_diag = np.where(self.is_tubular_mask, ri_cc_tubular, ri_cc_solid)

        # DC Internal Inductance
        li_cc_solid = self.model.mu / (8 * np.pi)
        
        # 1. Crie um array de raios internos seguro, substituindo 0 por 1 nos condutores sólidos.
        ri_safe_for_division = np.where(self.is_tubular_mask, self.ri, 1.0)
        
        # 2. Calcule a razão com o divisor seguro. Nenhuma divisão por zero ocorre aqui.
        ratio = self.ro / ri_safe_for_division
        
        # 3. Use o 'where' principal para escolher entre a razão (para tubos) e 1 (para sólidos).
        log_argument = np.where(self.is_tubular_mask, ratio, 1.0)
        
        # 4. O log é agora calculado sem nenhum aviso.
        log_term = np.log(log_argument)

        ro2, ri2 = self.ro**2, self.ri**2
        ro2ri2 = ro2 - ri2
        
        # Evita divisão por zero se ro2ri2 for zero (tubo defeituoso)
        term1_denom = 4 * ro2ri2
        term1 = np.divide((ro2 - 3 * ri2), term1_denom, where=term1_denom!=0)

        term2_denom = ro2ri2**2
        term2 = np.divide((ri2**2), term2_denom, where=term2_denom!=0) * log_term

        li_cc_tubular = (self.model.mu / (2 * np.pi)) * (term1 + term2)        
        li_cc_diag = np.where(self.is_tubular_mask, li_cc_tubular, li_cc_solid)

        return {
            'Ri_cc': self._build_diag_matrix(ri_cc_diag),
            'Li_cc': self._build_diag_matrix(li_cc_diag)
        }

    def matrix_impedance_bessel(self):
        """Calculates exact internal impedance using Bessel functions."""
        # (Este é o método 'impedance_matrix' da versão anterior, renomeado e limpo)
        dc_params = self.dc_parameters()
        ri_cc_diag = np.diag(dc_params['Ri_cc'])        
        jw_mu = self.jw[:, np.newaxis] * self.model.mu
        
        # Solid Wire
        with np.errstate(divide='ignore', invalid='ignore'):
            zi_hf_solid_diag = (1.0 / (2 * np.pi * self.ro)) * np.sqrt(jw_mu / self.model.sigma)
            bessel_arg_solid = np.sqrt(jw_mu * self.model.sigma) * self.ro
            zi_solid_diag = zi_hf_solid_diag * ss.iv(0, bessel_arg_solid) / ss.iv(1, bessel_arg_solid)
        
        zi_solid_diag[np.abs(bessel_arg_solid) < 1e-9] = ri_cc_diag

        # Tubular Wire
        with np.errstate(divide='ignore', invalid='ignore'):
            kappa = np.sqrt(jw_mu * self.model.sigma)
            kro, kri = kappa * self.ro, kappa * self.ri
            num = ss.iv(0, kro) * ss.kv(1, kri) + ss.iv(1, kri) * ss.kv(0, kro)
            den = ss.kv(1, kri) * ss.iv(1, kro) - ss.iv(1, kri) * ss.kv(1, kro)
            zi_tubular_diag = zi_hf_solid_diag * (num / den)
        
        zi_tubular_diag[np.abs(kro) < 1e-9] = ri_cc_diag        
        zi_diag = np.where(self.is_tubular_mask, zi_tubular_diag, zi_solid_diag)
        
        return {'Zi_bessel': self._build_diag_matrix(zi_diag)}
        
    def kelvin_impedance_solid_wires(self):
        """Calculates exact internal impedance for solid wires using Kelvin functions."""
        dc_params = self.dc_parameters()
        ri_cc_diag = np.diag(dc_params['Ri_cc'])

        # Skin depth: shape (F, N)
        skin_depth = 1 / np.sqrt(self.model.mu * np.pi * self.f[:, np.newaxis] * self.model.sigma)
        q = np.sqrt(2) * self.ro / skin_depth
        
        # Kelvin functions work on the (F, N) matrix q
        Be, Ke, Bep, Kep = ss.kelvin(q)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            scaling_factor = ri_cc_diag * (q / 2) / (Bep.imag**2 + Bep.real**2)
            # Resistance
            numerador_R = Be.real * Bep.imag - Be.imag * Bep.real
            Ri_val = scaling_factor * numerador_R
            
            # Reactance
            numerador_wL = Be.real * Bep.real + Be.imag * Bep.imag
            wLi_val = scaling_factor * numerador_wL

        zi_kelvin_diag = Ri_val + 1j * wLi_val
        
        # Kelvin functions are for solid wires, result for tubular is meaningless/inf
        zi_kelvin_diag = np.where(self.is_tubular_mask, np.inf, zi_kelvin_diag)
        zi_kelvin_diag[q < 1e-6] = ri_cc_diag

        return {'Zi_kelvin': self._build_diag_matrix(zi_kelvin_diag)}

    def approximations(self):
        """Calculates various common impedance approximations."""
        dc_params = self.dc_parameters()
        ri_cc_diag = np.diag(dc_params['Ri_cc'])        
        jw_mu = self.jw[:, np.newaxis] * self.model.mu

        # High-frequency impedance is the limit of Bessel
        zi_hf_diag = 1.0 / (2 * np.pi * self.ro) * np.sqrt(jw_mu / self.model.sigma)

        Ri_cc = self._build_diag_matrix(ri_cc_diag)
        Zi_hf = self._build_diag_matrix(zi_hf_diag)
        
        # Nahman's approximation
        Zi_nahman = Ri_cc + Zi_hf
        
        # Simple approximation
        Zi_approx = np.sqrt(Ri_cc**2 + Zi_hf**2)

        return {
            'Zi_hf': Zi_hf,
            'Zi_nahman': Zi_nahman,
            'Zi_approx': Zi_approx
        }

    def all_terms(self):
        """
        Calculates all internal impedance formulations and returns them in a single dictionary.
        This is the main method to be called from outside.
        """
        dc_params = self.dc_parameters()
        bessel_params = self.matrix_impedance_bessel()
        kelvin_params = self.kelvin_impedance_solid_wires()        
        approx_params = self.approximations()
        return {**dc_params, **bessel_params, **kelvin_params, **approx_params}

# =============================================================================
# CLASSE PARA PARÂMETROS EXTERNOS E MONTAGEM FINAL
# =============================================================================
class PerUnitParameters:
    """ Calculates external and final PUL parameters for overhead lines. """
    def __init__(self, model: MulticonductorTransmissionLine, f: np.ndarray):
        self.model = model
        self.f = np.asarray(f)
        self.num_freq = len(self.f)
        self.num_conductors = len(model.surfaces)
        self.jw = 1j * 2 * np.pi * self.f
        
        # Ground-return related constants
        er_1 = model.mtl_ref[0]['relative_permittivity']
        sigma_1 = model.mtl_ref[0]['conductivity']
        mur_1 = model.mtl_ref[0]['relative_permeability']
        
        self.jw_mu0_2pi = self.jw * sc.mu_0 / (2 * np.pi)
        self.jw_2pi_e0 = self.jw * 2 * np.pi * sc.epsilon_0
        self.k_air2 = -self.jw * sc.mu_0 * self.jw * sc.epsilon_0
        self.k_earth2 = -self.jw * mur_1 * sc.mu_0 * (sigma_1 + self.jw * er_1 * sc.epsilon_0)

    def _calc_quasi_tem(self, N, integral_form, **kwargs):
        """ Handles the full quasi-TEM formulation (S1, S2, T). """
        S1 = np.zeros((self.num_freq, N, N), dtype=complex)
        S2 = np.zeros((self.num_freq, N, N), dtype=complex)
        T = np.zeros((self.num_freq, N, N), dtype=complex)

        if integral_form == 'quad':
            # Loop over each frequency since adaptive quad is scalar
            for i in range(self.num_freq):
                ke2_i, ka2_i = self.k_earth2[i], self.k_air2[i]
                for n in range(N):
                    for m in range(N):
                        hnm = self.model.images_vertical_distance_matrix[n, m]
                        dnm = self.model.horizontal_separation_matrix[n, m]
                        S1[i, n, m] = 2 * sommerfeld_adaptive(hnm, dnm, ke2_i, ka2_i)
                        S2[i, n, m] = 2 * sommerfeld_adaptive(hnm, dnm, ke2_i, ka2_i, s_form='s2')
                        T[i, n, m] = 2 * t_sommerfeld_adaptive(hnm, dnm, ke2_i, ka2_i)
        else:
            # Vectorized Gauss-Legendre integration
            for n in range(N):
                for m in range(N):
                    hnm = self.model.images_vertical_distance_matrix[n, m]
                    dnm = self.model.horizontal_separation_matrix[n, m]
                    S1[:, n, m] = 2 * sommerfeld(hnm, dnm, self.k_earth2, self.k_air2)
                    S2[:, n, m] = 2 * sommerfeld(hnm, dnm, self.k_earth2, self.k_air2, s_form='s2')
                    T[:, n, m] = 2 * t_sommerfeld(hnm, dnm, self.k_earth2, self.k_air2)
        return S1, S2, T

    def _calc_quasi_tem_log(self, N, **kwargs):
        """ Handles the logarithmic approximation of the quasi-TEM formulation. """
        eta = np.sqrt(self.k_air2 - self.k_earth2)
        eta_3d = eta[:, np.newaxis, np.newaxis]
        
        eta_sqrt_matrix = eta_3d * self.model.D_matrix_ground_return
        n2 = self.k_earth2 / self.k_air2
        n2_3d = n2[:, np.newaxis, np.newaxis]

        log_s1 = 1 + 2 / eta_sqrt_matrix
        log_s2 = 1 + (n2_3d + 1) / eta_sqrt_matrix
        log_t = 1 + 2 * (n2_3d + 1) / eta_sqrt_matrix

        S1 = np.log(log_s1)
        S2 = 2 / (n2_3d + 1) * np.log(log_s2)
        T = 2 * np.log(2) + (2 * n2_3d / (n2_3d + 1)) * np.log(log_s2 / log_t)
        return S1, S2, T

    def _calc_simplified_approx(self, N, zg_form, integral_form):
        """ Handles Carson, Sunde, and Nakagawa approximations. """
        S1 = np.zeros((self.num_freq, N, N), dtype=complex)
        S2 = np.zeros((self.num_freq, N, N), dtype=complex)
        T = np.zeros((self.num_freq, N, N), dtype=complex) # Not used, but returned for consistency

        if zg_form == 'nakagawa':
            er_1 = self.model.mtl_ref[0]['relative_permittivity']
            sigma_1 = self.model.mtl_ref[0]['conductivity']
            mur_1 = self.model.mtl_ref[0]['relative_permeability']
            ke2_form = -self.jw * mur_1*sc.mu_0 * (sigma_1 + self.jw * (er_1 - 1)*sc.epsilon_0)
            n2_form = self.k_earth2 / self.k_air2
        elif zg_form == 'sunde':
            ke2_form = self.k_earth2
            n2_form = 1
        elif zg_form == 'carson':
            sigma_1 = self.model.mtl_ref[0]['conductivity']
            mur_1 = self.model.mtl_ref[0]['relative_permeability']
            ke2_form = -self.jw * mur_1 * sc.mu_0 * sigma_1
            n2_form = 1

        if integral_form == 'quad':
            for i in range(self.num_freq):
                ke2_i = ke2_form[i]
                n2_i = n2_form[i] if isinstance(n2_form, np.ndarray) else n2_form
                for n in range(N):
                    for m in range(N):
                        hnm = self.model.images_vertical_distance_matrix[n, m]
                        dnm = self.model.horizontal_separation_matrix[n, m]
                        S1[i, n, m] = 2 * sommerfeld_quasi_tem_approx_adaptive(hnm, dnm, ke2_i)
                        if zg_form == 'nakagawa':
                            S2[i, n, m] = 2 * sommerfeld_quasi_tem_approx_adaptive(hnm, dnm, ke2_i, n2=n2_i)
        else: 
            for n in range(N):
                for m in range(N):
                    hnm = self.model.images_vertical_distance_matrix[n, m]
                    dnm = self.model.horizontal_separation_matrix[n, m]
                    S1[:, n, m] = 2 * sommerfeld_quasi_tem_approx(hnm, dnm, ke2_form)                    
                    if zg_form == 'nakagawa':
                        S2[:, n, m] = 2 * sommerfeld_quasi_tem_approx(hnm, dnm, ke2_form, n2=n2_form)
        return S1, S2, T

    def _calc_log_approx(self, N, zg_form, **kwargs):
        """ Handles Deri and Sunde (log) approximations. """
        S2 = np.zeros((self.num_freq, N, N), dtype=complex)
        T = np.zeros((self.num_freq, N, N), dtype=complex)

        if zg_form == 'deri':
            sigma_1 = self.model.mtl_ref[0]['conductivity']
            mur_1 = self.model.mtl_ref[0]['relative_permeability']
            k_e2_form = -self.jw * mur_1 * sc.mu_0 * sigma_1
        else: # sunde_log
            k_e2_form = self.k_earth2
        
        p_dot = 1 / np.sqrt(-k_e2_form)
        
        h_array = np.array([c['center_point'][1] for c in self.model.surfaces])
        S1_diag = np.log((h_array + p_dot[:, np.newaxis]) / h_array)
        
        num_sq = (self.model.images_vertical_distance_matrix + 2 * p_dot[:, np.newaxis, np.newaxis])**2 + self.model.horizontal_separation_matrix**2
        num = np.sqrt(num_sq)
        den = self.model.d_matrix_ground_return
        
        S1 = np.log(num / den)
        S1[:, np.arange(N), np.arange(N)] = S1_diag
        return S1, S2, T

    def ground_return_parameters(self, zg_form='quasi_tem', integral_form='quad'):
        """
        Calculates external and ground-return intermediate matrices (S1, S2, T)
        by dispatching to the appropriate private calculation method.
        """
        N = self.num_conductors
        
        # External impedance term (common to all formulations)
        M = np.log(self.model.D_matrix_ground_return) - np.log(self.model.d_matrix_ground_return)
        ze = self.jw_mu0_2pi[:, np.newaxis, np.newaxis] * M
        
        # Map formulation strings to their corresponding calculation methods
        formulation_map = {
            'quasi_tem': self._calc_quasi_tem,
            'quasi_tem_log': self._calc_quasi_tem_log,
            'nakagawa': self._calc_simplified_approx,
            'sunde': self._calc_simplified_approx,
            'carson': self._calc_simplified_approx,
            'sunde_log': self._calc_log_approx,
            'deri': self._calc_log_approx
        }

        # Get the appropriate calculation function from the map
        calculation_func = formulation_map.get(zg_form)
        if not calculation_func:
            raise ValueError(f"Unknown zg_form: '{zg_form}'")
            
        # Call the selected method to get S1, S2, and T matrices
        S1, S2, T = calculation_func(N=N, zg_form=zg_form, integral_form=integral_form)
        
        return {
            'external_impedance_matrix': ze,
            'M_matrix': M,
            'S1': S1,
            'S2': S2,
            'T': T
        }

    def pul_matrices(self, zi, zg_form='quasi_tem'):
        """ Assembles the final PUL matrices and propagation parameters. """
        N = self.num_conductors
        
        # Get the intermediate matrices from the chosen ground formulation
        ground_params = self.ground_return_parameters(zg_form=zg_form)
        ze = ground_params['external_impedance_matrix']
        M, S1, S2, T = ground_params['M_matrix'], ground_params['S1'], ground_params['S2'], ground_params['T']
        
        # --- Assemble Zs and Ysh according to the formulation ---
        if zg_form in ['quasi_tem', 'quasi_tem_log']:
            zg = self.jw_mu0_2pi[:, np.newaxis, np.newaxis] * (S1 - (T + S2))
            M_minus_T = M[np.newaxis, :, :] - T
            # Inversion requires frequency loop
            inv_M_minus_T = np.zeros_like(M_minus_T, dtype=complex)
            for i in range(self.num_freq):
                inv_M_minus_T[i] = linalg.inv(M_minus_T[i])
            ysh = self.jw_2pi_e0[:, np.newaxis, np.newaxis] * inv_M_minus_T
        
        else: # Default for Nakagawa, Carson, Sunde, Deri, etc.
            zg = self.jw_mu0_2pi[:, np.newaxis, np.newaxis] * S1
            M_plus_S2 = M[np.newaxis, :, :] + S2
            # Inversion requires frequency loop
            inv_M_plus_S2 = np.zeros_like(M_plus_S2, dtype=complex)
            for i in range(self.num_freq):
                inv_M_plus_S2[i] = linalg.inv(M_plus_S2[i])
            ysh = self.jw_2pi_e0[:, np.newaxis, np.newaxis] * inv_M_plus_S2

        zs = zi + ze + zg

        # --- Calculate propagation parameters (requires loops) ---
        gamma_i = np.zeros((self.num_freq, N, N), dtype=complex)
        gamma_v = np.zeros((self.num_freq, N, N), dtype=complex)
        zc = np.zeros((self.num_freq, N, N), dtype=complex)
        yc = np.zeros((self.num_freq, N, N), dtype=complex)
        
        for i in range(self.num_freq):
            zs_i = zs[i]
            ysh_i = ysh[i]

            # Inverse of Y using LU decomposition
            # Solve the system Y * Y_inv = I to find Y_inv
            lu, piv = linalg.lu_factor(ysh_i)
            I = np.identity(ysh_i.shape[0])
            ysh_inv_i = linalg.lu_solve((lu, piv), I)

            # Inverse of Z using LU decomposition
            # Solve the system Z * Z_inv = I to find Z_inv
            lu, piv = linalg.lu_factor(zs_i)
            I = np.identity(zs_i.shape[0])
            zs_inv_i = linalg.lu_solve((lu, piv), I)

            # Currents and Voltages Propagation Constants Matrix
            gamma_i[i] = linalg.sqrtm(ysh_i @ zs_i)
            gamma_v[i] = linalg.sqrtm(zs_i @ ysh_i)
            
            # Characteristic impedance Zc = Y_inv * sqrt(Y*Z)
            zc[i] = ysh_inv_i @ gamma_i[i]

            # Characteristic admittance Yc = Z_inv * sqrt(Z*Y)
            yc[i] = zs_inv_i @ gamma_v[i]

        return {
            'series_impedance_matrix': zs,
            'shunt_admittance_matrix': ysh, 
            'internal_impedance_matrix': zi,
            'free-space_impedance_matrix': ze,
            'earth-return_impedance_matrix': zg,
            'propagation_voltage_matrix': gamma_v,
            'propagation_current_matrix': gamma_i,
            'characteristic_impedance_matrix': zc,
            'characteristic_admittance_matrix': yc
        }
    
    