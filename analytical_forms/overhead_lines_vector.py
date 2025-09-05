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

def sommerfeld(hnm, dnm, ke2, ka2, s_form='s1', pts=150):
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

def t_sommerfeld(hnm, dnm, ke2, ka2, pts=150):
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

def sommerfeld_quasi_tem_approx(hnm, dnm, ke2, n2=1, pts=150):
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
        
        # --- MODIFICAÇÃO PARA EVITAR DIVISÃO POR ZERO ---
        # 1. Crie um array de raios internos seguro, substituindo 0 por 1 nos condutores sólidos.
        ri_safe_for_division = np.where(self.is_tubular_mask, self.ri, 1.0)
        
        # 2. Calcule a razão com o divisor seguro. Nenhuma divisão por zero ocorre aqui.
        ratio = self.ro / ri_safe_for_division
        
        # 3. Use o 'where' principal para escolher entre a razão (para tubos) e 1 (para sólidos).
        log_argument = np.where(self.is_tubular_mask, ratio, 1.0)
        
        # 4. O log é agora calculado sem nenhum aviso.
        log_term = np.log(log_argument)
        # --- FIM DA MODIFICAÇÃO ---

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

    def bessel_impedance(self):
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
            den = ss.kv(1, kri) * ss.iv(1, kro) - ss.iv(1, kri) * ss.kv(1, kro)
            num_ext = ss.iv(0, kro) * ss.kv(1, kri) + ss.iv(1, kri) * ss.kv(0, kro)
            zi_tubular_diag = zi_hf_solid_diag * (num_ext / den)
        zi_tubular_diag[np.abs(kro) < 1e-9] = ri_cc_diag
        
        zi_diag = np.where(self.is_tubular_mask, zi_tubular_diag, zi_solid_diag)
        
        return {'Zi_bessel': self._build_diag_matrix(zi_diag)}
        
    def kelvin_impedance(self):
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
        bessel_params = self.bessel_impedance()
        kelvin_params = self.kelvin_impedance()        
        approx_params = self.approximations()
        return {**dc_params, **bessel_params, **kelvin_params, **approx_params}

# =============================================================================
# CLASSE PARA PARÂMETROS EXTERNOS E MONTAGEM FINAL
# =============================================================================
class PerUnitParameters:
    """
    Calculates external and final PUL parameters for overhead lines.
    """
    def __init__(self, model: MulticonductorTransmissionLine, f: np.ndarray):
        self.model = model
        self.f = np.asarray(f)
        self.num_freq = len(self.f)
        self.num_conductors = len(model.surfaces)
        self.jw = 1j * 2 * np.pi * self.f
        
        # Parâmetros do solo
        er_1 = model.mtl_ref[0]['relative_permittivity']
        sigma_1 = model.mtl_ref[0]['conductivity']
        mur_1 = model.mtl_ref[0]['relative_permeability']
        
        self.jw_mu0_2pi = self.jw * sc.mu_0 / (2 * np.pi)
        self.k_air2 = -self.jw * sc.mu_0 * self.jw * sc.epsilon_0
        self.k_earth2 = -self.jw * mur_1 * sc.mu_0 * (sigma_1 + self.jw * er_1 * sc.epsilon_0)

    def ground_return_parameters(self, form='quasi_tem'):
        """ Calculates external and ground-return impedance matrices. """
        N = self.num_conductors
        
        # Termo de impedância externa (independente da frequência)
        M = np.log(self.model.D_matrix_ground_return) - np.log(self.model.d_matrix_ground_return)
        # Transforma em 3D usando broadcasting com o vetor de frequência
        ze = self.jw_mu0_2pi[:, np.newaxis, np.newaxis] * M
        
        # Termo de retorno pelo solo (dependente da formulação)
        S1 = np.zeros((self.num_freq, N, N), dtype=complex)
        # S2 e T podem ser necessários para outras formulações
        
        # O laço (n, m) é mantido, mas as chamadas internas são vetorizadas
        for n in range(N):
            for m in range(N):
                hnm = self.model.vertical_separation_matrix[n, m]
                dnm = self.model.horizontal_separation_matrix[n, m]
                
                if form == 'quasi_tem':
                     S1[:, n, m] = 2 * sommerfeld(hnm, dnm, self.k_earth2, self.k_air2)
                # Adicionar outras formulações (elif) aqui se necessário
        
        zg = self.jw_mu0_2pi[:, np.newaxis, np.newaxis] * S1

        return {'external_impedance_matrix': ze, 'earth_return_impedance_matrix': zg}

    def calculate_pul_parameters(self, pul_internal, form='quasi_tem'):
        """ Assembles the final PUL matrices. """
        N = self.num_conductors
        
        # Pega a impedância interna pré-calculada
        zi = pul_internal['internal_impedance_matrix']
        
        # Calcula as impedâncias externas e de retorno pelo solo
        ground_params = self.ground_return_parameters(form=form)
        ze = ground_params['external_impedance_matrix']
        zg = ground_params['earth_return_impedance_matrix']
        
        # Matriz de Impedância Série Final
        zs = zi + ze + zg

        # Matriz de Admitância Shunt (baseada nos Coeficientes de Potencial de Maxwell)
        # Nota: Esta formulação não inclui o termo de correção do solo S2
        M = np.log(self.model.D_matrix_ground_return) - np.log(self.model.d_matrix_ground_return)
        Pe = (1 / (2 * np.pi * sc.epsilon_0)) * M
        
        # A inversão deve ser feita para cada frequência
        ysh = np.zeros((self.num_freq, N, N), dtype=complex)
        jw_2pi_e0 = self.jw * 2 * np.pi * sc.epsilon_0
        identity_N = np.identity(N)
        
        # Inversão da matriz de potencial (que é constante)
        inv_Pe = linalg.inv(Pe)
        
        for i in range(self.num_freq):
            ysh[i] = jw_2pi_e0[i] * inv_Pe
            
        return {'series_impedance_matrix': zs, 'shunt_admittance_matrix': ysh, 'zi': zi, 'ze': ze, 'zg': zg}