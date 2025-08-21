import numpy as np
from scipy.special import iv
from scipy.constants import mu_0, epsilon_0
from scipy.integrate import quad
from mtl_main.source import MulticonductorTransmissionLine


class PerUnitParameters:
    """ This class calculates PUL parameters using an MTL geometry model. """

    def __init__(self, model: MulticonductorTransmissionLine, f, sigma_1, er_1=1, mur_1=1, ge=0):
        # MTL Geometry Model
        self.model = model

        # Soil Permittivity [np.array]
        self.er_1 = er_1

        # Soil conductivity [np.array]
        self.sigma_1 = sigma_1

        # Soil Permeability [np.array]
        self.mur_1 = mur_1

        # External Conductance [np.array]
        self.ge = ge

        # Angular frequency (rad/s)
        self.jw = 1j * 2 * np.pi * f

        # Constant vacuum terms
        self.jw_mu0_2pi = self.jw * mu_0 / 2 / np.pi
        self.jw_2pi_e0 = self.jw * 2 * np.pi * epsilon_0

        # Air wave number (rad/m)
        self.ka2 = - self.jw * mu_0 * self.jw * epsilon_0

        # Earth wave number (rad/m)
        self.ke2 = - self.jw * self.mur_1 * mu_0 * (self.sigma_1 + self.jw * self.er_1 * epsilon_0)

    def internal_impedance(self, type_form='approx'):
        """ This method calculates the internal impedance of solid wires. """
        N = len(self.model.surfaces)
        Zi = np.zeros((N, N), dtype=complex)
        Ri_cc = np.zeros_like(Zi)
        Zi_hf = np.zeros_like(Zi)

        for p, conductor in enumerate(self.model.surfaces):
            # Angular frequency (rad/s)
            jw_mu = self.jw * self.model.mu[p]

            # Outer wire radius (m)
            ro = conductor['radius']

            # Wire conductivity (S/m)
            sigma = self.model.sigma[p]

            # High frequency impedance (ohm/m)
            # Derived assuming current conduction in a ring with
            # thickness equal to the penetration depth (ohm/m).
            zi_hf = 1.0/(2 * np.pi * ro) * np.sqrt(jw_mu/sigma)

            # Approximation Closed-Form Expression
            if type_form == 'approx':
                # Continuous current resistance (ohm/m)
                ri_cc = 1.0/(sigma * np.pi * ro**2)

                Ri_cc[p, p] = ri_cc
                Zi_hf[p, p] = zi_hf
                Zi[p, p] = np.sqrt(ri_cc ** 2 + zi_hf ** 2)

            else:
                arg = np.sqrt(jw_mu * sigma) * ro
                io = iv(0, arg)
                i1 = iv(1, arg)
                Zi[p, p] = zi_hf * (io/i1)

        return Zi, Ri_cc, Zi_hf

    def external_impedance(self):
        """ This method calculates the impedance matrix of the earth return path. """
        N = len(self.model.surfaces)
        Ze = np.zeros((N, N), dtype=complex)

        # Busca as matrizes pré-calculadas do modelo
        d_matrix = self.model.d_matrix_ground_return
        D_matrix = self.model.D_matrix_ground_return

        # for n, conductor_n in enumerate(self.mtl.surfaces):
        #     for m, conductor_m in enumerate(self.mtl.surfaces):
        #         # Distance between the conductors (m)
        #         dnm, Dnm, _, _ = self.conductor_distances(conductor_n['tag'], conductor_m['tag'])

        #         # External impedance (ohm/m)
        #         Ze[n, m] = self.jw_mu0_2pi * np.log(Dnm/dnm)

        # O cálculo agora é uma operação de matriz, muito mais eficiente
        # Nota: np.log(0) na diagonal de log(d_matrix) pode dar -inf. Tratar se necessário.
        # A diagonal de Ze geralmente é tratada de forma diferente (impedância interna).
        Ze = self.jw_mu0_2pi * (np.log(D_matrix) - np.log(d_matrix))

        return Ze

    def external_admittance(self, type_form='potentials'):
        """
        This method calculates the external admittance matrix of the system using a 
        vectorized approach for the Maxwell Potential Coefficients matrix (Pe).
        """

        N = len(self.model.surfaces)
        Pe = np.zeros((N, N))
        Ge = np.zeros((N, N))

        # Busca as matrizes pré-calculadas do modelo
        d_matrix = self.model.d_matrix_ground_return
        D_matrix = self.model.D_matrix_ground_return

        # Loop over the conductors
        # for n, _ in enumerate(self.mtl.surfaces):
        #     for m, _ in enumerate(self.mtl.surfaces):

        #         # Distance between the conductors (m)
        #         dnm, Dnm, _, _ = self.conductor_distances(n, m)

        #         # External impedance (ohm/m)
        #         Pe[n, m] = 1 / (2 * np.pi * epsilon_0) * np.log(Dnm/dnm)

        #         # External admittance (mho/m)
        #         Ge[n, n] = self.ge[n] * 1E-6

        # 3. Calculate the complete Pe matrix
        Pe = (1 / (2 * np.pi * epsilon_0)) * (np.log(D_matrix) - np.log(d_matrix))

        # External capacitance matrix
        if type_form == 'potentials':
            Ce = np.linalg.inv(Pe)

        elif type_form == 'indirect':
            Le = self.external_impedance() / self.jw
            Ce = np.linalg.inv(Le) * mu_0 * epsilon_0

        return Ge + self.jw * Ce

    def simplified_soil_admittance(self):
        """ This method calculates the external admittance matrix of the system. """

        # Earth return impedance matrix
        Zg = self.earth_return_impedance(type_form='carson')

        return -self.ke2 * np.linalg.inv(Zg)

    def sn_sommerfeld(self, hn_hm, dn_dm, ke2, ka2, s_form='s1'):
        """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""

        # soil refractive index
        if s_form == 's1':
            n2 = 1
        else:
            n2 = ke2 / ka2

        # Define the real and imaginary parts of the integrand
        def int_real(x, hn_hm, dn_dm, ke2, ka2, n2):
            num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
            den = np.sqrt(x**2 + ka2 - ke2) + n2 * x
            return (num/den).real

        def int_imag(x, hn_hm, dn_dm, ke2, ka2, n2):
            num = np.exp(-hn_hm * x) * np.cos(dn_dm * x)
            den = np.sqrt(x**2 + ka2 - ke2) + n2 * x
            return (num/den).imag

        # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
        quad_real, _ = quad(int_real, 0, np.inf, args=(hn_hm, dn_dm, ke2, ka2, n2))
        quad_imag, _ = quad(int_imag, 0, np.inf, args=(hn_hm, dn_dm, ke2, ka2, n2))

        return quad_real + 1j * quad_imag

    def sn_sommerfeld_gauss(self, hn_hm, dn_dm, ke2, ka2, s_form='s1', type_form='quad', pts=150):
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
                int_transformed_real, 0, np.pi/2, pts, hn_hm, dn_dm, ke2, ka2, n)
            gauss_legendre_imag = int_gauss_legendre(
                int_transformed_imag, 0, np.pi/2, pts, hn_hm, dn_dm, ke2, ka2, n)
            S_n = gauss_legendre_real + 1j * gauss_legendre_imag

        else:
            # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
            quad_real, _ = quad(int_quad_real, 0, np.inf,
                                args=(hn_hm, dn_dm, ke2, ka2, n))
            quad_imag, _ = quad(int_quad_imag, 0, np.inf,
                                args=(hn_hm, dn_dm, ke2, ka2, n))
            S_n = quad_real + 1j * quad_imag

        return S_n

    def t_sommerfeld(self, hn_hm, dn_dm, k_e2, k_a2, type_form='gauss_legendre', pts=150):
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
            gauss_legendre_real = _int_gauss_legendre(_int_transformed_real, 0, np.pi/2, pts, hn_hm, dn_dm, k_e2, k_a2)
            gauss_legendre_imag = _int_gauss_legendre(_int_transformed_imag, 0, np.pi/2, pts, hn_hm, dn_dm, k_e2, k_a2)
            T = gauss_legendre_real + 1j * gauss_legendre_imag

        else:
            # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
            quad_real, _ = quad(_int_quad_real, 0, np.inf, args=(hn_hm, dn_dm, k_e2, k_a2))
            quad_imag, _ = quad(_int_quad_imag, 0, np.inf, args=(hn_hm, dn_dm, k_e2, k_a2))
            T = quad_real + 1j * quad_imag

        return T

    def sn_sommerfeld_simplified(self, hn_hm, dn_dm, k_e2, n2=1, type_form='gauss_legendre', pts=150):
        """ Calculates the Carson integral using Gauss-Legendre quadrature and scipy.integrate.quad."""  

        # Define the real and imaginary parts of the integrand
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
            gauss_legendre_real = int_gauss_legendre(
                int_transformed_real, 0, np.pi/2, pts, hn_hm, dn_dm, k_e2, n2)
            gauss_legendre_imag = int_gauss_legendre(
                int_transformed_imag, 0, np.pi/2, pts, hn_hm, dn_dm, k_e2, n2)
            Sn = gauss_legendre_real + 1j * gauss_legendre_imag

        else:
            # Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
            quad_real, _ = quad(int_quad_real, 0, np.inf, args=(hn_hm, dn_dm, k_e2, n2))
            quad_imag, _ = quad(int_quad_imag, 0, np.inf, args=(hn_hm, dn_dm, k_e2, n2))
            Sn = quad_real + 1j * quad_imag

        return Sn

    def earth_return_impedance(self, type_form='quasi_tem'):
        """ This method calculates the impedance matrix of the earth return path. """

        # Impedance matrix
        N = len(self.model.surfaces)
        Zg = np.zeros((N, N), dtype=complex)

        # Complex depth (m)
        if type_form == 'approx_log':
            p_dot = 1 / np.sqrt(-self.ke2)

        if type_form == 'carson' or type_form == 'deri':
            # Soil wave number (rad/m)
            k_e2 = - self.jw * self.mur_1 * mu_0 * self.sigma_1

            # A. Deri Complex depth (m)
            p_dot = 1 / np.sqrt(-k_e2)

        # Loop over the conductors
        for n, conductor_n in enumerate(self.model.surfaces):
            for m, conductor_m in enumerate(self.model.surfaces):

                # Distance between the conductors (m)
                hn_hm = self.model.vertical_separation_matrix[n, m]
                dn_dm = self.model.horizontal_separation_matrix[n, m]

                # Quasi-TEM Integral Equation
                if type_form == 'quasi_tem':
                    S1 = 2 * self.sn_sommerfeld(hn_hm, dn_dm, self.ke2, self.ka2)

                # Quasi-TEM Logarithmic Approximation
                elif type_form == 'quasitem_log':
                    eta = np.sqrt(self.ka2 - self.ke2)
                    eta_sqrt = eta * np.sqrt(hn_hm**2 + dn_dm**2)
                    S1 = np.log(1 + 2 / eta_sqrt)

                # Sunde Integral Equation
                elif type_form == 'sunde':
                    if n == m:
                        dn_dm = 0

                    # Sommerfeld Integral
                    S1 = 2 * self.sn_sommerfeld(hn_hm, dn_dm, self.ke2, ka2=0)

                # Carson Integral Equation
                elif type_form == 'carson':
                    if n == m:
                        dn_dm = 0

                    # Sommerfeld Integral
                    S1 = 2 * self.sn_sommerfeld(hn_hm, dn_dm, k_e2, ka2=0)

                # Log. Approximation Closed-Form Expression
                else:
                    if n == m:
                        hn = conductor_n['center_point'][1]
                        S1 = np.log((hn + p_dot) / hn)

                    else:
                        num = np.sqrt((hn_hm + 2*p_dot)**2 + dn_dm**2)
                        den = np.sqrt(hn_hm**2 + dn_dm**2)
                        S1 = np.log(num / den)

                # Earth return impedance (ohm/m)
                Zg[n, m] = self.jw_mu0_2pi * S1

        return Zg

    def pul_extended_theory(self, type_form='nakagawa'):
        """ This method calculates the impedance matrix of the earth return path. """

        # Impedance matrix
        N = len(self.model.surfaces)
        M = np.zeros((N, N))
        S1 = np.zeros((N, N), dtype=complex)
        S2 = np.zeros((N, N), dtype=complex)
        T = np.zeros((N, N), dtype=complex)
        Zi = self.internal_impedance()[0]

        # Loop over the conductors
        for n, _ in enumerate(self.model.surfaces):
            for m, _ in enumerate(self.model.surfaces):

                # Distance between the conductors (m)
                # dnm, Dnm, hn_hm, dn_dm = self.conductor_distances(conductor_n['tag'], conductor_m['tag'])
                # Distance between the conductors (m)
                dnm = self.model.d_matrix_ground_return[n, m]
                Dnm = self.model.D_matrix_ground_return[n, m]
                hn_hm = self.model.vertical_separation_matrix[n, m]
                dn_dm = self.model.horizontal_separation_matrix[n, m]

                # External Impedance term
                M[n, m] = np.log(Dnm/dnm)

                # Quasi-TEM Integral Equation
                if type_form == 'quasitem':
                    S1[n, m] = 2 * self.sn_sommerfeld_gauss(hn_hm, dn_dm, self.ke2, self.ka2)
                    S2[n, m] = 2 * self.sn_sommerfeld_gauss(hn_hm, dn_dm, self.ke2, self.ka2, s_form='s2')
                    T[n, m] = 2 * self.t_sommerfeld(hn_hm, dn_dm, self.ke2, self.ka2)

                # Quasi-TEM Logarithmic Approximation
                elif type_form == 'quasitem_log':
                    eta = np.sqrt(self.ka2 - self.ke2)
                    eta_sqrt = eta * np.sqrt(hn_hm**2 + dn_dm**2)
                    n2 = self.ke2 / self.ka2

                    log_s1 = 1 + 2 / eta_sqrt
                    log_s2 = 1 + (n2 + 1) / eta_sqrt
                    log_t = 1 + 2*(n2 + 1) / eta_sqrt

                    S1[n, m] = np.log(log_s1)
                    S2[n, m] = 2 / (n2 + 1) * np.log(log_s2)
                    T[n, m] = 2 * np.log(2) + (2 * n2 / (n2 + 1)) * np.log(log_s2/log_t)

                # Nakagawa Integral Equation
                elif type_form == 'nakagawa':
                    # soil refractive index
                    n2_naka = self.ke2 / self.ka2

                    # Earth wave number (rad/m)
                    ke2_naka = - self.jw * self.mur_1 * mu_0 * \
                        (self.sigma_1 + self.jw * (self.er_1 - 1) * epsilon_0)

                    S1[n, m] = 2 * self.sn_sommerfeld_simplified(hn_hm, dn_dm, ke2_naka)
                    S2[n, m] = 2 * self.sn_sommerfeld_simplified(hn_hm, dn_dm, ke2_naka, n2=n2_naka)

                # Sunde Integral Equation
                elif type_form == 'sunde':
                    S1[n, m] = 2 * self.sn_sommerfeld_simplified(hn_hm, dn_dm, self.ke2)

                # Carson Integral Equation
                elif type_form == 'carson':
                    k_e2_carson = - self.jw * self.mur_1 * mu_0 * self.sigma_1
                    S1[n, m] = 2 * self.sn_sommerfeld_simplified(hn_hm, dn_dm, k_e2_carson)

        # Earth return impedance and admittance (ohm/m)
        if type_form == 'quasitem' or type_form == 'quasitem_log':
            Zg = self.jw_mu0_2pi * (M + S1 - (T + S2))
            Y = self.jw_2pi_e0 * np.linalg.inv(M - T)

        else:
            Zg = self.jw_mu0_2pi * (M + S1)
            Y = self.jw_2pi_e0 * np.linalg.inv(M + S2)

        # Propagation constant (rad/m)
        gamma = np.sqrt((Zi+Zg) @ Y)

        return {'Zi': Zi, 'Zg': Zg, 'Zs': Zi+Zg, 'Y': Y, 'γ': gamma}
