import numpy as np
from scipy.special import jv, jvp
from scipy.constants import mu_0
import matplotlib.pyplot as plt

from data.mtl import MulticonductorTransmissionLine


class Bifilar(MulticonductorTransmissionLine):
    """ This class contains the analytical formulation of the system. """

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

    def series_impedance(self, f):
        """
        This function calculates the series resistance of the system using 
        the high frequency approximation.

        Returns:
        tuple: A tuple containing the high frequency resistance, external inductance, 
        and matrix impedance.
        """
        # Angular frequency, rad/s [float]
        w = 2 * np.pi * f

        # Skin Depth [np.array]
        delta = np.sqrt(1 / (w / 2 * mu_0 * self.sigma))

        # Surface Resistance [np.array]
        Rs = 1 / (self.sigma * delta) # pylint: disable=invalid-name

        # Outer Radii of the conductors [np.array]
        ap = np.array([cp['radius'][1] for cp in self.mtl])

        # Matrix Distance [np.array]
        D = self.D_pq  # pylint: disable=invalid-name

        # High Frequency Resistance and External Inductance [np.array]
        N = len(self.mtl)-1  # pylint: disable=invalid-name
        Rhf = np.zeros((N, N))  # pylint: disable=invalid-name
        Lext = np.zeros_like(Rhf)  # pylint: disable=invalid-name
        Zi = np.zeros_like(Rhf, dtype=complex)  # pylint: disable=invalid-name

        # Constant Term and Bessel argument
        Xi = np.sqrt(2) * ap / delta # pylint: disable=invalid-name
        constant_term = 1 / (np.sqrt(2) * np.pi * ap * self.sigma * delta)

        for p in range(N):
            # 1st solution: High Frequency Approximation
            # These formulas account for proximity effect
            # only at high frequencies.These formulas
            # account for proximity effect only at high
            # frequencies.

            # Common fraction term
            D_2a = D[p][p+1] / 2 / ap[p] # pylint: disable=invalid-name

            # Surface resistance
            Rs_pia = Rs[p] / np.pi / ap[p] # pylint: disable=invalid-name

            # High Frequency Resistance (Ω/m)
            # Equation (2.64) [1]
            Rhf[p] = Rs_pia * D_2a / np.sqrt(D_2a ** 2 - 1)

            # External Inductance (H/m)
            # Equation (2.65) [1]
            Lext[p] = mu_0 / np.pi * np.arccosh(D_2a)

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
        Zs = 2 * Zi + 1j * w * Lext  # pylint: disable=invalid-name

        return Zs, Rhf, Lext

    def capacitance(self):
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
        # Para uma linha bifilar (condutores 0 e 1), s é D[0, 1].
        s = self.D_pq[0, 1]

        # r_w1, r_w2: raios dos dois condutores.
        radii = np.array([conductor['radius'][1] for conductor in self.mtl])
        if len(radii) != 2:
            raise ValueError("Esta fórmula de capacitância é válida apenas para uma linha de dois condutores (bifilar).")
        r_w1 = radii[0]
        r_w2 = radii[1]

        two_pi_epsilon = 2 * np.pi * epsilon
        
        # --- Cálculo da Capacitância Aproximada (Eq. 4.21) ---
        den_approx = np.log((s**2) / (r_w1 * r_w2))
        capacitance_approx = two_pi_epsilon / den_approx

        # --- Cálculo da Capacitância Exata (Eq. 4.39) ---
        arg_arccosh = (s**2 - r_w1**2 - r_w2**2) / (2 * r_w1 * r_w2) #
        den_exact = np.arccosh(arg_arccosh) #
        capacitance_exact = two_pi_epsilon / den_exact #

        # Retorna um dicionário com ambos os resultados
        return {
            'approximate': capacitance_approx,
            'exact': capacitance_exact
        }

    def plot_series_resistance(self, f, z, rhf, data):
        """
        This function plots the series resistance as a function of frequency.

        Parameters:
        freq (array): Frequency array.
        zi_matrix (list of matrices): Matrix containing impedance values.
        p (int): Row index in the impedance matrix.
        q (int): Column index in the impedance matrix.
        """

        # Extracting the data from the list_data
        p = data[0]
        Np = data[1] # pylint: disable=invalid-name
        D = data[2][0,1] # pylint: disable=invalid-name

        # Extracting the impedance elements
        zs = np.array([item[p] for item in z[0]])
        rhf = np.array([item[p] for item in rhf])

        # Asymptotic Series Resistance
        plt.plot(f[0], 1E3 * rhf, label='Asymptotic',
                color='red', linestyle='--')

        # Closed-Form Approximation Series Resistance
        plt.plot(f[0], 1E3 * np.real(zs), label='Skin Effect Only',
                color='black', linestyle='-')

        # MoM Series Resistance
        plt.scatter(f[1], 1E3 * np.real(z[1]),
                    label=fr'MoM-SO: $Np=Nq={Np}$ [1]',
                    color='blue', marker='x', s=50)

        # Optional: Additional plotting configurations like labels, grid, etc.
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1E0, 1E7)
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Resistance p.u.l. (Ω/km)')
        plt.title('Figure 2.4: P.u.l. series resistance, $R_{int}$, of the bifilar overhead line\n'
                r'$r_1=r_2=0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},' fr'D={D}\,m,'
                r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
        plt.grid(True)
        plt.show()

    def plot_series_inductance(self, f, z, lext_hf, data):
        """ This function plots the series inductance as a function of frequency. """

        # Extracting the data from the list_data
        p = data[0]
        Np = data[1]  # pylint: disable=invalid-name
        D = data[2][0, 1]  # pylint: disable=invalid-name

        # Extracting the impedance elements
        zs = np.array([item[p] for item in z[0]])
        lext = np.array([np.imag(zs) / (2 * np.pi * f) for zs, f in zip(zs, f[0])])
        lext_hf = np.array([item[p] for item in lext_hf])
        lext_mom = np.array([np.imag(z) / (2 * np.pi * f)
                            for z, f in zip(z[1], f[1])])

        # Closed-Form Approximation Series Inductance
        plt.plot(f[0], 1E6 * lext, label='Analytical (Skin Effect Only)',
                color='black', linestyle='-')

        # Asymptotic Series Inductance
        plt.plot(f[0], 1E6 * lext_hf, label='Analytical (Asymptotic)',
                color='red', linestyle='--')

        # MoM Series Inductance
        plt.scatter(f[1], 1E6 * lext_mom, label=fr'MoM-SO ($Np=Nq={Np}$) [1]',
                    color='blue', marker='x', s=40)

        plt.xscale('log')
        plt.legend()
        plt.xlim(1E0, 1E7)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Series Inductance p.u.l. (mH/km)')
        plt.title('Figure 2.5: P.u.l. inductance of the bifilar overhead line\n'
                r'$r_1=r_2=0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m},' fr'D={D}\,m,'
                r'\sigma = 5.952 \times 10^7\,\mathrm{S/m}$ [1]')
        plt.grid(False)
        plt.show()


