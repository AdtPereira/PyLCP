import numpy as np
from scipy.integrate import quad
from scipy.constants import mu_0

# Definir constantes (substitua pelos valores reais)
hn_hm = 20.0    # valor para h_i
dn_dm = 0.01   # valor para x_ij
w = 2 * np.pi * 1E3  # frequência angular
sigma_g = 1/200 # condutividade (S/m)
Nq = 100 # pylint: disable=invalid-name
k_e2 = -1j * (2 * np.pi * 1E3) * mu_0 * sigma_g
k_a2 = 0

def integrand_quad_real(λ, hn_hm, dn_dm, k_e2, k_a2):
    num = np.exp(-hn_hm * λ) * np.cos(dn_dm * λ)
    den = np.sqrt(λ**2 + k_a2 - k_e2) + λ
    return (num/den).real

def integrand_quad_imag(λ, hn_hm, dn_dm, k_e2, k_a2):
    num = np.exp(-hn_hm * λ) * np.cos(dn_dm * λ)
    den = np.sqrt(λ**2 + k_a2 - k_e2) + λ
    return (num/den).imag

def integrand_transformed_real(t, hn_hm, dn_dm, k_e2, k_a2):
    λ = np.tan(t)
    num = np.exp(-hn_hm * λ) * np.cos(dn_dm * λ)
    den = np.sqrt(λ**2 + k_a2 - k_e2) + λ
    return (num/den).real * (1 / np.cos(t)**2)

def integrand_transformed_imag(t, hn_hm, dn_dm, k_e2, k_a2):
    λ = np.tan(t)
    num = np.exp(-hn_hm * λ) * np.cos(dn_dm * λ)
    den = np.sqrt(λ**2 + k_a2 - k_e2) + λ
    return (num/den).imag * (1 / np.cos(t)**2)

def gauss_legendre_integrate(func, a, b, n, *args):
    """
    Integra a função `func` no intervalo [a, b] usando o método de Gauss-Legendre.
    - func: função a ser integrada.
    - a, b: limites de integração.
    - n: número de pontos de quadratura.
    - args: argumentos adicionais a serem passados para `func`.
    """
    [x, w] = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (x + 1) * (b - a) + a
    return np.sum(w * func(t, *args)) * 0.5 * (b - a)

# Calculate the real and imaginary parts of the integral using Gauss-Legendre
gauss_legendre_real = gauss_legendre_integrate(integrand_transformed_real, 0, np.pi/2, Nq, hn_hm, dn_dm, k_e2, k_a2)
gauss_legendre_imag = gauss_legendre_integrate(integrand_transformed_imag, 0, np.pi/2, Nq, hn_hm, dn_dm, k_e2, k_a2)
S1_gauss_legendre = gauss_legendre_real + 1j * gauss_legendre_imag # pylint: disable=invalid-name

# Calculate the real and imaginary parts of the integral using scipy.integrate.Quad
quad_real, _ = quad(integrand_quad_real, 0, np.inf, args=(hn_hm, dn_dm, k_e2, k_a2))
quad_imag, _ = quad(integrand_quad_imag, 0, np.inf, args=(hn_hm, dn_dm, k_e2, k_a2))
S1_quad = quad_real + 1j * quad_imag # pylint: disable=invalid-name

# Print the results
print(f"Resultado com Gauss-Legendre: {2 * S1_gauss_legendre}")
print(f"Resultado com quad: {2 * S1_quad}")
