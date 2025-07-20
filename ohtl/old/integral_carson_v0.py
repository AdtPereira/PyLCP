import numpy as np
from scipy.integrate import quad

# Definir constantes (substitua pelos valores reais)
h_i = 1.0    # valor para h_i
h_j = 1.0    # valor para h_j
x_ij = 1.0   # valor para x_ij
omega = 1.0  # frequência angular
mu_0 = 4 * np.pi * 1e-7  # permeabilidade magnética do vácuo (H/m)
sigma = 1.0 # condutividade (S/m)
g = 1.0     # constante g

# Função integranda original para usar no quad - Parte Real
def integrand_quad_real(λ, h_i, h_j, x_ij, omega, mu_0, sigma, g):
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma * g) + λ
    numerador = np.exp(-(h_i + h_j) * λ) * np.cos(x_ij * λ)
    return (numerador / complex_term).real

# Função integranda original para usar no quad - Parte Imaginária
def integrand_quad_imag(λ, h_i, h_j, x_ij, omega, mu_0, sigma, g):
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma * g) + λ
    numerador = np.exp(-(h_i + h_j) * λ) * np.cos(x_ij * λ)
    return (numerador / complex_term).imag

# Função integranda transformada com o cosseno para usar com gauss_legendre_integrate - Parte Real
def integrand_transformed_real(t, h_i, h_j, x_ij, omega, mu_0, sigma, g):
    λ = np.tan(t)
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma * g) + λ
    numerador = np.exp(-(h_i + h_j) * λ) * np.cos(x_ij * λ)
    return (numerador / complex_term).real * (1 / np.cos(t)**2)

# Função integranda transformada com o cosseno para usar com gauss_legendre_integrate - Parte Imaginária
def integrand_transformed_imag(t, h_i, h_j, x_ij, omega, mu_0, sigma, g):
    λ = np.tan(t)
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma * g) + λ
    numerador = np.exp(-(h_i + h_j) * λ) * np.cos(x_ij * λ)
    return (numerador / complex_term).imag * (1 / np.cos(t)**2)

# Implementação da integração usando Gauss-Legendre
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
    integral = np.sum(w * func(t, *args)) * 0.5 * (b - a)
    return integral

# Número de pontos de quadratura
n_points = 100

# Calcula a parte real da integral usando Gauss-Legendre
result_gauss_legendre_real = gauss_legendre_integrate(integrand_transformed_real, 0, np.pi/2, n_points, h_i, h_j, x_ij, omega, mu_0, sigma, g)

# Calcula a parte imaginária da integral usando Gauss-Legendre
result_gauss_legendre_imag = gauss_legendre_integrate(integrand_transformed_imag, 0, np.pi/2, n_points, h_i, h_j, x_ij, omega, mu_0, sigma, g)

# Calcula a parte real da integral usando quad
result_quad_real, _ = quad(integrand_quad_real, 0, np.inf, args=(h_i, h_j, x_ij, omega, mu_0, sigma, g))

# Calcula a parte imaginária da integral usando quad
result_quad_imag, _ = quad(integrand_quad_imag, 0, np.inf, args=(h_i, h_j, x_ij, omega, mu_0, sigma, g))

# Combina as partes real e imaginária para obter o resultado complexo
result_gauss_legendre = result_gauss_legendre_real + 1j * result_gauss_legendre_imag
result_quad = result_quad_real + 1j * result_quad_imag

# Imprime os resultados
print(f"Resultado com Gauss-Legendre: {result_gauss_legendre}")
print(f"Resultado com quad: {result_quad}")
