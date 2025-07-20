import numpy as np
from scipy.constants import mu_0
from scipy.integrate import quad

# Definir constantes (substitua pelos valores reais)
h_i = 10.0    # valor para h_i
h_j = 10.0    # valor para h_j
x_ij = 0.01   # valor para x_ij
omega = 2 * np.pi * 1E5  # frequência angular
sigma_g = 1/200 # condutividade (S/m)

# Função integranda original para usar no quad - Parte Real
def integrand_quad_real(λ, h_i, h_j, x_ij, omega, mu_0, sigma_g):
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma_g) + λ
    numerador = np.exp(-(h_i + h_j) * λ) * np.cos(x_ij * λ)
    return (numerador / complex_term).real

# Função integranda original para usar no quad - Parte Imaginária
def integrand_quad_imag(λ, h_i, h_j, x_ij, omega, mu_0, sigma_g):
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma_g) + λ
    numerador = np.exp(-(h_i + h_j) * λ) * np.cos(x_ij * λ)
    return (numerador / complex_term).imag

# Função integranda transformada com o termo do cosseno
def integrand_transformed_cos(t, h_i, h_j, x_ij, omega, mu_0, sigma_g):
    λ = np.tan(t)
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma_g) + λ
    numerador = np.exp(-(h_i + h_j) * λ) * np.cos(x_ij * λ)
    return (numerador / complex_term).real * (1 / np.cos(t)**2)

# Implementação da integração usando Gauss-Legendre
def gauss_legendre_integrate(func, a, b, n, *args):
    """
    Integra a função `func` no intervalo [a, b] usando o método de Gauss-Legendre.
    - func: função a ser integrada.
    - a, b: limites de integração.
    - n: número de pontos de quadratura.
    - args: argumentos adicionais a serem passados para `func`.
    """
    # Obtém os pontos de Legendre e pesos
    [x, w] = np.polynomial.legendre.leggauss(n)
    
    # Transforma os pontos para o intervalo [a, b]
    t = 0.5 * (x + 1) * (b - a) + a
    # Calcula a função nos pontos transformados e aplica os pesos
    integral = np.sum(w * func(t, *args)) * 0.5 * (b - a)
    return integral

# Número de pontos de quadratura
n_points = 100  # Escolha um valor maior para mais precisão

# Calcula a parte real da integral no intervalo [0, pi/2)
result_real = gauss_legendre_integrate(integrand_transformed_cos, 0, np.pi/2, n_points, h_i, h_j, x_ij, omega, mu_0, sigma_g)

print(f"Resultado da parte real da integral com cosseno: {result_real}")

# Calculando a parte imaginária da integral (similarmente)
def integrand_imag_transformed_cos(t, h_i, h_j, x_ij, omega, mu_0, sigma_g):
    λ = np.tan(t)
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma_g) + λ
    numerador = np.exp(-(h_i + h_j) * λ) * np.cos(x_ij * λ)
    return (numerador / complex_term).imag * (1 / np.cos(t)**2)

result_imag = gauss_legendre_integrate(integrand_imag_transformed_cos, 0, np.pi/2, n_points, h_i, h_j, x_ij, omega, mu_0, sigma_g)

print(f"Resultado da parte imaginária da integral com cosseno: {result_imag}")

# Combinar parte real e imaginária para obter o resultado completo
S_complex = result_real + 1j * result_imag

print(f"Resultado completo da integral com cosseno: {S_complex}")
print(f"2*Sij (gauss_legendre): {2*S_complex}")

# Calcula a parte real da integral usando quad
result_quad_real, _ = quad(integrand_quad_real, 0, np.inf, args=(h_i, h_j, x_ij, omega, mu_0, sigma_g))
result_quad_imag, _ = quad(integrand_quad_imag, 0, np.inf, args=(h_i, h_j, x_ij, omega, mu_0, sigma_g))
print(f"2*Sij (quad): {2*(result_quad_real+1j*result_quad_imag)}")