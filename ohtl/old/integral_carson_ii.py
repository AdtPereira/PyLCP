import numpy as np
from scipy.constants import mu_0

# Definir constantes (substitua pelos valores reais)
h_i = 10    # valor para h_i
omega = 2 * np.pi * 1E3  # frequência angular
sigma_g = 1/200 # condutividade (S/m)

# Função integranda transformada
def integrand_transformed(t, h_i, omega, mu_0, sigma_g):
    λ = np.tan(t)
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma_g) + λ
    numerador = np.exp(-2 * h_i * λ)
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

# Calcula a integral no intervalo [0, pi/2)
result_real = gauss_legendre_integrate(integrand_transformed, 0, np.pi/2, n_points, h_i, omega, mu_0, sigma_g)

print(f"Resultado da parte real da integral usando Gauss-Legendre: {result_real}")

# Calculando a parte imaginária da integral (similarmente)
def integrand_imag_transformed(t, h_i, omega, mu_0, sigma_g):
    λ = np.tan(t)
    complex_term = np.sqrt(λ**2 + 1j * omega * mu_0 * sigma_g) + λ
    numerador = np.exp(-2 * h_i * λ)
    return (numerador / complex_term).imag * (1 / np.cos(t)**2)

result_imag = gauss_legendre_integrate(integrand_imag_transformed, 0, np.pi/2, n_points, h_i, omega, mu_0, sigma_g)

print(f"Resultado da parte imaginária da integral usando Gauss-Legendre: {result_imag}")

# Combinar parte real e imaginária para obter o resultado completo
S_complex = result_real + 1j * result_imag

print(f"2*S_ii: {2*S_complex}")
