import numpy as np
import matplotlib.pyplot as plt

def solve_bifilar_mom_and_exact(D_over_R, NF, verbose=False):
    """
    Calcula a capacitância em uma linha bifilar usando o MoM com pontos de colocação
    rotacionados para garantir estabilidade.
    """
    epsilon_0 = 8.854187817e-12
    R = 1.0
    D = D_over_R * R
    h = D / 2.0
    centers = [np.array([-h, 0.0]), np.array([h, 0.0])]

    N_coeffs_per_wire = NF
    N_unknowns = 2 * N_coeffs_per_wire
    D_matrix = np.zeros((N_unknowns, N_unknowns))
    V_vector = np.zeros(N_unknowns)

    initial_angles = np.linspace(0, 2 * np.pi, N_coeffs_per_wire, endpoint=False)
    rotation_angle_rad = np.pi / (N_coeffs_per_wire + 1)
    match_angles = initial_angles + rotation_angle_rad

    for i in range(N_unknowns):
        field_wire_idx = i // N_coeffs_per_wire
        field_angle = match_angles[i % N_coeffs_per_wire]
        field_center = centers[field_wire_idx]
        V_vector[i] = 1.0 if field_wire_idx == 0 else -1.0
        field_point = field_center + R * np.array([np.cos(field_angle), np.sin(field_angle)])

        for j in range(N_unknowns):
            source_wire_idx = j // N_coeffs_per_wire
            source_harmonic_order = j % N_coeffs_per_wire
            source_center = centers[source_wire_idx]

            if source_wire_idx == field_wire_idx:
                if source_harmonic_order == 0:
                    D_matrix[i, j] = (-R / epsilon_0) * np.log(R)
                else:
                    k = source_harmonic_order
                    D_matrix[i, j] = (R / (2 * k * epsilon_0)) * np.cos(k * field_angle)
            else:
                dist_vec = field_point - source_center
                dist = np.linalg.norm(dist_vec)
                if source_harmonic_order == 0:
                    D_matrix[i, j] = (-R / epsilon_0) * np.log(dist)
                else:
                    k = source_harmonic_order
                    phi = np.arctan2(dist_vec[1], dist_vec[0])
                    D_matrix[i, j] = (R / (2 * k * epsilon_0)) * ((R / dist)**k) * np.cos(k * phi)

    sigma_coeffs = np.linalg.solve(D_matrix, V_vector)
    sigma_0_wire1 = sigma_coeffs[0]
    C_mom = np.pi * R * sigma_0_wire1
    C_exact = (np.pi * epsilon_0) / np.arccosh(D_over_R / 2.0)
    
    return C_mom, C_exact

if __name__ == '__main__':
    # --- Parâmetros da Simulação ---
    D_R_ratio = 2.1
    nf_range = range(1, 21)
    
    # --- ALTERAÇÃO APLICADA AQUI ---
    # Em vez de um valor fixo, calculamos o valor exato de referência uma vez.
    # O valor de NF usado aqui (ex: 1) não afeta o cálculo de C_exact.
    _, c_exact = solve_bifilar_mom_and_exact(D_R_ratio, NF=20, verbose=False)

    # --- Coleta de Dados ---
    nf_results = []
    capacitance_results = []
    C_FACTOR = 1e12  # Fator de conversão para pF/m
    
    for nf in nf_range:
        c_mom, _ = solve_bifilar_mom_and_exact(D_R_ratio, nf)
        nf_results.append(nf)
        capacitance_results.append(c_mom)

    # Separar dados em par e ímpar para plotar com marcadores diferentes
    nf_odd = [nf for nf in nf_results if nf % 2 != 0]
    cap_odd = [capacitance_results[i] * C_FACTOR for i, nf in enumerate(nf_results) if nf % 2 != 0]
    nf_even = [nf for nf in nf_results if nf % 2 == 0]
    cap_even = [capacitance_results[i] * C_FACTOR for i, nf in enumerate(nf_results) if nf % 2 == 0]

    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axhline(y=c_exact*C_FACTOR, color='k', linestyle='-', label=f'Valor Exato = {c_exact*C_FACTOR:.2f} pF/m')

    # Triângulos (Δ) para valores ímpares
    # Asteriscos (*) para valores pares
    ax.plot(nf_odd, cap_odd, linestyle='none', marker='^', markersize=8,  fillstyle='none', markeredgecolor='black', label='NF Ímpar (Δ)')
    ax.plot(nf_even, cap_even, linestyle='none', marker='*', markersize=9, color='black', label='NF Par (*)')

    # 3. Formatação e legendas do gráfico
    ax.set_title(f'Convergência da Capacitância para D/R = {D_R_ratio}')
    ax.set_xlabel('NF (Número de Coeficientes de Fourier)')
    ax.set_ylabel('Capacitância (pF/m)')
    ax.set_xlim(0, 21)
    ax.set_ylim(0, 100)
    # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.grid(True, linestyle=':', alpha=0.7)
    ax.legend()
    plt.show()