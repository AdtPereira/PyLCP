import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0
from matplotlib.patches import Circle

def calculate_collocation_points(R, D, NF):
    """
    Calcula as coordenadas dos pontos de colocação da fonte e da observação.

    Retorna tanto as coordenadas cartesianas quanto os ângulos em radianos.

    Args:
        R (float): Raio dos fios.
        D (float): Distância entre os centros dos fios.
        NF (int): Número de harmônicos (pontos por fio).

    Returns:
        dict: Dicionário contendo coordenadas cartesianas e ângulos.
    """
    centers = [np.array([-D / 2.0, 0.0]), np.array([D / 2.0, 0.0])]

    source_angles = np.linspace(0, 2 * np.pi, NF, endpoint=False)
    observation_angles = source_angles + np.pi / (NF + 1)

    points = {'source': [], 'observation': []}
    for i in range(2):
        center = centers[i]
        source_pts = center + R * np.array([np.cos(source_angles), np.sin(source_angles)]).T
        obs_pts = center + R * np.array([np.cos(observation_angles), np.sin(observation_angles)]).T
        points['source'].append(source_pts)
        points['observation'].append(obs_pts)

    return {
        'cartesian': points,
        'angles_rad': {
            'source': source_angles,
            'observation': observation_angles
        }
    }

def plot_collocation_points(points_data, R, D, coord_mode='Cartesian'):
    """
    Plota os condutores e os pontos de colocação.

    A anotação de cada ponto mostra as coordenadas cartesianas globais (x,y)
    e as coordenadas polares locais (r, θ) com o ângulo em radianos.

    Args:
        points_data (dict): Dicionário com dados dos pontos.
        R (float): Raio dos fios.
        D (float): Distância entre os centros dos fios.
    """
    fig, ax = plt.subplots(figsize=(12, 9))

    cartesian_points = points_data['cartesian']
    angles_rad = points_data['angles_rad']

    centers = [np.array([-D / 2.0, 0.0]), np.array([D / 2.0, 0.0])]
    wire1 = Circle(centers[0], R, facecolor='none', edgecolor='k', lw=1, ls='--')
    wire2 = Circle(centers[1], R, facecolor='none', edgecolor='k', lw=1, ls='--')
    ax.add_patch(wire1)
    ax.add_patch(wire2)

    source_w1, source_w2 = cartesian_points['source']
    obs_w1, obs_w2 = cartesian_points['observation']

    ax.plot(source_w1[:, 0], source_w1[:, 1], 'o', c='blue', label='Source')
    ax.plot(obs_w1[:, 0], obs_w1[:, 1], 'x', c='red', label='Observation')
    ax.plot(source_w2[:, 0], source_w2[:, 1], 'o', c='blue')
    ax.plot(obs_w2[:, 0], obs_w2[:, 1], 'x', c='red')

    # --- REVISÃO DA ANOTAÇÃO PARA APRESENTAR AMBOS OS SISTEMAS DE COORDENADAS ---
    point_sets = [
        (source_w1, angles_rad['source']),
        (obs_w1, angles_rad['observation']),
        (source_w2, angles_rad['source']),
        (obs_w2, angles_rad['observation'])
    ]

    for points, angles in point_sets:
        for i, pt in enumerate(points):
            # Coordenadas cartesianas globais
            global_x, global_y = pt[0], pt[1]
            
            # Coordenadas polares locais
            local_r = R
            local_theta_rad = angles[i]
            
            # Formatação da anotação em duas linhas
            if coord_mode == 'Polar':
                # Coordenadas polares locais
                annotation = (f'({local_r:.3f}, {local_theta_rad:.3f} rad)')
            elif coord_mode == 'Cartesian':
                # Coordenadas cartesianas globais
                annotation = (f'({global_x:.3f}, {global_y:.3f})')
            else:
                raise ValueError("coord_mode must be 'Cartesian' or 'Polar'.")
            
            ax.text(pt[0], pt[1] + R * 0.15, annotation, fontsize=8, ha='center', va='bottom', rotation=15)

    d_fator = 1.5
    ax.set_title('Mapa de Pontos de Colocação')
    ax.set_xlabel('Coordenada X (m)')
    ax.set_ylabel('Coordenada Y (m)')
    ax.set_xlim(-D * d_fator, D * d_fator)
    ax.set_ylim(-D, D)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    ax.legend()

def plot_charge_density(theta_deg, charge_density_mom, charge_density_exact, D_R_ratio, NF_order):
    """
    Plota o gráfico de comparação da densidade de carga (MoM vs. Exata).

    Args:
        theta_deg (array): Ângulos em graus para o eixo X.
        charge_density_mom (array): Densidade de carga calculada pelo MoM.
        charge_density_exact (array): Densidade de carga da solução analítica.
        D_R_ratio (float): Razão D/R para o título do gráfico.
        NF_order (int): Ordem de harmônicos para a legenda e título.
    """
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(theta_deg, charge_density_exact, color='r', linestyle='-', label='Solução Exata')
    ax.plot(theta_deg, charge_density_mom, color='k', linestyle='-.', label=f'MoM (NF={NF_order})')
    ax.set_title(f'Fig. 5. Charge distributions with \n D/R = {D_R_ratio} and {NF_order} harmonic expansion functions per wire.', fontsize=12)
    ax.set_xlabel('Ângulo (Graus)', fontsize=12)
    ax.set_ylabel('Densidade de Carga (C/m²)', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_xticks(np.arange(0, 361, 90))
    ax.set_xlim(0, 360)
    plt.tight_layout()
    ax.legend()

def generalized_capacitance(D_matrix: np.ndarray, 
                                 radii: list[float], 
                                 NF: int) -> np.ndarray:
    """
    Calcula a matriz de capacitância generalizada C a partir da matriz T (D^-1).

    Esta função implementa a Equação (22) do artigo de Clements et al. (1975):
    c_nm = 2 * pi * a_n * sum(T_1i^nm), que significa somar os elementos da
    primeira linha da submatriz T_nm.

    Args:
        D_matrix (np.ndarray): Suas dimensões são (N_cond * NF) x (N_cond * NF).
        radii (list[float]): Uma lista contendo o raio de cada condutor [a_0, a_1, ...].
        NF (int): O número de funções de base (harmônicos) por condutor.

    Returns:
        np.ndarray: A matriz de capacitância generalizada (N_cond x N_cond).
    """
    # Validação das dimensões
    total_functions = D_matrix.shape[0]
    if total_functions % NF != 0:
        raise ValueError("A dimensão da matriz D não é um múltiplo de NF.")

    N = total_functions // NF
    
    if len(radii) != N:
        raise ValueError("O número de raios fornecidos não corresponde ao número de condutores.")

    # Inicializa a matriz de capacitância C com zeros
    T_matrix = np.linalg.inv(D_matrix)
    C_matrix = np.zeros((N, N))

    # Itera sobre os condutores para preencher a matriz C
    # n: índice do condutor da carga (afeta a_n e a linha de T)
    # m: índice do condutor do potencial (afeta a coluna de T)
    for n in range(N):
        for m in range(N):
            # A "primeira linha" da submatriz T_nm corresponde à linha da matriz T
            # completa que está associada à função de base constante (k=0) do condutor 'n'.
            # Em um arranjo 0-indexado, este é o índice 'n * NF'.

            # As colunas da submatriz T_nm correspondem aos pontos de casamento
            # no condutor 'm'. Estes são os índices de 'm * NF' a '(m+1)*NF - 1'.
            
            # Soma dos elementos relevantes de T, conforme a equação
            sum_of_T_elements = np.sum(T_matrix[(n * NF), (m * NF):((m + 1) * NF)])

            # Calcula o elemento c_nm da matriz de capacitância
            C_matrix[n, m] = 2 * np.pi * radii[n] * sum_of_T_elements

    return C_matrix

def maxwellian_capacitance_matrix(generalized_capacitance_matrix):
    """
    Calcula a matriz de capacitância física (n x n) a partir da matriz de
    capacitância generalizada ((n+1) x (n+1)), seguindo a Eq. 5.21 de Clayton Paul.

    A fórmula implementada é:
    C_ij = c_ij - ( (soma da linha i de c) * (soma da coluna j de c) ) / (soma total de c)

    Onde 'c' é a matriz generalizada e 'C' é a matriz física resultante.
    Assume-se que o condutor de índice 0 da matriz generalizada é o de referência
    e está sendo eliminado.

    Args:
        matriz_generalizada (np.ndarray): A matriz de capacitância generalizada
                                        simétrica de ordem (n+1) x (n+1).

    Returns:
        np.ndarray: A matriz de capacitância física de ordem n x n.
        
    Raises:
        ValueError: Se a matriz de entrada não for quadrada ou se a soma de
                    seus elementos for zero.
    """

    gc = generalized_capacitance_matrix

    # --- Validação da entrada com assert ---
    assert isinstance(gc, np.ndarray), "A entrada deve ser um array NumPy."
    assert np.sum(gc) != 0, "A soma total dos elementos da matriz generalizada não pode ser zero."
    assert gc.ndim == 2, "A entrada deve ser uma matriz 2D (array de 2 dimensões)."
    assert gc.shape[0] == gc.shape[1], "A entrada deve ser uma matriz quadrada."
    assert gc.shape[0] >= 2, "A matriz generalizada deve ser de ordem mínima 2x2."

    # Ordem da matriz generalizada (N = n+1)
    N = gc.shape[0]

    # 2. Numerador: Soma de cada linha e de cada coluna
    # Para uma matriz simétrica, as somas das linhas e colunas são iguais.
    row_sum = np.sum(gc, axis=1)     # axis=1 soma ao longo das colunas
    column_sum = np.sum(gc, axis=0)  # axis=0 soma ao longo das linhas

    # assert np.equal(row_sum, column_sum).all(), "As somas das linhas e colunas devem ser iguais."

    # Inicializa a matriz de capacitância física n x n com zeros
    matrix_c = np.zeros((N - 1, N - 1), dtype=gc.dtype)

    # Itera sobre os índices da matriz física (de 1 a n na matriz original)
    # Condutor de índice 0 é o de referência e não é incluído na matriz física
    for i in range(1, N):
        for j in range(1, N):
            c_ij = gc[i, j]
            row_i_sum = row_sum[i]
            column_j_sum = column_sum[j]
            
            # Eq. 5.21 [2]
            matrix_c[i - 1, j - 1] = c_ij - (row_i_sum * column_j_sum) / np.sum(gc)

    return matrix_c

def bifilar_mom(R, D, NF, plot_data=True):
    """
    Calcula a densidade de carga e a capacitância em uma linha bifilar usando o MoM
    e a solução analítica exata para comparação.
    """

    collocation_data = calculate_collocation_points(R, D, NF)
    observation_points = np.vstack(collocation_data['cartesian']['observation'])

    N = 2 * NF
    D_matrix = np.zeros((N, N))
    V_vector = np.zeros(N)
    centers = [np.array([-D / 2.0, 0.0]), np.array([D / 2.0, 0.0])]

    for i in range(N):
        field_wire_idx = i // NF
        V_vector[i] = 1.0 if field_wire_idx == 0 else -1.0
        field_point = observation_points[i]

        for j in range(N):
            source_wire_idx = j // NF
            source_harmonic_order = j % NF

            field_angle_rad = np.arctan2(field_point[1] - centers[field_wire_idx][1],
                                         field_point[0] - centers[field_wire_idx][0])

            if source_wire_idx == field_wire_idx:
                if source_harmonic_order == 0:
                    D_matrix[i, j] = (-R / epsilon_0) * np.log(R)
                else:
                    k = source_harmonic_order
                    D_matrix[i, j] = (R / (2 * k * epsilon_0)) * np.cos(k * field_angle_rad)
            else:
                dist_vec = field_point - centers[source_wire_idx]
                dist = np.linalg.norm(dist_vec)
                if source_harmonic_order == 0:
                    D_matrix[i, j] = (-R / epsilon_0) * np.log(dist)
                else:
                    k = source_harmonic_order
                    phi = np.arctan2(dist_vec[1], dist_vec[0])
                    D_matrix[i, j] = (R / (2 * k * epsilon_0)) * ((R / dist)**k) * np.cos(k * phi)

    sigma_coeffs = np.linalg.solve(D_matrix, V_vector)

    delta_v = 2.0
    C_exact = (np.pi * epsilon_0) / np.arccosh(D / R / 2.0)
    theta_plot = np.linspace(0, 2 * np.pi, 2 * 360)
    numerator = (D**2 / (4 * R**2)) - 1
    denominator = (D / R) - 2 * np.cos(theta_plot)
    charge_density_exact = (C_exact * delta_v / R) * (numerator / denominator)

    charge_density_mom = np.zeros_like(theta_plot)
    for k in range(NF):
        charge_density_mom += sigma_coeffs[k] * np.cos(k * theta_plot)

    if plot_data:
        plot_charge_density(np.rad2deg(theta_plot), charge_density_mom, charge_density_exact, D/R, NF)
        plot_collocation_points(collocation_data, R, D)

    return {
        'NF': NF,
        'sigma_coeffs': sigma_coeffs,
        'D_matrix': D_matrix,
        'C_exact': C_exact,
    }

def convergence_rates(R, D):
    # Iterando sobre o intervalo de NF
    mom_results = [bifilar_mom(R, D, nf, plot_data=False) for nf in range(1, 21)]

    # Separar dados em par e ímpar para plotar com marcadores diferentes
    C_FACTOR = 1e12  # Fator de conversão para pF/m
    c_exact = (np.pi * epsilon_0) / np.arccosh(D / R / 2.0)

    nf_odd, nf_even, cap_odd, cap_even = [], [], [], []

    for mom in mom_results:
        general_cap_matrix = generalized_capacitance(mom['D_matrix'], [R, R], mom['NF'])
        mom_capacitance = maxwellian_capacitance_matrix(general_cap_matrix).item()

        if mom['NF'] % 2 != 0:
            nf_odd.append(mom['NF'])
            cap_odd.append(mom_capacitance * C_FACTOR)
        else:
            nf_even.append(mom['NF'])
            cap_even.append(mom_capacitance * C_FACTOR)
    
    # Triângulos (Δ) para valores ímpares e Asteriscos (*) para valores pares
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axhline(y=c_exact*C_FACTOR, color='k', linestyle='-', label=f'Exactly Value = {c_exact*C_FACTOR:.2f} pF/m')
    ax.plot(nf_odd, cap_odd, linestyle='none', marker='^', markersize=8,  fillstyle='none', markeredgecolor='black', label='Odd NF')
    ax.plot(nf_even, cap_even, linestyle='none', marker='*', markersize=9, color='black', label='Even NF')
    ax.set_title(f'Fourier Series approximation for D/R = {D_R_ratio}')
    ax.set_xlabel('NF - Fourier coefficients Number for each wire')
    ax.set_ylabel('Capacitance (pF/m)')
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 100)
    ax.grid(False)
    ax.legend()

if __name__ == "__main__":
    os.system('cls' if os.name == 'nt' else 'clear')

    NF_order = 2
    D_R_ratio = 2.1
    R_val = 0.01
    D_val = D_R_ratio * R_val
    C_exact = (np.pi * epsilon_0) / np.arccosh(D_val / R_val / 2.0)

    mom_data = bifilar_mom(R_val, D_val, NF_order)
    general_cap_matrix = generalized_capacitance(mom_data['D_matrix'], [R_val, R_val], NF_order)
    maxwellian_cap_matrix = maxwellian_capacitance_matrix(general_cap_matrix).item()
    convergence_rates(R_val, D_val)

    print(f"\n--- Two Bare Wires Problem ---")
    print(f"--- D/R Ratio: {D_R_ratio} and NF={NF_order} -----")    
    print(f"\nSigma Coefficients (Dim.: {mom_data['sigma_coeffs'].shape}):\n{mom_data['sigma_coeffs']}")
    print(f"\nD matrix (Dim.: {mom_data['D_matrix'].shape}):\n{mom_data['D_matrix']}")
    print(f"\nT matrix (Dim.: {np.linalg.inv(mom_data['D_matrix']).shape}):\n{np.linalg.inv(mom_data['D_matrix'])}")
    print(f"\nMoM Generalized Capacitance (F/m): \n{general_cap_matrix}")
    print(f"\nCapacitância Exata: {C_exact * 1E12:.4f} pF/m")
    print(f"Maxwellian Capacitance (MoM): {maxwellian_cap_matrix * 1E12:.4f} pF/m.")
    plt.show()