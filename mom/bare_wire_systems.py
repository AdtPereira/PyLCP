import numpy as np
import scipy.integrate as spi
import scipy.constants as spc

from utils.case_utils import *
from mtl_main.source import MulticonductorTransmissionLine

@staticmethod
def galerkin_integrand(phi_p, p_idx, q_idx, test_func_idx, basis_func_idx, self_instance):
    """
    Método estático para calcular o integrando do Método de Galerkin.
    
    Calcula f_pa(phi_p) * g_qb(phi_p), onde 'f' é a função de teste e 'g' é o potencial
    da função de base.

    Args:
        phi_p (float): Ângulo no condutor de observação 'p' (variável de integração).
        p_idx (int): Índice da superfície de observação.
        q_idx (int): Índice da superfície da fonte.
        test_func_idx (int): Índice da função de teste 'a' no condutor 'p'.
        basis_func_idx (int): Índice da função de base 'b' no condutor 'q'.
        self_instance (object): A instância da classe para acessar os dados do modelo.

    Returns:
        float: O valor do integrando.
    """
    # Obter dados das superfícies
    field_surface = self_instance.model.surfaces[p_idx]
    source_surface = self_instance.model.surfaces[q_idx]
    epsilon = self_instance.model.epsilon_out[source_surface['tag']]

    # --- 1. Calcular o valor da função de teste f_pa(phi_p) ---
    is_test_cos = (test_func_idx % 2 != 0)
    k_test = (test_func_idx + 1) // 2 if is_test_cos else test_func_idx // 2

    if k_test == 0:
        f_pa = 1.0
    elif is_test_cos:
        f_pa = np.cos(k_test * phi_p)
    else:
        f_pa = np.sin(k_test * phi_p)

    # --- 2. Calcular o potencial g_qb(phi_p) ---
    # Coordenadas do ponto de observação no condutor 'p'
    obs_point = np.array(field_surface['center_point']) + \
                field_surface['radius'] * np.array([np.cos(phi_p), np.sin(phi_p)])

    # Vetor do centro da fonte 'q' ao ponto de observação
    rho_b_vector = obs_point - np.array(source_surface['center_point'])
    rho_b = np.linalg.norm(rho_b_vector)
    theta_b = np.arctan2(rho_b_vector[1], rho_b_vector[0])

    is_basis_cos = (basis_func_idx % 2 != 0)
    k_basis = (basis_func_idx + 1) // 2 if is_basis_cos else basis_func_idx // 2
    
    g_qb = 0.0
    if k_basis == 0:  # Termo constante da fonte
        if p_idx == q_idx:
            g_qb = (-source_surface['radius'] / epsilon) * np.log(field_surface['radius'])
        else:
            g_qb = (-source_surface['radius'] / epsilon) * np.log(rho_b)
    else:  # Termos harmônicos da fonte
        if p_idx == q_idx:
            term = np.cos(k_basis * phi_p) if is_basis_cos else np.sin(k_basis * phi_p)
            g_qb = (source_surface['radius'] / (2 * k_basis * epsilon)) * term
        else:
            term = np.cos(k_basis * theta_b) if is_basis_cos else np.sin(k_basis * theta_b)
            g_qb = (source_surface['radius'] / (2 * k_basis * epsilon)) * ((source_surface['radius'] / rho_b)**k_basis) * term

    return f_pa * g_qb

@staticmethod
def maxwellian_capacitance(model, mom_data):
    """
    Calcula a matriz de capacitância física (Maxwelliana) de dimensão (N-1)x(N-1)
    a partir da matriz de capacitância generalizada de dimensão NxN.

    Este processo ocorre em duas etapas:
    1.  Primeiro, uma matriz Maxwelliana completa (NxN) é calculada usando a
        Equação 5.21, que é dada por:
        C_completa_ij = c_ij - (soma_linha_i * soma_coluna_j) / soma_total
    2.  Em seguida, a matriz é reduzida para (N-1)x(N-1) ao remover a linha e a
        coluna correspondentes ao condutor de referência, cujo índice é
        especificado pelo atributo da classe `self.idx_ref`.
    """
    cgen = mom_data['generalized_capacitance']
    idx_ref = model.mtl_idx_ref

    # --- Validações ---
    assert isinstance(cgen, np.ndarray), "A matriz de capacitância generalizada deve ser um array NumPy."
    assert cgen.ndim == 2 and cgen.shape[0] == cgen.shape[1], "A matriz de capacitância generalizada deve ser quadrada."
    assert cgen.shape[0] > 1, "O cálculo da capacitância Maxwelliana requer pelo menos 2 condutores."
    assert 0 <= idx_ref < cgen.shape[0], f"O índice de referência self.idx_ref ({idx_ref}) está fora do intervalo válido [0, {cgen.shape[0]-1}]."

    # --- Etapa 1: Calcular a matriz Maxwelliana completa (NxN) ---
    total_sum = np.sum(cgen)

    # Evita a divisão por zero
    assert np.abs(total_sum) > 1e-15, "A soma dos elementos da matriz de capacitância generalizada é zero, resultando em divisão por zero."

    correction_matrix = np.outer(np.sum(cgen, axis=1), np.sum(cgen, axis=0)) / total_sum
    C_full = cgen - correction_matrix

    # --- Etapa 2: Reduzir a matriz para (N-1)x(N-1) ---
    # Usa np.delete para remover a linha (axis=0) e a coluna (axis=1)
    # correspondentes ao índice do condutor de referência `self.idx_ref`.
    return np.delete(np.delete(C_full, idx_ref, axis=0), idx_ref, axis=1)

class BareWireMoMSolver:
    """
    Calcula a capacitância e distribuição de carga para sistemas de fios nus
    usando o Método dos Momentos (MoM) com expansão em séries harmônicas.

    Esta classe utiliza uma instância de MulticonductorTransmissionLine (MTL)
    para obter os parâmetros geométricos e elétricos do sistema. Ela então
    executa a simulação completa do MoM, preenchendo seus próprios atributos 
    de resultado.

    Executa a simulação completa do MoM, preenchendo todos os atributos de 
    resultado. A construção da matriz D agora inclui os termos de expansão
    constante, cossenoidal e senoidal, conforme as expressões (20a), (20b)
    e (20c) de Clements (1975).    

    Nesta classe, a ordem máxima da harmônica é definida por 'k',
    enquanto NF (número de coeficientes) é derivado como 2*k + 1.
    """
    def __init__(self, model: MulticonductorTransmissionLine):
        # MTL Geometry Model
        self.model = model
        
        # Número de coeficientes harmônicos de Fourier por condutor
        self.NF = [2*surface['fourier_order']+1 for surface in self.model.surfaces][0]
        
        # Atributos de resultado
        self.mom_data = {'collocation': {}, 'galerkin': {}}

        self.DR_ratio = (model.D_pq[0, 1]) / (model.surfaces[0]['radius'])
        assert self.DR_ratio > 2, "A razão D/R deve ser maior que 2 para garantir a convergência da solução."
        self.C_exact_bare_wires = np.pi * spc.epsilon_0 / np.arccosh(0.5 * self.DR_ratio)

    def _collocation_points(self):
        """
        Calcula e armazena os pontos de colocação, classificando-os em um dicionário
        aninhado pela 'tag' do condutor e pelo tipo de superfície ('conductor', 'sheath').
        """
        # Inicializa o dicionário principal que será o atributo da classe.
        collocation_data = {}

        # Equação (A.4b): Ângulo de rotação para o conjunto de pontos.
        delta = np.pi / (2 * self.NF)

        # Calcula os ângulos base, que são rotacionados por delta para obter
        # os ângulos dos pontos de observação (match points).
        base_angles = np.linspace(0, 2 * np.pi, self.NF, endpoint=False)
        match_angles = base_angles + delta

        # Itera sobre cada superfície definida na classe base MTL.
        for surface in self.model.surfaces:
            if surface['tag'] not in collocation_data:
                collocation_data[surface['tag']] = {}

            # Calcula as coordenadas cartesianas para os pontos de fonte e observação.
            match_points = np.array(surface['center_point']) + surface['radius'] * np.array([np.cos(match_angles), np.sin(match_angles)]).T
            
            # Preenche o dicionário para a superfície específica com seus dados.
            collocation_data[surface['tag']][surface['type']] = {
                'observation': {
                    'cartesian': match_points,
                    'angles_rad': match_angles
                }
            }
        
        self.mom_data['collocation']['data'] = collocation_data

    def _generalized_capacitance_clements(self):
        """
        Calcula a matriz de capacitância generalizada C a partir da matriz T (D^-1).
        """
        moment_matrix = self.mom_data['collocation']['moment_matrix']
        T_matrix = np.linalg.inv(moment_matrix)
        C_matrix = np.zeros((2, 2))

        for n in range(2):      # Índice do condutor da carga
            r_i = self.model.surfaces[n]['radius']
            for m in range(2):  # Índice do condutor do potencial
                sum_of_T_elements = np.sum(T_matrix[(n * self.NF), (m * self.NF):((m + 1) * self.NF)])
                C_matrix[n, m] = 2 * np.pi * r_i * sum_of_T_elements

        self.mom_data['collocation']['generalized_capacitance'] = C_matrix

    def _generalized_capacitance_savage(self):
        """
        Calcula a matriz de capacitância generalizada C a partir da matriz T (D^-1)
        seguindo a formulação de Savage (1993) para o Método de Galerkin.

        A formulação é dada por: C_ij = (2*pi)^2 * r_i * T_ij[0,0], onde T_ij[0,0]
        é o elemento superior esquerdo da submatriz correspondente da matriz inversa T.
        """
        # 1. Inverter a matriz D para obter a matriz T
        moment_matrix = self.mom_data['galerkin']['moment_matrix']
        T_matrix = np.linalg.inv(moment_matrix)

        # 2. Obter o número de condutores (superfícies)
        num_conductors = len(self.model.surfaces)
        C_matrix = np.zeros((num_conductors, num_conductors))

        # 3. Iterar sobre cada elemento da matriz de capacitância a ser calculada
        for i in range(num_conductors):      # Índice 'i' para o condutor da carga (linha)
            for j in range(num_conductors):  # Índice 'j' para o condutor do potencial (coluna)

                # 4. Obter o raio do condutor da carga 'i'
                # Isso é mais robusto que usar self.R, pois considera raios diferentes.
                r_i = self.model.surfaces[i]['radius']

                # 5. Localizar o elemento (0,0) da submatriz T_ij
                # Este é o elemento superior esquerdo do bloco que relaciona a observação
                # no condutor 'i' com a fonte no condutor 'j'.
                T_ij_00 = T_matrix[i * self.NF, j * self.NF]

                # 6. Calcular o elemento da capacitância C_ij conforme Equação 4.37
                C_matrix[i, j] = (2 * np.pi)**2 * r_i * T_ij_00

        self.mom_data['galerkin']['generalized_capacitance'] = C_matrix

    def run_collocation_method(self):
        """
        Executa a simulação completa do MoM, montando o sistema de equações para
        todas as superfícies (condutoras e dielétricas) com base nas novas
        estruturas de dados.
        """
        self._collocation_points()
        collocation_data = self.mom_data['collocation']['data']
        moment_matrix = np.zeros((self.model.N, self.model.N))
        V_vector = np.zeros(self.model.N)

        # 1. Preparar os índices e vetores do sistema
        # Pré-calcula o número de coeficientes (NF) para cada superfície
        nfs_per_surface = [2 * surface['fourier_order'] + 1 for surface in self.model.surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        # 2. Montar a Matriz [D] e o Vetor [V]
        # Loop sobre as superfícies de OBSERVAÇÃO p (linhas da matriz)
        for p, field_surface in enumerate(self.model.surfaces):
            tag_p = field_surface['tag']
            type_p = field_surface['type']
            radius_p = field_surface['radius']

            # Obtém os pontos de observação para a superfície p
            match_points = collocation_data[tag_p][type_p]['observation']['cartesian']

            # Preenche o vetor de potencial V para o bloco de linhas da superfície p
            if type_p == 'conductor':
                V_vector[offsets[p] : offsets[p] + nfs_per_surface[p]] = self.model.mtl[tag_p]['potential_to_infinity']
            
            # A condição de fronteira na bainha dielétrica resulta em 0 no lado direito da equação
            elif type_p == 'primary_insulation':
                V_vector[offsets[p] : offsets[p] + nfs_per_surface[p]] = 0.0

            # Loop sobre as superfícies de FONTE q (colunas da matriz)
            for q, source_surface in enumerate(self.model.surfaces):
                radius_q = source_surface['radius']
                epsilon = self.model.epsilon_out[source_surface['tag']]

                # Loop sobre cada ponto de observação m na superfície p
                for m in range(nfs_per_surface[p]):
                    row_idx = offsets[p] + m
                    
                    # Ângulo do ponto de observação relativo ao centro da sua PRÓPRIA superfície
                    rho_i_vector = match_points[m] - np.array(field_surface['center_point'])
                    theta_i = np.arctan2(rho_i_vector[1], rho_i_vector[0])

                    # Loop sobre cada função de base n na superfície q
                    for n in range(nfs_per_surface[q]):
                        col_idx = offsets[q] + n
                        
                        # Índice harmônico local da fonte
                        harmonic_idx = n
                        is_cosine_term = (harmonic_idx % 2 != 0)
                        k = (harmonic_idx + 1) // 2 if is_cosine_term else harmonic_idx // 2
                        
                        # Ângulo e vetor fonte 'b' relativo ao centro da superfície FONTE 'q'
                        rho_b_vector = match_points[m] - np.array(source_surface['center_point'])
                        rho_b = np.linalg.norm(rho_b_vector)
                        theta_b = np.arctan2(rho_b_vector[1], rho_b_vector[0])

                        # ========================================================================
                        # ==== INÍCIO DA LÓGICA DE CÁLCULO DO ELEMENTO DA MATRIZ D ===============
                        # ========================================================================

                        # === BLOCO 1: CÁLCULO DE POTENCIAL (φ) ==================================
                        # === Aplica a condição de contorno V = Vm nas superfícies condutoras. ===

                        # Auto-interação (Observador NA fronteira da fonte)
                        # Termo constante (k=0)
                        if harmonic_idx == 0:
                            if p == q:
                                moment_matrix[row_idx, col_idx] = (-radius_q / epsilon) * np.log(radius_p)
                            
                            # Interação mútua
                            else: 
                                moment_matrix[row_idx, col_idx] = (-radius_q / epsilon) * np.log(rho_b)
                        
                        # Termos harmônicos (k>0)
                        else:  
                            # Auto-interação
                            if p == q:
                                term = np.cos(k * theta_i) if is_cosine_term else np.sin(k * theta_i)
                                moment_matrix[row_idx, col_idx] = (radius_q / (2 * k * epsilon)) * term
                            
                            # Interação mútua
                            else:
                                term = np.cos(k * theta_b) if is_cosine_term else np.sin(k * theta_b)
                                moment_matrix[row_idx, col_idx] = (radius_q / (2 * k * epsilon)) * ((radius_q / rho_b)**k) * term
                        
                        # ========================================================================
                        # ==== FIM DA LÓGICA DE CÁLCULO DO ELEMENTO DA MATRIZ D ==================
                        # ========================================================================

        # 4. Armazenar os resultados no dicionário mom_data
        self.mom_data['collocation']['V_vector'] = V_vector
        self.mom_data['collocation']['moment_matrix'] = moment_matrix
        self.mom_data['collocation']['sigma_coeffs'] = np.linalg.solve(moment_matrix, V_vector)
        
        self._generalized_capacitance_clements()
        cap_matrix = maxwellian_capacitance(self.model, self.mom_data['collocation'])
        self.mom_data['collocation']['maxwellian_capacitance'] = cap_matrix

    def run_galerkin_method(self): 
        """
        Executa a simulação completa do MoM usando o Método de Galerkin.
        """
        moment_matrix = np.zeros((self.model.N, self.model.N))
        V_vector = np.zeros(self.model.N)
        
        nfs_per_surface = [2 * surface['fourier_order'] + 1 for surface in self.model.surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        # Loop sobre as superfícies de OBSERVAÇÃO p (linhas da matriz)
        for p, field_surface in enumerate(self.model.surfaces):
            # Loop sobre as superfícies de FONTE q (colunas da matriz)
            for q, source_surface in enumerate(self.model.surfaces):
                
                # Loop sobre as FUNÇÕES DE TESTE 'm' na superfície 'p'
                for m in range(nfs_per_surface[p]):
                    row_idx = offsets[p] + m
                    
                    # Loop sobre as FUNÇÕES DE BASE 'n' na superfície 'q'
                    for n in range(nfs_per_surface[q]):
                        col_idx = offsets[q] + n
                        
                        # --- Integração Numérica com scipy.integrate.quad ---
                        integral_value, _ = spi.quad(
                            galerkin_integrand, 0, 2 * np.pi,
                            args=(p, q, m, n, self),
                            limit=350
                        )
                        moment_matrix[row_idx, col_idx] = integral_value

            # Preenchimento do Vetor V conforme a formulação de Galerkin
            if field_surface['type'] == 'conductor':
                potential = self.model.mtl[field_surface['tag']]['potential_to_infinity']
                # Apenas o termo constante (m=0) da integral do lado direito é não-nulo
                V_vector[offsets[p]] = 2 * np.pi * potential
            
        # Resolver o sistema e obter os resultados
        self.mom_data['galerkin']['V_vector'] = V_vector
        self.mom_data['galerkin']['moment_matrix'] = moment_matrix
        self.mom_data['galerkin']['sigma_coeffs'] = np.linalg.solve(moment_matrix, V_vector)
        
        self._generalized_capacitance_savage()
        cap_matrix = maxwellian_capacitance(self.model, self.mom_data['galerkin'])
        self.mom_data['galerkin']['maxwellian_capacitance'] = cap_matrix

