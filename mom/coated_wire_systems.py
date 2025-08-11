import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from scipy.constants import epsilon_0

from mtl_data.mtl import MulticonductorTransmissionLine as MTL
from mtl_data.utils import *


class TwoCoatedWireSystem(MTL):
    """
    Calcula a capacitância e distribuição de carga para sistemas de fios nus
    usando o Método dos Momentos (MoM) com expansão em séries harmônicas.

    Esta classe herda de MulticonductorTransmissionLine (MTL) e a especializa
    para o caso de dois fios, derivando seus parâmetros de uma configuração 'mtl'.

    Executa a simulação completa do MoM, preenchendo todos os atributos de resultado.
    A construção da matriz D agora inclui os termos de expansão constante, cossenoidal e senoidal,
    conforme as expressões (20a), (20b) e (20c) de Clements (1975).    

    Nesta classe, a ordem máxima da harmônica é definida por 'k',
    enquanto NF (número de coeficientes) é derivado como 2*k + 1.
    """
    def __init__(self, mtl: dict):
        super().__init__(mtl)
        assert isinstance(mtl, dict), "O parâmetro mtl deve ser um dicionário com a configuração da linha."
        assert len(self.surfaces) > 1, "O sistema deve ter mais de dois condutores."

        # Raio do condutor 'p' (primeiro condutor)
        self.R = self.surfaces[0]['radius']
        self.D = self.D_pq[0, 1]
        self.DR_ratio = self.D / self.R

        assert self.DR_ratio > 2, "A razão D/R deve ser maior que 2 para garantir a convergência da solução."

        # Atributos de resultado
        self.collocation_data = None
        self.D_matrix = None
        self.T_matrix = None
        self.V_vector = None
        self.sigma_coeffs = None
        self.C_generalized = None
        self.C_maxwellian = None
        self.C_exact_bare_wires = None

    def _calculate_collocation_points(self):
        """
        Calcula e armazena os pontos de colocação, classificando-os em um dicionário
        aninhado pela 'tag' do condutor e pelo tipo de superfície ('conductor', 'sheath').
        """
        # Inicializa o dicionário principal que será o atributo da classe.
        self.collocation_data = {}

        # Os ângulos são os mesmos para todas as superfícies, pois NF é constante.
        source_angles = np.linspace(0, 2 * np.pi, self.NF, endpoint=False) + (np.pi / 2)
        field_angles = source_angles + np.pi / self.NF

        # Itera sobre cada superfície definida na classe base MTL.
        for surface in self.surfaces:
            tag = surface['tag']
            surface_type = surface['type']
            center = np.array(surface['center_point'])
            radius = surface['radius']

            # Cria o dicionário para a 'tag' do condutor, se ainda não existir.
            if tag not in self.collocation_data:
                self.collocation_data[tag] = {}

            # Calcula as coordenadas cartesianas para os pontos de fonte e observação.
            source_points = center + radius * np.array([np.cos(source_angles), np.sin(source_angles)]).T
            field_points = center + radius * np.array([np.cos(field_angles), np.sin(field_angles)]).T
            
            # Preenche o dicionário para a superfície específica com seus dados.
            self.collocation_data[tag][surface_type] = {
                'source': {
                    'cartesian': source_points,
                    'angles_rad': source_angles,
                    'angles_deg': np.degrees(source_angles)
                },
                'observation': {
                    'cartesian': field_points,
                    'angles_rad': field_angles,
                    'angles_deg': np.degrees(field_angles)
                }
            }

    def _calculate_generalized_capacitance(self):
        """
        Calcula a matriz de capacitância generalizada C a partir da matriz T (D^-1).
        """
        self.T_matrix = np.linalg.inv(self.D_matrix)
        C_matrix = np.zeros((2, 2))

        for n in range(2):  # Índice do condutor da carga
            for m in range(2):  # Índice do condutor do potencial
                sum_of_T_elements = np.sum(self.T_matrix[(n * self.NF), (m * self.NF):((m + 1) * self.NF)])
                C_matrix[n, m] = 2 * np.pi * self.R * sum_of_T_elements
        
        self.C_generalized = C_matrix

    def _calculate_maxwellian_capacitance(self):
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
        gc = self.C_generalized

        # --- Validações ---
        assert isinstance(gc, np.ndarray), "A matriz de capacitância generalizada deve ser um array NumPy."
        assert gc.ndim == 2 and gc.shape[0] == gc.shape[1], "A matriz de capacitância generalizada deve ser quadrada."
        num_conductors = gc.shape[0]
        assert num_conductors > 1, "O cálculo da capacitância Maxwelliana requer pelo menos 2 condutores."
        assert hasattr(self, 'idx_ref'), "O atributo 'idx_ref' (índice do condutor de referência) não foi encontrado."
        assert 0 <= self.idx_ref < num_conductors, f"O índice de referência self.idx_ref ({self.idx_ref}) está fora do intervalo válido [0, {num_conductors-1}]."

        # --- Etapa 1: Calcular a matriz Maxwelliana completa (NxN) ---
        total_sum = np.sum(gc)

        # Evita a divisão por zero
        if np.abs(total_sum) < 1e-15:
            raise ValueError("A soma dos elementos da matriz de capacitância generalizada é zero, resultando em divisão por zero.")

        row_sums = np.sum(gc, axis=1)
        col_sums = np.sum(gc, axis=0)

        correction_matrix = np.outer(row_sums, col_sums) / total_sum
        C_full = gc - correction_matrix

        # --- Etapa 2: Reduzir a matriz para (N-1)x(N-1) ---
        # Usa np.delete para remover a linha (axis=0) e a coluna (axis=1)
        # correspondentes ao índice do condutor de referência `self.idx_ref`.
        self.C_maxwellian = np.delete(np.delete(C_full, self.idx_ref, axis=0), self.idx_ref, axis=1)

    def run_simulation(self):
        """
        Executa a simulação completa do MoM, implementando a física para
        as fronteiras condutoras e dielétricas.

        Esta versão revisada distingue entre pontos de observação internos e
        externos a uma fronteira de fonte, implementando as fórmulas de potencial
        das Tabelas II.a e II.b de Clements (1975).

        NOTA: A condição de contorno do vetor deslocamento elétrico na
        superfície da bainha ainda precisa ser implementada. Esta versão calcula
        o potencial em todas as fronteiras.
        """
        self._calculate_collocation_points()

        # 1. Separar as superfícies por tipo para garantir a ordem de bloco correta.
        conductor_surfaces = [s for s in self.surfaces if s['type'] == 'conductor']
        sheath_surfaces = [s for s in self.surfaces if s['type'] == 'sheath']
        ordered_surfaces = conductor_surfaces + sheath_surfaces

        # 1. Preparar os índices e vetores do sistema
        nfs_per_surface = [2 * surface['fourier_order'] + 1 for surface in ordered_surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        self.D_matrix = np.zeros((self.N, self.N))
        self.V_vector = np.zeros(self.N)

        # 2. Montar a Matriz [D] e o Vetor [V]
        # Loop sobre as superfícies de OBSERVAÇÃO p (linhas da matriz)
        for p, field_surface in enumerate(ordered_surfaces):
            tag_p = field_surface['tag']
            type_p = field_surface['type']
            nf_p = nfs_per_surface[p]
            offset_p = offsets[p]

            # Obtém os pontos de observação para a superfície p
            observation_points = self.collocation_data[tag_p][type_p]['observation']['cartesian']

            # Preenche o vetor de potencial V para o bloco de linhas da superfície p
            if type_p == 'conductor':
                self.V_vector[offset_p : offset_p + nf_p] = self.mtl[tag_p]['potential_to_infinity']
            
            # A condição de fronteira na bainha dielétrica resulta em 0 no lado direito [cite: 222]
            elif type_p == 'sheath':
                self.V_vector[offset_p : offset_p + nf_p] = 0.0

            # Loop sobre as superfícies de FONTE q (colunas da matriz)
            for q, source_surface in enumerate(ordered_surfaces):
                tag_q = source_surface['tag']
                type_q = source_surface['type']
                center_q = np.array(source_surface['center_point'])
                nf_q = nfs_per_surface[q]
                offset_q = offsets[q]

                # Obtém os pontos de fonte para a superfície q
                source_points = self.collocation_data[tag_q][type_q]['source']['cartesian']

                # Loop sobre cada ponto de observação 'm' na superfície 'p'
                for m in range(nf_p):
                    row_idx = offset_p + m
                    
                    # Ângulo e vetor de observação 'i' relativo ao centro da superfície FONTE 'q'
                    rho_i_vector = observation_points[m] - center_q
                    rho_i = np.linalg.norm(rho_i_vector)
                    theta_i = np.arctan2(rho_i_vector[1], rho_i_vector[0])

                    # Loop sobre cada função de base 'n' na superfície 'q'
                    for n in range(nf_q):
                        col_idx = offset_q + n
                        source_harmonic_idx = n
                        
                        # Trigonometric Term at observation point 
                        is_cosine_term = (source_harmonic_idx % 2 != 0)
                        k = (source_harmonic_idx + 1) // 2 if is_cosine_term else source_harmonic_idx // 2                  
                        harmonic_term = np.cos(k * theta_i) if is_cosine_term else np.sin(k * theta_i)                                    
                        
                        # Vetor fonte 'rho_b' relativo ao centro da superfície FONTE 'q'
                        rho_b = np.linalg.norm(source_points[m] - center_q)

                        # ========================================================================
                        # ==== INÍCIO DA LÓGICA DE CÁLCULO DO ELEMENTO DA MATRIZ D ===============
                        # ========================================================================
                        is_observer_inside = (rho_i < rho_b) and not np.isclose(rho_i, rho_b)
                        
                        # === BLOCO 1: CÁLCULO DE POTENCIAL (φ) ==================================
                        # === Aplica a condição de contorno V = Vm nas superfícies condutoras. ===

                        if type_p == 'conductor':  
                            # --- TABELA II.b: rho_i < rho_b (Interação para Observador DENTRO da fronteira da fonte) --- 
                            if is_observer_inside:
                                pass
                            
                            # --- TABELA II.a: rho_i >= rho_b (Interação para Observador FORA ou SOBRE a fronteira da fonte) --- 
                            else:                         
                                if source_harmonic_idx == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = - rho_b * np.log(rho_i) / epsilon_0                                
                                
                                else: # Harmonic Terms (k>0)
                                    self.D_matrix[row_idx, col_idx] = rho_b**(k+1) / (2 * k * epsilon_0 * rho_i**k) * harmonic_term

                        # === BLOCO 2: CONDIÇÃO DE CONTORNO DO VETOR DESLOCAMENTO (εE) =======================
                        # === Aplica (1 - εr) * Er = 0 nas superfícies da bainha. ============================
                        
                        elif type_p == 'sheath':
                            er = field_surface['relative_permittivity']
                            if source_harmonic_idx == 0: # Constant Term (k=0)
                                self.D_matrix[row_idx, col_idx] = (1 - er) * (rho_b / rho_i)

                            else: # Harmonic Terms (k>0)
                                self.D_matrix[row_idx, col_idx] = 0.5 * (1 - er) * (rho_b / rho_i)**(k + 1) * harmonic_term

                        # ===============================================================
                        # ==== FIM DA LÓGICA DE CÁLCULO DO ELEMENTO DA MATRIZ D =========
                        # ===============================================================

        # 3. Resolver o sistema e obter os resultados
        self.sigma_coeffs = np.linalg.solve(self.D_matrix, self.V_vector)
        self.C_exact_bare_wires = (np.pi * epsilon_0) / np.arccosh(self.DR_ratio / 2.0)
        self._calculate_generalized_capacitance()
        self._calculate_maxwellian_capacitance()

    def print_results(self):
        """Imprime um resumo dos resultados da simulação."""
        if self.C_maxwellian is None:
            print("Executando simulação primeiro...")
            self.run_simulation()

        print(f"\nSurfaces Dim: {len(self.surfaces)}.")
        print(f"\nD Matrix Shape: {self.D_matrix.shape}.")
        if self.NF < 4: 
            matrix_viewer(self.D_matrix, "D Matrix")
            print(f"\nSurfaces (len: {len(self.surfaces)}): \n{self.surfaces}")
            print(f"\nSigma Coefficients (Shape: {self.sigma_coeffs.shape}): \n{self.sigma_coeffs}")
        
        matrix_viewer(self.C_generalized, "MoM Generalized Capacitance Matrix (F/m)")
        matrix_viewer(self.C_maxwellian, "Maxwellian Bifilar Capacitance (MoM) (F/m)")

    def plot_collocation_points(self):
        """
        Gera um gráfico interativo dos pontos de colocação usando Plotly,
        refletindo a nova estrutura de dicionário de self.collocation_data.
        """
        if self.collocation_data is None:
            self._calculate_collocation_points()

        # 1. Preparar os dados para o Plotly a partir da nova estrutura aninhada
        plot_data = []
        # Itera sobre cada 'tag' de condutor no dicionário (ex: 0, 1)
        for tag, conductor_surfaces in self.collocation_data.items():
            # Itera sobre cada superfície desse condutor (ex: 'conductor', 'sheath')
            for surface_type, surface_data in conductor_surfaces.items():
                
                # Busca o raio correspondente na lista self.surfaces,
                # pois ele não está em self.collocation_data.
                matching_surface = next(s for s in self.surfaces if s['tag'] == tag and s['type'] == surface_type)
                radius = matching_surface['radius']

                # Extrai e adiciona os dados de pontos de fonte
                for i, pt in enumerate(surface_data['source']['cartesian']):
                    plot_data.append({
                        'x': pt[0], 'y': pt[1],
                        'type': 'Source',
                        'tag': tag,
                        'surface': surface_type,
                        'radius': radius,
                        'angle_rad': surface_data['source']['angles_rad'][i]
                    })

                # Extrai e adiciona os dados de pontos de observação
                for i, pt in enumerate(surface_data['observation']['cartesian']):
                    plot_data.append({
                        'x': pt[0], 'y': pt[1],
                        'type': 'Observation',
                        'tag': tag,
                        'surface': surface_type,
                        'radius': radius,
                        'angle_rad': surface_data['observation']['angles_rad'][i]
                    })

        df = pd.DataFrame(plot_data)
        fig = go.Figure()

        # 3. Adicionar as formas dos círculos (esta parte não muda, pois já itera sobre self.surfaces)
        for surface in self.surfaces:
            fig.add_shape(type="circle",
                        xref="x", yref="y",
                        x0=surface['center_point'][0] - surface['radius'], y0=surface['center_point'][1] - surface['radius'],
                        x1=surface['center_point'][0] + surface['radius'], y1=surface['center_point'][1] + surface['radius'],
                        line_color="Black", fillcolor="LightGray", opacity=0.7)

        # 4. Adicionar os pontos de colocação a partir do DataFrame
        for pt_type, color, symbol in [('Source', 'blue', 'circle'), ('Observation', 'red', 'x-thin')]:
            df_subset = df[df['type'] == pt_type]
            fig.add_trace(go.Scatter(
                x=df_subset['x'], y=df_subset['y'],
                mode='markers',
                marker=dict(color=color, symbol=symbol, size=8, line=dict(width=1, color='DarkSlateGrey')),
                name=pt_type,
                # Atualiza o customdata e o hovertemplate para exibir as novas informações
                customdata=df_subset[['tag', 'surface', 'radius', 'angle_rad']],
                hovertemplate=(
                    f"<b>{pt_type}</b><br>"
                    "Condutor (tag): %{customdata[0]}<br>"
                    "Superfície: %{customdata[1]}<br>"
                    "Coord X: %{x:.4f} m<br>"
                    "Coord Y: %{y:.4f} m<br>"
                    "Ângulo: %{customdata[3]:.3f} rad<br>"
                    "Raio: %{customdata[2]:.4f} m"
                    "<extra></extra>"
                )
            ))

        # 5. Configurar o layout do gráfico (não muda)
        fig.update_layout(
            title='Mapa Interativo de Pontos de Colocação',
            xaxis_title='Coordenada X (m)',
            yaxis_title='Coordenada Y (m)',
            yaxis_scaleanchor="x",
            yaxis_scaleratio=1,
            legend_title_text='Tipo de Ponto',
            template='plotly_white'
        )
        fig.show()

    def plot_charge_density(self, tag_to_plot=1):
        """
        Plota a densidade de carga para um condutor específico, alinhando
        dinamicamente a solução exata com a geometria real do sistema.

        Args:
            tag_to_plot (int): A 'tag' do condutor para o qual a densidade de
                            carga será plotada.
        """
        if self.sigma_coeffs is None:
            self.run_simulation()

        # 1. Obter dados do condutor a ser plotado e de seu par
        all_tags = list(self.mtl.keys())
        if len(all_tags) != 2:
            print("Erro: plot_charge_density foi projetado para sistemas de 2 condutores.")
            return
        other_tag = next(tag for tag in all_tags if tag != tag_to_plot)

        center_plot = np.array(self.mtl[tag_to_plot]['center_point'])
        center_other = np.array(self.mtl[other_tag]['center_point'])

        conductor_surface = next(s for s in self.surfaces if s['tag'] == tag_to_plot and s['type'] == 'conductor')
        R = conductor_surface['radius']
        D = np.linalg.norm(center_plot - center_other)
        DR_ratio = D / R

        # 2. Calcular a Solução Analítica com Alinhamento e Sinal Corretos
        C_exact = (np.pi * epsilon_0) / np.arccosh(DR_ratio / 2.0)
        theta_plot = np.linspace(0, 2 * np.pi, 360)

        vec_to_other = center_other - center_plot
        angle_of_max_charge = np.arctan2(vec_to_other[1], vec_to_other[0])
        denominator = DR_ratio - 2 * np.cos(theta_plot - angle_of_max_charge)

        delta_v = self.mtl[tag_to_plot]['potential_to_infinity'] - self.mtl[other_tag]['potential_to_infinity']
        numerator = (DR_ratio**2 / 4) - 1
        
        # CORREÇÃO FINAL: Remover o abs() para preservar o sinal da carga
        charge_density_exact = (C_exact * delta_v / R) * (numerator / denominator)

        # 3. Reconstruir a Solução MoM para o Condutor Correto
        nfs_per_surface = [2 * s['fourier_order'] + 1 for s in self.surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)
        
        surface_index = next(i for i, s in enumerate(self.surfaces) if s['tag'] == tag_to_plot and s['type'] == 'conductor')
        
        offset = offsets[surface_index]
        nf = nfs_per_surface[surface_index]
        coeffs_to_plot = self.sigma_coeffs[offset : offset + nf]

        charge_density_mom = np.zeros_like(theta_plot)
        charge_density_mom += coeffs_to_plot[0]
        max_k = (nf - 1) // 2
        for k in range(1, max_k + 1):
            cos_coeff = coeffs_to_plot[2 * k - 1]
            sin_coeff = coeffs_to_plot[2 * k]
            charge_density_mom += cos_coeff * np.cos(k * theta_plot) + sin_coeff * np.sin(k * theta_plot)

        # 4. Geração do Gráfico
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(np.rad2deg(theta_plot), charge_density_exact, 'r-', label='Solução Exata')
        ax.plot(np.rad2deg(theta_plot), charge_density_mom, 'k-.', label=f'MoM (tag={tag_to_plot})')
        ax.set_title(f'Distribuição de Carga (Condutor {tag_to_plot}) com D/R = {DR_ratio:.2f}')
        ax.set_xlabel('Ângulo (Graus)'); ax.set_ylabel('Densidade de Carga (C/m²)')
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.set_xticks(np.arange(0, 361, 90)); ax.set_xlim(0, 360)
        ax.legend()
        plt.tight_layout()
             
    def plot_harmonic_coefficients(self):
        """
        Reproduz e expande a Figura 4(c) de Clements (1975), mostrando a magnitude
        de todos os coeficientes da série harmônica com indexação ajustada.

        O gráfico mostra a razão entre a magnitude de cada coeficiente harmônico
        e a magnitude do coeficiente constante. A plotagem segue a convenção:
        - j=1: Termo Constante
        - j=2, 4, 6,...: Coeficientes Cossenoidais
        - j=3, 5, 7,...: Coeficientes Senoidais
        """
        if self.sigma_coeffs is None:
            print("Executando a simulação para obter os coeficientes...")
            self.run_simulation()

        # Isola os coeficientes do primeiro condutor
        coeffs_condutor1 = self.sigma_coeffs[:self.NF]

        # O coeficiente constante (alpha_n1) está no índice 0 do código
        constant_term = coeffs_condutor1[0]
        if np.abs(constant_term) < 1e-15: # Evita divisão por zero
            print("Coeficiente constante é próximo de zero. Não é possível normalizar.")
            return

        # --- Mapeamento de Índices para Plotagem ---
        # j=1: Termo Constante. Sua razão normalizada é 1.0.
        plot_j_const = [1]
        ratio_const = [1.0]

        # j=2, 4, 6,...: Coeficientes Cossenoidais (índices 1, 3, 5,... no código)
        code_indices_cos = np.arange(1, self.NF, 2)
        plot_j_cos = code_indices_cos + 1 # Mapeia [1, 3, 5] para [2, 4, 6]
        coeffs_cos = coeffs_condutor1[code_indices_cos]
        ratio_cos = np.abs(coeffs_cos) / np.abs(constant_term)

        # j=3, 5, 7,...: Coeficientes Senoidais (índices 2, 4, 6,... no código)
        code_indices_sin = np.arange(2, self.NF, 2)
        plot_j_sin = code_indices_sin + 1 # Mapeia [2, 4, 6] para [3, 5, 7]
        coeffs_sin = coeffs_condutor1[code_indices_sin]
        ratio_sin = np.abs(coeffs_sin) / np.abs(constant_term)

        # --- Geração do Gráfico ---
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(12, 8))

        # Plota cada série com um marcador distinto
        ax.plot(plot_j_const, ratio_const, marker='s', markersize=6, linestyle='none',
                color='blue', label='Termo Constante (j=1)')

        ax.plot(plot_j_cos, ratio_cos, marker='^', markersize=6, linestyle='none',
                fillstyle='none', markeredgecolor='black', label='Coeficientes Cossenoidais (j=2, 4, ...)')

        ax.plot(plot_j_sin, ratio_sin, marker='o', markersize=6, linestyle='none',
                color='red', label='Coeficientes Senoidais (j=3, 5, ...)')

        ax.set_title(f'Magnitude Normalizada dos Coeficientes Harmônicos (d/a = {self.DR_ratio:.1f})', fontsize=14)
        ax.set_ylabel(r'$|\alpha_{nj} / \alpha_{n1}|$', fontsize=12)
        ax.set_xlabel('Índice do Coeficiente (j)', fontsize=12)
        ax.legend()
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=-0.05)
        ax.set_xticks(np.arange(1, self.NF + 1))
        ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        