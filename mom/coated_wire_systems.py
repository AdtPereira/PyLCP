import numpy as np
import pandas as pd
import scipy.constants as sc
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from mtl_main.utils import *
from mtl_main.source import MulticonductorTransmissionLine

class MulticonductorCoatedWireSystems:
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
    def __init__(self, model: MulticonductorTransmissionLine):
        # MTL Geometry Model
        self.model = model
        
        # Raio do condutor 'p' (primeiro condutor)
        self.R = self.model.surfaces[0]['radius']
        self.D = self.model.D_pq[0, 1]
        self.DR_ratio = self.D / self.R
        assert self.DR_ratio > 2, "A razão D/R deve ser maior que 2 para garantir a convergência da solução."
        self.C_exact_bare_wires = np.pi * sc.epsilon_0 / np.arccosh(0.5 * self.DR_ratio)

        # Superfícies de condutores e isolamento
        self.conductor_surfaces = [s for s in self.model.surfaces if s['type'] == 'conductor']
        self.insulation_surfaces = [s for s in self.model.surfaces if s['type'] == 'primary_insulation']
        self.ordered_surfaces = self.conductor_surfaces + self.insulation_surfaces

        # Atributos de resultado
        self.collocation_data = None
        self.D_matrix = None
        self.T_matrix = None
        self.V_vector = None
        self.sigma_coeffs = None
        self.C_generalized = None
        self.C_maxwellian = None

    def _calculate_collocation_points(self):
        """
        Calcula e armazena os pontos de colocação, classificando-os em um dicionário
        aninhado pela 'tag' do condutor e pelo tipo de superfície ('conductor', 'primary_insulation').
        """
        # Inicializa o dicionário principal que será o atributo da classe.
        self.collocation_data = {}

        # Equação (A.4a): Ângulo de separação entre os pontos de colocação.
        theta = 2 * np.pi / self.model.NF

        # Equação (A.4b): Ângulo de rotação para o conjunto de pontos.
        delta = np.pi / (2 * self.model.NF)

        # Calcula os ângulos base, que são rotacionados por delta para obter
        # os ângulos dos pontos de observação (match points).
        base_angles = np.linspace(0, 2 * np.pi, self.model.NF, endpoint=False)
        field_angles = base_angles + delta

        # Os pontos de fonte são posicionados na metade do caminho entre os
        # pontos de observação para garantir a estabilidade numérica.
        source_angles = field_angles - (theta / 2)

        # Itera sobre cada superfície definida na classe base MTL.
        for surface in self.model.surfaces:
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
        Calcula a matriz de capacitância generalizada C de forma robusta.

        Esta versão revisada corrige a falha da implementação anterior, garantindo
        que a matriz C seja construída corretamente independentemente dos valores ou
        da ordem das 'tags' dos condutores.
        """
        self.T_matrix = np.linalg.inv(self.D_matrix)

        num_conductors = len(self.conductor_surfaces)
        nfs_per_surface = [2 * s['fourier_order'] + 1 for s in self.ordered_surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        # 2. Criar um mapa para fácil acesso às propriedades e offsets de cada superfície
        surface_map = {}
        for i, surface in enumerate(self.ordered_surfaces):
            tag = surface['tag']
            if tag not in surface_map:
                surface_map[tag] = {}
            surface_map[tag][surface['type']] = {
                'radius': surface['radius'],
                'offset': offsets[i],
                'nf': nfs_per_surface[i]
            }
            
        # --- INÍCIO DA LÓGICA REVISADA ---
        
        # 3. Garantir uma ordem consistente para a matriz de capacitância
        ### Obter uma lista ordenada das tags dos condutores. Essencial para consistência.
        sorted_conductor_tags = sorted([s['tag'] for s in self.conductor_surfaces])
        
        ### Criar um mapa de 'tag' para o índice da matriz (0, 1, 2...).
        tag_to_idx = {tag: i for i, tag in enumerate(sorted_conductor_tags)}

        # 4. Calcular a matriz de capacitância
        C_matrix = np.zeros((num_conductors, num_conductors))

        ### Loop sobre os TAGS dos condutores, não sobre índices genéricos.
        for i_tag in sorted_conductor_tags:
            for j_tag in sorted_conductor_tags:
                
                # Obter os índices corretos da matriz a partir das tags
                row = tag_to_idx[i_tag]
                col = tag_to_idx[j_tag]
                
                # Informações do bloco de colunas do condutor 'j_tag'
                info_cond_j = surface_map[j_tag]['conductor']
                col_start_j = info_cond_j['offset']
                col_end_j = col_start_j + info_cond_j['nf']

                # Termo 1 (Eq. 5.48): Contribuição da superfície do condutor 'i_tag'.
                info_cond_i = surface_map[i_tag]['conductor']
                row_idx_cond_i = info_cond_i['offset']
                radius_cond_i = info_cond_i['radius']
                sum_bij = np.sum(self.T_matrix[row_idx_cond_i, col_start_j:col_end_j])
                term1 = 2 * np.pi * radius_cond_i * sum_bij

                # Termo 2 (Eq. 5.48): Contribuição da superfície da bainha 'i_tag'.
                term2 = 0.0
                if 'primary_insulation' in surface_map[i_tag]:
                    surface_i = surface_map[i_tag]['primary_insulation']
                    row_idx_i = surface_i['offset']
                    sum_b_prime_ij = np.sum(self.T_matrix[row_idx_i, col_start_j:col_end_j])
                    term2 = 2 * np.pi * surface_i['radius'] * sum_b_prime_ij

                ### Atribuir o valor à posição correta na matriz usando os índices mapeados.
                C_matrix[row, col] = term1 + term2

        self.C_generalized = C_matrix

    def _calculate_maxwellian_capacitance(self):
        """
        Calcula a matriz de capacitância física (Maxwelliana) de dimensão (N-1)x(N-1),
        replicando fielmente a lógica e a ordenação do código RIBBON.FOR.

        A ordenação da matriz final é baseada na sequência original dos condutores,
        simplesmente removendo a linha/coluna do condutor de referência,
        conforme implementado no código-fonte de referência.
        """
        gc = self.C_generalized
        
        # Supondo que self.model.idx_ref já foi corrigido para 0-base no __init__
        ref_idx = self.model.mtl_idx_ref 

        # --- Validações ---
        num_conductors = gc.shape[0]
        assert 0 <= ref_idx < num_conductors, f"Índice de referência ({ref_idx}) inválido."

        # --- Etapa 1: Calcular as somas necessárias, como em RIBBON.FOR ---
        total_sum = np.sum(gc)
        if np.abs(total_sum) < 1e-15:
            raise ValueError("A soma dos elementos da matriz de capacitância generalizada é próxima de zero.")
        
        row_sums = np.sum(gc, axis=1)
        col_sums = np.sum(gc, axis=0)

        # --- Etapa 2: Calcular a matriz (N-1)x(N-1) com a ordenação simples de RIBBON.FOR ---
        
        # Obter os índices originais dos condutores, exceto o de referência.
        # A ordem é a natural dos índices (0, 1, 2, ... N-1), que corresponde a I=1,N do FORTRAN.
        final_indices = [i for i in range(num_conductors) if i != ref_idx]
        
        # Inicializar a matriz final (N-1)x(N-1)
        C_maxwellian = np.zeros((num_conductors - 1, num_conductors - 1))

        # Preencher a matriz final iterando sobre os índices preservando a ordem original
        for i_new, i_orig in enumerate(final_indices):
            for j_new, j_orig in enumerate(final_indices):
                
                # Aplica a fórmula de RIBBON.FOR / Eq. (5.21)
                # Os índices i_orig e j_orig correspondem diretamente às linhas/colunas de gc, row_sums e col_sums
                correction_term = (row_sums[i_orig] * col_sums[j_orig]) / total_sum
                C_maxwellian[i_new, j_new] = gc[i_orig, j_orig] - correction_term
                
        self.C_maxwellian = C_maxwellian

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

        # 1. Preparar os índices e vetores do sistema
        nfs_per_surface = [2 * surface['fourier_order'] + 1 for surface in self.ordered_surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        self.D_matrix = np.zeros((self.model.N, self.model.N))
        self.V_vector = np.zeros(self.model.N)

        # 2. Montar a Matriz [D] e o Vetor [V]
        # Loop sobre as superfícies de OBSERVAÇÃO p (linhas da matriz)
        for p, field_surface in enumerate(self.ordered_surfaces):
            tag_p = field_surface['tag']
            type_p = field_surface['type']
            radius_p = field_surface['radius']
            center_p = np.array(field_surface['center_point'])
            nf_p = nfs_per_surface[p]
            offset_p = offsets[p]

            # Obtém os pontos de observação para a superfície p
            observation_points = self.collocation_data[tag_p][type_p]['observation']['cartesian']

            # Preenche o vetor de potencial V para o bloco de linhas da superfície p
            if type_p == 'conductor':
                self.V_vector[offset_p : offset_p + nf_p] = self.model.mtl[tag_p]['potential_to_infinity']
            
            # A condição de fronteira na bainha dielétrica resulta em 0 no lado direito [cite: 222]
            elif type_p == 'primary_insulation':
                self.V_vector[offset_p : offset_p + nf_p] = 0.0

            # Loop sobre as superfícies de FONTE q (colunas da matriz)
            for q, source_surface in enumerate(self.ordered_surfaces):
                tag_q = source_surface['tag']
                type_q = source_surface['type']
                epsilon = self.model.epsilon_out[source_surface['tag']]
                center_q = np.array(source_surface['center_point'])
                nf_q = nfs_per_surface[q]
                offset_q = offsets[q]

                # Obtém os pontos de fonte para a superfície q
                source_points = self.collocation_data[tag_q][type_q]['source']['cartesian']

                # Loop sobre cada ponto de observação 'm' na superfície 'p'
                for m in range(nf_p):
                    row_idx = offset_p + m
                    
                    # Vetor aponta do centro da superfície FONTE 'q' para o ponto de OBSERVAÇÃO 'm'.
                    rho_i_vector = observation_points[m] - center_q
                    rho_i = np.linalg.norm(rho_i_vector)
                    theta_i = np.arctan2(rho_i_vector[1], rho_i_vector[0])

                    # Vetor unitário (un_rho_i) do centro da superfície 'q' até o ponto de OBSERVAÇÃO 'm'.
                    un_rho_i = rho_i_vector / rho_i

                    # Vetor normal unitário (un_p) do centro da superfície 'p' até o ponto de OBSERVAÇÃO 'm'.
                    un_p = (observation_points[m] - center_p) / radius_p

                    # Loop sobre cada função de base 'n' na superfície 'q'
                    for n in range(nf_q):
                        col_idx = offset_q + n
                        harmonic_ord = n
                        
                        # Trigonometric Term at observation point 
                        is_cosine_term = (harmonic_ord % 2 != 0)
                        k = (harmonic_ord + 1) // 2 if is_cosine_term else harmonic_ord // 2                  
                        harmonic_term = np.cos(k * theta_i) if is_cosine_term else np.sin(k * theta_i)
                        k2epsilon = k * 2 * epsilon
                        
                        # Vetor fonte 'rho_b' relativo ao centro da superfície FONTE 'q'
                        rho_b = np.linalg.norm(source_points[m] - center_q)

                        # ====================================================================================
                        # ==== INÍCIO DA LÓGICA DE CÁLCULO DO ELEMENTO DA MATRIZ D ===========================
                        # ====================================================================================
                        is_observer_inside = (rho_i < rho_b) and not np.isclose(rho_i, rho_b)
                        
                        # === BLOCO 1: CÁLCULO DE POTENCIAL (φ) ==============================================
                        # === Aplica a condição de contorno V = Vm nas superfícies condutoras. ===============

                        if type_p == 'conductor':  
                            # --- TABELA II.b: rho_i < rho_b (Interação para Observador DENTRO da fronteira da fonte) --- 
                            if is_observer_inside:
                                if harmonic_ord == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = - rho_b * np.log(rho_b) / epsilon                                
                                
                                else: # Harmonic Terms (k>0)
                                    self.D_matrix[row_idx, col_idx] = rho_i**k / k2epsilon / rho_b**(k-1) * harmonic_term
                            
                            # --- TABELA II.a: rho_i >= rho_b (Interação para Observador FORA ou SOBRE a fronteira da fonte) --- 
                            else:                         
                                if harmonic_ord == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = - rho_b * np.log(rho_i) / epsilon                                
                                
                                else: # Harmonic Terms (k>0)
                                    self.D_matrix[row_idx, col_idx] = rho_b**(k+1) / k2epsilon / rho_i**k * harmonic_term

                        # === BLOCO 2: CONDIÇÃO DE CONTORNO DO VETOR DESLOCAMENTO (εE) =======================
                        # === Aplica continuidade da componente normal de D sobre a bainha dielétrica ========
                        
                        elif type_p == 'primary_insulation':
                            er = field_surface['relative_permittivity']

                            # Produto escalar dos vetores unitários em RIBBON.FOR: COS(TH - ANG)
                            RDN = np.dot(un_p, un_rho_i)

                            # 4. RIBBON.FOR: TDN = -sin(TH-ANG) = sin(ANG-TH)
                            TDN = np.cross(un_p, un_rho_i)

                            # --- TABELA II.a: rho_i = rho_b (Interação para Observador SOBRE a fronteira dielétrica) --- 
                            if type_q == 'primary_insulation' and tag_p == tag_q:
                                if harmonic_ord == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = (0 - 1) * (rho_b / rho_i) * RDN

                                else: # Harmonic Terms (k>0)
                                    self.D_matrix[row_idx, col_idx] = - 0.5 * (er + 1) * (rho_b / rho_i)**(k-1) * RDN * harmonic_term
                                
                            # --- TABELA II.a: rho_i >= rho_b (Interação para Observador FORA da fronteira dielétrica) --- 
                            else:   
                                if harmonic_ord == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = (er - 1) * (rho_b / rho_i) * RDN
                                
                                else: # Harmonic Terms (k>0)
                                    if is_cosine_term:
                                        harmonic_term = np.cos(k * theta_i) * RDN - np.sin(k * theta_i) * TDN
                                    else:
                                        harmonic_term = np.sin(k * theta_i) * RDN + np.cos(k * theta_i) * TDN
                                    self.D_matrix[row_idx, col_idx] = 0.5 * (er - 1) * (rho_b / rho_i)**(k+1) * harmonic_term

                        # ====================================================================================
                        # ==== FIM DA LÓGICA DE CÁLCULO DO ELEMENTO DA MATRIZ D ==============================
                        # ====================================================================================

        # 3. Resolver o sistema e obter os resultados
        self.sigma_coeffs = np.linalg.solve(self.D_matrix, self.V_vector)
        self._calculate_generalized_capacitance()
        self._calculate_maxwellian_capacitance()

    def print_results(self):
        """Imprime um resumo dos resultados da simulação."""
        if self.C_maxwellian is None:
            print("Executando simulação primeiro...")
            self.run_simulation()

        if self.model.NF < 4: 
            matrix_viewer(self.D_matrix, "D Matrix")
            # matrix_viewer(self.sigma_coeffs, "Sigma Coefficients")
        else:
            print(f"\nD Matrix Shape: {self.D_matrix.shape}.")
        
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
            # Itera sobre cada superfície desse condutor (ex: 'conductor', 'primary_insulation')
            for surface_type, surface_data in conductor_surfaces.items():
                
                # Busca o raio correspondente na lista self.model.surfaces,
                # pois ele não está em self.collocation_data.
                matching_surface = next(s for s in self.model.surfaces if s['tag'] == tag and s['type'] == surface_type)
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

        # 3. Adicionar as formas dos círculos (esta parte não muda, pois já itera sobre self.model.surfaces)
        for surface in self.model.surfaces:
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
        all_tags = list(self.model.mtl.keys())
        if len(all_tags) != 2:
            print("Erro: plot_charge_density foi projetado para sistemas de 2 condutores.")
            return
        other_tag = next(tag for tag in all_tags if tag != tag_to_plot)

        center_plot = np.array(self.model.mtl[tag_to_plot]['center_point'])
        center_other = np.array(self.model.mtl[other_tag]['center_point'])

        conductor_surface = next(s for s in self.model.surfaces if s['tag'] == tag_to_plot and s['type'] == 'conductor')
        R = conductor_surface['radius']
        D = np.linalg.norm(center_plot - center_other)
        DR_ratio = D / R

        # 2. Calcular a Solução Analítica com Alinhamento e Sinal Corretos
        theta_plot = np.linspace(0, 2 * np.pi, 360)

        vec_to_other = center_other - center_plot
        angle_of_max_charge = np.arctan2(vec_to_other[1], vec_to_other[0])
        denominator = DR_ratio - 2 * np.cos(theta_plot - angle_of_max_charge)

        delta_v = self.model.mtl[tag_to_plot]['potential_to_infinity'] - self.model.mtl[other_tag]['potential_to_infinity']
        numerator = (DR_ratio**2 / 4) - 1
        
        # CORREÇÃO FINAL: Remover o abs() para preservar o sinal da carga
        charge_density_exact = (self.C_exact_bare_wires * delta_v / R) * (numerator / denominator)

        # 3. Reconstruir a Solução MoM para o Condutor Correto
        nfs_per_surface = [2 * s['fourier_order'] + 1 for s in self.model.surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)
        
        surface_index = next(i for i, s in enumerate(self.model.surfaces) if s['tag'] == tag_to_plot and s['type'] == 'conductor')
        
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
        coeffs_condutor1 = self.sigma_coeffs[:self.model.NF]

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
        code_indices_cos = np.arange(1, self.model.NF, 2)
        plot_j_cos = code_indices_cos + 1 # Mapeia [1, 3, 5] para [2, 4, 6]
        coeffs_cos = coeffs_condutor1[code_indices_cos]
        ratio_cos = np.abs(coeffs_cos) / np.abs(constant_term)

        # j=3, 5, 7,...: Coeficientes Senoidais (índices 2, 4, 6,... no código)
        code_indices_sin = np.arange(2, self.model.NF, 2)
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
        ax.set_xticks(np.arange(1, self.model.NF + 1))
        ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        