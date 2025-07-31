import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0
from matplotlib.patches import Circle

from .mtl import MulticonductorTransmissionLine as MTL


class TwoBareWireSystem(MTL):
    """
    Calcula a capacitância e distribuição de carga de uma linha bifilar (2 fios)
    usando o Método dos Momentos (MoM).

    Esta classe herda de MulticonductorTransmissionLine (MTL) e a especializa
    para o caso de dois fios, derivando seus parâmetros de uma configuração 'mtl'.

    Executa a simulação completa do MoM, preenchendo todos os atributos de resultado.
    A construção da matriz D agora inclui os termos de expansão constante, cossenoidal e senoidal,
    conforme as expressões (20a), (20b) e (20c) de Clements (1975).
    """
    def __init__(self, mtl: dict):
        # Chama o construtor da classe pai para processar a configuração.
        super().__init__(mtl)

        assert isinstance(mtl, dict), "O parâmetro mtl deve ser um dicionário com a configuração da linha."
        assert len(mtl) == 2, "A classe BifilarMoM foi projetada para modelos com 2 condutores."

        # Parâmetros de entrada
        c_idx = 0
        self.R = self.surfaces[c_idx]['radius']
        self.D = self.D_pq[c_idx, c_idx + 1]
        self.NF = self.surfaces[c_idx]['fourier_order'] + 1
        self.DR_ratio = self.D / self.R

        assert self.DR_ratio > 2, "A razão D/R deve ser maior que 2 para garantir a convergência da solução."

        # Atributos de resultado (inicializados como None)
        self.collocation_data = None
        self.D_matrix = None
        self.T_matrix = None
        self.V_vector = None
        self.sigma_coeffs = None
        self.C_generalized = None
        self.C_maxwellian = None
        self.C_exact = None

    def _calculate_collocation_points(self):
        """
        Calcula e armazena os pontos de colocação (fonte e observação).
        """
        centers = [np.array([-self.D / 2.0, 0.0]), np.array([self.D / 2.0, 0.0])]
        source_angles = np.linspace(0, 2 * np.pi, self.NF, endpoint=False)
        observation_angles = source_angles + np.pi / (self.NF + 1)

        points = {'source': [], 'observation': []}
        for i in range(2):
            center = centers[i]
            source_pts = center + self.R * np.array([np.cos(source_angles), np.sin(source_angles)]).T
            obs_pts = center + self.R * np.array([np.cos(observation_angles), np.sin(observation_angles)]).T
            points['source'].append(source_pts)
            points['observation'].append(obs_pts)

        self.collocation_data = {
            'cartesian': points,
            'angles_rad': {'source': source_angles, 'observation': observation_angles}
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
        Calcula a matriz de capacitância física (Maxwelliana) a partir da generalizada.
        """
        gc = self.C_generalized
        assert isinstance(gc, np.ndarray) and np.sum(gc) != 0 and gc.shape == (2, 2)

        row_sum = np.sum(gc, axis=1)
        column_sum = np.sum(gc, axis=0)
        
        # Para um sistema de 2 condutores, a matriz maxwelliana é 1x1
        # C11 = c11 - (sum(row1) * sum(col1)) / sum(total)
        c11 = gc[1, 1]
        row_1_sum = row_sum[1]
        column_1_sum = column_sum[1]
        
        self.C_maxwellian = c11 - (row_1_sum * column_1_sum) / np.sum(gc)

    def run_simulation(self):
        """
        Executa a simulação completa do MoM, preenchendo todos os atributos de resultado.
        """
        self._calculate_collocation_points()
        
        observation_points = np.vstack(self.collocation_data['cartesian']['observation'])
        centers = [np.array([-self.D / 2.0, 0.0]), np.array([self.D / 2.0, 0.0])]
        
        N = 2 * self.NF
        self.D_matrix = np.zeros((N, N))
        self.V_vector = np.zeros(N)

        for i in range(N):
            field_idx = i // self.NF
            self.V_vector[i] = 1.0 if field_idx == 0 else -1.0
            field_point = observation_points[i]
            field_angle = np.arctan2(field_point[1] - centers[field_idx][1], field_point[0] - centers[field_idx][0])

            for j in range(N):
                source_idx = j // self.NF

                # Índice local da função de base no condutor de origem (0 a NF-1)
                source_harmonic_order = j % self.NF

                # Termos de auto-interação
                if source_idx == field_idx:  
                    if source_harmonic_order == 0:
                        self.D_matrix[i, j] = (-self.R / epsilon_0) * np.log(self.R)
                    else:
                        k = source_harmonic_order
                        self.D_matrix[i, j] = (self.R / (2 * k * epsilon_0)) * np.cos(k * field_angle)
                
                # Termos de interação mútua
                else:
                    dist_vec = field_point - centers[source_idx]
                    dist = np.linalg.norm(dist_vec)
                    if source_harmonic_order == 0:
                        self.D_matrix[i, j] = (-self.R / epsilon_0) * np.log(dist)
                    else:
                        k = source_harmonic_order
                        phi = np.arctan2(dist_vec[1], dist_vec[0])
                        self.D_matrix[i, j] = (self.R / (2 * k * epsilon_0)) * ((self.R / dist)**k) * np.cos(k * phi)

                # # Termo constante (k = 0)
                # # Eq. (20a) [1]
                # if source_harmonic_order == 0:
                #     # Termo de auto-interação
                #     if source_idx == field_idx:
                #         self.D_matrix[i, j] = (-self.R / epsilon_0) * np.log(self.R)
                    
                #     # Termo de interação mútua
                #     else:
                #         dist_vec = field_point - centers[source_idx]
                #         dist = np.linalg.norm(dist_vec)
                #         self.D_matrix[i, j] = (-self.R / epsilon_0) * np.log(dist)

                # # Termos Harmônicos (k > 0)
                # else:
                #     # Índices ímpares -> cosseno; Índices pares -> seno
                #     is_cosine_term = (source_harmonic_order % 2 != 0)


        self.sigma_coeffs = np.linalg.solve(self.D_matrix, self.V_vector)
        self.C_exact = (np.pi * epsilon_0) / np.arccosh(self.DR_ratio / 2.0)

        # Calcula as matrizes de capacitância
        self._calculate_generalized_capacitance()
        self._calculate_maxwellian_capacitance()
    
    def print_results(self):
        """Imprime um resumo dos resultados da simulação."""
        if self.C_maxwellian is None:
            print("Executando simulação primeiro...")
            self.run_simulation()

        if self.NF < 3: 
            print(f"\nD Matrix (Shape: {self.D_matrix.shape}):\n{self.D_matrix}")
            print(f"\nT Matrix (Inverse of D) (Shape: {self.T_matrix.shape}):\n{self.T_matrix}")
        
        print(f"\n--- Results for D/R = {self.DR_ratio} and NF={self.NF} ---")
        print(f"\nSigma Coefficients (Shape: {self.sigma_coeffs.shape}):\n{self.sigma_coeffs}")
        print(f"\nMoM Generalized Capacitance Matrix (F/m): \n{self.C_generalized}")
        print(f"\nExact Capacitance: {self.C_exact * 1E12:.4f} pF/m")
        print(f"\nMaxwellian Bifilar Capacitance (MoM): {self.C_maxwellian * 1E12:.4f} pF/m.")

    def plot_collocation_points(self, coord_mode='Cartesian'):
        """
        Plota os condutores e os pontos de colocação.

        Args:
            coord_mode (str): 'Cartesian' ou 'Polar' para o tipo de anotação.
        """
        if self.collocation_data is None:
            self.run_simulation()
            
        fig, ax = plt.subplots(figsize=(12, 9))
        centers = [np.array([-self.D / 2.0, 0.0]), np.array([self.D / 2.0, 0.0])]
        wire1 = Circle(centers[0], self.R, facecolor='none', edgecolor='k', lw=1, ls='--')
        wire2 = Circle(centers[1], self.R, facecolor='none', edgecolor='k', lw=1, ls='--')
        ax.add_patch(wire1)
        ax.add_patch(wire2)

        source_w1, source_w2 = self.collocation_data['cartesian']['source']
        obs_w1, obs_w2 = self.collocation_data['cartesian']['observation']

        ax.plot(source_w1[:, 0], source_w1[:, 1], 'o', c='blue', label='Fonte')
        ax.plot(obs_w1[:, 0], obs_w1[:, 1], 'x', c='red', label='Observação')
        ax.plot(source_w2[:, 0], source_w2[:, 1], 'o', c='blue')
        ax.plot(obs_w2[:, 0], obs_w2[:, 1], 'x', c='red')

        point_sets = [(source_w1, self.collocation_data['angles_rad']['source']),
                      (obs_w1, self.collocation_data['angles_rad']['observation']),
                      (source_w2, self.collocation_data['angles_rad']['source']),
                      (obs_w2, self.collocation_data['angles_rad']['observation'])]

        for points, angles in point_sets:
            for i, pt in enumerate(points):
                if coord_mode == 'Polar':
                    annotation = f'({self.R:.3f}, {angles[i]:.3f} rad)'
                elif coord_mode == 'Cartesian':
                    annotation = f'({pt[0]:.3f}, {pt[1]:.3f})'
                else:
                    raise ValueError("coord_mode deve ser 'Cartesian' ou 'Polar'.")
                ax.text(pt[0], pt[1] + self.R * 0.15, annotation, fontsize=8, ha='center', va='bottom', rotation=15)
        
        ax.set_title('Mapa de Pontos de Colocação')
        ax.set_xlabel('Coordenada X (m)'); ax.set_ylabel('Coordenada Y (m)')
        ax.set_xlim(-self.D * 1.5, self.D * 1.5); ax.set_ylim(-self.D, self.D)
        ax.grid(True, linestyle='--', alpha=0.6); ax.set_aspect('equal', adjustable='box')
        plt.tight_layout(); ax.legend()

    def plot_charge_density(self):
        """
        Plota o gráfico de comparação da densidade de carga (MoM vs. Exata).
        """
        if self.sigma_coeffs is None:
            self.run_simulation()

        theta_plot = np.linspace(0, 2 * np.pi, 720)
        
        # Solução exata
        delta_v = 2.0
        numerator = (self.DR_ratio**2 / 4) - 1
        denominator = self.DR_ratio - 2 * np.cos(theta_plot)
        charge_density_exact = (self.C_exact * delta_v / self.R) * (numerator / denominator)

        # Solução MoM
        charge_density_mom = np.zeros_like(theta_plot)
        for k in range(self.NF):
            charge_density_mom += self.sigma_coeffs[k] * np.cos(k * theta_plot)

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(np.rad2deg(theta_plot), charge_density_exact, color='r', linestyle='-', label='Solução Exata')
        ax.plot(np.rad2deg(theta_plot), charge_density_mom, color='k', linestyle='-.', label=f'MoM (NF={self.NF})')
        ax.set_title(f'Distribuição de Carga com D/R = {self.DR_ratio} e NF={self.NF}')
        ax.set_xlabel('Ângulo (Graus)'); ax.set_ylabel('Densidade de Carga (C/m²)')
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.set_xticks(np.arange(0, 361, 90)); ax.set_xlim(0, 360)
        plt.tight_layout(); ax.legend()
        
    @staticmethod
    def plot_convergence_rates(MTL, nf_max=20):
        """ Plota a convergência da capacitância em função de NF, usando um modelo base. """
        print(f"\nGerando gráfico de convergência até NF={nf_max}...")
        
        C_FACTOR = 1e12  # Fator de conversão para pF/m
        
        # Extrai R e D da configuração base para calcular o valor exato.
        R = MTL['data'][0]['radius'][1]
        center1 = np.array(MTL['data'][0]['center_point'])
        center2 = np.array(MTL['data'][1]['center_point'])
        D = np.linalg.norm(center1 - center2)

        c_exact = (np.pi * epsilon_0) / np.arccosh(D / R / 2.0)
        nf_range = range(1, nf_max + 1)

        nf_odd, nf_even, cap_odd, cap_even = [], [], [], []

        for nf in nf_range:
            # Cria uma cópia temporária do modelo para modificar NF sem alterar o original.
            temp_config = copy.deepcopy(MTL)
            temp_config['data'][0]['fourier_order'] = nf
            temp_config['data'][1]['fourier_order'] = nf

            sim = TwoBareWireSystem(temp_config)
            sim.run_simulation()
            
            if nf % 2 != 0:
                nf_odd.append(nf)
                cap_odd.append(sim.C_maxwellian * C_FACTOR)
            else:
                nf_even.append(nf)
                cap_even.append(sim.C_maxwellian * C_FACTOR)
        
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.axhline(y=c_exact * C_FACTOR, color='k', linestyle='-', label=f'Valor Exato = {c_exact*C_FACTOR:.2f} pF/m')
        ax.plot(nf_odd, cap_odd, linestyle='none', marker='^', markersize=8, fillstyle='none', markeredgecolor='black', label='NF Ímpar')
        ax.plot(nf_even, cap_even, linestyle='none', marker='*', markersize=9, color='black', label='NF Par')
        ax.set_title(f'Convergência da Capacitância para D/R = {D/R:.2f}')
        ax.set_xlabel('NF - Número de Coeficientes de Fourier por Fio')
        ax.set_ylabel('Capacitância (pF/m)')
        ax.set_xticks(np.arange(0, nf_max + 1, 2))
        ax.set_xlim(0, nf_max); ax.set_ylim(bottom=0)
        ax.set_ylim(0, max(cap_odd + cap_even) * 1.1)
        ax.grid(False)
        ax.legend()
        plt.tight_layout()


class MulticonductorBareWireSystems(MTL):
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
        assert len(self.surfaces) > 1, "A classe MulticonductorBareWireSystems foi projetada para modelos com mais de 1 superfície."
        
        # Número de coeficientes harmônicos de Fourier por condutor
        self.NF = [2*surface['fourier_order']+1 for surface in self.surfaces][0]
        
        # Raio das superfícies
        self.R = [surface['radius'] for surface in self.surfaces][0]       

        # Coordenadas dos centros dos condutores
        self.center_points = [np.array(surface['center']) for surface in self.surfaces] 
        
        self.D = self.D_pq[0, 1]
        self.DR_ratio = self.D / self.R

        assert self.DR_ratio > 2, "A razão D/R deve ser maior que 2 para garantir a convergência da solução."

        # Atributos de resultado (inicializados como None)
        self.collocation_data = None
        self.D_matrix = None
        self.T_matrix = None
        self.V_vector = None
        self.sigma_coeffs = None
        self.C_generalized = None
        self.C_maxwellian = None
        self.C_exact = None

    def _calculate_collocation_points(self):
        """
        Calcula e armazena os pontos de colocação (fonte e observação).
        Esta versão corrigida utiliza os centros dos condutores herdados da classe base MTL.
        """
        points = {'source': [], 'observation': []}
        source_angles = np.linspace(0, 2 * np.pi, self.NF, endpoint=False)
        observation_angles = source_angles + np.pi / self.NF

        for i in range(len(self.surfaces)):
            center = self.center_points[i]
            radius = self.surfaces[i]['radius'] 
            source_pts = center + radius * np.array([np.cos(source_angles), np.sin(source_angles)]).T
            obs_pts = center + radius * np.array([np.cos(observation_angles), np.sin(observation_angles)]).T
            points['source'].append(source_pts)
            points['observation'].append(obs_pts)

        self.collocation_data = {
            'cartesian': points,
            'angles_rad': {'source': source_angles, 'observation': observation_angles}
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
        Calcula a matriz de capacitância física (Maxwelliana) a partir da generalizada.
        """
        gc = self.C_generalized
        assert isinstance(gc, np.ndarray) and np.sum(gc) != 0 and gc.shape == (2, 2)

        row_sum = np.sum(gc, axis=1)
        column_sum = np.sum(gc, axis=0)
        
        # Para um sistema de 2 condutores, a matriz maxwelliana é 1x1
        # C11 = c11 - (sum(row1) * sum(col1)) / sum(total)
        c11 = gc[1, 1]
        row_1_sum = row_sum[1]
        column_1_sum = column_sum[1]
        
        self.C_maxwellian = c11 - (row_1_sum * column_1_sum) / np.sum(gc)

    def run_simulation(self):
        """ Executa a simulação completa do MoM, preenchendo todos os atributos de resultado. """
        self._calculate_collocation_points()
        
        observation_points = np.vstack(self.collocation_data['cartesian']['observation'])
        centers = [np.array([-self.D / 2.0, 0.0]), np.array([self.D / 2.0, 0.0])]
        
        N = self.N
        self.D_matrix = np.zeros((N, N))
        self.V_vector = np.zeros(N)

        for i in range(N):
            field_idx = i // self.NF
            self.V_vector[i] = 1.0 if field_idx == 0 else -1.0
            field_point = observation_points[i]
            field_angle = np.arctan2(field_point[1] - centers[field_idx][1], field_point[0] - centers[field_idx][0])

            for j in range(N):
                source_idx = j // self.NF

                # Índice local da função de base no condutor de origem (0 a NF-1)
                source_harmonic_order = j % self.NF

                # # Termos de auto-interação
                # if source_idx == field_idx:  
                #     if source_harmonic_order == 0:
                #         self.D_matrix[i, j] = (-self.R / epsilon_0) * np.log(self.R)
                #     else:
                #         k = source_harmonic_order
                #         self.D_matrix[i, j] = (self.R / (2 * k * epsilon_0)) * np.cos(k * field_angle)
                
                # # Termos de interação mútua
                # else:
                #     dist_vec = field_point - centers[source_idx]
                #     dist = np.linalg.norm(dist_vec)
                #     if source_harmonic_order == 0:
                #         self.D_matrix[i, j] = (-self.R / epsilon_0) * np.log(dist)
                #     else:
                #         k = source_harmonic_order
                #         phi = np.arctan2(dist_vec[1], dist_vec[0])
                #         self.D_matrix[i, j] = (self.R / (2 * k * epsilon_0)) * ((self.R / dist)**k) * np.cos(k * phi)

                # Termo constante (k = 0)
                # Eq. (20a) [1]
                if source_harmonic_order == 0:
                    # Termo de auto-interação
                    if source_idx == field_idx:
                        self.D_matrix[i, j] = (-self.R / epsilon_0) * np.log(self.R)
                    
                    # Termo de interação mútua
                    else:
                        dist_vec = field_point - centers[source_idx]
                        dist = np.linalg.norm(dist_vec)
                        self.D_matrix[i, j] = (-self.R / epsilon_0) * np.log(dist)

                # Termos Harmônicos (k > 0)
                else:
                    # Índices ímpares -> cosseno; Índices pares -> seno
                    is_cosine_term = (source_harmonic_order % 2 != 0)
                    k = (source_harmonic_order + 1) // 2 if is_cosine_term else source_harmonic_order // 2

                    # Termos de auto-interação
                    if source_idx == field_idx:
                        term = np.cos(k * field_angle) if is_cosine_term else np.sin(k * field_angle)
                        self.D_matrix[i, j] = (self.R / (2 * k * epsilon_0)) * term

                    # Termos de interação mútua
                    else:
                        dist_vec = field_point - centers[source_idx]
                        dist = np.linalg.norm(dist_vec)
                        phi = np.arctan2(dist_vec[1], dist_vec[0])
                        term = np.cos(k * phi) if is_cosine_term else np.sin(k * phi)
                        self.D_matrix[i, j] = (self.R / (2 * k * epsilon_0)) * ((self.R / dist)**k) * term

        self.sigma_coeffs = np.linalg.solve(self.D_matrix, self.V_vector)
        self.C_exact = (np.pi * epsilon_0) / np.arccosh(self.DR_ratio / 2.0)

        # Calcula as matrizes de capacitância
        self._calculate_generalized_capacitance()
        self._calculate_maxwellian_capacitance()
    
    def print_results(self):
        """Imprime um resumo dos resultados da simulação."""
        if self.C_maxwellian is None:
            print("Executando simulação primeiro...")
            self.run_simulation()

        if self.NF < 3: 
            print(f"\nD Matrix (Shape: {self.D_matrix.shape}):\n{self.D_matrix}")
            print(f"\nT Matrix (Inverse of D) (Shape: {self.T_matrix.shape}):\n{self.T_matrix}")
        
        print(f"\n--- Results for D/R = {self.DR_ratio} and NF={self.NF} ---")
        print(f"\nSigma Coefficients (Shape: {self.sigma_coeffs.shape}):\n{self.sigma_coeffs}")
        print(f"\nMoM Generalized Capacitance Matrix (F/m): \n{self.C_generalized}")
        print(f"\nExact Capacitance: {self.C_exact * 1E12:.4f} pF/m")
        print(f"\nMaxwellian Bifilar Capacitance (MoM): {self.C_maxwellian * 1E12:.4f} pF/m.")

    def plot_collocation_points(self):
        """
        Alternativa de plotagem que gera um gráfico interativo usando Plotly.
        As coordenadas e outros detalhes aparecem ao passar o mouse sobre os pontos.
        """
        try:
            import plotly.graph_objects as go
            import pandas as pd
        except ImportError:
            print("Para usar este método, instale as bibliotecas necessárias: pip install plotly pandas")
            return

        if self.collocation_data is None:
            self._calculate_collocation_points()

        # 1. Preparar os dados para o Plotly
        plot_data = []
        for i in range(len(self.surfaces)):
            for pt_type, marker_symbol in [('source', 'circle'), ('observation', 'cross')]:
                points = self.collocation_data['cartesian'][pt_type][i]
                angles = self.collocation_data['angles_rad'][pt_type]
                radius = self.surfaces[i]['radius']
                for j, pt in enumerate(points):
                    plot_data.append({
                        'x': pt[0], 'y': pt[1],
                        'type': pt_type.capitalize(),
                        'conductor': f'C{i+1}',
                        'radius': radius,
                        'angle_rad': angles[j]
                    })
        df = pd.DataFrame(plot_data)

        # 2. Criar a figura
        fig = go.Figure()

        # 3. Adicionar os condutores como formas
        for surface in self.surfaces:
            fig.add_shape(type="circle",
                          xref="x", yref="y",
                          x0=surface['center'][0] - surface['radius'], y0=surface['center'][1] - surface['radius'],
                          x1=surface['center'][0] + surface['radius'], y1=surface['center'][1] + surface['radius'],
                          line_color="Black", fillcolor="LightGray", opacity=0.7)

        # 4. Adicionar os pontos de colocação (Fonte e Observação)
        for pt_type, color, symbol in [('Source', 'blue', 'circle'), ('Observation', 'red', 'x-thin')]:
            df_subset = df[df['type'] == pt_type]
            fig.add_trace(go.Scatter(
                x=df_subset['x'], y=df_subset['y'],
                mode='markers',
                marker=dict(color=color, symbol=symbol, size=8, line=dict(width=1, color='DarkSlateGrey')),
                name=pt_type,
                customdata=df_subset[['conductor', 'radius', 'angle_rad']],
                hovertemplate=(
                    f"<b>{pt_type}</b><br>"
                    "Condutor: %{customdata[0]}<br>"
                    "Coord X: %{x:.4f} m<br>"
                    "Coord Y: %{y:.4f} m<br>"
                    "Ângulo: %{customdata[2]:.3f} rad<br>"
                    "Raio: %{customdata[1]:.4f} m"
                    "<extra></extra>" # Remove o 'trace' box
                )
            ))

        # 5. Configurar o layout
        fig.update_layout(
            title='Mapa Interativo de Pontos de Colocação',
            xaxis_title='Coordenada X (m)',
            yaxis_title='Coordenada Y (m)',
            yaxis_scaleanchor="x", # Garante a proporção 1:1 (aspecto 'equal')
            yaxis_scaleratio=1,
            legend_title_text='Tipo de Ponto',
            template='plotly_white'
        )
        fig.show()
                    
    def plot_charge_density(self):
        """
        Plota o gráfico de comparação da densidade de carga, reconstruindo a partir da série completa.
        """
        if self.sigma_coeffs is None:
            self.run_simulation()

        theta_plot = np.linspace(0, 2 * np.pi, 360)
        
        # Solução exata para o caso bifilar
        delta_v = 2.0
        numerator = (self.DR_ratio**2 / 4) - 1
        denominator = self.DR_ratio - 2 * np.cos(theta_plot)
        charge_density_exact = (self.C_exact * delta_v / self.R) * (numerator / denominator)

        # Reconstrução da solução MoM com a série completa
        charge_density_mom = np.zeros_like(theta_plot)
        coeffs_condutor1 = self.sigma_coeffs[:self.NF]
        
        # Termo constante
        charge_density_mom += coeffs_condutor1[0]
        
        # Termos harmônicos
        max_k = (self.NF - 1) // 2
        for k in range(1, max_k + 1):
            cos_coeff = coeffs_condutor1[2 * k - 1]
            sin_coeff = coeffs_condutor1[2 * k]
            charge_density_mom += cos_coeff * np.cos(k * theta_plot) + sin_coeff * np.sin(k * theta_plot)

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(np.rad2deg(theta_plot), charge_density_exact, color='r', linestyle='-', label='Solução Exata')
        ax.plot(np.rad2deg(theta_plot), charge_density_mom, color='k', linestyle='-.', label=f'MoM (k_max={max_k})')
        ax.set_title(f'Distribuição de Carga com D/R = {self.DR_ratio:.2f}')
        ax.set_xlabel('Ângulo (Graus)'); ax.set_ylabel('Densidade de Carga (C/m²)')
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.set_xticks(np.arange(0, 361, 90)); ax.set_xlim(0, 360)
        plt.tight_layout()
        ax.legend()

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
        
    @staticmethod
    def plot_convergence_rates(MTL, nf_max=20):
        """ Plota a convergência da capacitância em função de NF, usando um modelo base. """
        print(f"\nGerando gráfico de convergência até NF={nf_max}...")
        
        C_FACTOR = 1e12  # Fator de conversão para pF/m
        
        # Extrai R e D da configuração base para calcular o valor exato.
        R = MTL['data'][0]['radius'][1]
        center1 = np.array(MTL['data'][0]['center_point'])
        center2 = np.array(MTL['data'][1]['center_point'])
        D = np.linalg.norm(center1 - center2)

        c_exact = (np.pi * epsilon_0) / np.arccosh(D / R / 2.0)
        nf_range = range(1, nf_max + 1)

        nf_odd, nf_even, cap_odd, cap_even = [], [], [], []

        for nf in nf_range:
            # Cria uma cópia temporária do modelo para modificar NF sem alterar o original.
            temp_config = copy.deepcopy(MTL)
            temp_config['data'][0]['fourier_order'] = nf
            temp_config['data'][1]['fourier_order'] = nf

            sim = TwoBareWireSystem(temp_config)
            sim.run_simulation()
            
            if nf % 2 != 0:
                nf_odd.append(nf)
                cap_odd.append(sim.C_maxwellian * C_FACTOR)
            else:
                nf_even.append(nf)
                cap_even.append(sim.C_maxwellian * C_FACTOR)
        
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.axhline(y=c_exact * C_FACTOR, color='k', linestyle='-', label=f'Valor Exato = {c_exact*C_FACTOR:.2f} pF/m')
        ax.plot(nf_odd, cap_odd, linestyle='none', marker='^', markersize=8, fillstyle='none', markeredgecolor='black', label='NF Ímpar')
        ax.plot(nf_even, cap_even, linestyle='none', marker='*', markersize=9, color='black', label='NF Par')
        ax.set_title(f'Convergência da Capacitância para D/R = {D/R:.2f}')
        ax.set_xlabel('NF - Número de Coeficientes de Fourier por Fio')
        ax.set_ylabel('Capacitância (pF/m)')
        ax.set_xticks(np.arange(0, nf_max + 1, 2))
        ax.set_xlim(0, nf_max); ax.set_ylim(bottom=0)
        ax.set_ylim(0, max(cap_odd + cap_even) * 1.1)
        ax.grid(False)
        ax.legend()
        plt.tight_layout()

