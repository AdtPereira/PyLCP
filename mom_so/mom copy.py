import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0
from matplotlib.patches import Circle

from .mtl import MulticonductorTransmissionLine as MTL

class BifilarMoM(MTL):
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

            sim = BifilarMoM(temp_config)
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

    
