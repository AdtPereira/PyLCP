import copy
import numpy as np
import pandas as pd
import scipy.constants as spc
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from utils.case_utils import *
from mtl_main.source import MulticonductorTransmissionLine
from mom.bare_wire_systems import BareWireMoMSolver

class MoMVisualizer:
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
    def __init__(self, solver: BareWireMoMSolver, case_name: str = ""):
        if solver.mom_data['collocation'] is {} and solver.mom_data['galerkin'] is {}:
            raise ValueError("O objeto solver fornecido não foi executado. Chame solver.run_...() primeiro.")
        
        self.solver = solver
        self.case_name = case_name
        
        self.sigma_coeffs_col = solver.mom_data['collocation'].get('sigma_coeffs', None)
        self.sigma_coeffs_gal = solver.mom_data['galerkin'].get('sigma_coeffs', None)

        self.figsize = (12, 5)
        self.conductor_tag = 1

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.case_name, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

    def plot_convergence_rates(self, mtl: dict, nf_max=20):
        """ Plota a convergência da capacitância em função de NF, usando um modelo base. """
        print(f"\nGerando gráfico de convergência até NF={nf_max}...")
        
        C_FACTOR = 1e12  # Fator de conversão para pF/m
        
        # Extrai R e D da configuração base para calcular o valor exato.
        R = mtl[0]['radius'][1]
        center1 = np.array(mtl[0]['center_point'])
        center2 = np.array(mtl[1]['center_point'])
        D = np.linalg.norm(center1 - center2)
        c_exact = (np.pi * spc.epsilon_0) / np.arccosh(D / R / 2.0)
        nf_range = range(1, nf_max + 1)

        nf_odd, nf_even = [], [],
        cap_odd_col, cap_even_col = [], []
        cap_odd_gal, cap_even_gal = [], []

        for nf in nf_range:
            print(f"  Calculando NF={nf}...", end='\r')
            temp_mtl = copy.deepcopy(mtl)
            temp_mtl[0]['fourier_order'] = nf
            temp_mtl[1]['fourier_order'] = nf

            mtl_model = MulticonductorTransmissionLine(temp_mtl)
            solver = BareWireMoMSolver(mtl_model)
            solver.run_collocation_method()
            solver.run_galerkin_method()
            
            if nf % 2 != 0:
                nf_odd.append(nf)
                cap_odd_col.append(solver.mom_data['collocation']['maxwellian_capacitance'].item() * C_FACTOR)
                cap_odd_gal.append(solver.mom_data['galerkin']['maxwellian_capacitance'].item() * C_FACTOR)
            else:
                nf_even.append(nf)
                cap_even_col.append(solver.mom_data['collocation']['maxwellian_capacitance'].item() * C_FACTOR)
                cap_even_gal.append(solver.mom_data['galerkin']['maxwellian_capacitance'].item() * C_FACTOR)

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.axhline(y=c_exact * C_FACTOR, color='k', linestyle=':', label=f'Exactly Value = {c_exact*C_FACTOR:.2f} pF/m')
        
        ax.plot(nf_odd, cap_odd_col, linestyle='none', marker='x', markersize=4, fillstyle='none', markeredgecolor='black', label='MoM Collocation')
        ax.plot(nf_even, cap_even_col, linestyle='none', marker='x', markersize=4, fillstyle='none', markeredgecolor='black')
        
        ax.plot(nf_odd, cap_odd_gal, linestyle='none', marker='o', markersize=4, color='red', label='MoM Galerkin')
        ax.plot(nf_even, cap_even_gal, linestyle='none', marker='o', markersize=4, color='red')
        
        ax.set_title(f'Convergência da Capacitância para D/R = {D/R:.2f}')
        ax.set_xlabel('NF - Número de Coeficientes de Fourier por Fio')
        ax.set_ylabel('Capacitância (pF/m)')
        ax.set_xticks(np.arange(0, nf_max + 1, 2))
        ax.set_xlim(0, nf_max); ax.set_ylim(bottom=0)
        ax.set_ylim(0, max(cap_odd_col + cap_even_col) * 1.1)
        ax.grid(False)
        ax.legend()
        save_figure(fig, self.results_dir, base_filename='convergence_rates')
        plt.tight_layout()

    def plot_collocation_points(self):
        """
        Gera um gráfico interativo dos pontos de colocação usando Plotly,
        refletindo a nova estrutura de dicionário de self.collocation_data.
        """
        plot_data = []
        surfaces = self.solver.model.surfaces        
        collocation_data = self.solver.mom_data['collocation']['data']
        
        # Itera sobre cada 'tag' de condutor no dicionário (ex: 0, 1)
        for tag, conductor_surfaces in collocation_data.items():
            for surface_type, surface_data in conductor_surfaces.items():                
                matching_surface = next(s for s in surfaces if s['tag'] == tag and s['type'] == surface_type)
                for i, pt in enumerate(surface_data['observation']['cartesian']):
                    plot_data.append({
                        'x': pt[0], 'y': pt[1],
                        'type': 'Observation',
                        'tag': tag,
                        'surface': surface_type,
                        'radius': matching_surface['radius'],
                        'angle_rad': surface_data['observation']['angles_rad'][i]
                    })

        df = pd.DataFrame(plot_data)
        fig = go.Figure()

        for surface in surfaces:
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
 
    def plot_surface_charge_density(self, tag_to_plot=1, comsol_data: dict = None):
        """
        Plota a densidade de carga para um condutor específico, alinhando
        dinamicamente a solução exata com a geometria real do sistema.

        Args:
            tag_to_plot (int): A 'tag' do condutor para o qual a densidade de
                            carga será plotada.
        """
        mtl = self.solver.model.mtl
        surfaces = self.solver.model.surfaces
        tag_to_plot = self.conductor_tag

        # 1. Obter dados do condutor a ser plotado e de seu par
        all_tags = list(mtl.keys())
        if len(all_tags) != 2:
            print("Erro: plot_charge_density foi projetado para sistemas de 2 condutores.")
            return
        other_tag = next(tag for tag in all_tags if tag != tag_to_plot)

        center_plot = np.array(mtl[tag_to_plot]['center_point'])
        center_other = np.array(mtl[other_tag]['center_point'])
        conductor_surface = next(s for s in surfaces if s['tag'] == tag_to_plot and s['type'] == 'conductor')
        
        R = conductor_surface['radius']
        D = np.linalg.norm(center_plot - center_other)
        DR_ratio = D / R

        # 2. Calcular a Solução Analítica com Alinhamento e Sinal Corretos
        theta_plot = np.linspace(0, 2 * np.pi, 360)
        vec_to_other = center_other - center_plot
        angle_of_max_charge = np.arctan2(vec_to_other[1], vec_to_other[0])
        denominator = DR_ratio - 2 * np.cos(theta_plot - angle_of_max_charge)

        delta_v = mtl[tag_to_plot]['potential_to_infinity'] - mtl[other_tag]['potential_to_infinity']
        numerator = (DR_ratio**2 / 4) - 1
        charge_density_exact = (self.solver.C_exact_bare_wires * delta_v / R) * (numerator / denominator)

        # 3. Reconstruir a Solução MoM para o Condutor Correto
        nfs_per_surface = [2 * s['fourier_order'] + 1 for s in surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)        
        surface_index = next(i for i, s in enumerate(surfaces) if s['tag'] == tag_to_plot and s['type'] == 'conductor')
        
        offset = offsets[surface_index]
        nf = nfs_per_surface[surface_index]
        coeffs_col = self.sigma_coeffs_col[offset : offset + nf]
        coeffs_gal = self.sigma_coeffs_gal[offset : offset + nf]

        charge_density_col = np.zeros_like(theta_plot)
        charge_density_gal = np.zeros_like(theta_plot)
        charge_density_col += coeffs_col[0]
        charge_density_gal += coeffs_gal[0]
        
        max_k = (nf - 1) // 2
        for k in range(1, max_k + 1):
            cos_kt = np.cos(k * theta_plot)
            sin_kt = np.sin(k * theta_plot)
            charge_density_col += coeffs_col[2*k-1] * cos_kt + coeffs_col[2*k] * sin_kt
            charge_density_gal += coeffs_gal[2*k-1] * cos_kt + coeffs_gal[2*k] * sin_kt

        # 4. Geração do Gráfico
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.plot(np.rad2deg(theta_plot), charge_density_exact * 1e9, 'k-', linewidth=1, label='Solução Exata')
        ax.plot(np.rad2deg(theta_plot), charge_density_gal * 1e9, 'b--', linewidth=1, label=f'MoM Galerkin')
        ax.plot(np.rad2deg(theta_plot), charge_density_col * 1e9, 'b:', linewidth=1, label=f'MoM Colocação')

        # --- Plot COMSOL Data if provided ---
        if comsol_data is not None and isinstance(comsol_data, dict):
            comsol_df = comsol_data.get('curve_1')

            if comsol_df is not None and not comsol_df.empty:
                angle_deg = (comsol_df['arc_length'] / R) * (180 / np.pi)
                ax.plot(angle_deg, comsol_df['surface_charge_density'], 'r.', markersize=2, label='COMSOL')
        
        ax.set_xlim(0, 360)
        ax.set_ylim(bottom=0)
        ax.set_title(f'Distribuição de Carga (Condutor {tag_to_plot}) com D/R = {DR_ratio:.2f}')
        ax.set_xlabel('Ângulo (Graus)'); ax.set_ylabel('Densidade de Carga (nC/m²)')
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.set_xticks(np.arange(0, 361, 90)); ax.set_xlim(0, 360)
        ax.legend()
        save_figure(fig, self.results_dir, base_filename='surface_charge_density')
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
        NF = self.solver.NF
        sigma_coeffs_col = self.solver.mom_data['collocation'].get('sigma_coeffs', None)
        sigma_coeffs_gal = self.solver.mom_data['galerkin'].get('sigma_coeffs', None)
        
        # Isola os coeficientes do primeiro condutor
        c1_gal = sigma_coeffs_gal[:NF]
        c1_col = sigma_coeffs_col[:NF]

        # j=2, 4, 6,...: Coeficientes Cossenoidais (índices 1, 3, 5,... no código)
        plot_j_cos = np.arange(1, NF, 2) + 1 
        ratio_cos_gal = np.abs(c1_gal[np.arange(1, NF, 2)]) / np.abs(c1_gal[0])
        ratio_cos_col = np.abs(c1_col[np.arange(1, NF, 2)]) / np.abs(c1_col[0])

        # j=3, 5, 7,...: Coeficientes Senoidais (índices 2, 4, 6,... no código)
        plot_j_sin = np.arange(2, NF, 2) + 1 
        ratio_sin_gal = np.abs(c1_gal[np.arange(2, NF, 2)]) / np.abs(c1_gal[0])
        ratio_sin_col = np.abs(c1_col[np.arange(2, NF, 2)]) / np.abs(c1_col[0])

        # --- Geração do Gráfico ---
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)

        # Plota cada série com um marcador distinto
        ax.plot([1], [1.0], marker='s', markersize=6, linestyle='none',
                fillstyle='none', markeredgecolor='black', label='Termo Constante (j=1)')

        ax.plot(plot_j_cos, ratio_cos_gal, marker='^', markersize=6, linestyle='none',
                fillstyle='none', markeredgecolor='black', label='Coeficientes Cossenoidais (Galerkin)')
        
        ax.plot(plot_j_cos, ratio_cos_col, marker='^', markersize=6, linestyle='none',
                fillstyle='none', markeredgecolor='red', label='Coeficientes Cossenoidais (Colocação)')

        ax.plot(plot_j_sin, ratio_sin_gal, marker='o', markersize=4, linestyle='none',
                fillstyle='none', markeredgecolor='black', label='Coeficientes Senoidais (Galerkin)')
        
        ax.plot(plot_j_sin, ratio_sin_col, marker='o', markersize=4, linestyle='none',
                fillstyle='none', markeredgecolor='red', label='Coeficientes Senoidais (Colocação)')

        ax.set_title(f'Magnitude Normalizada dos Coeficientes Harmônicos (d/a = {self.solver.DR_ratio:.1f})', fontsize=14)
        ax.set_ylabel(r'$|\alpha_{nj} / \alpha_{n1}|$', fontsize=12)
        ax.set_xlabel('Índice do Coeficiente (j)', fontsize=12)
        ax.legend()
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=-0.05)
        ax.set_xticks(np.arange(1, NF + 1))
        ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.7)
        save_figure(fig, self.results_dir, base_filename='harmonic_coefficients')
        plt.tight_layout()

    def print_terminal_results(self):
        """Imprime um resumo dos resultados da simulação."""

        for method in ['collocation', 'galerkin']:
            mom_data = self.solver.mom_data[method]
            
            if mom_data is not {}:
                print(f"\n============== pyMoM {method.capitalize()} Method =============")
                if self.solver.NF < 4:
                    matrix_viewer(mom_data['moment_matrix'], "Moment Matrix")
                    matrix_viewer(mom_data['sigma_coeffs'], "Sigma Coefficients Vector")
                
                print(f"\nMoment Matrix Shape: {mom_data['moment_matrix'].shape}.")
                matrix_viewer(mom_data['generalized_capacitance'], "MoM Generalized Capacitance Matrix (F/m)")
                matrix_viewer(mom_data['maxwellian_capacitance'], "Maxwellian Bifilar Capacitance (MoM) (F/m)")
                
            else:
                print(f"\nNenhum dado disponível para o método {method.capitalize()}. Execute o solver primeiro.")
            
        print("\n")
    