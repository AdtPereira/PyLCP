import os
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from utils.case_utils import *
from mtl_main.source import MulticonductorTransmissionLine

class DeContiModels:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, file_path: str, pul_data: dict, autoSave: bool = True):
        self.script_path = Path(file_path)
        self.autoSave = autoSave
        self.pul_data = pul_data
        
        self.f = pul_data['frequencies']
        self.w = 2 * np.pi * self.f
        
        self.figsize = (12, 5)
        self.xlim = tuple(pul_data['frequencies'][[0, -1]])

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.script_path.stem, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

        self.comparison = [
            {'key': 'magalhaes_xue', 'label': 'Magalhães/Xue',  'color': 'black', 'linestyle': '-', 'linewidth': 2},
            {'key': 'sunde', 'label': 'Sunde', 'color': 'gray', 'linestyle': ':', 'linewidth': 2},
            {'key': 'pollaczek', 'label': 'Pollaczek', 'color': 'blue', 'linestyle': '-', 'linewidth': 2},
            {'key': 'ametani', 'label': 'Ametani', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 2},
            {'key': 'deconti', 'label': 'De Conti', 'color': 'red', 'linestyle': ':', 'linewidth': 1},
            {'key': 'saad', 'label': 'Saad', 'color': 'red', 'linestyle': ':', 'linewidth': 1},
            {'key': 'wedepohl', 'label': 'Wedepohl and Wilcox', 'color': 'darkgray', 'linestyle': '-', 'linewidth': 1},
        ]
        
        self.paper_2023 = [
            {'key': 'p100',   'label': r'Integral Form. $\rho_e=100 \;\Omega m$',  'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
            {'key': 'p1000',  'label': r'Integral Form. $\rho_e=1000 \;\Omega m$', 'color': 'blue', 'linestyle': '-', 'linewidth': 2.0},
            {'key': 'p10000', 'label': r'Integral Form. $\rho_e=10000 \;\Omega m$','color': 'darkgreen', 'linestyle': '-', 'linewidth': 2.0},
            {'key': 'p100_deConti', 'label': 'De Conti Approx.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
            {'key': 'p1000_deConti', 'label': '', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
            {'key': 'p10000_deConti', 'label': '', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
        ]

        self.plot_configs = {
            'ground_return_impedance': {
                'suptitle': r'P.u.l. ground-return impedance for single buried bare-wire for $\rho_e=100 \;\Omega$ m and $\epsilon_r=10$',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'p': 0, 'q': 0,
                'x_lim': {'resistance': (1E4, 1E8), 'inductance': (1E-1, 1E8)},
                'y_lim': {'resistance': (0, 100), 'inductance': (0, 3.0)},
                'series_to_plot': self.comparison
            },
            'impedance': {
                'suptitle': fr"Flat arrangement's mutual ground-return impedance for $\varepsilon_{{r1}} = {10}$",
                'norm_title': 'Absolute Impedance',
                'angle_title': 'Angle of Impedance',
                'p': 2, 'q': 0,
                'y_lim': {'norm': (0, 25), 'angle': (20, 90)},
                'series_to_plot': self.paper_2023
            },
            'potential': {
                'suptitle': r"Flat arrangement's mutual ground-return potential coefficients for $\varepsilon_{r1} = 10$",
                'norm_title': 'Absolute Value',
                'angle_title': 'Angle of Potential Coefficient',
                'p': 2, 'q': 0,
                'y_lim': {'norm': (0, 14), 'angle': (-90, 90)},
                'series_to_plot': self.paper_2023
            }, 
        }

    def ground_return_impedance(self, graph_key='ground_return_impedance'):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            Zg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['impedance_matrix']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': series['linewidth']}
            ax1.plot(self.f, np.real(Zg[:, p, q]) , **style)
            ax2.plot(self.f, np.imag(Zg[:, p, q]) / self.w * 1e6, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['resistance'])
        ax1.set_ylim(config['y_lim']['resistance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_s \, (\Omega/m)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['inductance'])
        ax2.set_ylim(config['y_lim']['inductance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_s \, (\mu H/m)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=f'{graph_key}')

    def _impedance_subplots(self, graph_key, graph_form='norm_and_angle'):
        config = self.plot_configs[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            Zg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['impedance_matrix']
            style = {
                'label': series.get('label', ''),
                'color': series.get('color', 'black'),
                'linestyle': series.get('linestyle', '-'),
                'linewidth': series.get('linewidth', 1.0)}            

            # Configure left subplot (Resistance)
            if graph_form == 'norm_and_angle':
                ax1.plot(self.f, np.abs(Zg[:, p, q]), **style)
                ax2.plot(self.f, np.angle(Zg[:, p, q], deg=True), **style)

                ax1.set_xscale('log')
                ax1.set_xlim(self.xlim)
                ax1.set_ylim(config['y_lim']['norm'])
                ax1.legend(fontsize='small')
                ax1.set_xlabel('Frequency (Hz)')
                ax1.set_ylabel(fr'$|Zg_{{{p+1}{q+1}}}| \, (\Omega/m)$')
                ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
                ax1.set_title(config['norm_title'])

                # Configure right subplot (Inductance)
                ax2.set_xscale('log')
                ax2.set_xlim(self.xlim)
                ax2.set_ylim(config['y_lim']['angle'])
                ax2.legend(fontsize='small')
                ax2.set_xlabel('Frequency (Hz)')
                ax2.set_ylabel(fr'Angle of $Zg_{{{p+1}{q+1}}}$ (Degrees)')
                ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
                ax2.set_title(config['angle_title'])
                plt.tight_layout(rect=[0, 0, 1, 0.96])

            elif graph_form=='resistance_and_inductance':
                ax1.plot(self.f, np.real(Zg[:, p, q]) * 1e3, **style)
                ax2.plot(self.f, np.imag(Zg[:, p, q]) / self.w * 1e6, **style)

                ax1.set_xscale('log')
                ax1.set_yscale('log')
                ax1.set_xlim(self.xlim)
                # ax1.set_ylim(config['y_lim']['norm'])
                ax1.legend(fontsize='small')
                ax1.set_xlabel('Frequency (Hz)')
                # ax1.set_ylabel(fr'$|Zg_{{{p+1}{q+1}}}| \, (\Omega/m)$')
                ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
                ax1.set_title(config['norm_title'])

                # Configure right subplot (Inductance)
                ax2.set_xscale('log')
                ax2.set_xlim(self.xlim)
                # ax2.set_ylim(config['y_lim']['angle'])
                ax2.legend(fontsize='small')
                ax2.set_xlabel('Frequency (Hz)')
                # ax2.set_ylabel(fr'Angle of $Zg_{{{p+1}{q+1}}}$ (Degrees)')
                ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
                ax2.set_title(config['angle_title'])
                plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=f'{graph_key}_{graph_form}')

    def _potential_subplots(self, graph_key):
        config = self.plot_configs[graph_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            Pg = self.pul_data['scenarios'][series['key']]['earth_return_parameters']['potential_coefficient']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'markersize': 3, 'linewidth': 1.0}
            ax1.plot(self.f, np.abs(Pg[:, p, q]) * 1e-9, **style)
            ax2.plot(self.f, np.angle(Pg[:, p, q], deg=True), **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_xlim(self.xlim)
        ax1.set_ylim(config['y_lim']['norm'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$|Pg_{{{p+1}{q+1}}}| \times 10^9 \, (\Omega m s^{{-1}})$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['norm_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(self.xlim)
        ax2.set_ylim(config['y_lim']['angle'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'Angle of $Pg_{{{p+1}{q+1}}}$ (Degrees)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['angle_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename=f'{graph_key}')

    def fig_3(self, graph_form='norm_and_angle'):
        """Plots the data corresponding to Figure 3 from the reference."""
        self._impedance_subplots('impedance', graph_form)

    def fig_6(self):
        """Plots the data corresponding to Figure 6 from the reference."""
        self._potential_subplots('potential')

    def log_matricial_pul_parameters(self, discrete_pul_data, scale_units=True):
        """ Generates a terminal log report by slicing the discrete vectorized results. """
        def _matrix_to_string(matrix: np.ndarray) -> str:
            # ... (função auxiliar sem alteração)
            lines = []
            s_rows = [[f"{val:11.4e}" for val in row] for row in matrix]
            for row in s_rows:
                lines.append("  ".join(row))
            return "\n".join(lines)

        def _print_real_matrix(name: str, matrix_data: np.ndarray, unit: str):
            # ... (função auxiliar sem alteração)
            lines = _matrix_to_string(matrix_data).split('\n')
            num_rows = len(lines)
            middle_row_idx = num_rows // 2
            
            print("") 
            for i in range(num_rows):
                name_part = f"{name} = " if i == middle_row_idx else " " * (len(name) + 3)
                unit_part = f" {unit}" if i == middle_row_idx else ""
                print(f"{name_part}{lines[i]}{unit_part}")
            print("") 

        param_mapping = {
            '[Zs]': 'series_impedance_matrix',
            '[Ysh]': 'shunt_admittance_matrix',
            '[Zi]': 'internal_impedance_matrix',
            '[Z0]': 'earth-return_impedance_matrix',
        }

        print("\n--- Per-Unit-Length (PUL) Parameters Report (Discrete Frequencies) ---")
        
        # Itera sobre as frequências discretas e seus índices
        for i, freq in enumerate(discrete_pul_data.get('frequencies', [])):
            print("\n" + "="*80)
            print(f"Frequency: {freq:,.0f} Hz")
            print("="*80)
            w = 2 * np.pi * freq

            # Itera sobre as matrizes (Zs, Ysh, etc.)
            for name, matrix_key in param_mapping.items():
                if matrix_key in discrete_pul_data:
                    matrix_3d = discrete_pul_data[matrix_key]
                    matrix = matrix_3d[i, :, :]
                    
                    # (Lógica de impressão das matrizes, sem alteração)
                    if 'Z' in name: unit = "[Ohm/m]"
                    elif 'Y' in name: unit = "[S/m]"
                    else: unit = ""

                    r_lines = _matrix_to_string(matrix.real).split('\n')
                    x_lines = _matrix_to_string(matrix.imag).split('\n')
                    num_rows = len(r_lines)
                    middle_row_idx = num_rows // 2
                    
                    print("") 
                    for i_row in range(num_rows):
                        name_part = f"{name} = " if i_row == middle_row_idx else " " * (len(name) + 3)
                        op_part = " + j ".center(7) if i_row == middle_row_idx else " " * 7
                        unit_part = f" {unit}" if i_row == middle_row_idx else ""
                        print(f"{name_part}{r_lines[i_row]}{op_part}{x_lines[i_row]}{unit_part}")
                    print("")
                    
                    if name == '[Zs]':
                        l_matrix = (matrix.imag / w)
                        # ... (resto da lógica de L e C sem alteração)
                        l_unit = '[H/m]'
                        if scale_units:
                            l_matrix *= 1e6
                            l_unit = '[uH/m]'
                        _print_real_matrix('[L]', l_matrix, l_unit)

                    if name == '[Ysh]':
                        c_matrix = (matrix.imag / w)
                        c_unit = '[F/m]'
                        if scale_units:
                            c_matrix *= 1e12
                            c_unit = '[pF/m]'
                        _print_real_matrix('[C]', c_matrix, c_unit)
        print("\n" + "="*80)

class InternalLinesModels:
    """
    Classe para gerar gráficos a partir dos resultados vetorizados
    do módulo de linhas aéreas.
    """
    def __init__(self, pul_data: dict, pul_data_tubular: dict, model: MulticonductorTransmissionLine):
        """
        Inicializa o plotter com os dados das simulações tubular e sólida.

        Args:
            pul_data_tubular (dict): Dicionário com os resultados do modelo tubular.
            pul_data_solid (dict): Dicionário com os resultados do modelo sólido equivalente.
            model_tubular (MulticonductorTransmissionLine): Objeto do modelo MTL tubular.
        """
        self.pul_data = pul_data
        self.pul_data_tubular = pul_data_tubular
        self.model = model
        self.f = pul_data_tubular['frequencies']
        self.w = 2 * np.pi * self.f

        # Skin Depth (m)
        self.skin_depth = 1 / np.sqrt(self.model.mu * np.pi * self.f * self.model.sigma)

        # Plotter Parameters
        self.figsize = (12, 5)

    def internal_solid_conductors(self, p=0):
        """
        Gera o gráfico das características de impedância interna de um condutor
        cilíndrico sólido com base nos dados de simulação vetorizados.

        Args:
            p (int): O índice do condutor a ser analisado (padrão é 0).
        """
        conductor = self.model.mtl[p+1]
        ro = conductor['radius'][1] if isinstance(conductor['radius'], list) else conductor['radius']

        # Extrai os parâmetros DC (matrizes 2D) e pega o valor diagonal para o condutor 'p'
        Ri_cc = self.pul_data['internal']['Ri_cc'][p, p].real
        Li_cc = self.pul_data['internal']['Li_cc'][p, p].real

        # Extrai a impedância de Bessel (matriz 3D) e fatia para obter o vetor do condutor 'p'
        # A fatia [:, p, p] pega o elemento da diagonal (p, p) para todas as frequências (:)
        Zi = self.pul_data['internal']['Zi_bessel'][:, p, p]
        
        # Calcula a indutância interna a partir da reatância
        Li = np.divide(Zi.imag, self.w, out=np.zeros_like(self.w), where=self.w != 0)

        # Os cálculos das razões já são vetorizados
        R_ratio = Zi.real / Ri_cc
        wL_R_ratio = Zi.imag / Ri_cc
        L_ratio = Li / Li_cc
        wL_div_R = np.divide(Zi.imag, Zi.real, out=np.zeros_like(Zi.imag), where=Zi.real != 0)

        # --- Configuração do Gráfico (sem alteração na lógica) ---
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Internal parameters for bare-wire conductor with $\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm', fontsize=13, y=0.97)

        x_axis = ro / self.skin_depth
        ax.plot(x_axis, R_ratio,    'k-', lw=1, label=r"$R_i / R_{i(cc)}$")
        ax.plot(x_axis, wL_R_ratio, 'b-', lw=1, label=r"$\omega L_i / R_{i(cc)}$")
        ax.plot(x_axis, wL_div_R,   'g-', lw=1, label=r"$\omega L_i / R_i$")
        ax.plot(x_axis, L_ratio,    'r-', lw=1, label=r"$L_i / L_{i(cc)}$")

        ax.set_xscale('log')
        ax.set_xlabel(r'$(r_o / \delta)$', fontsize=12)
        ax.set_xlim(left=1e-1, right=max(x_axis))
        ax.set_ylim(0, max(np.max(R_ratio), np.max(wL_R_ratio)) * 0.15)
        ax.grid(True, which="both", ls="--", color='0.7')
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(fontsize=12, frameon=True)
        plt.tight_layout(rect=[0, 0, 1, 1])

    def internal_impedance(self, p=0):
        """
        Este método gera gráficos comparativos da impedância interna de um condutor
        usando a formulação exata, a aproximação de Nahman e Holt, e uma terceira
        aproximação.

        Os gráficos gerados são:
        1. Módulo da impedância interna |Z'_i| vs. Frequência.
        2. Ângulo da impedância interna arg(Z'_i) vs. Frequência (em graus).
        """
        # Pega as propriedades do condutor especificado 'p'
        conductor = self.model.mtl[p+1]
        ro = conductor['radius'][1] if isinstance(conductor['radius'], list) else conductor['radius']
        
        # Extrai os parâmetros DC (matrizes 2D) e pega o valor diagonal para o condutor 'p'
        Zi_bessel = self.pul_data['internal']['Zi_bessel'][:, p, p]
        Zi_kelvin = self.pul_data['internal']['Zi_kelvin'][:, p, p]    

        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize)
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Comparação de Modelos de Impedância Interna ($\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm)', fontsize=14, y=0.98)

        # --- Gráfico 1: Módulo ---
        ax1.set_title('Módulo da Impedância')
        ax1.plot(self.f, np.abs(Zi_bessel), 'k-', label='Bessel')
        ax1.plot(self.f, np.abs(Zi_kelvin), 'r--', lw=1, label='Kelvin')
        ax1.set_xscale('log')
        ax1.set_xlim(left=min(self.f), right=max(self.f))
        # ax1.set_ylim(0, np.max(mod_bessel) * 1.1)
        ax1.set_xlabel('Frequência (Hz)', fontsize=12)
        ax1.set_ylabel(r"Módulo $|Z_i|$ ($\Omega / m$)", fontsize=12)
        ax1.grid(True, which="both", ls=":", color='0.7')
        ax1.legend(loc='upper left')

        # --- Gráfico 2: Ângulo ---
        ax2.set_title('Ângulo da Impedância')
        ax2.plot(self.f, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Bessel')
        ax2.plot(self.f, np.angle(Zi_kelvin, deg=True), 'r--', lw=1, label='Kelvin')
        ax2.set_xscale('log')
        ax2.set_xlim(left=min(self.f), right=max(self.f))
        # ax2.set_ylim(0, np.max(angle_bessel) * 1.1)
        ax2.set_xlabel('Frequência (Hz)', fontsize=12)
        ax2.set_ylabel(r"Ângulo $Z_i$ (Graus)", fontsize=12)
        ax2.grid(True, which="both", ls=":", color='0.7')
        ax2.legend(loc='upper left')
        plt.tight_layout(rect=[0, 0, 1, 0.95])

    def nahman_holt_comparison(self, p=0):
        """
        Este método gera gráficos comparativos da impedância interna de um condutor
        usando a formulação exata, a aproximação de Nahman e Holt, e uma terceira
        aproximação.

        Os gráficos gerados são:
        1. Módulo da impedância interna |Z'_i| vs. Frequência.
        2. Ângulo da impedância interna arg(Z'_i) vs. Frequência (em graus).
        """
        # Pega as propriedades do condutor especificado 'p'
        conductor = self.model.mtl[p+1]
        ro = conductor['radius'][1] if isinstance(conductor['radius'], list) else conductor['radius']
        
        # Extrai os parâmetros DC (matrizes 2D) e pega o valor diagonal para o condutor 'p'
        Zi_bessel = self.pul_data['internal']['Zi_bessel'][:, p, p]
        Zi_nahman = self.pul_data['internal']['Zi_nahman'][:, p, p]  
        Zi_approx = self.pul_data['internal']['Zi_approx'][:, p, p]  
        
        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize)
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Comparação de Modelos de Impedância Interna ($\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm)', fontsize=14, y=0.98)

        # --- Gráfico 1: Módulo ---
        ax1.set_title('Módulo da Impedância')
        ax1.plot(self.f, np.abs(Zi_bessel), 'k-', label='Exactly')
        ax1.plot(self.f, np.abs(Zi_nahman), 'r-', lw=1, label='Nahman e Holt')
        ax1.plot(self.f, np.abs(Zi_approx), 'g--', lw=1, label='Nahman e Holt (Modified)')
        ax1.set_xscale('log')
        ax1.set_xlim(left=min(self.f), right=max(self.f))
        ax1.set_ylim(0, np.max(np.abs(Zi_bessel)) * 1.1)
        ax1.set_xlabel('Frequência (Hz)', fontsize=12)
        ax1.set_ylabel(r"Módulo $|Z_i|$ ($\Omega / m$)", fontsize=12)
        ax1.grid(True, which="both", ls=":", color='0.7')
        ax1.legend(loc='upper left')

        # --- Gráfico 2: Ângulo ---
        ax2.set_title('Ângulo da Impedância')
        ax2.plot(self.f, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Exactly')
        ax2.plot(self.f, np.angle(Zi_nahman, deg=True), 'r-', lw=1, label='Nahman e Holt')
        ax2.plot(self.f, np.angle(Zi_approx, deg=True), 'g--', lw=1, label='Nahman e Holt (Modified)')
        ax2.set_xscale('log')
        ax2.set_xlim(left=min(self.f), right=max(self.f))
        ax2.set_ylim(0, np.max(np.angle(Zi_bessel, deg=True)) * 1.1)
        ax2.set_xlabel('Frequência (Hz)', fontsize=12)
        ax2.set_ylabel(r"Ângulo $Z_i$ (Graus)", fontsize=12)
        ax2.grid(True, which="both", ls=":", color='0.7')
        ax2.legend(loc='upper left')
        plt.tight_layout(rect=[0, 0, 1, 0.95])

    def internal_tubular_characteristics(self, p=0):
        """
        Gera um gráfico comparativo da impedância interna de um condutor tubular
        com a de um condutor sólido, usando dados pré-calculados.

        Args:
            p (int): O índice do condutor a ser analisado.
        """
        conductor = self.model.mtl[p+1]
        ri, ro = conductor['radius'] if isinstance(conductor['radius'], list) else (0, conductor['radius'])

        # Acessa a matriz 3D 'Zi_bessel' diretamente, sem a chave aninhada.
        Zi_tubular = self.pul_data_tubular['internal']['Zi_bessel'][:, p, p]
        Zi_solid = self.pul_data['internal']['Zi_bessel'][:, p, p]

        # --- 2. Cálculo das Razões e Plotagem ---

        R_ratio = np.divide(Zi_tubular.real, Zi_solid.real, 
                            out=np.ones_like(self.f), where=Zi_solid.real != 0)
        L_ratio = np.divide(Zi_tubular.imag, Zi_solid.imag,
                            out=np.ones_like(self.f), where=Zi_solid.imag != 0)
        
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(f'Internal parameters for tubular bare-wire conductor with\n'
                     fr'$\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.2f}$ mm, $r_i = {ri*1e3:.2f}$ mm', fontsize=12)

        x_axis = ro / self.skin_depth
        ax.plot(x_axis, R_ratio, 'k-', lw=1.5, label=r"$R_{i(\text{tubular})} / R_{i(\text{solid})}$")
        ax.plot(x_axis, L_ratio, 'r--', lw=1.5, label=r"$L_{i(\text{tubular})} / L_{i(\text{solid})}$")

        ax.set_xlabel(r'$(r_o / \delta)$', fontsize=12)
        ax.set_xlim(0, 8)
        ax.set_ylim(0, 2)
        ax.grid(True, which="both", ls="--", color='0.7')
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(fontsize=12, frameon=True)
        plt.tight_layout(rect=[0, 0.02, 1, 0.95])

    