""" 
Este script executa simulações de impedância de linhas de transmissão coaxiais
usando tanto uma abordagem analítica (formulação de Patel) quanto uma abordagem numérica
(Método dos Momentos - MoM-SO). Ele gera gráficos comparativos dos resultados
obtidos por ambas as metodologias.

Este arquivo é parte do projeto PyLCP, que é um pacote Python para análise de linhas de transmissão.

REFERENCES:
[1] PATEL, Utkarsh R. A Surface Admittance Approach For Fast Calculation of the 
    Series Impedance of Cables Including Skin, Proximity, and Ground Return Effects.
    2014. University of Toronto, Graduate Department of The Edward S. Rogers Sr. 
    Department of Electrical & Computer Engineering. 

[2] U. R. Patel, B. Gustavsen and P. Triverio, "An Equivalent Surface Current Approach
    for the Computation of the Series Impedance of Power Cables with Inclusion of Skin
    and Proximity Effects," in IEEE Transactions on Power Delivery, vol. 28, no. 4, pp.
    2474-2482, Oct. 2013, doi: 10.1109/TPWRD.2013.2267098.

[3] U. R. Patel, B. Gustavsen and P. Triverio, "Application of the MoM-SO Method for 
    Accurate Impedance Calculation of Single-Core Cables Enclosed by a Conducting Pipe," 
    Proc. International Conference on Power Systems Transients (IPST 2013), Vancouver, 
    Canada July 18-20, 2013. https://www.ipstconf.org/Proc_IPST2013.php

[4] A. Ametani, "A General Formulation of Impedance and Admittance of Cables," in IEEE
    Transactions on Power Apparatus and Systems, vol. PAS-99, no. 3, pp. 902-910, May
    1980, doi: 10.1109/TPAS.1980.319718.

[5] A. Ametani, "Wave Propagation Characteristics of Cables," in IEEE Transactions on
    Power Apparatus and Systems, vol. PAS-99, no. 2, pp. 499-505, March 1980, 
    doi: 10.1109/TPAS.1980.319685.
"""
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# RAIZ DO PROJETO E DIRETÓRIOS
try:
    os.system('cls' if os.name == 'nt' else 'clear')
    script_dir = Path(__file__).resolve().parent
    print(f"Script directory: {script_dir}")
    project_root = script_dir.parents[1]
    print(f"Project root: {project_root}")
    sys.path.append(str(project_root))
    print("Caminhos do projeto configurados com sucesso.")
except IndexError:
    raise FileNotFoundError(
        "Não foi possível encontrar a raiz do projeto. "
        "Certifique-se de que o script está em 'examples/coated_wires'."
    )

# IMPORTAÇÕES DOS MÓDULOS E MODELO DE DADOS
try:
    from mtl_main.models_scc import PRYSMIAN_138kV_CORE_SHEATH as MODEL
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.scc import InternalPerUnitParameters, PerUnitParameters
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

class ModelPlotter:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, mtl_model, freq, pul, discrete_pul):
        self.freq_data = freq
        self.pul_data = pul
        self.discrete_pul = discrete_pul
        self.f = self.freq_data['Analytically']
        self.w = 2 * np.pi * self.f
        self.ro = mtl_model.surfaces[0]['radius']
        self.h1 = mtl_model.surfaces[0]['center_point'][1]
        self.rho_1 = 1/mtl_model.mtl_ref[0]['conductivity']
        self.epsr_1 = mtl_model.mtl_ref[0]['relative_permittivity']

        self.nakagawa_series = [
            {'key': 'magalhaes_xue', 'label': 'Magalhães/Xue',  'color': 'black', 'linestyle': '-', 'linewidth': 2},
            {'key': 'sunde', 'label': 'Sunde', 'color': 'gray', 'linestyle': ':', 'linewidth': 2},
            {'key': 'pollaczek', 'label': 'Pollaczek', 'color': 'darkblue', 'linestyle': '-', 'linewidth': 2},
            {'key': 'ametani', 'label': 'Ametani', 'color': 'cyan', 'linestyle': '-', 'linewidth': 1},
            {'key': 'deconti', 'label': 'De Conti', 'color': 'red', 'linestyle': '--', 'linewidth': 1},
            {'key': 'saad', 'label': 'Saad', 'color': 'red', 'linestyle': '--', 'linewidth': 1},
        ]
        
        self.plot_configs = {
            'core': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. parameters of core conductor, $z_{{cs}}$ [2]',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'series_to_plot': [
                    {'key': 'zcs', 'label': r'$z_{cs} = z_{11} + z_{12} + z_{2i}$', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
                    {'key': 'z11', 'label': r'$z_{11}$: Internal impedance of core outer surface', 'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.0},
                    {'key': 'z12', 'label': r'$z_{12}$: Core outer insulator impedance', 'color': 'darkblue',  'linestyle': '--', 'linewidth': 1.0},
                    {'key': 'z2i', 'label': r'$z_{2i}$: Internal impedance of sheath inner surface', 'color': 'red',   'linestyle': '-.', 'linewidth': 1.0},
                ]
            },
            'sheath': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. parameters of sheath conductor, $z_{{s3}} [2]$',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'series_to_plot': [
                    {'key': 'zs3', 'label': r'$z_{s3} = z_{20} + z_{23}$', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
                    {'key': 'z20', 'label': r'$z_{20}$: Internal impedance of sheath outer surface', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
                    {'key': 'z2m', 'label': r'$z_{2m}$: Sheath mutual impedance', 'color': 'darkblue',  'linestyle': '-.', 'linewidth': 1.0},
                    {'key': 'z23', 'label': r'$z_{23}$: Sheath outer insulator impedance', 'color': 'darkgreen',  'linestyle': ':', 'linewidth': 1.0},
                ]
            },
            'core-sheath': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. Internal Impedance Matrix, $[z_{{i}}]$ [2]',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'series_to_plot': [
                    {'key': 'Zcc', 'label': r'$z_{cc} = z_{cs} + z_{s3} - 2z_{2m}$: Core self-impedance ', 'color': 'black', 'linestyle': '-', 'linewidth': 1.0},
                    {'key': 'Zcs', 'label': r'$z_{cs} = z_{20} + z_{23} - z_{2m}$: Mutual impedance between the core and sheath ', 'color': 'darkblue',  'linestyle': '--', 'linewidth': 1.0},
                    {'key': 'Zss', 'label': r'$z_{ss} = z_{20} + z_{23}$: Sheath self-impedance', 'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.0},
                ]
            },
            'internal-parameters': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. parameters of SCC [2]',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'series_to_plot': [
                    {'key': 'z11', 'label': r'$z_{11}$: Internal impedance of core outer surface', 'color': 'black', 'linestyle': '-', 'linewidth': 1.0},
                    {'key': 'z2m', 'label': r'$z_{2m}$: Sheath mutual impedance', 'color': 'darkblue', 'linestyle': '-', 'linewidth': 1.0},
                    {'key': 'z2i', 'label': r'$z_{2i}$: Internal impedance of sheath inner surface', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 1.0},
                    {'key': 'z20', 'label': r'$z_{20}$: Internal impedance of sheath outer surface', 'color': 'darkgray', 'linestyle': '-', 'linewidth': 1.0},
                ]
            },
            'ground_return_impedance': {
                'suptitle': fr'Prysmian 138 kV SCC P.u.l. ground-return impedance for $\rho_1={self.rho_1:.0f} \;\Omega$ m and $\epsilon_{{r1}}={self.epsr_1:.0f}$',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'x_lim': {'resistance': (1E3, 1E7), 'inductance': (1E3, 1E7)},
                # 'y_lim': {'resistance': (0, 200), 'inductance': (0.5, 3.0)},
                'p': 0,
                'q': 0,
                'series_to_plot': self.nakagawa_series
            },
        }

    def _ground_return_impedance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zg_raw = np.array([item['earth-return_impedance_matrix'][p, q] for item in self.pul_data[series['key']]])
            rg = np.real(zg_raw) 
            lg = np.imag(zg_raw) / (2 * np.pi * self.f) * 1e6  # Inductance in mH/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': series['linewidth']}
            ax1.plot(self.f, rg, **style)
            ax2.plot(self.f, lg, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_s \, (\Omega/m)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_xlim(config['x_lim']['resistance'])
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_s \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_xlim(config['x_lim']['inductance'])
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def ground_return_impedance(self):
        self._ground_return_impedance_subplots('ground_return_impedance')

    def core_parameters(self, config_key='core'):
        """
        This method now pre-calculates the impedance components, including
        the summed 'zcs', and stores them in a dictionary for easy plotting.
        """
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        # Step 1: Pre-calculate base impedance arrays from simulation results.
        z11 = np.array([data['zcs']['z11'] for data in self.pul_data['internal']])
        z12 = np.array([data['zcs']['z12'] for data in self.pul_data['internal']])
        z2i = np.array([data['zcs']['z2i'] for data in self.pul_data['internal']])

        # Step 2: Store all impedance data in a dictionary for easy access.
        impedance_data = {
            'zcs': z11 + z12 + z2i,
            'z11': z11,
            'z12': z12,
            'z2i': z2i,
        }

        # Step 3: Iterate through the plotting configuration and plot from pre-calculated data.
        for series_config in config['series_to_plot']:
            z_values = impedance_data[series_config['key']]
            r_values = np.real(z_values) 
            l_values = np.imag(z_values) / self.w * 1e6  # Convert to uH/m
            
            # Define the plot style from the configuration.
            style = {
                'label': series_config['label'], 
                'color': series_config['color'], 
                'linestyle': series_config['linestyle'], 
                'linewidth': series_config['linewidth']
            }
            
            # Plot resistance and inductance on their respective subplots.
            ax1.plot(self.f, r_values, **style)
            ax2.plot(self.f, l_values, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'Resistance ($\Omega$/m)')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Inductance ($\mu$H/m)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def sheath_parameters(self, config_key='sheath'):
        """
        This method now pre-calculates the impedance components, including
        the summed 'zcs', and stores them in a dictionary for easy plotting.
        """
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        # Step 1: Pre-calculate base impedance arrays from simulation results.
        z20 = np.array([data['zs3']['z20'] for data in self.pul_data['internal']])
        z23 = np.array([data['zs3']['z23'] for data in self.pul_data['internal']])

        # Step 2: Store all impedance data in a dictionary for easy access.
        impedance_data = {
            'zs3': z20 + z23,
            'z20': z20,
            'z23': z23,
            'z2m': np.array([data['z2m'] for data in self.pul_data['internal']]),
        }

        # Step 3: Iterate through the plotting configuration and plot from pre-calculated data.
        for series_config in config['series_to_plot']:
            z_values = impedance_data[series_config['key']]
            r_values = np.real(z_values) 
            l_values = np.imag(z_values) / self.w * 1e6  # Convert to uH/m
            
            # Define the plot style from the configuration.
            style = {
                'label': series_config['label'], 
                'color': series_config['color'], 
                'linestyle': series_config['linestyle'], 
                'linewidth': series_config['linewidth']
            }
            
            # Plot resistance and inductance on their respective subplots.
            ax1.plot(self.f, r_values, **style)
            ax2.plot(self.f, l_values, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'Resistance ($\Omega$/m)')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Inductance ($\mu$H/m)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def core_sheath_internal_impedance_matrix(self, config_key='core-sheath'):
        """
        This method now pre-calculates the impedance components, including
        the summed 'zcs', and stores them in a dictionary for easy plotting.
        """
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        # Step 1: Pre-calculate base impedance arrays from simulation results.
        z11 = np.array([data['zcs']['z11'] for data in self.pul_data['internal']])
        z12 = np.array([data['zcs']['z12'] for data in self.pul_data['internal']])
        z2i = np.array([data['zcs']['z2i'] for data in self.pul_data['internal']])
        z20 = np.array([data['zs3']['z20'] for data in self.pul_data['internal']])
        z23 = np.array([data['zs3']['z23'] for data in self.pul_data['internal']])
        z2m = np.array([data['z2m'] for data in self.pul_data['internal']])
        zcs = z11 + z12 + z2i
        zs3 = z20 + z23

        # Step 2: Store all impedance data in a dictionary for easy access.
        impedance_data = {
            'Zcc': zcs + zs3 - 2 * z2m,
            'Zss': zs3,
            'Zcs': zs3 - z2m,
        }

        # Step 3: Iterate through the plotting configuration and plot from pre-calculated data.
        for series_config in config['series_to_plot']:
            z_values = impedance_data[series_config['key']]
            r_values = np.real(z_values) 
            l_values = np.imag(z_values) / self.w * 1e6  # Convert to uH/m
            
            # Define the plot style from the configuration.
            style = {
                'label': series_config['label'], 
                'color': series_config['color'], 
                'linestyle': series_config['linestyle'], 
                'linewidth': series_config['linewidth']
            }
            
            # Plot resistance and inductance on their respective subplots.
            ax1.plot(self.f, r_values, **style)
            ax2.plot(self.f, l_values, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'Resistance ($\Omega$/m)')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Inductance ($\mu$H/m)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def core_sheath_internal_parameters(self, config_key='internal-parameters'):
        """
        This method now pre-calculates the impedance components, including
        the summed 'zcs', and stores them in a dictionary for easy plotting.
        """
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        # Step 2: Store all impedance data in a dictionary for easy access.
        impedance_data = {
            'z11': np.array([data['zcs']['z11'] for data in self.pul_data['internal']]),
            'z2m': np.array([data['z2m'] for data in self.pul_data['internal']]),
            'z2i': np.array([data['zcs']['z2i'] for data in self.pul_data['internal']]),
            'z20': np.array([data['zs3']['z20'] for data in self.pul_data['internal']]),
        }

        # Step 3: Iterate through the plotting configuration and plot from pre-calculated data.
        for series_config in config['series_to_plot']:
            z_values = impedance_data[series_config['key']]
            r_values = np.real(z_values) 
            l_values = np.imag(z_values) / self.w * 1e6  # Convert to uH/m
            
            # Define the plot style from the configuration.
            style = {
                'label': series_config['label'], 
                'color': series_config['color'], 
                'linestyle': series_config['linestyle'], 
                'linewidth': series_config['linewidth']
            }
            
            # Plot resistance and inductance on their respective subplots.
            ax1.plot(self.f, r_values, **style)
            ax2.plot(self.f, l_values, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'Resistance ($\Omega$/m)')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Inductance ($\mu$H/m)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def log_matricial_pul_parameters(self, scale_units=True):
        """
        Generates and prints a terminal log report of the calculated PUL parameters.
        It displays complex matrices for impedance [Z] and admittance [Y].
        Additionally, it calculates and displays the real matrices for 
        inductance [L] (from [Zs]) and capacitance [C] (from [Ysh]).

        Args:
            scale_units (bool): If True, scales inductance to microhenries (µH/m)
                                and capacitance to picofarads (pF/m) for better
                                readability. Defaults to True.
        """
        
        def _matrix_to_string(matrix: np.ndarray) -> str:
            """
            Formats a NumPy matrix into a multi-line string for display.
            """
            lines = []
            s_rows = [[f"{val:11.4e}" for val in row] for row in matrix]
            for row in s_rows:
                lines.append("  ".join(row))
            return "\n".join(lines)

        def _print_real_matrix(name: str, matrix_data: np.ndarray, unit: str):
            """
            Helper function to format and print a real-valued matrix.
            """
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
            '[Zi]': 'internal_impedance_matrix',
            '[Z0]': 'earth-return_impedance_matrix',
            '[Yg]': 'earth-return_admittance_matrix',
            '[Zs]': 'series_impedance_matrix',
            '[Ysh]': 'shunt_admittance_matrix',
        }

        print("\n--- Per-Unit-Length (PUL) Parameters Report ---")
        
        # Iterate through the frequencies and their corresponding data
        for freq, params in sorted(self.discrete_pul.items()):
            print("\n" + "="*80)
            print(f"Frequency: {freq:,.0f} Hz")
            print("="*80)
            
            # Iterate over the defined matrices
            w = 2 * np.pi * freq            
            for name, key in param_mapping.items():
                if key in params:
                    matrix = params[key]
                    
                    # --- Print the primary complex matrix (Z or Y) ---
                    if 'Z' in name: unit = "[Ohm/m]"
                    elif 'Y' in name: unit = "[S/m]"
                    else: unit = ""

                    r_lines = _matrix_to_string(matrix.real).split('\n')
                    x_lines = _matrix_to_string(matrix.imag).split('\n')
                    num_rows = len(r_lines)
                    middle_row_idx = num_rows // 2
                    
                    print("") 
                    for i in range(num_rows):
                        name_part = f"{name} = " if i == middle_row_idx else " " * (len(name) + 3)
                        op_part = " + j ".center(7) if i == middle_row_idx else " " * 7
                        unit_part = f" {unit}" if i == middle_row_idx else ""
                        print(f"{name_part}{r_lines[i]}{op_part}{x_lines[i]}{unit_part}")
                    print("") 
                    
                    # --- Conditionally print the derived real matrix (L or C) ---
                    if name == '[Zs]':
                        if scale_units:
                            l_matrix = (matrix.imag / w) 
                            l_unit = '[H/m]'
                        else:
                            l_matrix = matrix.imag / w
                            l_unit = '[H/m]'
                        _print_real_matrix('[L]', l_matrix, l_unit)

                    if name == '[Ysh]':
                        if scale_units:
                            c_matrix = (matrix.imag / w) 
                            c_unit = '[F/m]'
                        else:
                            c_matrix = matrix.imag / w
                            c_unit = '[F/m]'
                        _print_real_matrix('[C]', c_matrix, c_unit)

        print("\n" + "="*80)

if __name__ == "__main__":
    st = time.time()
    frequency = {'Analytically': np.logspace(1, 7, num=80)}

    # Define models for different physical scenarios
    mtl_model = MulticonductorTransmissionLine(MODEL)

    # Define the calculation scenarios
    scenarios = {
        'internal': {'mtl': mtl_model},
        'magalhaes_xue': {'mtl': mtl_model, 'zg_form': 'magalhaes_xue'},
        'sunde': {'mtl': mtl_model, 'zg_form': 'sunde'},
        'pollaczek': {'mtl': mtl_model, 'zg_form': 'pollaczek'},
        'ametani': {'mtl': mtl_model, 'zg_form': 'ametani'},
        'deconti': {'mtl': mtl_model, 'zg_form': 'deconti'},
        'saad': {'mtl': mtl_model, 'zg_form': 'saad'},
    }
    
    pul_parameters = {key: [] for key in scenarios}
    for f in frequency['Analytically']:
        for key, params in scenarios.items():
            # Internal Parameters
            if key == 'internal':
                pul = InternalPerUnitParameters(params['mtl'], f)
                pul_parameters[key].append(pul.parameters_by_bessel())
            
            # Ground Return Parameters
            else:
                pul = PerUnitParameters(params['mtl'], f)
                pul_parameters[key].append(
                    pul.ground_return_parameters(zg_form=params['zg_form'])
                )

    discrete_pul_parameters = {}
    for f in [1e2, 1e4, 1e5]:
        pul = PerUnitParameters(mtl_model, f)
        internal = InternalPerUnitParameters(mtl_model, f)
        discrete_pul_parameters[f] = pul.quasi_tem_pul(
            internal.internal_matrices(), zg_form='ametani', yg_form='ametani')

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = ModelPlotter(mtl_model, frequency, pul_parameters, discrete_pul_parameters)
    plotter.core_parameters()
    plotter.sheath_parameters()   
    plotter.core_sheath_internal_impedance_matrix() 
    plotter.core_sheath_internal_parameters()
    plotter.ground_return_impedance()
    plotter.log_matricial_pul_parameters()
    MTLRepresentation(mtl_model, units='millimeter').ground_return_systems()
    plt.show()