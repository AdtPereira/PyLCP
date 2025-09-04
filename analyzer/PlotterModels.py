import numpy as np
import matplotlib.pyplot as plt

class PrysmianModels:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, pul_parameters, discrete_pul_parameters):
        self.pul_data = pul_parameters
        self.discrete_pul = discrete_pul_parameters
        self.f = pul_parameters['frequencies']
        self.w = 2 * np.pi * self.f
        self.rho_1 = 1000
        self.epsr_1 = 1.0

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

    def _plot_generic_impedance(self, config_key, impedance_data):
        """Generic plotting function for impedance-like data."""
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series_config in config['series_to_plot']:
            z = impedance_data[series_config['key']]
            r = np.real(z)
            l = np.imag(z) / self.w * 1e6  # Inductance in µH/m
            style = {k: series_config[k] for k in ['label', 'color', 'linestyle', 'linewidth']}
            
            ax1.plot(self.f, r, **style)
            ax2.plot(self.f, l, **style)

        # Configure axes
        ax1.set(xscale='log', xlabel='Frequency (Hz)', ylabel=r'Resistance ($\Omega$/m)', title=config['resistance_title'])
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.legend(fontsize='small')
        
        ax2.set(xscale='log', xlabel='Frequency (Hz)', ylabel=r'Inductance ($\mu$H/m)', title=config['inductance_title'])
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.legend(fontsize='small')

        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def core_parameters(self, config_key='core'):
        # Access the dictionary of internal parameter arrays directly.
        internal_params = self.pul_data['internal']
        z11 = internal_params['zcs']['z11']
        z12 = internal_params['zcs']['z12']
        z2i = internal_params['zcs']['z2i']

        impedance_data = {
            'zcs': z11 + z12 + z2i,
            'z11': z11,
            'z12': z12,
            'z2i': z2i,
        }
        self._plot_generic_impedance(config_key, impedance_data)

    def sheath_parameters(self, config_key='sheath'):
        internal_params = self.pul_data['internal']
        z20 = internal_params['zs3']['z20']
        z23 = internal_params['zs3']['z23']
        z2m = internal_params['z2m']
        
        impedance_data = {
            'zs3': z20 + z23,
            'z20': z20,
            'z23': z23,
            'z2m': z2m,
        }
        self._plot_generic_impedance(config_key, impedance_data)

    def core_sheath_internal_impedance_matrix(self, config_key='core-sheath'):
        internal_params = self.pul_data['internal']
        z11 = internal_params['zcs']['z11']
        z12 = internal_params['zcs']['z12']
        z2i = internal_params['zcs']['z2i']
        z20 = internal_params['zs3']['z20']
        z23 = internal_params['zs3']['z23']
        z2m = internal_params['z2m']
        
        zcs = z11 + z12 + z2i
        zs3 = z20 + z23

        impedance_data = {
            'Zcc': zcs + zs3 - 2 * z2m,
            'Zss': zs3,
            'Zcs': zs3 - z2m,
        }
        self._plot_generic_impedance(config_key, impedance_data)

    def core_sheath_internal_parameters(self, config_key='internal-parameters'):
        internal_params = self.pul_data['internal']
        impedance_data = {
            'z11': internal_params['zcs']['z11'],
            'z2m': internal_params['z2m'],
            'z2i': internal_params['zcs']['z2i'],
            'z20': internal_params['zs3']['z20'],
        }
        self._plot_generic_impedance(config_key, impedance_data)

    def _ground_return_impedance_subplots(self, config_key):
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zg_3d = self.pul_data[series['key']]['earth-return_impedance_matrix']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': series['linewidth']}
            ax1.plot(self.f, np.real(zg_3d[:, p, q]), **style)
            ax2.plot(self.f, np.imag(zg_3d[:, p, q]) / self.w * 1e6, **style) # Inductance in µH/m

        ax1.set(xscale='log', xlim=config['x_lim']['resistance'], xlabel='Frequency (Hz)', ylabel=r'$R_s \, (\Omega/m)$', title=config['resistance_title'])
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.legend(fontsize='small')
        
        ax2.set(xscale='log', xlim=config['x_lim']['inductance'], xlabel='Frequency (Hz)', ylabel=r'$L_s \, (\mu H/m)$', title=config['inductance_title'])
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.legend(fontsize='small')
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
    def ground_return_impedance(self):
        self._ground_return_impedance_subplots('ground_return_impedance')

    # def log_matricial_pul_parameters(self, scale_units=True):
    #     """
    #     Generates and prints a terminal log report of the calculated PUL parameters.
    #     It displays complex matrices for impedance [Z] and admittance [Y].
    #     Additionally, it calculates and displays the real matrices for 
    #     inductance [L] (from [Zs]) and capacitance [C] (from [Ysh]).

    #     Args:
    #         scale_units (bool): If True, scales inductance to microhenries (µH/m)
    #                             and capacitance to picofarads (pF/m) for better
    #                             readability. Defaults to True.
    #     """
        
    #     def _matrix_to_string(matrix: np.ndarray) -> str:
    #         """
    #         Formats a NumPy matrix into a multi-line string for display.
    #         """
    #         lines = []
    #         s_rows = [[f"{val:11.4e}" for val in row] for row in matrix]
    #         for row in s_rows:
    #             lines.append("  ".join(row))
    #         return "\n".join(lines)

    #     def _print_real_matrix(name: str, matrix_data: np.ndarray, unit: str):
    #         """
    #         Helper function to format and print a real-valued matrix.
    #         """
    #         lines = _matrix_to_string(matrix_data).split('\n')
    #         num_rows = len(lines)
    #         middle_row_idx = num_rows // 2
            
    #         print("") 
    #         for i in range(num_rows):
    #             name_part = f"{name} = " if i == middle_row_idx else " " * (len(name) + 3)
    #             unit_part = f" {unit}" if i == middle_row_idx else ""
    #             print(f"{name_part}{lines[i]}{unit_part}")
    #         print("") 

    #     param_mapping = {
    #         '[Zi]': 'internal_impedance_matrix',
    #         '[Z0]': 'earth-return_impedance_matrix',
    #         '[Yg]': 'earth-return_admittance_matrix',
    #         '[Zs]': 'series_impedance_matrix',
    #         '[Ysh]': 'shunt_admittance_matrix',
    #     }

    #     print("\n--- Per-Unit-Length (PUL) Parameters Report ---")
        
    #     # Iterate through the frequencies and their corresponding data
    #     for freq, params in sorted(self.discrete_pul.items()):
    #         print("\n" + "="*80)
    #         print(f"Frequency: {freq:,.0f} Hz")
    #         print("="*80)
            
    #         # Iterate over the defined matrices
    #         w = 2 * np.pi * freq            
    #         for name, key in param_mapping.items():
    #             if key in params:
    #                 matrix = params[key]
                    
    #                 # --- Print the primary complex matrix (Z or Y) ---
    #                 if 'Z' in name: unit = "[Ohm/m]"
    #                 elif 'Y' in name: unit = "[S/m]"
    #                 else: unit = ""

    #                 r_lines = _matrix_to_string(matrix.real).split('\n')
    #                 x_lines = _matrix_to_string(matrix.imag).split('\n')
    #                 num_rows = len(r_lines)
    #                 middle_row_idx = num_rows // 2
                    
    #                 print("") 
    #                 for i in range(num_rows):
    #                     name_part = f"{name} = " if i == middle_row_idx else " " * (len(name) + 3)
    #                     op_part = " + j ".center(7) if i == middle_row_idx else " " * 7
    #                     unit_part = f" {unit}" if i == middle_row_idx else ""
    #                     print(f"{name_part}{r_lines[i]}{op_part}{x_lines[i]}{unit_part}")
    #                 print("") 
                    
    #                 # --- Conditionally print the derived real matrix (L or C) ---
    #                 if name == '[Zs]':
    #                     if scale_units:
    #                         l_matrix = (matrix.imag / w) 
    #                         l_unit = '[H/m]'
    #                     else:
    #                         l_matrix = matrix.imag / w
    #                         l_unit = '[H/m]'
    #                     _print_real_matrix('[L]', l_matrix, l_unit)

    #                 if name == '[Ysh]':
    #                     if scale_units:
    #                         c_matrix = (matrix.imag / w) 
    #                         c_unit = '[F/m]'
    #                     else:
    #                         c_matrix = matrix.imag / w
    #                         c_unit = '[F/m]'
    #                     _print_real_matrix('[C]', c_matrix, c_unit)

    #     print("\n" + "="*80)

    def log_matricial_pul_parameters(self, scale_units=True):
        """
        Generates and prints a terminal log report from the discrete PUL parameters.
        It iterates through the specified frequencies and scenarios.
        """
        def _matrix_to_string(matrix: np.ndarray) -> str:
            # ... (função auxiliar sem alteração) ...
            lines = []
            s_rows = [[f"{val:11.4e}" for val in row] for row in matrix]
            for row in s_rows:
                lines.append("  ".join(row))
            return "\n".join(lines)

        def _print_real_matrix(name: str, matrix_data: np.ndarray, unit: str):
            # ... (função auxiliar sem alteração) ...
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

        print("\n--- Per-Unit-Length (PUL) Parameters Report (Discrete Frequencies) ---")
        
        if not self.discrete_pul:
            print("No discrete data provided for logging.")
            return

        # Itera sobre as frequências discretas (100, 10k, 100k)
        for freq, scenarios_at_freq in sorted(self.discrete_pul.items()):
            print("\n" + "="*80)
            print(f"Frequency: {freq:,.0f} Hz")
            print("="*80)
            w = 2 * np.pi * freq

            # Itera sobre os cenários (p100, deconti, etc.) para essa frequência
            for scenario_key, params in scenarios_at_freq.items():
                print(f"\n--- Scenario: {scenario_key} ---")
                
                # Itera sobre as matrizes (Zs, Ysh, etc.) para esse cenário
                for name, matrix_key in param_mapping.items():
                    if matrix_key in params:
                        matrix = params[matrix_key]
                        
                        # (Lógica de impressão das matrizes, sem alteração)
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
                        
                        if name == '[Zs]':
                            l_matrix = (matrix.imag / w)
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
        
class DeContiModels:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, pul_parameters):
        self.pul_data = pul_parameters
        self.f = pul_parameters['frequencies']
        self.w = 2 * np.pi * self.f

        self.comparison = [
            {'key': 'magalhaes_xue', 'label': 'Magalhães/Xue',  'color': 'black', 'linestyle': '-', 'linewidth': 2},
            {'key': 'sunde', 'label': 'Sunde', 'color': 'gray', 'linestyle': ':', 'linewidth': 2},
            {'key': 'pollaczek', 'label': 'Pollaczek', 'color': 'darkblue', 'linestyle': '-', 'linewidth': 2},
            {'key': 'ametani', 'label': 'Ametani', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 2},
            {'key': 'deconti', 'label': 'De Conti', 'color': 'red', 'linestyle': ':', 'linewidth': 1},
            {'key': 'saad', 'label': 'Saad', 'color': 'red', 'linestyle': ':', 'linewidth': 1},
            {'key': 'wedepohl', 'label': 'Wedepohl and Wilcox', 'color': 'darkgray', 'linestyle': '-', 'linewidth': 1},
        ]
        
        self.paper_2023 = [
            {'key': 'p100',   'label': r'$\rho_e=100 \;\Omega m$',  'color': 'black', 'linestyle': '-', 'marker': 'o'},
            {'key': 'p1000',  'label': r'$\rho_e=1000 \;\Omega m$', 'color': 'black', 'linestyle': '-', 'marker': 's'},
            {'key': 'p10000', 'label': r'$\rho_e=10000 \;\Omega m$','color': 'black', 'linestyle': '-', 'marker': '^'},
            {'key': 'p100_deConti', 'label': '', 'color': 'red', 'linestyle': '--', 'marker': None},
            {'key': 'p1000_deConti', 'label': '', 'color': 'red', 'linestyle': '--', 'marker': None},
            {'key': 'p10000_deConti', 'label': '', 'color': 'red', 'linestyle': '--', 'marker': None},
        ]

        self.plot_configs = {
            'ground_return_impedance': {
                'suptitle': r'P.u.l. ground-return impedance for single buried bare-wire for $\rho_e=100 \;\Omega$ m and $\epsilon_r=10$',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'p': 0, 'q': 0,
                'x_lim': {'resistance': (1E5, 1E8), 'inductance': (1E-1, 1E8)},
                'y_lim': {'resistance': (0, 200), 'inductance': (0.5, 3.0)},
                'series_to_plot': self.comparison
            },
            'impedance': {
                'suptitle': r"Flat arrangement's mutual ground-return impedance for $\varepsilon_{r1} = 10$",
                'norm_title': 'Absolute Impedance',
                'angle_title': 'Angle of Impedance',
                'p': 2,
                'q': 0,
                'x_lim': {'norm': (1E4, 1E7), 'angle': (1E4, 1E7)},
                'y_lim': {'norm': (0, 25), 'angle': (20, 90)},
                'series_to_plot': self.paper_2023
            },
            'potential': {
                'suptitle': r"Flat arrangement's mutual ground-return potential coefficients for $\varepsilon_{r1} = 10$",
                'norm_title': 'Absolute Value',
                'angle_title': 'Angle of Potential Coefficient',
                'p': 2,
                'q': 0,
                'x_lim': {'norm': (1E4, 1E7), 'angle': (1E4, 1E7)},
                'y_lim': {'norm': (0, 14), 'angle': (-90, 90)},
                'series_to_plot': self.paper_2023
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
            zg_3d = self.pul_data[series['key']]['earth-return_impedance_matrix']
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': series['linewidth']}
            ax1.plot(self.f, np.real(zg_3d[:, p, q]) , **style)
            ax2.plot(self.f, np.imag(zg_3d[:, p, q]) / (2 * np.pi * self.f) * 1e6, **style) # Inductance in mH/km

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
        # ax2.set_ylim(config['y_lim']['inductance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_s \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _impedance_subplots(self, config_key):
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zg_3d = self.pul_data[series['key']]['earth-return_impedance_matrix']
            zg = zg_3d[:, p, q]
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'],
                     'marker': series['marker'], 'markersize': 3, 'linewidth': 1.0}
            ax1.plot(self.f, np.abs(zg), **style)
            ax2.plot(self.f, np.angle(zg, deg=True), **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['norm'])
        ax1.set_ylim(config['y_lim']['norm'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$|Zg_{{{p+1}{q+1}}}| \, (\Omega/m)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['norm_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['angle'])
        ax2.set_ylim(config['y_lim']['angle'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'Angle of $Zg_{{{p+1}{q+1}}}$ (Degrees)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['angle_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _potential_subplots(self, config_key):
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            pg_3d = self.pul_data[series['key']]['earth-return_potential_coefficient']
            pg = pg_3d[:, p, q]
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'],
                     'marker': series['marker'], 'markersize': 3, 'linewidth': 1.0}
            ax1.plot(self.f, np.abs(pg) * 1e-9, **style)
            ax2.plot(self.f, np.angle(pg, deg=True), **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['norm'])
        ax1.set_ylim(config['y_lim']['norm'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$|Pg_{{{p+1}{q+1}}}| \times 10^9 \, (\Omega m s^{{-1}})$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['norm_title'])

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['angle'])
        ax2.set_ylim(config['y_lim']['angle'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'Angle of $Pg_{{{p+1}{q+1}}}$ (Degrees)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['angle_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def impedance_comparison(self):
        """Plots the data corresponding to Figure 3 from the reference."""
        self._impedance_subplots('impedance')

    def potential_comparison(self):
        """Plots the data corresponding to Figure 6 from the reference."""
        self._potential_subplots('potential')

    def ground_return_impedance(self):
        self._ground_return_impedance_subplots('ground_return_impedance')

class XueModels:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    This version is adapted for the vectorized data structure.
    """
    def __init__(self, pul_parameters):
        self.pul_data = pul_parameters
        self.f = pul_parameters['frequencies']
        self.w = 2 * np.pi * self.f

        self.xue_series = [
            {'key': 'p100_deconti',      'label': 'De Conti Approx.', 'color': 'red', 'linestyle': ':'},
            {'key': 'p100_er20_deconti', 'label': '', 'color': 'red', 'linestyle': ':'},
            {'key': 'p500_deconti',      'label': '', 'color': 'red', 'linestyle': ':'},
            {'key': 'p100',      'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20', 'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--'},
            {'key': 'p500',      'label': r'$\rho_e=500 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-.'},
        ]

        self.nakagawa_carson_series = [
            {'key': 'p100',         'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$ (Nakagawa)',   'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20',    'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$ (Nakagawa)',  'color': 'black', 'linestyle': '-.'},
            {'key': 'p100_carson',  'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$ (Carson)',     'color': 'red',   'linestyle': '--'}
        ]

        self.attenuation_constant_series = [
            {'key': 'p100',         'label': r'$\rho_e=100 \;\Omega m$ (Nakagawa)', 'color': 'black',   'linestyle': '-'},
            {'key': 'p100_carson',  'label': r'$\rho_e=100 \;\Omega m$ (Carson)',   'color': 'red',     'linestyle': '-'},
            {'key': 'p2000',        'label': r'$\rho_e=2000 \;\Omega m$ (Nakagawa)','color': 'black',   'linestyle': '--'},
            {'key': 'p2000_carson', 'label': r'$\rho_e=2000 \;\Omega m$ (Carson)',  'color': 'red',     'linestyle': '--'},
        ]

        self.plot_configs = {
            'fig419': {
                'suptitle': 'Figure 4.19: P.u.l. Self-impedance of phase - a sheath with Magalhães/Xue formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, # Sheath of phase-a
                'q': 1, # Sheath of phase-a
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (0.5, 2.0)},
                'series_to_plot': self.xue_series
            },
            'fig421': {
                'suptitle': 'Figure 4.21: P.u.l. Mutual impedance between phase - a and phase - b sheaths with Magalhães/Xue formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1, # Sheath of phase-a
                'q': 3, # Sheath of phase-b
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (0.0, 1.5)},
                'series_to_plot': self.xue_series
            },
            'fig423': {
                'suptitle': 'Figure 4.23: P.u.l. Self-admittance of phase - a sheath with Magalhães/Xue formulation [1]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 1, # Sheath of phase-a
                'q': 1, # Sheath of phase-a
                'x_lim': {'conductance': (1E3, 1E7), 'capacitance': (1E3, 1E7)},
                'y_lim': {'conductance': (0.0, 20), 'capacitance': (0.0, 3.0)},
                'series_to_plot': self.xue_series
            },
        }

    def _plot_impedance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            # --- MODIFICATION START ---
            # Instead of a list comprehension, we now use direct NumPy slicing.
            # self.pul_data[series['key']] is a dictionary containing the 3D matrices.
            # We get the 3D matrix and slice it: [all_frequencies, row_p, col_q]
            zs_matrix_3d = self.pul_data[series['key']]['series_impedance_matrix']
            zs = zs_matrix_3d[:, p, q]
            # --- MODIFICATION END ---
            
            rs = np.real(zs) * 1e3  # Resistance in Ohm/km
            ls = np.imag(zs) / self.w * 1e6  # Inductance in mH/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': 1.0}
            ax1.plot(self.f, rs, **style)
            ax2.plot(self.f, ls, **style)

        # ... (the rest of the method for plotting axes remains unchanged) ...
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(config['x_lim']['resistance'])
        ax1.set_ylim(config['y_lim']['resistance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$Rs_{{{p+1}{q+1}}} \, (\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['inductance'])
        ax2.set_ylim(config['y_lim']['inductance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$Ls_{{{p+1}{q+1}}} \, (mH/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['inductance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _plot_admittance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            # --- MODIFICATION START ---
            # Using direct NumPy slicing for admittance as well.
            ysh_matrix_3d = self.pul_data[series['key']]['shunt_admittance_matrix']
            ysh = ysh_matrix_3d[:, p, q]
            # --- MODIFICATION END ---

            g_pq = np.real(ysh) * 1e3  # Conductance in S/km
            c_pq = np.imag(ysh) / self.w * 1e9  # Capacitance in uF/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, g_pq, **style)
            ax2.plot(self.f, c_pq, **style)

        # ... (the rest of the method for plotting axes remains unchanged) ...
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['conductance'])
        ax1.set_ylim(config['y_lim']['conductance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$G_{{{p+1}{q+1}}} \, (S/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['conductance_title'])
        
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['capacitance'])
        ax2.set_ylim(config['y_lim']['capacitance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$C_{{{p+1}{q+1}}} \, (\mu F/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['capacitance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def plot_fig419(self):
        """Plots the data corresponding to Figure 4.19 from the reference."""
        self._plot_impedance_subplots('fig419')

    def plot_fig421(self):
        """Plots the data corresponding to Figure 4.21 from the reference."""
        self._plot_impedance_subplots('fig421')

    def plot_fig423(self):
        """Plots the data corresponding to Figure 4.23 from the reference."""
        self._plot_admittance_subplots('fig423')