import numpy as np
import matplotlib.pyplot as plt
from utils.case_utils import *

class PrysmianModels:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, pul_parameters, discrete_pul_data):
        self.pul_data = pul_parameters
        self.discrete_pul_data = discrete_pul_data
        self.discrete_frequencies = discrete_pul_data.get('frequencies', [])
        
        # --- Lógica para dados de plotagem (vetorizados) ---
        self.f = pul_parameters['frequencies']
        self.w = 2 * np.pi * self.f
        if 'internal' in pul_parameters:
            self.internal_data = self.pul_data['internal']

        self.rho_1 = 1000
        self.epsr_1 = 1.0
        self.figsize = (12, 5)

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
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
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
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
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

    def log_matricial_pul_parameters(self, scale_units=True):
        """
        Generates a terminal log report by slicing the discrete vectorized results.
        """
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
            '[Zi]': 'internal_impedance_matrix',
            '[Z0]': 'earth-return_impedance_matrix',
            '[Yg]': 'earth-return_admittance_matrix',
            '[Zs]': 'series_impedance_matrix',
            '[Ysh]': 'shunt_admittance_matrix',
        }

        print("\n--- Per-Unit-Length (PUL) Parameters Report (Discrete Frequencies) ---")
        
        if not self.discrete_frequencies:
            print("No discrete frequencies found for logging.")
            return

        # Itera sobre as frequências discretas e seus índices
        for i, freq in enumerate(self.discrete_frequencies):
            print("\n" + "="*80)
            print(f"Frequency: {freq:,.0f} Hz")
            print("="*80)
            w = 2 * np.pi * freq

            # Como os dados discretos representam um único cenário, não há laço de cenário
            # Itera sobre as matrizes (Zs, Ysh, etc.)
            for name, matrix_key in param_mapping.items():
                if matrix_key in self.discrete_pul_data:
                    # Fatia a matriz 3D para obter a matriz 2D da frequência atual
                    matrix_3d = self.discrete_pul_data[matrix_key]
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

    def compare_internal_impedance_matrix(self, config_key='core-sheath'):
        internal_params = self.pul_data['internal']
        Ri = self.pul_data['internal_matrices']['resistance_matrix']
        Li = self.pul_data['internal_matrices']['inductance_matrix']

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
        
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        ax1.plot(self.f, Ri[:, 0, 0], color='red', linestyle='--', linewidth=1.5, label='Ri')
        ax2.plot(self.f, Li[:, 0, 0], color='red', linestyle='--', linewidth=1.5, label='Li')
        ax1.plot(self.f, Ri[:, 1, 1], color='red', linestyle='--', linewidth=1.5, label='')
        ax2.plot(self.f, Li[:, 1, 1], color='red', linestyle='--', linewidth=1.5, label='')
        ax1.plot(self.f, Ri[:, 0, 1], color='red', linestyle='--', linewidth=1.5, label='')
        ax2.plot(self.f, Li[:, 0, 1], color='red', linestyle='--', linewidth=1.5, label='')
        
        for series_config in config['series_to_plot']:
            z = impedance_data[series_config['key']]
            style = {k: series_config[k] for k in ['label', 'color', 'linestyle', 'linewidth']}
            
            ax1.plot(self.f, np.real(z), **style)
            ax2.plot(self.f, np.imag(z) / self.w, **style)

        # Configure axes
        ax1.set(xscale='log', xlabel='Frequency (Hz)', ylabel=r'Resistance ($\Omega$/m)', title=config['resistance_title'])
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.legend(fontsize='small')
        
        ax2.set(xscale='log', xlabel='Frequency (Hz)', ylabel=r'Inductance (H/m)', title=config['inductance_title'])
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.legend(fontsize='small')

        plt.tight_layout(rect=[0, 0, 1, 0.96])
