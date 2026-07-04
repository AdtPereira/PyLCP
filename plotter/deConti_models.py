import os
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from utils.case_utils import *
from mtl_main.source import MulticonductorTransmissionLine
from .models_base import BasePlotter

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
            # ... (helper function, unchanged)
            lines = []
            s_rows = [[f"{val:11.4e}" for val in row] for row in matrix]
            for row in s_rows:
                lines.append("  ".join(row))
            return "\n".join(lines)

        def _print_real_matrix(name: str, matrix_data: np.ndarray, unit: str):
            # ... (helper function, unchanged)
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
        
        # Iterate over the discrete frequencies and their indices
        for i, freq in enumerate(discrete_pul_data.get('frequencies', [])):
            print("\n" + "="*80)
            print(f"Frequency: {freq:,.0f} Hz")
            print("="*80)
            w = 2 * np.pi * freq

            # Iterate over the matrices (Zs, Ysh, etc.)
            for name, matrix_key in param_mapping.items():
                if matrix_key in discrete_pul_data:
                    matrix_3d = discrete_pul_data[matrix_key]
                    matrix = matrix_3d[i, :, :]

                    # (Matrix printing logic, unchanged)
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
                        # ... (rest of the L and C logic, unchanged)
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

class InternalLinesModels(BasePlotter):
    """
    Class to generate plots from the vectorized results of the
    overhead-line module.
    """
    def __init__(self, file_path: str, pul_data: dict, pul_data_tubular: dict,
                 model: MulticonductorTransmissionLine, model_tubular: MulticonductorTransmissionLine = None,
                 comsol: dict = None, autoSave: bool = True):
        """
        Initializes the plotter with the tubular and solid simulation data.

        Args:
            file_path (str): Path of the calling script (used to locate Results/).
            pul_data (dict): Dictionary with the solid model's results.
            pul_data_tubular (dict): Dictionary with the tubular model's results.
            model (MulticonductorTransmissionLine): Solid MTL model object.
            model_tubular (MulticonductorTransmissionLine, optional): Tubular/hollow
                MTL model object, used to extract the correct (ri, ro) in the
                hollow-conductor-specific plots. If omitted, uses `model`'s geometry.
            comsol (dict, optional): Internal impedance measured via COMSOL
                (see ComsolPostProcessor.get_bare_and_hollow_wire_internal_impedance).
                'Zi_measured' refers to the solid conductor (overlay in
                internal_impedance()/internal_solid_conductors()) and 'Zi_hollow' to
                the hollow conductor (overlay in hollow_conductor_impedance()/
                internal_hollow_conductors()).
        """
        super().__init__(file_path, pul_data, plot_config={}, autoSave=autoSave)
        self.pul_data_tubular = pul_data_tubular
        self.model = model
        self.model_tubular = model_tubular
        self.comsol = comsol
        self.f = pul_data_tubular['frequencies']
        self.w = 2 * np.pi * self.f

        # Skin depth (m) — depends only on the material's mu/sigma, which are
        # the same for the solid and tubular models (only the radius geometry differs).
        self.skin_depth = 1 / np.sqrt(self.model.mu * np.pi * self.f * self.model.sigma)

        # Plotter Parameters
        self.figsize = (12, 5)

    def _conductor_radius(self, p, tubular=False):
        """Returns (ri, ro) of conductor 'p' (ri=0 if solid)."""
        model = self.model_tubular if (tubular and self.model_tubular is not None) else self.model
        conductor = model.mtl[p + 1]
        radius = conductor['radius']
        return (radius[0], radius[1]) if isinstance(radius, list) else (0.0, radius)

    def internal_solid_conductors(self, p=0):
        """
        Generates the internal-impedance characteristics plot for a solid
        cylindrical conductor, based on the vectorized simulation data.

        Args:
            p (int): Index of the conductor to analyze (default is 0).
        """
        _, ro = self._conductor_radius(p)

        # Extract the DC parameters (2D matrices) and take the diagonal value for conductor 'p'
        Ri_cc = self.pul_data['internal']['Ri_cc'][p, p].real
        Li_cc = self.pul_data['internal']['Li_cc'][p, p].real

        # Extract the Bessel impedance (3D matrix) and slice to get the vector for conductor 'p'
        # The [:, p, p] slice takes the diagonal (p, p) element for all frequencies (:)
        Zi = self.pul_data['internal']['Zi_bessel'][:, p, p]

        # Compute the internal inductance from the reactance
        Li = np.divide(Zi.imag, self.w, out=np.zeros_like(self.w), where=self.w != 0)

        # The ratio calculations are already vectorized
        R_ratio = Zi.real / Ri_cc
        wL_R_ratio = Zi.imag / Ri_cc
        L_ratio = Li / Li_cc
        wL_div_R = np.divide(Zi.imag, Zi.real, out=np.zeros_like(Zi.imag), where=Zi.real != 0)

        # --- Plot Configuration (logic unchanged) ---
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Internal parameters for bare-wire conductor with $\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm', fontsize=13, y=0.97)

        x_axis = ro / self.skin_depth
        ax.plot(x_axis, R_ratio,    'k-', lw=1, label=r"$R_i / R_{i(cc)}$")
        ax.plot(x_axis, wL_R_ratio, 'b-', lw=1, label=r"$\omega L_i / R_{i(cc)}$")
        ax.plot(x_axis, wL_div_R,   'g-', lw=1, label=r"$\omega L_i / R_i$")
        ax.plot(x_axis, L_ratio,    'r-', lw=1, label=r"$L_i / L_{i(cc)}$")

        if self.comsol is not None:
            f_cmsl = self.comsol['frequencies']
            w_cmsl = 2 * np.pi * f_cmsl
            skin_depth_cmsl = 1 / np.sqrt(self.model.mu * np.pi * f_cmsl * self.model.sigma)
            x_cmsl = ro / skin_depth_cmsl
            Zi_cmsl = self.comsol['Zi_measured']
            Li_cmsl = np.divide(Zi_cmsl.imag, w_cmsl, out=np.zeros_like(w_cmsl), where=w_cmsl != 0)
            wL_div_R_cmsl = np.divide(Zi_cmsl.imag, Zi_cmsl.real, out=np.zeros_like(Zi_cmsl.imag), where=Zi_cmsl.real != 0)

            scatter_kwargs = dict(marker='o', s=25, facecolors='none', zorder=10)
            ax.scatter(x_cmsl, Zi_cmsl.real / Ri_cc,   edgecolors='k', label='COMSOL', **scatter_kwargs)
            ax.scatter(x_cmsl, Zi_cmsl.imag / Ri_cc,   edgecolors='b', **scatter_kwargs)
            ax.scatter(x_cmsl, wL_div_R_cmsl,          edgecolors='g', **scatter_kwargs)
            ax.scatter(x_cmsl, Li_cmsl / Li_cc,        edgecolors='r', **scatter_kwargs)

        ax.set_xscale('log')
        ax.set_xlabel(r'$(r_o / \delta)$', fontsize=12)
        ax.set_xlim(left=1e-1, right=max(x_axis))
        ax.set_ylim(0, max(np.max(R_ratio), np.max(wL_R_ratio)) * 0.15)
        ax.grid(True, which="both", ls="--", color='0.7')
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(fontsize=12, frameon=True)
        plt.tight_layout(rect=[0, 0, 1, 1])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='internal_solid_conductors')

    def internal_impedance(self, p=0):
        """
        This method generates comparison plots of a conductor's internal
        impedance using the exact formulation, the Nahman and Holt
        approximation, and a third approximation.

        The generated plots are:
        1. Internal impedance magnitude |Z'_i| vs. Frequency.
        2. Internal impedance angle arg(Z'_i) vs. Frequency (in degrees).
        """
        # Get the properties of the specified conductor 'p'
        _, ro = self._conductor_radius(p)

        # Extract the DC parameters (2D matrices) and take the diagonal value for conductor 'p'
        Zi_bessel = self.pul_data['internal']['Zi_bessel'][:, p, p]
        Zi_kelvin = self.pul_data['internal']['Zi_kelvin'][:, p, p]

        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize)
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Internal Impedance Model Comparison ($\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm)', fontsize=14, y=0.98)

        # --- Plot 1: Magnitude ---
        ax1.set_title('Impedance Magnitude')
        ax1.plot(self.f, np.abs(Zi_bessel), 'k-', label='Bessel')
        ax1.plot(self.f, np.abs(Zi_kelvin), 'r--', lw=1, label='Kelvin')
        if self.comsol is not None:
            ax1.scatter(self.comsol['frequencies'], np.abs(self.comsol['Zi_measured']),
                        marker='o', s=12, facecolors='none', edgecolors='black', zorder=10, label='COMSOL')
        ax1.set_xscale('log')
        ax1.set_xlim(left=min(self.f), right=max(self.f))
        # ax1.set_ylim(0, np.max(mod_bessel) * 1.1)
        ax1.set_xlabel('Frequency (Hz)', fontsize=12)
        ax1.set_ylabel(r"Magnitude $|Z_i|$ ($\Omega / m$)", fontsize=12)
        ax1.grid(True, which="both", ls=":", color='0.7')
        ax1.legend(loc='upper left')

        # --- Plot 2: Angle ---
        ax2.set_title('Impedance Angle')
        ax2.plot(self.f, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Bessel')
        ax2.plot(self.f, np.angle(Zi_kelvin, deg=True), 'r--', lw=1, label='Kelvin')
        if self.comsol is not None:
            ax2.scatter(self.comsol['frequencies'], np.angle(self.comsol['Zi_measured'], deg=True),
                        marker='o', s=12, facecolors='none', edgecolors='black', zorder=10, label='COMSOL')
        ax2.set_xscale('log')
        ax2.set_xlim(left=min(self.f), right=max(self.f))
        # ax2.set_ylim(0, np.max(angle_bessel) * 1.1)
        ax2.set_xlabel('Frequency (Hz)', fontsize=12)
        ax2.set_ylabel(r"Angle $Z_i$ (Degrees)", fontsize=12)
        ax2.grid(True, which="both", ls=":", color='0.7')
        ax2.legend(loc='upper left')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='internal_impedance')

    def internal_hollow_conductors(self, p=0):
        """
        Generates the internal-impedance characteristics plot for a hollow
        (hollow/tubular) conductor, analogous to internal_solid_conductors(),
        based on the tubular model's vectorized simulation data.

        Args:
            p (int): Index of the conductor to analyze (default is 0).
        """
        ri, ro = self._conductor_radius(p, tubular=True)
        sigma_model = self.model_tubular if self.model_tubular is not None else self.model

        # Extract the DC parameters (2D matrices) and take the diagonal value for conductor 'p'
        Ri_cc = self.pul_data_tubular['internal']['Ri_cc'][p, p].real
        Li_cc = self.pul_data_tubular['internal']['Li_cc'][p, p].real

        # Extract the Bessel impedance (3D matrix) and slice to get the vector for conductor 'p'
        Zi = self.pul_data_tubular['internal']['Zi_bessel'][:, p, p]

        # Compute the internal inductance from the reactance
        Li = np.divide(Zi.imag, self.w, out=np.zeros_like(self.w), where=self.w != 0)

        # The ratio calculations are already vectorized
        R_ratio = Zi.real / Ri_cc
        wL_R_ratio = Zi.imag / Ri_cc
        L_ratio = Li / Li_cc
        wL_div_R = np.divide(Zi.imag, Zi.real, out=np.zeros_like(Zi.imag), where=Zi.real != 0)

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)
        sigma_str = format_scientific_notation(sigma_model.sigma[p])
        fig.suptitle(fr'Internal parameters for hollow bare-wire conductor with $\sigma={sigma_str}$ S/m, '
                     fr'$r_o={ro*1e3:.2f}$ mm, $r_i={ri*1e3:.2f}$ mm', fontsize=13, y=0.97)

        x_axis = ro / self.skin_depth
        ax.plot(x_axis, R_ratio,    'k-', lw=1, label=r"$R_i / R_{i(cc)}$")
        ax.plot(x_axis, wL_R_ratio, 'b-', lw=1, label=r"$\omega L_i / R_{i(cc)}$")
        ax.plot(x_axis, wL_div_R,   'g-', lw=1, label=r"$\omega L_i / R_i$")
        ax.plot(x_axis, L_ratio,    'r-', lw=1, label=r"$L_i / L_{i(cc)}$")

        if self.comsol is not None and self.comsol.get('Zi_hollow') is not None:
            f_cmsl = self.comsol['frequencies']
            w_cmsl = 2 * np.pi * f_cmsl
            skin_depth_cmsl = 1 / np.sqrt(self.model.mu * np.pi * f_cmsl * self.model.sigma)
            x_cmsl = ro / skin_depth_cmsl
            Zi_cmsl = self.comsol['Zi_hollow']
            Li_cmsl = np.divide(Zi_cmsl.imag, w_cmsl, out=np.zeros_like(w_cmsl), where=w_cmsl != 0)
            wL_div_R_cmsl = np.divide(Zi_cmsl.imag, Zi_cmsl.real, out=np.zeros_like(Zi_cmsl.imag), where=Zi_cmsl.real != 0)

            scatter_kwargs = dict(marker='o', s=25, facecolors='none', zorder=10)
            ax.scatter(x_cmsl, Zi_cmsl.real / Ri_cc,   edgecolors='k', label='COMSOL', **scatter_kwargs)
            ax.scatter(x_cmsl, Zi_cmsl.imag / Ri_cc,   edgecolors='b', **scatter_kwargs)
            ax.scatter(x_cmsl, wL_div_R_cmsl,          edgecolors='g', **scatter_kwargs)
            ax.scatter(x_cmsl, Li_cmsl / Li_cc,        edgecolors='r', **scatter_kwargs)

        ax.set_xscale('log')
        ax.set_xlabel(r'$(r_o / \delta)$', fontsize=12)
        ax.set_xlim(left=1e-1, right=max(x_axis))
        ax.set_ylim(0, max(np.max(R_ratio), np.max(wL_R_ratio)) * 0.15)
        ax.grid(True, which="both", ls="--", color='0.7')
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(fontsize=12, frameon=True)
        plt.tight_layout(rect=[0, 0, 1, 1])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='internal_hollow_conductors')

    def hollow_conductor_impedance(self, p=0):
        """
        This method generates comparison plots of the internal impedance of a
        hollow (hollow/tubular) conductor, analogous to internal_impedance(),
        using the exact (Bessel) formulation vs. COMSOL (Kelvin does not
        apply to tubular conductors).

        The generated plots are:
        1. Internal impedance magnitude |Z'_i| vs. Frequency.
        2. Internal impedance angle arg(Z'_i) vs. Frequency (in degrees).
        """
        ri, ro = self._conductor_radius(p, tubular=True)
        sigma_model = self.model_tubular if self.model_tubular is not None else self.model

        # Zi_kelvin is not plotted here: kelvin_impedance_solid_wires() is only
        # valid for solid conductors and returns 'inf' for tubular conductors (see
        # analytical_forms/overhead_lines.py), so the comparison is restricted to
        # Bessel (valid for both solid and tubular) vs. COMSOL.
        Zi_bessel = self.pul_data_tubular['internal']['Zi_bessel'][:, p, p]

        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize)
        sigma_str = format_scientific_notation(sigma_model.sigma[p])
        fig.suptitle(fr'Internal Impedance Model Comparison — Hollow Conductor '
                     fr'($\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.2f}$ mm, $r_i={ri*1e3:.2f}$ mm)', fontsize=14, y=0.98)

        comsol_hollow = self.comsol.get('Zi_hollow') if self.comsol is not None else None

        # --- Plot 1: Magnitude ---
        ax1.set_title('Impedance Magnitude')
        ax1.plot(self.f, np.abs(Zi_bessel), 'k-', label='Bessel')
        if comsol_hollow is not None:
            ax1.scatter(self.comsol['frequencies'], np.abs(comsol_hollow),
                        marker='o', s=12, facecolors='none', edgecolors='black', zorder=10, label='COMSOL')
        ax1.set_xscale('log')
        ax1.set_xlim(left=min(self.f), right=max(self.f))
        ax1.set_xlabel('Frequency (Hz)', fontsize=12)
        ax1.set_ylabel(r"Magnitude $|Z_i|$ ($\Omega / m$)", fontsize=12)
        ax1.grid(True, which="both", ls=":", color='0.7')
        ax1.legend(loc='upper left')

        # --- Plot 2: Angle ---
        ax2.set_title('Impedance Angle')
        ax2.plot(self.f, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Bessel')
        if comsol_hollow is not None:
            ax2.scatter(self.comsol['frequencies'], np.angle(comsol_hollow, deg=True),
                        marker='o', s=12, facecolors='none', edgecolors='black', zorder=10, label='COMSOL')
        ax2.set_xscale('log')
        ax2.set_xlim(left=min(self.f), right=max(self.f))
        ax2.set_xlabel('Frequency (Hz)', fontsize=12)
        ax2.set_ylabel(r"Angle $Z_i$ (Degrees)", fontsize=12)
        ax2.grid(True, which="both", ls=":", color='0.7')
        ax2.legend(loc='upper left')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='hollow_conductor_impedance')

    def nahman_holt_comparison(self, p=0):
        """
        This method generates comparison plots of a conductor's internal
        impedance using the exact formulation, the Nahman and Holt
        approximation, and a third approximation.

        The generated plots are:
        1. Internal impedance magnitude |Z'_i| vs. Frequency.
        2. Internal impedance angle arg(Z'_i) vs. Frequency (in degrees).
        """
        # Get the properties of the specified conductor 'p'
        _, ro = self._conductor_radius(p)

        # Extract the DC parameters (2D matrices) and take the diagonal value for conductor 'p'
        Zi_bessel = self.pul_data['internal']['Zi_bessel'][:, p, p]
        Zi_nahman = self.pul_data['internal']['Zi_nahman'][:, p, p]
        Zi_approx = self.pul_data['internal']['Zi_approx'][:, p, p]

        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize)
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Internal Impedance Model Comparison ($\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm)', fontsize=14, y=0.98)

        # --- Plot 1: Magnitude ---
        ax1.set_title('Impedance Magnitude')
        ax1.plot(self.f, np.abs(Zi_bessel), 'k-', label='Exactly')
        ax1.plot(self.f, np.abs(Zi_nahman), 'r-', lw=1, label='Nahman and Holt')
        ax1.plot(self.f, np.abs(Zi_approx), 'g--', lw=1, label='Nahman and Holt (Modified)')
        ax1.set_xscale('log')
        ax1.set_xlim(left=min(self.f), right=max(self.f))
        ax1.set_ylim(0, np.max(np.abs(Zi_bessel)) * 1.1)
        ax1.set_xlabel('Frequency (Hz)', fontsize=12)
        ax1.set_ylabel(r"Magnitude $|Z_i|$ ($\Omega / m$)", fontsize=12)
        ax1.grid(True, which="both", ls=":", color='0.7')
        ax1.legend(loc='upper left')

        # --- Plot 2: Angle ---
        ax2.set_title('Impedance Angle')
        ax2.plot(self.f, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Exactly')
        ax2.plot(self.f, np.angle(Zi_nahman, deg=True), 'r-', lw=1, label='Nahman and Holt')
        ax2.plot(self.f, np.angle(Zi_approx, deg=True), 'g--', lw=1, label='Nahman and Holt (Modified)')
        ax2.set_xscale('log')
        ax2.set_xlim(left=min(self.f), right=max(self.f))
        ax2.set_ylim(0, np.max(np.angle(Zi_bessel, deg=True)) * 1.1)
        ax2.set_xlabel('Frequency (Hz)', fontsize=12)
        ax2.set_ylabel(r"Angle $Z_i$ (Degrees)", fontsize=12)
        ax2.grid(True, which="both", ls=":", color='0.7')
        ax2.legend(loc='upper left')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='nahman_holt_comparison')

    def internal_tubular_characteristics(self, p=0):
        """
        Generates a plot comparing a tubular conductor's internal impedance
        with that of a solid conductor, using pre-calculated data.

        Args:
            p (int): Index of the conductor to analyze.
        """
        ri, ro = self._conductor_radius(p, tubular=True)
        sigma_model = self.model_tubular if self.model_tubular is not None else self.model

        # Access the 3D 'Zi_bessel' matrix directly, without the nested key.
        Zi_tubular = self.pul_data_tubular['internal']['Zi_bessel'][:, p, p]
        Zi_solid = self.pul_data['internal']['Zi_bessel'][:, p, p]

        # --- Ratio Calculation and Plotting ---

        R_ratio = np.divide(Zi_tubular.real, Zi_solid.real,
                            out=np.ones_like(self.f), where=Zi_solid.real != 0)
        L_ratio = np.divide(Zi_tubular.imag, Zi_solid.imag,
                            out=np.ones_like(self.f), where=Zi_solid.imag != 0)

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=self.figsize)
        sigma_str = format_scientific_notation(sigma_model.sigma[p])
        fig.suptitle(f'Internal parameters for tubular bare-wire conductor with\n'
                     fr'$\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.2f}$ mm, $r_i = {ri*1e3:.2f}$ mm', fontsize=12)

        x_axis = ro / self.skin_depth
        ax.plot(x_axis, R_ratio, 'k-', lw=1.5, label=r"$R_{i(\text{tubular})} / R_{i(\text{solid})}$")
        ax.plot(x_axis, L_ratio, 'r--', lw=1.5, label=r"$L_{i(\text{tubular})} / L_{i(\text{solid})}$")

        if (self.comsol is not None and self.comsol.get('Zi_hollow') is not None
                and self.comsol.get('Zi_measured') is not None):
            f_cmsl = self.comsol['frequencies']
            skin_depth_cmsl = 1 / np.sqrt(self.model.mu * np.pi * f_cmsl * self.model.sigma)
            x_cmsl = ro / skin_depth_cmsl
            Zi_hollow_cmsl = self.comsol['Zi_hollow']
            Zi_solid_cmsl = self.comsol['Zi_measured']

            R_ratio_cmsl = np.divide(Zi_hollow_cmsl.real, Zi_solid_cmsl.real,
                                      out=np.ones_like(f_cmsl), where=Zi_solid_cmsl.real != 0)
            L_ratio_cmsl = np.divide(Zi_hollow_cmsl.imag, Zi_solid_cmsl.imag,
                                      out=np.ones_like(f_cmsl), where=Zi_solid_cmsl.imag != 0)

            scatter_kwargs = dict(marker='o', s=25, facecolors='none', zorder=10)
            ax.scatter(x_cmsl, R_ratio_cmsl, edgecolors='k', label='COMSOL', **scatter_kwargs)
            ax.scatter(x_cmsl, L_ratio_cmsl, edgecolors='r', **scatter_kwargs)

        ax.set_xlabel(r'$(r_o / \delta)$', fontsize=12)
        ax.set_xlim(0, 8)
        ax.set_ylim(0, 2)
        ax.grid(True, which="both", ls="--", color='0.7')
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(fontsize=12, frameon=True)
        plt.tight_layout(rect=[0, 0.02, 1, 0.95])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='internal_tubular_characteristics')

    