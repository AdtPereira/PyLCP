import os
import sys
import time
import copy
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import scipy.constants as sc

# --- PROJECT ROOT AND DIRECTORIES CONFIGURATION ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    script_dir = Path(__file__).resolve().parent
    print(f"Script directory: {script_dir}")
    project_root = script_dir.parents[1]
    print(f"Project root: {project_root}")
    sys.path.append(str(project_root))
    print("Project paths configured successfully.")
except IndexError:
    raise FileNotFoundError(
        "Could not find the project root. "
        "Make sure the script is in 'examples/coated_wires'."
    )

# --- MODULES AND DATA MODEL IMPORTS ---
try:
    from mtl_main.models_wires import SINGLE_WIRE_R10_LIMA as MODEL
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_formulation.overhead_lines import PerUnitParameters
    print("Modules and data model imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

class LimaModelPlotter:
    """
    A highly refactored class to handle plotting for the Lima model results.
    It uses a configuration-driven approach to minimize code duplication,
    with specific settings for each plot type.
    """
    def __init__(self, freq, pul, p=0, q=0):
        self.freq_data = freq
        self.pul_data = pul
        self.p = p
        self.q = q
        
        self.f = self.freq_data['Analytically']
        self.w = 2 * np.pi * self.f
        
        title_suffix = (
            'of the single overhead line\n'
            r'$r_1 = 0.01\,\mathrm{m}, h_1 = 10\,\mathrm{m}, \sigma = 6.496 \times 10^7\,\mathrm{S/m}$ [2]'
        )
        
        # --- Centralized Plotting Configuration with Specific Overrides ---
        self.plot_configs = {
            'attenuation': {
                'title': f'Attenuation Constant {title_suffix}',
                'ylabel': r'Attenuation Constant, $\alpha$ (Np/km)',
                'quantity_key': 'gamma_v',
                'transform': lambda x: 1E3 * np.real(x),
                'xscale': 'log',
                'xlim': (1E2, 1E10),
                'single_plot': {
                    'models': ['quasi_tem', 'quasi_tem_log', 'nakagawa', 'sunde', 'carson'],
                    'ylim': (0, 4),
                    'grid': False
                },
                'soil_plot': {
                    'ylim': (0, 1.6)
                }
            },
            'phase_velocity': {
                'title': f'Phase Velocity {title_suffix}',
                'ylabel': r'Phase Velocity ($v_p/c$)',
                'quantity_key': 'gamma_v',
                'transform': lambda x: self.w / np.imag(x) / sc.c,
                'xscale': 'log',
                'xlim': (1E0, 1E10),
                'single_plot': {
                    'models': ['quasi_tem', 'quasi_tem_log', 'nakagawa', 'sunde', 'carson'],
                    'ylim': (0.5, 1.05),
                    'grid': True
                },
                'soil_plot': {
                    'ylim': (0.5, 1.05)
                }
            },
            'impedance': {
                'title': f'Characteristic Impedance {title_suffix}',
                'ylabel': r'Absolute Characteristic Impedance, $|Z_c| \, (\Omega)$',
                'quantity_key': 'zc',
                'transform': np.abs,
                'xscale': 'log',
                'xlim': (1E0, 1E10),
                'single_plot': {
                    'models': ['quasi_tem', 'quasi_tem_log', 'nakagawa', 'sunde', 'carson'],
                    'ylim': (0, 1200),
                    'grid': True
                },
                'soil_plot': {
                    'ylim': (0, 1200)
                }
            }
        }

        self.model_styles = {
            'quasi_tem': {'color': 'black', 'linestyle': '-', 'label': 'Quasi-TEM (Integral Eq.)'},
            'quasi_tem_log': {'color': 'blue', 'linestyle': '--', 'label': 'Quasi-TEM (Approx. Log.)'},
            'nakagawa': {'color': 'red', 'linestyle': '-', 'label': 'Wise (1948) & Nakagawa (1981)'},
            'sunde': {'color': 'green', 'linestyle': '-', 'label': 'Sunde (1968)'},
            'carson': {'color': 'black', 'linestyle': ':', 'label': 'Carson (1926)'}
        }

        self.soil_configs = {
            'keys': ['1e+6', '5e-3', '2e-4'],
            'styles': {
                '1e+6': {'color': 'black', 'linestyle': '--', 'label': r'$\rho_g = 1 \, \mu \Omega m$'},
                '5e-3': {'color': 'black', 'linestyle': '-', 'label': r'$\rho_g = 200 \, \Omega m$'},
                '2e-4': {'color': 'red', 'linestyle': ':', 'label': r'$\rho_g = 5000 \, \Omega m$'}
            }
        }

    def _plot_single_comparison(self, config_key):
        """Generic method to create a single plot comparing different models."""
        config = self.plot_configs[config_key]
        plot_settings = config['single_plot']        
        
        plt.figure(figsize=(8, 5))
        for model_key in plot_settings['models']:
            data_source = self.pul_data[model_key]
            if isinstance(data_source, dict):
                data_source = data_source['5e-3']

            raw_data = [z[config['quantity_key']][self.p, self.q] for z in data_source]
            y_values = config['transform'](raw_data)
            plt.plot(self.f, y_values, **self.model_styles[model_key])
        
        plt.title(config['title'])
        plt.ylabel(config['ylabel'])
        plt.xlabel('Frequency (Hz)')
        plt.xscale(config['xscale'])
        plt.xlim(config['xlim'])
        plt.ylim(plot_settings['ylim'])
        plt.legend(loc='best', fontsize='small')
        if plot_settings.get('grid', False):
            plt.grid(True, linestyle='--', alpha=0.6)

    def _plot_soil_comparison(self, config_key):
        """Generic method for 1x2 subplots comparing soil resistivity effects."""
        # ### ALTERAÇÃO AQUI ###
        # Fetch common and plot-specific settings
        config = self.plot_configs[config_key]
        plot_settings = config['soil_plot']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5), sharey=True)
        fig.suptitle(config['title'], fontsize=11)
        
        main_models = {'nakagawa': ax1, 'quasi_tem': ax2}
        titles = {'nakagawa': 'Nakagawa Form.', 'quasi_tem': 'Quasi-TEM Integral Approx.'}
        
        for model_key, ax in main_models.items():
            for soil_key in self.soil_configs['keys']:
                data_source = self.pul_data[model_key][soil_key]
                raw_data = [z[config['quantity_key']][self.p, self.q] for z in data_source]
                y_values = config['transform'](raw_data)
                ax.plot(self.f, y_values, **self.soil_configs['styles'][soil_key])

            # Use the specific ylim from plot_settings
            ax.set(xscale=config['xscale'], xlim=config['xlim'], ylim=plot_settings['ylim'],
                   xlabel='Frequency (Hz)', title=titles[model_key])
            ax.legend(loc='best', fontsize='small')
            ax.grid(True, which='both', ls='--', linewidth=0.5)

        ax1.set_ylabel(config['ylabel'])
        plt.tight_layout(rect=[0, 0, 1, 0.95])

    def plot_attenuation_constant(self):
        self._plot_single_comparison('attenuation')

    def plot_attenuation_constant_soil(self):
        self._plot_soil_comparison('attenuation')

    def plot_phase_velocity(self):
        self._plot_single_comparison('phase_velocity')

    def plot_phase_velocity_soil(self):
        self._plot_soil_comparison('phase_velocity')

    def plot_characteristic_impedance(self):
        self._plot_single_comparison('impedance')

    def plot_characteristic_impedance_soil(self):
        self._plot_soil_comparison('impedance')

if __name__ == "__main__":
    st = time.time()
    frequency = {'Analytically': np.logspace(0, 10, num=140)}

    # --- Data-Driven Calculation Setup ---
    pul_parameters = {}
    soil_conductivities = {'5e-3': MODEL, '1e+6': None, '2e-4': None}
    
    # Create model variations
    model_b = copy.deepcopy(MODEL); model_b[0]['conductivity'] = 1e+6
    model_c = copy.deepcopy(MODEL); model_c[0]['conductivity'] = 2e-4
    soil_conductivities['1e+6'] = model_b
    soil_conductivities['2e-4'] = model_c
    
    models_to_run_soil_variant = ['nakagawa', 'quasi_tem']
    models_to_run_single = ['sunde', 'carson', 'quasi_tem_log', 'sunde_log', 'deri']

    # Initialize dictionary structure
    for model in models_to_run_soil_variant:
        pul_parameters[model] = {key: [] for key in soil_conductivities}
    for model in models_to_run_single:
        pul_parameters[model] = []

    # --- Refactored Calculation Loop ---
    for f in frequency['Analytically']:
        # Calculate models with soil variations
        for key, model_data in soil_conductivities.items():
            pul_instance = PerUnitParameters(MulticonductorTransmissionLine(model_data), f)
            for form_key in models_to_run_soil_variant:
                pul_parameters[form_key][key].append(pul_instance.pul_extended_theory(form=form_key))
        
        # Calculate models with default soil
        pul_default = PerUnitParameters(MulticonductorTransmissionLine(MODEL), f)
        for form_key in models_to_run_single:
            pul_parameters[form_key].append(pul_default.pul_extended_theory(form=form_key))

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")

    # --- Plotting ---
    plotter = LimaModelPlotter(frequency, pul_parameters)
    plotter.plot_attenuation_constant()
    plotter.plot_attenuation_constant_soil()
    plotter.plot_phase_velocity()
    plotter.plot_phase_velocity_soil()
    plotter.plot_characteristic_impedance()
    plotter.plot_characteristic_impedance_soil()

    MTLRepresentation(MulticonductorTransmissionLine(MODEL), units='meter').ground_return_systems()
    plt.show()