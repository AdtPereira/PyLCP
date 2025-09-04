""" 
Este script executa simulações de impedância de linhas de transmissão coaxiais
usando tanto uma abordagem analítica (formulação de Patel) quanto uma abordagem numérica
(Método dos Momentos - MoM-SO). Ele gera gráficos comparativos dos resultados
obtidos por ambas as metodologias.

Este arquivo é parte do projeto PyLCP, que é um pacote Python para análise de linhas de transmissão.

REFERENCES:
[1] XUE, Haoyan. General Formulation and Accurate Evaluation of Earth-Return Parameters
    for Overhead / Underground Cables. PhD thesis, Department of Electrical Engineering,
    École Polytechnique de Montréal, Université de Montréal, August 2018. 

[2] A. Ametani, T. Ohno and N. Nagaoka, Cable System Transients: Theory, Modeling and 
    Simulation, Wiley-IEEE Press, 2015.

[3] A. De Conti, N. Duarte and R. Alipio, "Closed-Form Expressions for the Calculation of the 
    Ground-Return Impedance and Admittance of Underground Cables," in IEEE Transactions on Power 
    Delivery, vol. 38, no. 4, pp. 2891-2900, Aug. 2023, doi: 10.1109/TPWRD.2023.3264614.

"""
import os
import sys
import time
import numpy as np
import copy
import scipy.constants as sc
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
    from mtl_main.models_scc import XUE_FLAT_ARRANGEMENT as MODEL
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
    def __init__(self, mtl_model, freq, pul):
        self.freq_data = freq
        self.pul_data = pul
        self.f = self.freq_data['Analytically']
        self.w = 2 * np.pi * self.f

        self.xue_series = [
            {'key': 'p100',      'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-'},
            {'key': 'p100_er20', 'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--'},
            {'key': 'p500',      'label': r'$\rho_e=500 \;\Omega m, \epsilon_r=1$',  'color': 'black', 'linestyle': '-.'},
            # {'key': 'p100_deconti',      'label': '', 'color': 'red', 'linestyle': ':'},
            # {'key': 'p100_er20_deconti', 'label': '', 'color': 'red', 'linestyle': ':'},
            # {'key': 'p500_deconti',      'label': '', 'color': 'red', 'linestyle': ':'}
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
                'p': 1,
                'q': 1,
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (0.5, 2.0)},
                'series_to_plot': self.xue_series
            },
            'fig421': {
                'suptitle': 'Figure 4.21: P.u.l. Mutual impedance between phase - a and phase - b sheaths with Magalhães/Xue formulation [1]',
                'resistance_title': 'P.u.l. series resistance',
                'inductance_title': 'P.u.l. series inductance',
                'p': 1,
                'q': 3,
                'x_lim': {'resistance': (1E4, 1E7), 'inductance': (1E4, 1E7)},
                'y_lim': {'resistance': (1E0, 1E5), 'inductance': (0.0, 1.5)},
                'series_to_plot': self.xue_series
            },
            'fig423': {
                'suptitle': 'Figure 4.23: P.u.l. Self-admittance of phase - a sheath with Magalhães/Xue formulation [1]',
                'conductance_title': 'P.u.l. shunt conductance',
                'capacitance_title': 'P.u.l. shunt capacitance',
                'p': 1,
                'q': 1,
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
            zs = np.array([item['series_impedance_matrix'][p, q] for item in self.pul_data[series['key']]])
            rs = np.real(zs) * 1e3  # Resistance in Ohm/km
            ls = np.imag(zs) / (2 * np.pi * self.f) * 1e6  # Inductance in mH/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, rs, **style)
            ax2.plot(self.f, ls, **style)

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(config['x_lim']['resistance'])
        ax1.set_ylim(config['y_lim']['resistance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$Rs_{{{p+1}{q+1}}} \, (\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['resistance_title'])

        # Configure right subplot (Inductance)
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
            ysh = np.array([item['shunt_admittance_matrix'][p, q] for item in self.pul_data[series['key']]])
            g_pq = np.real(ysh) * 1e3  # Conductance in S/km
            c_pq = np.imag(ysh) / (2 * np.pi * self.f) * 1e9  # Capacitance in uF/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, g_pq, **style)
            ax2.plot(self.f, c_pq, **style)

        # Configure left subplot (Conductance)
        ax1.set_xscale('log')
        ax1.set_xlim(config['x_lim']['conductance'])
        ax1.set_ylim(config['y_lim']['conductance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(fr'$G_{{{p+1}{q+1}}} \, (S/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['conductance_title'])

        # Configure right subplot (Capacitance)
        ax2.set_xscale('log')
        ax2.set_xlim(config['x_lim']['capacitance'])
        ax2.set_ylim(config['y_lim']['capacitance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(fr'$C_{{{p+1}{q+1}}} \, (\mu F/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['capacitance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _plot_propagation_constant_subplots(self, config_key):
        config = self.plot_configs[config_key]
        p, q = config['p'], config['q']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            gamma_v = np.array([item['gamma_v'][p, q] for item in self.pul_data[series['key']]])
            alfa = gamma_v.real * 1e3
            phase_vel = self.w / gamma_v.imag / sc.c
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, alfa, **style)
            ax2.plot(self.f, phase_vel, **style)

        # Configure left subplot (attenuation)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1E3, 1E9)
        ax1.set_ylim(config['y_lim']['attenuation'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel( r'Attenuation Constant, $\alpha_{\nu}$ (Np/km)')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['attenuation_title'])

        # Configure right subplot (phase velocity)
        ax2.set_xscale('log')
        ax2.set_xlim(1E3, 1E9)
        ax2.set_ylim(config['y_lim']['phase_velocity'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'Phase Velocity, $c_{\nu}/c_0$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['phase_velocity_title'])
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

if __name__ == "__main__":
    st = time.time()
    frequency = {'Analytically': np.logspace(3, 7, num=40)}

    # Define models for different physical scenarios
    mtl_model_a = MulticonductorTransmissionLine(MODEL)
    
    model_b_data = copy.deepcopy(MODEL)
    model_b_data[0]['relative_permittivity'] = 20
    mtl_model_b = MulticonductorTransmissionLine(model_b_data)

    model_c_data = copy.deepcopy(MODEL)
    model_c_data[0]['conductivity'] = 0.002 # rho = 500 Ohm.m
    mtl_model_c = MulticonductorTransmissionLine(model_c_data)

    # Define the calculation scenarios
    scenarios = {
        'p100':      {'mtl': mtl_model_a, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p100_er20': {'mtl': mtl_model_b, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p500':      {'mtl': mtl_model_c, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        # 'p100_deconti':      {'mtl': mtl_model_a, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        # 'p100_er20_deconti': {'mtl': mtl_model_b, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        # 'p500_deconti':      {'mtl': mtl_model_c, 'zg_form': 'deconti', 'yg_form': 'deconti'},
    }
    
    pul_parameters = {key: [] for key in scenarios}
    for f in frequency['Analytically']:
        # Internal Impedance elements
        internal = InternalPerUnitParameters(mtl_model_a, f)

        # Ground-return elements
        for key, params in scenarios.items():
            pul = PerUnitParameters(params['mtl'], f)
            pul_parameters[key].append(pul.quasi_tem_pul(
                internal.internal_matrices(), zg_form=params['zg_form'], yg_form=params['yg_form'])
            )

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = ModelPlotter(mtl_model_a, frequency, pul_parameters)
    plotter.plot_fig419()
    plotter.plot_fig421()
    plotter.plot_fig423()
    MTLRepresentation(mtl_model_a, units='millimeter').ground_return_systems()
    plt.show()