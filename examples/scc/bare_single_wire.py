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
import copy
import numpy as np
import pandas as pd
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
    from mtl_main.models_scc import BARE_SINGLE_WIRE as MODEL
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_formulation.scc import PerUnitParameters
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

class DeContiModelPlotter:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, mtl_model, freq, pul, p=0, q=0):
        self.freq_data = freq
        self.pul_data = pul
        self.p = p
        self.q = q
        self.f = self.freq_data['Analytically']
        self.w = 2 * np.pi * self.f
        self.ro = mtl_model.surfaces[0]['radius']
        self.h1 = mtl_model.surfaces[0]['center_point'][1]

        self.nakagawa_series = [
            {'key': 'magalhaes_xue', 'label': 'Magalhães/Xue',  'color': 'black', 'linestyle': '-', 'linewidth': 2},
            {'key': 'sunde', 'label': 'Sunde', 'color': 'gray', 'linestyle': ':', 'linewidth': 2},
            {'key': 'pollaczek', 'label': 'Pollaczek', 'color': 'darkblue', 'linestyle': '-', 'linewidth': 2},
            {'key': 'ametani', 'label': 'Ametani', 'color': 'cyan', 'linestyle': '-', 'linewidth': 1},
            {'key': 'deconti', 'label': 'De Conti', 'color': 'red', 'linestyle': '--', 'linewidth': 1},
            {'key': 'saad', 'label': 'Saad', 'color': 'red', 'linestyle': '--', 'linewidth': 1},
            {'key': 'wedepohl', 'label': 'Wedepohl and Wilcox', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 1},
        ]
        
        self.plot_configs = {
            'ground_return_impedance': {
                'suptitle': r'P.u.l. ground-return impedance for single buried bare-wire for $\rho_e=100 \;\Omega$ m and $\epsilon_r=10$',
                'resistance_title': 'P.u.l. resistance',
                'inductance_title': 'P.u.l. inductance',
                'x_lim': {'resistance': (1E5, 1E8), 'inductance': (1E-1, 1E8)},
                'y_lim': {'resistance': (0, 200), 'inductance': (0.5, 3.0)},
                'series_to_plot': self.nakagawa_series
            },
        }

    def _plot_impedance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            zg_raw = np.array([item['zg'][self.p, self.q] for item in self.pul_data[series['key']]])
            rg = np.real(zg_raw) 
            lg = np.imag(zg_raw) / (2 * np.pi * self.f) * 1e6  # Inductance in mH/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle'], 'linewidth': series['linewidth']}
            ax1.plot(self.f, rg, **style)
            ax2.plot(self.f, lg, **style)

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

    def _plot_admittance_subplots(self, config_key):
        """
        Generic method to create a 1x2 subplot for shunt conductance (left)
        and shunt capacitance (right) based on a configuration key.
        """
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            ysh_raw = np.array([item['ysh'][self.p, self.q] for item in self.pul_data[series['key']]])
            cond = np.real(ysh_raw) * 1e3  # Conductance in S/km
            cap = np.imag(ysh_raw) / (2 * np.pi * self.f) * 1e12  # Capacitance in nF/km
            style = {'label': series['label'], 'color': series['color'], 'linestyle': series['linestyle']}
            ax1.plot(self.f, cond, **style)
            ax2.plot(self.f, cap, **style)

        # Configure left subplot (Conductance)
        ax1.set_xscale('log')
        ax1.set_xlim(1E3, 1E9)
        ax1.set_ylim(config['y_lim']['conductance'])
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$G \, (S/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(config['conductance_title'])

        # Configure right subplot (Capacitance)
        ax2.set_xscale('log')
        ax2.set_xlim(1E3, 1E9)
        ax2.set_ylim(config['y_lim']['capacitance'])
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$C \, (nF/km)$')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(config['capacitance_title'])
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def _plot_propagation_constant_subplots(self, config_key):
        config = self.plot_configs[config_key]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        fig.suptitle(config['suptitle'], fontsize=12, y=0.98)

        for series in config['series_to_plot']:
            gamma_v = np.array([item['gamma_v'][self.p, self.q] for item in self.pul_data[series['key']]])
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

    def plot_ground_return_impedance(self):
        self._plot_impedance_subplots('ground_return_impedance')

if __name__ == "__main__":
    st = time.time()
    frequency = {'Analytically': np.logspace(-1, 8, num=200)}

    # Define models for different physical scenarios
    mtl_model = MulticonductorTransmissionLine(MODEL)
    
    # Define the calculation scenarios
    scenarios = {
        'magalhaes_xue': {'mtl': mtl_model, 'form': 'magalhaes_xue'},
        'sunde': {'mtl': mtl_model, 'form': 'sunde'},
        'pollaczek': {'mtl': mtl_model, 'form': 'pollaczek'},
        'ametani': {'mtl': mtl_model, 'form': 'ametani'},
        'deconti': {'mtl': mtl_model, 'form': 'deconti'},
        'saad': {'mtl': mtl_model, 'form': 'saad'},
        'wedepohl': {'mtl': mtl_model, 'form': 'wedepohl'},
    }
    
    pul_parameters = {key: [] for key in scenarios}
    for f in frequency['Analytically']:
        for key, params in scenarios.items():
            pul_parameters[key].append(
                PerUnitParameters(params['mtl'], f).ground_return_parameters(form=params['form'])
            )

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")
    plotter = DeContiModelPlotter(mtl_model, frequency, pul_parameters)
    plotter.plot_ground_return_impedance()
    MTLRepresentation(mtl_model, units='millimeter').ground_return_systems()
    plt.show()