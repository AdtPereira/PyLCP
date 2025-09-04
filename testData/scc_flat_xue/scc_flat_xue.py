# C:\git\PyLCP>
# python -m testData.scc_flat_xue.scc_flat_xue
# python -m cProfile -o profile.stats -m testData.scc_flat_xue.scc_flat_xue
# snakeviz profile.stats

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

import sys
import os
import time
import copy
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# --- Configure project root for module imports ---
try:
    os.system('cls' if os.name == 'nt' else 'clear')
    project_root = Path(__file__).resolve().parents[2]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root)) 
    print(f"Project root configured at: {project_root}")
except IndexError:
    raise RuntimeError("Could not find project root. Ensure the directory structure is correct.")

# --- Import custom modules ---
try:
    from utils.case_utils import load_json_parameters
    from .plotter import ModelPlotter
    from models import single_core_cables as scc 
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_formulation.scc import InternalPerUnitParameters, PerUnitParameters
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    model = scc.three_phase_flat_model(input_json, show_model=True)
    frequency = {'Analytically': np.logspace(3, 7, num=40)}

    # Define models for different physical scenarios
    mtl_model_a = MulticonductorTransmissionLine(model)
    model_b_data = copy.deepcopy(model)
    model_b_data[0]['relative_permittivity'] = 20
    mtl_model_b = MulticonductorTransmissionLine(model_b_data)
    model_c_data = copy.deepcopy(model)
    model_c_data[0]['conductivity'] = 0.002 # rho = 500 Ohm.m
    mtl_model_c = MulticonductorTransmissionLine(model_c_data)

    # Define the calculation scenarios
    scenarios = {
        'p100':      {'mtl': mtl_model_a, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p100_er20': {'mtl': mtl_model_b, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p500':      {'mtl': mtl_model_c, 'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue'},
        'p100_deconti':      {'mtl': mtl_model_a, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        'p100_er20_deconti': {'mtl': mtl_model_b, 'zg_form': 'deconti', 'yg_form': 'deconti'},
        'p500_deconti':      {'mtl': mtl_model_c, 'zg_form': 'deconti', 'yg_form': 'deconti'},
    }
    
    pul_parameters = {key: [] for key in scenarios}
    pul_parameters['frequencies'] = frequency['Analytically']
    for f in frequency['Analytically']:
        # Internal Impedance elements
        internal = InternalPerUnitParameters(mtl_model_a, f)

        # Ground-return elements
        for key, value in scenarios.items():
            pul = PerUnitParameters(value['mtl'], f)
            pul_parameters[key].append(pul.quasi_tem_pul(
                internal.internal_matrices(), zg_form=value['zg_form'], yg_form=value['yg_form'])
            )

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = ModelPlotter(pul_parameters)
    plotter.plot_fig419()
    plotter.plot_fig421()
    plotter.plot_fig423()
    MTLRepresentation(mtl_model_a, units='centimeter').ground_return_systems()
    plt.show()    

if __name__ == "__main__":
    main()