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
from pathlib import Path
import matplotlib.pyplot as plt

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
    from mtl_main.models_scc import DECONTI_SINGLE_PHASE as MODEL
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from scc.scc_data import scc_models_list
    from scc.scc_systm import SystemType
    from scc.scc_parameters import GroundReturnImpedance
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

def plot_ground_return_impedance(freq, pul_parameters, p, q):
    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    # fig.suptitle('')

    f = np.array(freq['Analytically'])
    zg = np.array([item[p, q] for item in pul_parameters['de_conti']])
    lg = np.imag(zg) / (2 * np.pi * f)

    # Extracting the impedance elements from zi_matrix
    ax1.plot(f, np.real(zg), label='Approx. De Conti et al.', color='red', linestyle='--')
    ax1.set_xscale('log')
    ax1.set_xlim(1E4, 1E7)
    ax1.set_ylim(0, 35)
    ax1.legend()
    ax1.set_xlabel('Frequency (Hz)')
    ax1.set_ylabel(r'$R_s \, (\Omega/m)$')
    ax1.grid(False)
    ax1.set_title('P.u.l. ground return resistance of the single buried bare-wire line')

    # Extracting the impedance elements from zi_matrix
    ax2.plot(f, 1E6 * lg, color='red', linestyle='--')
    ax2.set_xscale('log')
    ax2.set_xlim(1E-1, 1E7)
    ax2.set_ylim(0, 3.5)
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel(r'$L_s \, (mH/km)$')
    ax2.grid(False)
    ax2.set_title('P.u.l. ground return inductance of the single buried bare-wire line')
    plt.tight_layout(rect=[0, 0, 1, 0.96])

if __name__ == "__main__":
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    st = time.time()

    # Calculate the series impedance for each frequency
    frequency = {
        'Analytically': np.logspace(-1, 7, num=200),
        'Numerically': np.logspace(-1, 7, num=30)
    }

    # Dictionary to hold the series impedance calculations
    pul_parameters = {
        'de_conti': [],
    }

    mtl_model = MulticonductorTransmissionLine(MODEL)
    syst = SystemType(Model=scc_models_list[4], rhog=100, erg=10, Syst_id='#1')

    # Analytical Formulation
    for f in frequency['Analytically']:
        jw = 1j * 2 * np.pi * f
        tl = GroundReturnImpedance(jw, syst)
        pul_parameters['de_conti'].append(tl.DeConti(syst))

    plot_ground_return_impedance(frequency, pul_parameters, p=0, q=0)
    MTLRepresentation(mtl_model, units='millimeter').ground_return_systems()
    plt.show()