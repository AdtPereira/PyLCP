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
import scipy.io as sio

# RAIZ DO PROJETO E DIRETÓRIOS
os.system('cls' if os.name == 'nt' else 'clear')
try:
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
    from mtl_main.models_scc import PRYSMIAN_SINGLE_PHASE as MODEL
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from scc.scc_data import scc_models_list
    from scc.scc_systm import SystemType
    from scc.scc_nlt import MonoNetworkTopology, NumericalLaplaceTransform
    from scc.scc_parameters import TransmissionLineParameters
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

def load_data(file_path):
    """Loads MATLAB data from the .mat file."""
    m_data = sio.loadmat(file_path)    
    return m_data['Vktd'][0], m_data['Vmtd'][0], m_data['t'][0], \
            m_data['Vk'][0], m_data['Vm'][0], m_data['f'][0]

def TimeDomain_Graph(Vktd_m, Vmtd_m, t_m, Vktd, Vmtd, NP):
    """Plots the comparison between MATLAB and Python results."""
    fig, axs = plt.subplots(2, 1, figsize=(8, 8))

    axs[0].plot(1E3 * t_m[0:NP], Vktd_m[0:NP], color='black', label='MATLAB')
    axs[0].plot(1E3 * nlt.t[0:NP], Vktd[0:NP], color='red', linestyle='--', label='Python')   
    axs[0].set_xlim(0, 0.2) 
    axs[0].set_ylim(0, 1.2) 
    axs[0].legend()
    
    axs[1].plot(1E3 * t_m[0:NP], Vmtd_m[0:NP], color='black', label='MATLAB')
    axs[1].plot(1E3 * nlt.t[0:NP], Vmtd[0:NP], color='red', linestyle='--', label='Python')
    axs[1].set_xlim(0, 0.2) 
    axs[1].set_ylim(0, 2) 
    axs[1].legend()

    plt.tight_layout()

if __name__ == "__main__":
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    st = time.time()
    os.system('cls' if os.name == 'nt' else 'clear')
    print("Iniciando cálculos da impedância p.u.l. do cabo coaxial...")

    # Rotinas Analítica e MoM-SO
    mtl_model = MulticonductorTransmissionLine(MODEL)
    syst = SystemType(Model=scc_models_list[0], rhog=100, erg=1, Syst_id='#1')
    tl = TransmissionLineParameters(s=1e3, syst=syst)

    # # Load the MATLAB or ATP data
    # Vktd_mat, Vmtd_mat, t_mat, Vk_mat, Vm_mat, f_mat = load_data('C:\\Users\\adilt\\OneDrive\\1 ACADEMIA\\MODELOS\\2.PRYSMIAN\\PRY_M01_1SCC_1C.mat')

    # # Define system configuration and soil parameters
    # # Model Name: PRY_M01_1SCC_1C
    # scc = scc_models_list[0]
    # syst = SystemType(Model=scc, rhog=100, erg=1, Syst_id='#1')

    # ## Transmission Line and NLT Constants
    # LX = 200                        # Line distance
    # RS = 1                          # Load Resistance [Ohm]
    # N = 1024*8                      # Point numbers
    # T = 1E-3                        # Simulation maximum time [s]
    # display = int(np.floor(N * 0.20));   # Graph window display

    # # Transmission Line Topology
    # nlt = NumericalLaplaceTransform(N, T)
    # VkVm = []
    # for s_value in nlt.s[0:round(N/2)+1]:
    #     tl = MonoNetworkTopology(s_value, syst, RS, LX)
    #     tl.StepUnitSource(ksi=0, T=T)
    #     tl.TerminalVoltages(tl.VS)
    #     VkVm.append(tl.V)

    # # Calculate Vk and Vm
    # Vk = np.array([v[0] for v in VkVm]).flatten()
    # Vm = np.array([v[1] for v in VkVm]).flatten()

    # # Extend Vk and Vm to include the conjugate
    # for k in range(round(N/2)+1, N):
    #     Vk = np.append(Vk, np.conj(Vk[2 * round(N/2) - k]))
    #     Vm = np.append(Vm, np.conj(Vm[2 * round(N/2) - k]))

    # # Main NLT routine
    # Vktd, Vmtd = nlt.Main_NLT(Vk, Vm)

    # TimeDomain_Graph(Vktd_mat, Vmtd_mat, t_mat, Vktd, Vmtd, display)

    MTLRepresentation(mtl_model, units='millimeter').ground_return_systems()
    plt.show()