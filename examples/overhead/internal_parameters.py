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
    from mtl_main.models_wires import SINGLE_WIRE_R05 as MODEL
    from mtl_main.graphics import MTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from mtl_main.utils import *
    from analytical_formulation.overhead_lines import PerUnitParameters
    from analytical_formulation.overhead_lines import verify_kelvin_functions
    print("Módulos e modelo de dados importados com sucesso.")
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

def plot_internal_solid_characteristics(pul_itens, mtl_model, p=0):
    """
    Gera o gráfico das características de impedância interna de um condutor
    cilíndrico sólido com base nos dados de simulação fornecidos.

    Args:
        pul_parameters (dict): Dicionário onde as chaves são frequências e os valores
                               são dicionários contendo os dados calculados.
        mtl_model (MulticonductorTransmissionLine): Objeto do modelo MTL
                                                   contendo as propriedades físicas.
        frequencies (np.ndarray): Array de frequências usadas na simulação, na
                                  ordem correta para a plotagem.
        p (int): O índice do condutor a ser analisado (padrão é 0).
    """
    ro = mtl_model.surfaces[p]['radius']
    freq = np.array(sorted(pul_itens.keys()))
    omega = 2 * np.pi * freq
    skin_depths = np.array([pul_itens[f]['skin_depth'] for f in freq])
    
    Ri_cc = pul_itens[freq[0]]['itens']['Ri_cc'][p, p].real
    Li_cc = pul_itens[freq[0]]['itens']['Li_cc'][p, p].real
    Zi = np.array([pul_itens[f]['itens']['Zi_bessel'][p, p] for f in freq])
    Li = np.divide(Zi.imag, omega, out=np.zeros_like(Zi.imag), where=omega!=0)

    R_ratio = Zi.real / Ri_cc
    wL_R_ratio = Zi.imag / Ri_cc
    L_ratio = Li / Li_cc
    wL_div_R = np.divide(Zi.imag, Zi.real, out=np.zeros_like(Zi.imag), where=Zi.real!=0)

    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(8, 6))
    sigma = format_scientific_notation(mtl_model.sigma[p])
    fig.suptitle(fr'Internal parameters for bare-wire conductor with $\sigma={sigma}$ S/m, $r_o={ro*1e3:.1f}$ mm', fontsize=13, y=0.97)

    x_axis = ro / skin_depths
    ax.plot(x_axis, R_ratio,    'k-', lw=1, label=r"$R_i / R_{i(cc)}$")
    ax.plot(x_axis, wL_R_ratio, 'b-', lw=1, label=r"$\omega L_i / R_{i(cc)}$")
    ax.plot(x_axis, wL_div_R,   'g-', lw=1, label=r"$\omega L_i / R_i$")
    ax.plot(x_axis, L_ratio,    'r-', lw=1, label=r"$L_i / L_{i(cc)}$")

    ax.set_xscale('log')
    ax.set_xlabel(r'$(r_o / \delta)$', fontsize=12)
    ax.set_xlim(left=min(x_axis), right=max(x_axis))
    ax.set_ylim(0, max(np.max(R_ratio), np.max(wL_R_ratio)) * 0.15)
    ax.grid(True, which="both", ls="--", color='0.7')
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(fontsize=12, frameon=True)
    plt.tight_layout(rect=[0, 0, 1, 1])

def plot_internal_impedance(pul_itens, mtl_model, p=0):
    """
    Este método gera gráficos comparativos da impedância interna de um condutor
    usando a formulação exata, a aproximação de Nahman e Holt, e uma terceira
    aproximação.

    Os gráficos gerados são:
    1. Módulo da impedância interna |Z'_i| vs. Frequência.
    2. Ângulo da impedância interna arg(Z'_i) vs. Frequência (em graus).
    """
    freq = np.array(sorted(pul_itens.keys()))
    Zi_bessel = np.array([pul_itens[f]['itens']['Zi_bessel'][p, p] for f in freq])
    Zi_kelvin = np.array([pul_itens[f]['itens']['Zi_kelvin'][p, p] for f in freq])

    plt.style.use('default')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    sigma = format_scientific_notation(mtl_model.sigma[p])
    ro = mtl_model.surfaces[p]['radius']
    fig.suptitle(fr'Comparação de Modelos de Impedância Interna ($\sigma={sigma}$ S/m, $r_o={ro*1e3:.1f}$ mm)', fontsize=14, y=0.98)

    # --- Gráfico 1: Módulo ---
    ax1.set_title('Módulo da Impedância')
    ax1.plot(freq, np.abs(Zi_bessel), 'k-', label='Bessel')
    ax1.plot(freq, np.abs(Zi_kelvin), 'r--', lw=1, label='Kelvin')
    ax1.set_xscale('log')
    ax1.set_xlim(left=min(freq), right=max(freq))
    # ax1.set_ylim(0, np.max(mod_bessel) * 1.1)
    ax1.set_xlabel('Frequência (Hz)', fontsize=12)
    ax1.set_ylabel(r"Módulo $|Z_i|$ ($\Omega / m$)", fontsize=12)
    ax1.grid(True, which="both", ls=":", color='0.7')
    ax1.legend(loc='upper left')

    # --- Gráfico 2: Ângulo ---
    ax2.set_title('Ângulo da Impedância')
    ax2.plot(freq, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Bessel')
    ax2.plot(freq, np.angle(Zi_kelvin, deg=True), 'r--', lw=1, label='Kelvin')
    ax2.set_xscale('log')
    ax2.set_xlim(left=min(freq), right=max(freq))
    # ax2.set_ylim(0, np.max(angle_bessel) * 1.1)
    ax2.set_xlabel('Frequência (Hz)', fontsize=12)
    ax2.set_ylabel(r"Ângulo $Z_i$ (Graus)", fontsize=12)
    ax2.grid(True, which="both", ls=":", color='0.7')
    ax2.legend(loc='upper left')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

def plot_nahman_holt_comparison(pul_itens, mtl_model, p=0):
    """
    Este método gera gráficos comparativos da impedância interna de um condutor
    usando a formulação exata, a aproximação de Nahman e Holt, e uma terceira
    aproximação.

    Os gráficos gerados são:
    1. Módulo da impedância interna |Z'_i| vs. Frequência.
    2. Ângulo da impedância interna arg(Z'_i) vs. Frequência (em graus).
    """
    freq = np.array(sorted(pul_itens.keys()))
    Zi_bessel = np.array([pul_itens[f]['itens']['Zi_bessel'][p, p] for f in freq])
    Zi_approx = np.array([pul_itens[f]['itens']['Zi_approx'][p, p] for f in freq])
    Zi_nahman = np.array([pul_itens[f]['itens']['Zi_nahman'][p, p] for f in freq])

    plt.style.use('default')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    sigma = format_scientific_notation(mtl_model.sigma[p])
    ro = mtl_model.surfaces[p]['radius']
    fig.suptitle(fr'Comparação de Modelos de Impedância Interna ($\sigma={sigma}$ S/m, $r_o={ro*1e3:.1f}$ mm)', fontsize=14, y=0.98)

    # --- Gráfico 1: Módulo ---
    ax1.set_title('Módulo da Impedância')
    ax1.plot(freq, np.abs(Zi_bessel), 'k-', label='Exactly')
    ax1.plot(freq, np.abs(Zi_nahman), 'r-', lw=1, label='Nahman e Holt')
    ax1.plot(freq, np.abs(Zi_approx), 'g--', lw=1, label='Nahman e Holt (Modified)')
    ax1.set_xscale('log')
    ax1.set_xlim(left=min(freq), right=max(freq))
    ax1.set_ylim(0, np.max(np.abs(Zi_bessel)) * 1.1)
    ax1.set_xlabel('Frequência (Hz)', fontsize=12)
    ax1.set_ylabel(r"Módulo $|Z_i|$ ($\Omega / m$)", fontsize=12)
    ax1.grid(True, which="both", ls=":", color='0.7')
    ax1.legend(loc='upper left')

    # --- Gráfico 2: Ângulo ---
    ax2.set_title('Ângulo da Impedância')
    ax2.plot(freq, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Exactly')
    ax2.plot(freq, np.angle(Zi_nahman, deg=True), 'r-', lw=1, label='Nahman e Holt')
    ax2.plot(freq, np.angle(Zi_approx, deg=True), 'g--', lw=1, label='Nahman e Holt (Modified)')
    ax2.set_xscale('log')
    ax2.set_xlim(left=min(freq), right=max(freq))
    ax2.set_ylim(0, np.max(np.angle(Zi_bessel, deg=True)) * 1.1)
    ax2.set_xlabel('Frequência (Hz)', fontsize=12)
    ax2.set_ylabel(r"Ângulo $Z_i$ (Graus)", fontsize=12)
    ax2.grid(True, which="both", ls=":", color='0.7')
    ax2.legend(loc='upper left')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

def plot_internal_tubular_characteristics(solid_pul_itens, tubular_pul_itens, mtl_model, p=0):
    """
    Gera o gráfico das características de impedância interna de um condutor
    cilíndrico sólido com base nos dados de simulação fornecidos.

    Args:
        pul_parameters (dict): Dicionário onde as chaves são frequências e os valores
                               são dicionários contendo os dados calculados.
        mtl_model (MulticonductorTransmissionLine): Objeto do modelo MTL
                                                   contendo as propriedades físicas.
        frequencies (np.ndarray): Array de frequências usadas na simulação, na
                                  ordem correta para a plotagem.
        p (int): O índice do condutor a ser analisado (padrão é 0).
    """
    ri, ro = mtl_model.mtl[p+1]['radius']
    freq = np.array(sorted(solid_pul_itens.keys()))
    skin_depths = np.array([solid_pul_itens[f]['skin_depth'] for f in freq])
    Zi_solid = np.array([solid_pul_itens[f]['itens']['Zi_bessel'][p, p] for f in freq])
    Zi_tubular = np.array([tubular_pul_itens[f]['itens']['Zi_bessel']['external_return'][p, p] for f in freq])
    R_ratio = Zi_tubular.real / Zi_solid.real
    L_ratio = Zi_tubular.imag / Zi_solid.imag

    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(8, 6))
    sigma = format_scientific_notation(mtl_model.sigma[p])
    fig.suptitle(f'Internal parameters for tubular bare-wire conductor with\n'
             f'$\sigma={sigma}$ S/m, $r_o={ro*1e3:.2f}$ mm, $r_i = {ri*1e3:.2f}$ mm', 
             fontsize=12)

    x_axis = ro / skin_depths
    ax.plot(x_axis, R_ratio, 'k-', lw=1, label=r"$R_{i(tubular)} / R_{i(solid)}$")
    ax.plot(x_axis, L_ratio, 'r-', lw=1, label=r"$L_{i(tubular)} / L_{i(solid)}$")

    ax.set_xlabel(r'$(r_o / \delta)$', fontsize=12)
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 2)
    ax.grid(True, which="both", ls="--", color='0.7')
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(fontsize=12, frameon=True)
    plt.tight_layout(rect=[0, 0, 1, 1])

if __name__ == "__main__":
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    st = time.time()

    frequency = {
        'Analytically': np.logspace(0, 8, num=200),
        'Numerically': np.logspace(0, 8, num=30)
    }

    solid_pul_itens, tubular_pul_itens = {}, {}
    tubular_model = copy.deepcopy(MODEL)
    tubular_model[1]['radius'][0] = 0.5 * tubular_model[1]['radius'][1]

    solid_mtl_model = MulticonductorTransmissionLine(MODEL)
    tubular_mtl_model = MulticonductorTransmissionLine(tubular_model)
    for f in frequency['Analytically']:
        solid_pul = PerUnitParameters(solid_mtl_model, f)
        tubular_pul = PerUnitParameters(tubular_mtl_model, f)

        solid_pul_itens[f] = {
            'skin_depth': solid_pul.skin_depth,
            'itens': solid_pul.internal_impedance_elements_solid_wires()
        }

        tubular_pul_itens[f] = {
            'skin_depth': tubular_pul.skin_depth,
            'itens': tubular_pul.internal_impedance_elements_tubular_wires()
        }

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.2f} seconds.\n")
    plot_internal_solid_characteristics(solid_pul_itens, solid_mtl_model)
    plot_internal_impedance(solid_pul_itens, solid_mtl_model)
    plot_nahman_holt_comparison(solid_pul_itens, solid_mtl_model)
    plot_internal_tubular_characteristics(solid_pul_itens, tubular_pul_itens, tubular_mtl_model)
    MTLRepresentation(solid_mtl_model, units='millimeter').ground_return_systems()
    MTLRepresentation(tubular_mtl_model, units='millimeter').ground_return_systems()
    verify_kelvin_functions()
    plt.show()