# Adicione esta classe ao seu arquivo de plotters, como o PlotterModels.py

import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from mtl_main.source import MulticonductorTransmissionLine
from utils.case_utils import *

class InternalLinesModels:
    """
    Classe para gerar gráficos a partir dos resultados vetorizados
    do módulo de linhas aéreas.
    """
    def __init__(self, pul_data: dict, pul_data_tubular: dict, model: MulticonductorTransmissionLine):
        """
        Inicializa o plotter com os dados das simulações tubular e sólida.

        Args:
            pul_data_tubular (dict): Dicionário com os resultados do modelo tubular.
            pul_data_solid (dict): Dicionário com os resultados do modelo sólido equivalente.
            model_tubular (MulticonductorTransmissionLine): Objeto do modelo MTL tubular.
        """
        self.pul_data = pul_data
        self.pul_data_tubular = pul_data_tubular
        self.model = model
        self.f = pul_data_tubular['frequencies']
        self.w = 2 * np.pi * self.f

        # Skin Depth (m)
        self.skin_depth = 1 / np.sqrt(self.model.mu * np.pi * self.f * self.model.sigma)

    def solid_conductors_characteristics(self, p=0):
        """
        Gera o gráfico das características de impedância interna de um condutor
        cilíndrico sólido com base nos dados de simulação vetorizados.

        Args:
            p (int): O índice do condutor a ser analisado (padrão é 0).
        """
        # Pega as propriedades do condutor especificado 'p'
        conductor = self.model.mtl[p+1]
        ro = conductor['radius'][1] if isinstance(conductor['radius'], list) else conductor['radius']

        # Extrai os parâmetros DC (matrizes 2D) e pega o valor diagonal para o condutor 'p'
        Ri_cc = self.pul_data['internal']['Ri_cc'][p, p].real
        Li_cc = self.pul_data['internal']['Li_cc'][p, p].real

        # Extrai a impedância de Bessel (matriz 3D) e fatia para obter o vetor do condutor 'p'
        # A fatia [:, p, p] pega o elemento da diagonal (p, p) para todas as frequências (:)
        Zi = self.pul_data['internal']['Zi_bessel'][:, p, p]
        
        # Calcula a indutância interna a partir da reatância
        Li = np.divide(Zi.imag, self.w, out=np.zeros_like(self.w), where=self.w != 0)

        # Os cálculos das razões já são vetorizados
        R_ratio = Zi.real / Ri_cc
        wL_R_ratio = Zi.imag / Ri_cc
        L_ratio = Li / Li_cc
        wL_div_R = np.divide(Zi.imag, Zi.real, out=np.zeros_like(Zi.imag), where=Zi.real != 0)

        # --- Configuração do Gráfico (sem alteração na lógica) ---
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 6))
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Internal parameters for bare-wire conductor with $\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm', fontsize=13, y=0.97)

        x_axis = ro / self.skin_depth
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

    def internal_impedance(self, p=0):
        """
        Este método gera gráficos comparativos da impedância interna de um condutor
        usando a formulação exata, a aproximação de Nahman e Holt, e uma terceira
        aproximação.

        Os gráficos gerados são:
        1. Módulo da impedância interna |Z'_i| vs. Frequência.
        2. Ângulo da impedância interna arg(Z'_i) vs. Frequência (em graus).
        """
        # Pega as propriedades do condutor especificado 'p'
        conductor = self.model.mtl[p+1]
        ro = conductor['radius'][1] if isinstance(conductor['radius'], list) else conductor['radius']
        
        # Extrai os parâmetros DC (matrizes 2D) e pega o valor diagonal para o condutor 'p'
        Zi_bessel = self.pul_data['internal']['Zi_bessel'][:, p, p]
        Zi_kelvin = self.pul_data['internal']['Zi_kelvin'][:, p, p]    

        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Comparação de Modelos de Impedância Interna ($\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm)', fontsize=14, y=0.98)

        # --- Gráfico 1: Módulo ---
        ax1.set_title('Módulo da Impedância')
        ax1.plot(self.f, np.abs(Zi_bessel), 'k-', label='Bessel')
        ax1.plot(self.f, np.abs(Zi_kelvin), 'r--', lw=1, label='Kelvin')
        ax1.set_xscale('log')
        ax1.set_xlim(left=min(self.f), right=max(self.f))
        # ax1.set_ylim(0, np.max(mod_bessel) * 1.1)
        ax1.set_xlabel('Frequência (Hz)', fontsize=12)
        ax1.set_ylabel(r"Módulo $|Z_i|$ ($\Omega / m$)", fontsize=12)
        ax1.grid(True, which="both", ls=":", color='0.7')
        ax1.legend(loc='upper left')

        # --- Gráfico 2: Ângulo ---
        ax2.set_title('Ângulo da Impedância')
        ax2.plot(self.f, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Bessel')
        ax2.plot(self.f, np.angle(Zi_kelvin, deg=True), 'r--', lw=1, label='Kelvin')
        ax2.set_xscale('log')
        ax2.set_xlim(left=min(self.f), right=max(self.f))
        # ax2.set_ylim(0, np.max(angle_bessel) * 1.1)
        ax2.set_xlabel('Frequência (Hz)', fontsize=12)
        ax2.set_ylabel(r"Ângulo $Z_i$ (Graus)", fontsize=12)
        ax2.grid(True, which="both", ls=":", color='0.7')
        ax2.legend(loc='upper left')
        plt.tight_layout(rect=[0, 0, 1, 0.95])

    def nahman_holt_comparison(self, p=0):
        """
        Este método gera gráficos comparativos da impedância interna de um condutor
        usando a formulação exata, a aproximação de Nahman e Holt, e uma terceira
        aproximação.

        Os gráficos gerados são:
        1. Módulo da impedância interna |Z'_i| vs. Frequência.
        2. Ângulo da impedância interna arg(Z'_i) vs. Frequência (em graus).
        """
        # Pega as propriedades do condutor especificado 'p'
        conductor = self.model.mtl[p+1]
        ro = conductor['radius'][1] if isinstance(conductor['radius'], list) else conductor['radius']
        
        # Extrai os parâmetros DC (matrizes 2D) e pega o valor diagonal para o condutor 'p'
        Zi_bessel = self.pul_data['internal']['Zi_bessel'][:, p, p]
        Zi_nahman = self.pul_data['internal']['Zi_nahman'][:, p, p]  
        Zi_approx = self.pul_data['internal']['Zi_approx'][:, p, p]  
        
        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(fr'Comparação de Modelos de Impedância Interna ($\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.1f}$ mm)', fontsize=14, y=0.98)

        # --- Gráfico 1: Módulo ---
        ax1.set_title('Módulo da Impedância')
        ax1.plot(self.f, np.abs(Zi_bessel), 'k-', label='Exactly')
        ax1.plot(self.f, np.abs(Zi_nahman), 'r-', lw=1, label='Nahman e Holt')
        ax1.plot(self.f, np.abs(Zi_approx), 'g--', lw=1, label='Nahman e Holt (Modified)')
        ax1.set_xscale('log')
        ax1.set_xlim(left=min(self.f), right=max(self.f))
        ax1.set_ylim(0, np.max(np.abs(Zi_bessel)) * 1.1)
        ax1.set_xlabel('Frequência (Hz)', fontsize=12)
        ax1.set_ylabel(r"Módulo $|Z_i|$ ($\Omega / m$)", fontsize=12)
        ax1.grid(True, which="both", ls=":", color='0.7')
        ax1.legend(loc='upper left')

        # --- Gráfico 2: Ângulo ---
        ax2.set_title('Ângulo da Impedância')
        ax2.plot(self.f, np.angle(Zi_bessel, deg=True), 'k-', lw=1, label='Exactly')
        ax2.plot(self.f, np.angle(Zi_nahman, deg=True), 'r-', lw=1, label='Nahman e Holt')
        ax2.plot(self.f, np.angle(Zi_approx, deg=True), 'g--', lw=1, label='Nahman e Holt (Modified)')
        ax2.set_xscale('log')
        ax2.set_xlim(left=min(self.f), right=max(self.f))
        ax2.set_ylim(0, np.max(np.angle(Zi_bessel, deg=True)) * 1.1)
        ax2.set_xlabel('Frequência (Hz)', fontsize=12)
        ax2.set_ylabel(r"Ângulo $Z_i$ (Graus)", fontsize=12)
        ax2.grid(True, which="both", ls=":", color='0.7')
        ax2.legend(loc='upper left')
        plt.tight_layout(rect=[0, 0, 1, 0.95])

    def tubular_characteristics(self, p=0):
        """
        Gera um gráfico comparativo da impedância interna de um condutor tubular
        com a de um condutor sólido, usando dados pré-calculados.

        Args:
            p (int): O índice do condutor a ser analisado.
        """
        conductor = self.model.mtl[p+1]
        ri, ro = conductor['radius'] if isinstance(conductor['radius'], list) else (0, conductor['radius'])

        # Acessa a matriz 3D 'Zi_bessel' diretamente, sem a chave aninhada.
        Zi_tubular = self.pul_data_tubular['internal']['Zi_bessel'][:, p, p]
        Zi_solid = self.pul_data['internal']['Zi_bessel'][:, p, p]

        # --- 2. Cálculo das Razões e Plotagem ---

        R_ratio = np.divide(Zi_tubular.real, Zi_solid.real, 
                            out=np.ones_like(self.f), where=Zi_solid.real != 0)
        L_ratio = np.divide(Zi_tubular.imag, Zi_solid.imag,
                            out=np.ones_like(self.f), where=Zi_solid.imag != 0)
        
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 6))
        sigma_str = format_scientific_notation(self.model.sigma[p])
        fig.suptitle(f'Internal parameters for tubular bare-wire conductor with\n'
                     fr'$\sigma={sigma_str}$ S/m, $r_o={ro*1e3:.2f}$ mm, $r_i = {ri*1e3:.2f}$ mm', fontsize=12)

        x_axis = ro / self.skin_depth
        ax.plot(x_axis, R_ratio, 'k-', lw=1.5, label=r"$R_{i(\text{tubular})} / R_{i(\text{solid})}$")
        ax.plot(x_axis, L_ratio, 'r--', lw=1.5, label=r"$L_{i(\text{tubular})} / L_{i(\text{solid})}$")

        ax.set_xlabel(r'$(r_o / \delta)$', fontsize=12)
        ax.set_xlim(0, 8)
        ax.set_ylim(0, 2)
        ax.grid(True, which="both", ls="--", color='0.7')
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(fontsize=12, frameon=True)
        plt.tight_layout(rect=[0, 0.02, 1, 0.95])