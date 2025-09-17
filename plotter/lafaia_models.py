import os
import numpy as np
import matplotlib.pyplot as plt
from utils.case_utils import *

class LafaiaModels:
    """
    A highly refactored class to handle plotting for the Lafaia model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    It also saves the generated figures to a case-specific results directory.
    """
    def __init__(self, pul_data, case_name):
        """
        Initializes the PatelModels plotter.

        Args:
            pul_data (dict): A dictionary containing the per-unit-length data.
            case_name (str): The name of the case study, used to determine
                             the output directory for saving plots.
        """
        self.pul_data = pul_data
        self.case_name = case_name
        self.f = pul_data['analytical']['frequencies']
        self.f_mom = pul_data['numerical']['frequencies']
        self.w = 2 * np.pi * self.f
        self.figsize = (12, 5)
        self.r_factor = 1e3  # Convert Ohm/m to Ohm/km
        self.l_factor = 1e6  # Convert H/m to mH/km

        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.case_name, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

    def internal_impedance_elements(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle('Fig. 2.6: P.u.l. series impedance of a coaxial cable with sheath path return [Patel, 2014]', fontsize=12)

        # The slice [:, 0, 0] extracts the (0,0) element for all frequencies.        
        zcs = self.pul_data['analytical']['internal_parameters']['hybrid']['zcs']
        # rs = self.pul_data['numerical']['series_resistance_matrix'][:, 0, 0] 
        # ls = self.pul_data['numerical']['series_inductance_matrix'][:, 0, 0] 
        
        # Equivalent single conductor impedance. Eq. (2.10a) [AMETANI, 2015]
        Zcs = zcs['z11'] + zcs['z12'] + zcs['z2i'] 

        if 'comsol' in self.pul_data and self.pul_data['comsol'] is not None:
            cmsl_data = self.pul_data['comsol']['core_sheath_return'] 
            w = 2 * np.pi * cmsl_data['freq']
            ax1.scatter(cmsl_data['freq'], np.real(cmsl_data['coil_impedance']) * self.r_factor, label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)
            ax2.scatter(cmsl_data['freq'], np.imag(cmsl_data['coil_impedance']) / w * self.l_factor, label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)

            ax1.scatter(cmsl_data['freq'], cmsl_data['r11'] * self.r_factor, marker='o', facecolors='darkgray', s=10, zorder=3)
            ax2.scatter(cmsl_data['freq'], cmsl_data['l11'] * self.l_factor, marker='o', facecolors='darkgray', s=10, zorder=3)

            ax1.scatter(cmsl_data['freq'], cmsl_data['r2i'] * self.r_factor, marker='o', facecolors='darkblue', s=10, zorder=3)
            ax2.scatter(cmsl_data['freq'], cmsl_data['l2i'] * self.l_factor, marker='o', facecolors='darkblue', s=10, zorder=3)

            ax2.scatter(cmsl_data['freq'], cmsl_data['l12'] * self.l_factor, marker='o', facecolors='darkgreen', s=10, zorder=3)

        # ax1.scatter(self.f_mom, rs * self.r_factor, label='MoM-SO', marker='o', facecolors='none', edgecolors='k', s=50, zorder=3)
        # ax2.scatter(self.f_mom, ls * self.l_factor, label='MoM-SO', marker='o', facecolors='none', edgecolors='k', s=50, zorder=3)
        
        ax1.plot(self.f, np.real(zcs['z11']) * self.r_factor, label=r'$R_{11}$: internal resistance of core outer surface', linestyle='--', color='darkgray', linewidth=1.0) 
        ax1.plot(self.f, np.real(zcs['z12']) * self.r_factor, label=r'$R_{12}$: core outer insulator resistance', linestyle='--', color='green', linewidth=1.0) 
        ax1.plot(self.f, np.real(zcs['z2i']) * self.r_factor, label=r'$R_{2i}$: internal resistance of sheath inner surface', linestyle='--', color='blue', linewidth=1.0) 
        ax1.plot(self.f, np.real(Zcs) * self.r_factor, label=r'$R_{cs}=R_{11}+R_{12}+R_{2i}$ [1]', linestyle='-', color='black', linewidth=1.0)
        
        ax2.plot(self.f, np.imag(zcs['z11']) / self.w * self.l_factor, label=r'$L_{11}$: internal inductance of core outer surface', linestyle='--', color='darkgray', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zcs['z12']) / self.w * self.l_factor, label=r'$L_{12}$: core outer insulator inductance', linestyle='--', color='green', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zcs['z2i']) / self.w * self.l_factor, label=r'$L_{2i}$: internal inductance of sheath inner surface', linestyle='--', color='blue', linewidth=1.0) 
        ax2.plot(self.f,  np.imag(Zcs) / self.w * self.l_factor, label=r'$L_{cs}=L_{11}+L_{12}+L_{2i}$ [1]', linestyle='-', color='black', linewidth=1.0)  

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1e0, 1e6)
        ax1.set_ylim(1e-2, 1e1)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_{cs}$ $(\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(r'Internal Resistance, $R_{cs}$')

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1e0, 1e6)
        ax2.set_ylim(0, 0.2)
        # ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_{cs}$ (mH/km)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(r'Internal Inductance, $L_{cs}$')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        save_figure_multiformat(fig, self.results_dir, base_filename='patel_internal_impedance_elements')

    def internal_impedance_matrix_js_method(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(r'P.u.l. internal impedance matrix of a coaxial cable with $J_s$ Method [Yin, 1990]', fontsize=12)

        if 'comsol' in self.pul_data and self.pul_data['comsol'] is not None:
            core_exc = self.pul_data['comsol']['core']
            sheath_exc = self.pul_data['comsol']['sheath']
            Vs11 = core_exc['core_coil_voltage']       # Core Voltage source in the core excitation
            Vs12 = core_exc['sheath_coil_voltage']     # Sheath Voltage source in the core excitation
            Vs21 = sheath_exc['core_coil_voltage']     # Core Voltage source in the sheath excitation
            Vs22 = sheath_exc['sheath_coil_voltage']   # Sheath Voltage source in the sheath excitation
            w = 2 * np.pi * core_exc['freq']

            ax1.scatter(core_exc['freq'], np.real(Vs11) * self.r_factor,       label='mf.VCoil_Core (Core Exc.)', marker='o', edgecolor='black', facecolors='none', s=30)
            ax2.scatter(core_exc['freq'], np.imag(Vs11) / w * self.l_factor,   label='mf.VCoil_core (Core Exc.)', marker='o', edgecolor='black', facecolors='none', s=30)

            ax1.scatter(core_exc['freq'], np.real(Vs12) * self.r_factor,       label='mf.VCoil_Sheath (Core Exc.)', marker='o', edgecolor='darkgreen', facecolors='none', s=30)
            ax2.scatter(core_exc['freq'], np.imag(Vs12) / w * self.l_factor,   label='mf.VCoil_Sheath (Core Exc.)', marker='o', edgecolor='darkgreen', facecolors='none', s=30)

            ax1.scatter(sheath_exc['freq'], np.real(Vs21) * self.r_factor,     label='mf.VCoil_Core (Sheath Exc.)', marker='x', facecolors='darkgreen', s=12)
            ax2.scatter(sheath_exc['freq'], np.imag(Vs21) / w * self.l_factor, label='mf.VCoil_Core (Sheath Exc.)', marker='x', facecolors='darkgreen', s=12)

            ax1.scatter(sheath_exc['freq'], np.real(Vs22) * self.r_factor,     label='mf.VCoil_Sheath (Sheath Exc.)', marker='x', facecolors='darkblue', s=12)
            ax2.scatter(sheath_exc['freq'], np.imag(Vs22) / w * self.l_factor, label='mf.VCoil_Sheath (Sheath Exc.)', marker='x', facecolors='darkblue', s=12)

        # core self-impedance
        zi = self.pul_data['analytical']['internal_impedance_matrix']['hybrid']
        ax1.plot(self.f, np.real(zi[:, 0, 0]) * self.r_factor, label=r'$R_{cc}$: core self-resistance', linestyle='-', color='black', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zi[:, 0, 0]) / self.w * self.l_factor, label=r'$L_{cc}$: core self-inductance', linestyle='-', color='black', linewidth=1.0) 
        
        # mutual impedance between the core and sheath
        ax1.plot(self.f, np.real(zi[:, 0, 1]) * self.r_factor, label=r'$R_{cs}$: core-sheath mutual resistance', linestyle='-', color='darkgreen', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zi[:, 0, 1]) / self.w * self.l_factor, label=r'$L_{cs}$: core-sheath mutual inductance', linestyle='-', color='darkgreen', linewidth=1.0) 

        # sheath self-impedance
        ax1.plot(self.f, np.real(zi[:, 1, 1]) * self.r_factor, label=r'$R_{ss}$: sheath self-resistance', linestyle='-', color='darkblue', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zi[:, 1, 1]) / self.w * self.l_factor, label=r'$L_{ss}$: sheath self-inductance', linestyle='-', color='darkblue', linewidth=1.0) 

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1e0, 1e6)
        ax1.set_ylim(1e-3, 1e1)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$[R]$ $(\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(r'Internal Resistance, $R$')

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1e0, 1e6)
        ax2.set_ylim(0, 0.35)
        # ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$[L]$ (mH/km)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(r'Internal Inductance, $L$')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        save_figure_multiformat(fig, self.results_dir, base_filename='patel_internal_impedance_matrix_js_method')

    def internal_impedance_matrix_energy_method(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(r'P.u.l. internal impedance matrix of a coaxial cable with Loss-Energy Method [Yin, 1990]', fontsize=12)

        if 'comsol' in self.pul_data and self.pul_data['comsol'] is not None:
            core_exc = self.pul_data['comsol']['core']
            sheath_exc = self.pul_data['comsol']['sheath']

            ax1.scatter(core_exc['freq'],
                         (core_exc['r11']+core_exc['r2i']) * self.r_factor,
                         label='mf.r11+mf.r2i (Core Exc.)', marker='o', edgecolor='black', facecolors='none', s=30)
            ax2.scatter(core_exc['freq'],
                         (core_exc['l11']+core_exc['l2i']+core_exc['l12']+core_exc['l13']) * self.l_factor,
                         label='mf.L11+mf.L2i+mf.L12+mf.L13 (Core Exc.)', marker='o', edgecolor='black', facecolors='none', s=30)

            ax1.scatter(core_exc['freq'], 0.5 * core_exc['r2i'] * self.r_factor, label='mf.r2i (Core Exc.)', marker='o', edgecolor='darkgreen', facecolors='none', s=12)
            ax2.scatter(core_exc['freq'], 0.5 * core_exc['l2i'] * self.l_factor, label='mf.L2i (Core Exc.)', marker='o', edgecolor='darkgreen', facecolors='none', s=12)

            ax1.scatter(sheath_exc['freq'], sheath_exc['r2i'] * self.r_factor, label='mf.r2i (Sheath Exc.)', marker='x', facecolors='darkblue', s=12)
            ax2.scatter(sheath_exc['freq'], sheath_exc['l2i'] * self.l_factor, label='mf.L2i (Sheath Exc.)', marker='x', facecolors='darkblue', s=12)

        # core self-impedance
        zi = self.pul_data['analytical']['internal_impedance_matrix']['hybrid']
        ax1.plot(self.f, np.real(zi[:, 0, 0]) * self.r_factor, label=r'$R_{cc}$: core self-resistance', linestyle='-', color='black', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zi[:, 0, 0]) / self.w * self.l_factor, label=r'$L_{cc}$: core self-inductance', linestyle='-', color='black', linewidth=1.0) 
        
        # mutual impedance between the core and sheath
        ax1.plot(self.f, np.real(zi[:, 0, 1]) * self.r_factor, label=r'$R_{cs}$: mutual resistance between the core and sheath', linestyle='-', color='darkgreen', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zi[:, 0, 1]) / self.w * self.l_factor, label=r'$L_{cs}$: mutual inductance between the core and sheath', linestyle='-', color='darkgreen', linewidth=1.0) 

        # sheath self-impedance
        ax1.plot(self.f, np.real(zi[:, 1, 1]) * self.r_factor, label=r'$R_{ss}$: sheath self-resistance', linestyle='-', color='darkblue', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zi[:, 1, 1]) / self.w * self.l_factor, label=r'$L_{ss}$: sheath self-inductance', linestyle='-', color='darkblue', linewidth=1.0) 

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1e0, 1e6)
        ax1.set_ylim(1e-3, 1e1)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$[R]$ $(\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title(r'Internal Resistance, $R$')

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1e0, 1e6)
        ax2.set_ylim(0, 0.35)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$[L]$ (mH/km)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title(r'Internal Inductance, $L$')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        save_figure_multiformat(fig, self.results_dir, base_filename='patel_internal_impedance_matrix_energy_method')

