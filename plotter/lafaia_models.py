import os
import numpy as np
import matplotlib.pyplot as plt
from utils.case_utils import *
from analytical_forms.single_core_cable import capacitance_matrix_from_energy_method

class LafaiaModels:
    """
    A highly refactored class to handle plotting for the Lafaia model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    It also saves the generated figures to a case-specific results directory.
    """
    def __init__(self, pul_data: dict, case_name: str, autoSave: bool = True):
        """
        Initializes the PatelModels plotter.

        Args:
            pul_data (dict): A dictionary containing the per-unit-length data.
            case_name (str): The name of the case study, used to determine
                             the output directory for saving plots.
        """
        self.case_name = case_name
        self.autoSave = autoSave
        self.pul_0 = pul_data[0]
        self.pul_01 = pul_data[1]
        self.pul_02 = pul_data[2]
        self.pul_3 = pul_data[3]
        self.cmsl = pul_data[0]['comsol']
        self.f = pul_data[1]['analytical']['frequencies']
        self.w = 2 * np.pi * self.f

        self.r_factor = 1e3  # Convert Ohm/m to Ohm/km
        self.l_factor = 1e6  # Convert H/m to mH/km
        self.c_factor = 1e9  # Convert F/m to uF/km
        self.figsize = (12, 5)
        
        # Assumes the script is run from the project's root directory.
        self.results_dir = os.path.join('testData', self.case_name, 'Results')
        os.makedirs(self.results_dir, exist_ok=True)

        # Capacitance matrix from energy method
        print("Calculating capacitance matrix from energy method...")
        w11 = self.cmsl['cmsl_shunt_params']['wcc'][0]
        w22 = self.cmsl['cmsl_shunt_params']['wss'][0]
        w12 = self.cmsl['cmsl_shunt_params']['wcs'][0]
        self.cap_matrix_energy_method = capacitance_matrix_from_energy_method(np.array([w11, w22, w12]), v0=1.0)

    def internal_impedance_elements(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle('Fig. 2.6: P.u.l. series impedance of a single core cable with sheath path return [Patel, 2014]', fontsize=12)

        # The slice [:, 0, 0] extracts the (0,0) element for all frequencies.        
        zcs = self.pul_01['analytical']['internal_parameters']['zcs']
        
        # Equivalent single conductor impedance. Eq. (2.10a) [AMETANI, 2015]
        Zcs = zcs['z11'] + zcs['z12'] + zcs['z2i'] 

        if 'comsol' in self.pul_0 and self.cmsl is not None:
            data = self.cmsl['cmsl_core_sheath'] 
            
            ax1.scatter(self.f, np.real(data['coil_impedance']) * self.r_factor,
                        label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)
            
            ax2.scatter(self.f, np.imag(data['coil_impedance']) / self.w * self.l_factor,
                        label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)

            ax1.scatter(self.f, data['r11'] * self.r_factor, marker='o', facecolors='darkgray', s=10, zorder=3)
            
            ax2.scatter(self.f, data['l11'] * self.l_factor, marker='o', facecolors='darkgray', s=10, zorder=3)

            ax1.scatter(self.f, data['r2i'] * self.r_factor, marker='o', facecolors='darkblue', s=10, zorder=3)
            
            ax2.scatter(self.f, data['l2i'] * self.l_factor, marker='o', facecolors='darkblue', s=10, zorder=3)

            ax2.scatter(self.f, data['l12'] * self.l_factor, marker='o', facecolors='darkgreen', s=10, zorder=3)

        ax1.plot(self.f, np.real(zcs['z11']) * self.r_factor,
                 label=r'$R_{11}$: internal resistance of core outer surface', linestyle='--', color='darkgray', linewidth=1.0) 
        
        ax1.plot(self.f, np.real(zcs['z12']) * self.r_factor,
                 label=r'$R_{12}$: core outer insulator resistance', linestyle='--', color='green', linewidth=1.0) 
        
        ax1.plot(self.f, np.real(zcs['z2i']) * self.r_factor,
                 label=r'$R_{2i}$: internal resistance of sheath inner surface', linestyle='--', color='blue', linewidth=1.0) 
        
        ax1.plot(self.f, np.real(Zcs) * self.r_factor,
                 label=r'$R_{cs}=R_{11}+R_{12}+R_{2i}$ [1]', linestyle='-', color='black', linewidth=1.0)
        
        ax2.plot(self.f, np.imag(zcs['z11']) / self.w * self.l_factor,
                 label=r'$L_{11}$: internal inductance of core outer surface', linestyle='--', color='darkgray', linewidth=1.0) 
        
        ax2.plot(self.f, np.imag(zcs['z12']) / self.w * self.l_factor,
                 label=r'$L_{12}$: core outer insulator inductance', linestyle='--', color='green', linewidth=1.0) 
        
        ax2.plot(self.f, np.imag(zcs['z2i']) / self.w * self.l_factor,
                 label=r'$L_{2i}$: internal inductance of sheath inner surface', linestyle='--', color='blue', linewidth=1.0) 
        
        ax2.plot(self.f,  np.imag(Zcs) / self.w * self.l_factor,
                 label=r'$L_{cs}=L_{11}+L_{12}+L_{2i}$ [1]', linestyle='-', color='black', linewidth=1.0)  

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
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='internal_impedance_elements')

    def internal_impedance_matrix_js_method(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(r'Single core cable installed in HDPE tube [Lafaia, 2015] and $J_s$ Method [Yin, 1990]', fontsize=12)

        if 'comsol' in self.pul_0 and self.cmsl is not None:
            Vs11 = self.cmsl['cmsl_core_exc']['core_coil_voltage']       # Core Voltage source in the core excitation
            Vs12 = self.cmsl['cmsl_core_exc']['sheath_coil_voltage']     # Sheath Voltage source in the core excitation
            Vs21 = self.cmsl['cmsl_sheath_exc']['core_coil_voltage']     # Core Voltage source in the sheath excitation
            Vs22 = self.cmsl['cmsl_sheath_exc']['sheath_coil_voltage']   # Sheath Voltage source in the sheath excitation

            ax1.scatter(self.f, np.real(Vs11) * self.r_factor,
                        label='mf.VCoil_Core (Core Exc.)', marker='o', edgecolor='black', facecolors='none', s=30)
            
            ax2.scatter(self.f, np.imag(Vs11) / self.w * self.l_factor,
                        label='mf.VCoil_core (Core Exc.)', marker='o', edgecolor='black', facecolors='none', s=30)

            ax1.scatter(self.f, np.real(Vs12) * self.r_factor,
                        label='mf.VCoil_Sheath (Core Exc.)', marker='o', edgecolor='darkgreen', facecolors='none', s=30)
            
            ax2.scatter(self.f, np.imag(Vs12) / self.w * self.l_factor,
                        label='mf.VCoil_Sheath (Core Exc.)', marker='o', edgecolor='darkgreen', facecolors='none', s=30)

            ax1.scatter(self.f, np.real(Vs21) * self.r_factor,
                        label='mf.VCoil_Core (Sheath Exc.)', marker='x', facecolors='darkgreen', s=12)
            
            ax2.scatter(self.f, np.imag(Vs21) / self.w * self.l_factor,
                        label='mf.VCoil_Core (Sheath Exc.)', marker='x', facecolors='darkgreen', s=12)

            ax1.scatter(self.f, np.real(Vs22) * self.r_factor,
                        label='mf.VCoil_Sheath (Sheath Exc.)', marker='x', facecolors='darkblue', s=12)
            
            ax2.scatter(self.f, np.imag(Vs22) / self.w * self.l_factor,
                        label='mf.VCoil_Sheath (Sheath Exc.)', marker='x', facecolors='darkblue', s=12)

        # core self-impedance
        zi = self.pul_01['analytical']['internal_matrices']['impedance_matrix']        
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
        ax1.set_title(r'P.u.l. internal resistance matrix, $R$')

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
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='internal_impedance_matrix_js_method')

    def internal_impedance_matrix_energy_method(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(r'P.u.l. internal impedance matrix of a single core cable with Loss-Energy Method [Yin, 1990]', fontsize=12)

        if 'comsol' in self.pul_0 and self.cmsl is not None:
            core, sheath = self.cmsl['cmsl_core_exc'], self.cmsl['cmsl_sheath_exc']

            ax1.scatter(self.f, (core['r11']+core['r2i']) * self.r_factor,
                         label='mf.r11+mf.r2i (Core Exc.)', marker='o', edgecolor='black', facecolors='none', s=30)
            
            ax2.scatter(self.f, (core['l11']+core['l2i']+core['l12']+core['l13']) * self.l_factor,
                         label='mf.L11+mf.L2i+mf.L12+mf.L13 (Core Exc.)', marker='o', edgecolor='black', facecolors='none', s=30)

            ax1.scatter(self.f, 0.5 * core['r2i'] * self.r_factor,
                        label='mf.r2i (Core Exc.)', marker='o', edgecolor='darkgreen', facecolors='none', s=12)
            
            ax2.scatter(self.f, 0.5 * core['l2i'] * self.l_factor,
                        label='mf.L2i (Core Exc.)', marker='o', edgecolor='darkgreen', facecolors='none', s=12)

            ax1.scatter(self.f, sheath['r2i'] * self.r_factor,
                        label='mf.r2i (Sheath Exc.)', marker='x', facecolors='darkblue', s=12)
            
            ax2.scatter(self.f, sheath['l2i'] * self.l_factor,
                        label='mf.L2i (Sheath Exc.)', marker='x', facecolors='darkblue', s=12)

        # core self-impedance
        zi = self.pul_01['analytical']['internal_matrices']['impedance_matrix']
        ax1.plot(self.f, np.real(zi[:, 0, 0]) * self.r_factor, label=r'$R_{cc}$: core self-resistance', linestyle='-', color='black', linewidth=1.0) 
        ax2.plot(self.f, np.imag(zi[:, 0, 0]) / self.w * self.l_factor, label=r'$L_{cc}$: core self-inductance', linestyle='-', color='black', linewidth=1.0) 
        
        # mutual impedance between the core and sheath
        ax1.plot(self.f, np.real(zi[:, 0, 1]) * self.r_factor,
                 label=r'$R_{cs}$: mutual resistance between the core and sheath', linestyle='-', color='darkgreen', linewidth=1.0) 
        
        ax2.plot(self.f, np.imag(zi[:, 0, 1]) / self.w * self.l_factor,
                 label=r'$L_{cs}$: mutual inductance between the core and sheath', linestyle='-', color='darkgreen', linewidth=1.0) 

        # sheath self-impedance
        ax1.plot(self.f, np.real(zi[:, 1, 1]) * self.r_factor,
                 label=r'$R_{ss}$: sheath self-resistance', linestyle='-', color='darkblue', linewidth=1.0) 
        
        ax2.plot(self.f, np.imag(zi[:, 1, 1]) / self.w * self.l_factor,
                 label=r'$L_{ss}$: sheath self-inductance', linestyle='-', color='darkblue', linewidth=1.0) 

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
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='internal_impedance_matrix_energy_method')

    def internal_admittance_elements(self):
        """
        Plots the elements of the shunt internal capacitance matrix in a 1x3 subplot figure.
        Each subplot is dedicated to an element: Ccc, Ccs, and Css.
        """
        # Create a 1x3 subplot figure, with a wider size and a shared X-axis
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=self.figsize, sharex=True)
        fig.suptitle('P.u.l. Shunt Internal Capacitance Matrix \n Single core cable installed in HDPE tube [Lafaia, 2015]', fontsize=10)

        # --- Data Preparation ---
        # The capacitance matrix is frequency-independent, so we extract the 2D matrix
        # and then broadcast it across the frequency axis for plotting.
        cap_comsol = self.cap_matrix_energy_method * self.c_factor
        cap_case01 = self.pul_01['analytical']['internal_matrices']['capacitance_matrix'] * self.c_factor
        cap_case02 = self.pul_02['analytical']['internal_matrices']['capacitance_matrix'] * self.c_factor
        cap_case31 = self.pul_3['analytical']['internal_matrices']['capacitance_matrix'] * self.c_factor

        # --- Plotting Ccc (Core Self-Capacitance) on ax1 ---
        ax1.scatter(self.f, np.full_like(self.f, cap_comsol[0, 0]), label='Case 0', marker='o', c='black', s=5, zorder=3)
        ax1.scatter(self.f, np.full_like(self.f, cap_case01[0, 0]), label='Case 1', marker='o', facecolors='none', edgecolors='black', s=20)
        ax1.scatter(self.f, np.full_like(self.f, cap_case02[0, 0]), label='Case 2', marker='s', facecolors='none', edgecolors='darkgreen', s=25)
        ax1.scatter(self.f, np.full_like(self.f, cap_case31[0, 0]), label='Case 3', marker='^', facecolors='none', edgecolors='darkblue', s=30)
        ax1.text(1e3, 0.15 + cap_comsol[0, 0], f'$C_{{cc}} = {cap_comsol[0, 0]:.4f} \mu F/km$', 
                fontsize='small', color='black', va='center', ha='center', bbox=dict(facecolor='white', edgecolor='none', pad=2.0))
        
        # --- Plotting Ccs (Core-Sheath Mutual Capacitance) on ax2 ---
        ax2.scatter(self.f, np.full_like(self.f, cap_comsol[0, 1]), label='Case 0', marker='o', c='black', s=5, zorder=3)
        ax2.scatter(self.f, np.full_like(self.f, cap_case01[0, 1]), label='Case 1', marker='o', facecolors='none', edgecolors='black', s=20)
        ax2.scatter(self.f, np.full_like(self.f, cap_case02[0, 1]), label='Case 2', marker='s', facecolors='none', edgecolors='darkgreen', s=25)
        ax2.scatter(self.f, np.full_like(self.f, cap_case31[0, 1]), label='Case 3', marker='^', facecolors='none', edgecolors='darkblue', s=30)
        ax2.text(1e3, 0.15 + cap_comsol[0, 1], f'$C_{{cs}} = {cap_comsol[0, 1]:.4f} \mu F/km$', 
                fontsize='small', color='black', va='center', ha='center', bbox=dict(facecolor='white', edgecolor='none', pad=2.0))

        # --- Plotting Css (Sheath Self-Capacitance) on ax3 ---
        ax3.scatter(self.f, np.full_like(self.f, cap_comsol[1, 1]),
                     label=fr'Case 0: {cap_comsol[1, 1]:.4f} $\mu F/km$', marker='o', c='black', s=5, zorder=3)
        ax3.scatter(self.f, np.full_like(self.f, cap_case01[1, 1]),
                     label=fr'Case 1: {cap_case01[1, 1]:.4f} $\mu F/km$', marker='o', facecolors='none', edgecolors='black', s=20)
        ax3.scatter(self.f, np.full_like(self.f, cap_case02[1, 1]),
                     label=fr'Case 2: {cap_case02[1, 1]:.4f} $\mu F/km$', marker='s', facecolors='none', edgecolors='darkgreen', s=25)
        ax3.scatter(self.f, np.full_like(self.f, cap_case31[1, 1]),
                     label=fr'Case 3: {cap_case31[1, 1]:.4f} $\mu F/km$', marker='^', facecolors='none', edgecolors='darkblue', s=30)

        # --- Axis Configuration ---
        # Common settings for all subplots
        for ax in [ax1, ax2, ax3]:
            ax.set_xscale('log')
            ax.set_xlim(1e0, 1e6)
            ax.set_xlabel('Frequency (Hz)')
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax.legend(fontsize='small')

        # Specific settings for ax1 (Ccc)
        ax1.set_ylabel(r'$C_{cc} \quad (\mu F/km)$')
        ax1.set_title(r'Core Self-Capacitance ($C_{cc}$)', fontsize=10)
        ax1.set_ylim(-0.25, 2.0)

        # Specific settings for ax2 (Ccs)
        ax2.set_ylabel(r'$C_{cs} \quad (\mu F/km)$')
        ax2.set_title(r'Core-Sheath Mutual Capacitance ($C_{cs}$)', fontsize=10)
        ax2.set_ylim(-2.0, 0.25)

        # Specific settings for ax3 (Css)
        ax3.set_ylabel(r'$C_{ss} \quad (\mu F/km)$')
        ax3.set_title(r'Sheath Self-Capacitance ($C_{ss}$)', fontsize=10)
        ax3.set_ylim(-0.25, 2.0)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        if self.autoSave:
            save_figure(fig, self.results_dir, base_filename='internal_admittance_matrix')
