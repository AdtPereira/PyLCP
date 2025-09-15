import numpy as np
import matplotlib.pyplot as plt
from utils.case_utils import *

class PatelModels:
    """
    A highly refactored class to handle plotting for the Xue model results.
    It uses a configuration-driven approach to generate complex subplot figures.
    """
    def __init__(self, pul_data):
        self.pul_data = pul_data
        self.f = pul_data['analytical']['frequencies']
        self.f_mom = pul_data['numerical']['frequencies']
        self.w = 2 * np.pi * self.f
        self.figsize = (12, 5)
        self.r_factor = 1e3  # Convert Ohm/m to Ohm/km
        self.l_factor = 1e6  # Convert H/m to mH/km

    def internal_impedance_elements(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle('Fig. 2.6: P.u.l. series impedance of a coaxial cable with sheath path return [Patel, 2014]', fontsize=12)

        # The slice [:, 0, 0] extracts the (0,0) element for all frequencies.        
        zcs = self.pul_data['analytical']['internal_parameters']['hybrid']['zcs']
        rs = self.pul_data['numerical']['series_resistance_matrix'][:, 0, 0] 
        rs = self.pul_data['numerical']['series_resistance_matrix'][:, 0, 0] 
        ls = self.pul_data['numerical']['series_inductance_matrix'][:, 0, 0] 
        
        # Equivalent single conductor impedance. Eq. (2.10a) [AMETANI, 2015]
        Zcs = zcs['z11'] + zcs['z12'] + zcs['z2i'] 

        if 'comsol' in self.pul_data and self.pul_data['comsol'] is not None:
            cmsl_data = self.pul_data['comsol']['core_sheath_return']            
            ax1.scatter(cmsl_data['freq'], cmsl_data['coil_resistance'] * self.r_factor, label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)
            ax2.scatter(cmsl_data['freq'], cmsl_data['coil_inductance'] * self.l_factor, label='COMSOL (mf)', marker='o', facecolors='black', s=10, zorder=2)

            ax1.scatter(cmsl_data['freq'], cmsl_data['r11'] * self.r_factor, marker='o', facecolors='black', s=10, zorder=3)
            ax2.scatter(cmsl_data['freq'], cmsl_data['l11'] * self.l_factor, marker='o', facecolors='black', s=10, zorder=3)

            ax1.scatter(cmsl_data['freq'], cmsl_data['r2i'] * self.r_factor, marker='o', facecolors='darkblue', s=10, zorder=3)
            ax2.scatter(cmsl_data['freq'], cmsl_data['l2i'] * self.l_factor, marker='o', facecolors='darkblue', s=10, zorder=3)

            ax2.scatter(cmsl_data['freq'], cmsl_data['l12'] * self.l_factor, marker='o', facecolors='darkgreen', s=10, zorder=3)

        ax1.scatter(self.f_mom, rs * self.r_factor, label='MoM-SO', marker='o', facecolors='none', edgecolors='k', s=50, zorder=3)
        ax2.scatter(self.f_mom, ls * self.l_factor, label='MoM-SO', marker='o', facecolors='none', edgecolors='k', s=50, zorder=3)
        
        ax1.plot(self.f, np.real(zcs['z11']) * self.r_factor, label=r'$R_{11}$: internal resistance of core outer surface', linestyle='--', color='black', linewidth=1.0) 
        ax1.plot(self.f, np.real(zcs['z12']) * self.r_factor, label=r'$R_{12}$: core outer insulator resistance', linestyle='--', color='green', linewidth=1.0) 
        ax1.plot(self.f, np.real(zcs['z2i']) * self.r_factor, label=r'$R_{2i}$: internal resistance of sheath inner surface', linestyle='--', color='blue', linewidth=1.0) 
        ax1.plot(self.f, np.real(Zcs) * self.r_factor, label=r'$R_{cs}=R_{11}+R_{12}+R_{2i}$ [1]', linestyle='-', color='black', linewidth=1.0)
        
        ax2.plot(self.f, np.imag(zcs['z11']) / self.w * self.l_factor, label=r'$L_{11}$: internal inductance of core outer surface', linestyle='--', color='black', linewidth=1.0) 
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
        ax1.set_title('Internal Resistance, $R_{cs}$')

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1e0, 1e6)
        ax2.set_ylim(0, 0.2)
        # ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_{cs}$ (mH/km)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title('Internal Inductance, $L_{cs}$')
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def internal_impedance_matrix(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle('P.u.l. internal impedance matrix of a coaxial cable [Ametani, 2015]', fontsize=12)

        # The slice [:, 0, 0] extracts the (0,0) element for all frequencies.        
        zi = self.pul_data['analytical']['internal_impedance_matrix']['approximation']

        if 'comsol' in self.pul_data and self.pul_data['comsol'] is not None:
            core_exc = self.pul_data['comsol']['core']
            sheath_exc = self.pul_data['comsol']['sheath']

            ax1.scatter(core_exc['freq'], core_exc['core_coil_resistance'] * self.r_factor, label='mf.RCoil (Core Exc.)', marker='x', facecolors='black', s=15)
            ax2.scatter(core_exc['freq'], core_exc['core_coil_inductance'] * self.l_factor, label='mf.LCoil (Core Exc.)', marker='x', facecolors='black', s=15)

            ax1.scatter(sheath_exc['freq'], sheath_exc['r2i'] * self.r_factor, label='mf.r2i (Sheath Exc.)', marker='o', facecolors='darkblue', s=10)
            ax2.scatter(sheath_exc['freq'], sheath_exc['l2i'] * self.l_factor, label='mf.L2i (Sheath Exc.)', marker='o', facecolors='darkblue', s=10)

            ax1.scatter(core_exc['freq'], 0.5 * core_exc['r2i'] * self.r_factor, label='mf.r2i (Core Exc.)', marker='x', facecolors='darkgreen', s=15)
            ax2.scatter(core_exc['freq'], 0.5 * core_exc['l2i'] * self.l_factor, label='mf.L2i (Core Exc.)', marker='x', facecolors='darkgreen', s=15)

        # core self-impedance
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
        ax1.set_ylim(1e-2, 1e1)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$[R]$ $(\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title('Internal Resistance, $R_{cs}$')

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1e0, 1e6)
        ax2.set_ylim(0, 0.2)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$[L]$ (mH/km)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title('Internal Inductance, $L_{cs}$')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
