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

    def internal_impedance_matrix(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        r_factor = 1e3  # Convert Ohm/m to Ohm/km
        l_factor = 1e6  # Convert H/m to mH/km

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle('Fig. 2.6: P.u.l. internal impedance matrix of a coaxial cable [Ametani, 2015]', fontsize=12)

        # The slice [:, 0, 0] extracts the (0,0) element for all frequencies.        
        zi = self.pul_data['analytical']['internal_impedance_matrix']['approximation']
        
        # core self-impedance
        ax1.plot(self.f, np.real(zi[:, 0, 0]) * r_factor,
                    label=r'$R_{cc}$: core self-resistance', linestyle=':', color='black', linewidth=2) 
        ax2.plot(self.f, np.imag(zi[:, 0, 0]) / self.w * l_factor,
                    label=r'$L_{cc}$: core self-inductance', linestyle=':', color='black', linewidth=2) 
        
        # mutual impedance between the core and sheath
        ax1.plot(self.f, np.real(zi[:, 0, 1]) * r_factor,
                    label=r'$R_{cs}$: mutual resistance between the core and sheath', linestyle=':', color='darkgreen', linewidth=2) 
        ax2.plot(self.f, np.imag(zi[:, 0, 1]) / self.w * l_factor,
                    label=r'$L_{cs}$: mutual inductance between the core and sheath', linestyle=':', color='darkgreen', linewidth=2) 

        # sheath self-impedance
        ax1.plot(self.f, np.real(zi[:, 1, 1]) * r_factor,
                    label=r'$R_{ss}$: sheath self-resistance', linestyle=':', color='darkblue', linewidth=2) 
        ax2.plot(self.f, np.imag(zi[:, 1, 1]) / self.w * l_factor,
                    label=r'$L_{ss}$: sheath self-inductance', linestyle=':', color='darkblue', linewidth=2) 

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1e0, 1e6)
        ax1.set_ylim(1e-3, 1e1)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        # ax1.set_ylabel(r'$R_{cs}$ $(\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title('Internal Resistance, $R_{cs}$')

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1e0, 1e6)
        ax2.set_ylim(0, 0.2)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        # ax2.set_ylabel(r'$L_{cs}$ (mH/km)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title('Internal Inductance, $L_{cs}$')
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    def series_impedance_matrix(self):
        """
        Generic method to create a 1x2 subplot for series resistance (left)
        and series inductance (right) based on a configuration key.
        """
        r_factor = 1e3  # Convert Ohm/m to Ohm/km
        l_factor = 1e6  # Convert H/m to mH/km

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle('Fig. 2.6: P.u.l. series impedance of a coaxial cable [Patel, 2014]', fontsize=12)

        # The slice [:, 0, 0] extracts the (0,0) element for all frequencies.        
        zcs = self.pul_data['analytical']['internal_parameters']['approximation']['zcs']
        rs_mom = self.pul_data['numerical']['series_resistance_matrix'][:, 0, 0] 
        ls_mom = self.pul_data['numerical']['series_inductance_matrix'][:, 0, 0] 
        
        z11 = zcs['z11']   # internal impedance of core outer surface
        z12 = zcs['z12']   # core outer insulator impedance
        z2i = zcs['z2i']   # internal impedance of sheath inner surface
        zcs = z11 + z12 + z2i   # Equivalent single conductor impedance. Eq. (2.10a) [AMETANI, 2015]

        # Resistance in Ohm/km
        # Plotar dados do COMSOL se existirem no dicionário pul_data
        if 'comsol' in self.pul_data and self.pul_data['comsol'] is not None:
            comsol_df = self.pul_data['comsol']
            freq_comsol = comsol_df['Frequency (Hz)']
            zs_comsol = comsol_df['Zs (Ω/m)']

            ax1.scatter(freq_comsol, np.real(zs_comsol) * r_factor,
                        label='COMSOL (mf/ec)', marker='x', facecolors='k', s=12, zorder=2)

            ax2.scatter(freq_comsol, np.imag(zs_comsol) / (2*np.pi*freq_comsol) * l_factor,
                        label='COMSOL (mf/ec)', marker='x', facecolors='k', s=12, zorder=2)

        ax1.scatter(self.f_mom, rs_mom * r_factor,
                 label='MoM-SO', marker='o', facecolors='none', edgecolors='k', s=40, zorder=2)
        ax1.plot(self.f, np.real(z11) * r_factor,
                 label=r'$R_{11}$: internal resistance of core outer surface', linestyle=':', color='black', linewidth=2) 
        ax1.plot(self.f, np.real(z12) * r_factor,
                 label=r'$R_{12}$: core outer insulator resistance', linestyle=':', color='darkgreen', linewidth=2) 
        ax1.plot(self.f, np.real(z2i) * r_factor,
                 label=r'$R_{2i}$: internal resistance of sheath inner surface', linestyle=':', color='darkblue', linewidth=2) 
        ax1.plot(self.f, np.real(zcs) * r_factor,
                 label=r'$R_{cs}=R_{11}+R_{12}+R_{2i}$ [1]', linestyle='--', color='red', linewidth=1.0, zorder=1)

        # Inductance in mH/km
        ax2.scatter(self.f_mom, ls_mom * l_factor,
                 label='MoM-SO', marker='o', facecolors='none', edgecolors='k', s=40, zorder=2)
        ax2.plot(self.f, np.imag(z11) / self.w * l_factor,
                 label=r'$L_{11}$: internal inductance of core outer surface', linestyle=':', color='black', linewidth=2) 
        ax2.plot(self.f, np.imag(z12) / self.w * l_factor,
                 label=r'$L_{12}$: core outer insulator inductance', linestyle=':', color='darkgreen', linewidth=2) 
        ax2.plot(self.f, np.imag(z2i) / self.w * l_factor,
                 label=r'$L_{2i}$: internal inductance of sheath inner surface', linestyle=':', color='darkblue', linewidth=2) 
        ax2.plot(self.f,  np.imag(zcs) / self.w * l_factor,
                 label=r'$L_{cs}=L_{11}+L_{12}+L_{2i}$ [1]', linestyle='--', color='red', linewidth=1.0, zorder=1)    

        # Configure left subplot (Resistance)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(1e0, 1e6)
        ax1.set_ylim(1e-2, 1e1)
        ax1.legend(fontsize='small')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel(r'$R_{cs}$ $(\Omega/km)$')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.set_title('Series Resistance, $R_{cs}$')

        # Configure right subplot (Inductance)
        ax2.set_xscale('log')
        ax2.set_xlim(1e0, 1e6)
        ax2.set_ylim(0, 0.2)
        ax2.legend(fontsize='small')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel(r'$L_{cs}$ (mH/km)')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax2.set_title('Series Inductance, $L_{cs}$')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
