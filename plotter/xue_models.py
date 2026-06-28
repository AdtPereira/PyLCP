import numpy as np
import matplotlib.pyplot as plt


class XueModels:
    """
    Plotter for the Xue overhead-line model results.
    Expects pul_data with flat scenario keys (p100, p100_carson, etc.)
    and pul_data[key]['series_impedance_matrix'] of shape (n_freq, N, N).
    """

    def __init__(self, pul_data: dict):
        self.pul_data = pul_data
        self.f = pul_data['frequencies']
        self.w = 2 * np.pi * self.f
        self.figsize = (12, 5)
        self.xlim = (self.f[0], self.f[-1])

        self.overhead_line = [
            {'key': 'p100',         'label': r'Nakagawa — $\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=1$',   'color': 'black',  'linestyle': '-'},
            {'key': 'p100_er20',    'label': r'Nakagawa — $\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=20$',  'color': 'black',  'linestyle': '--'},
            {'key': 'p2000',        'label': r'Nakagawa — $\rho_e=2000\;\Omega\mathrm{m},\;\epsilon_r=1$',  'color': 'black',  'linestyle': '-.'},
            {'key': 'p100_carson',  'label': r'Carson — $\rho_e=100\;\Omega\mathrm{m}$',                    'color': 'red',    'linestyle': '-'},
            {'key': 'p2000_carson', 'label': r'Carson — $\rho_e=2000\;\Omega\mathrm{m}$',                   'color': 'red',    'linestyle': '-.'},
        ]

    def overhead_series_impedance_matrix(self, p: int = 0, q: int = 0):
        """Plots series resistance (left) and inductance (right) for all scenarios."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize, sharey=False)
        fig.suptitle(
            r'P.u.l. Series Impedance — Single overhead line [Xue, 2012]',
            fontsize=12, y=0.98
        )

        for series in self.overhead_line:
            key = series['key']
            if key not in self.pul_data:
                continue
            zs = self.pul_data[key]['series_impedance_matrix']
            style = {k: v for k, v in series.items() if k != 'key'}
            ax1.plot(self.f, zs[:, p, q].real * 1e3, **style)
            ax2.plot(self.f, zs[:, p, q].imag / self.w * 1e6, **style)

        for ax, ylabel in [(ax1, r'$R_s\;(\Omega/\mathrm{km})$'), (ax2, r'$L_s\;(\mathrm{mH/km})$')]:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(self.xlim)
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel(ylabel)
            ax.legend(fontsize='small')
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
