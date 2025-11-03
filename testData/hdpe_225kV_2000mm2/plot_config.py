import numpy as np

comsol_layout_template = [
    {
        'base_key': 'rho_g_100_epsr1_1_mf',
        'label': r'COMSOL ($\rho_e=100, \epsilon_r=1$)',
        'marker': 'o', 's': 35, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 10
    },
    {
        'base_key': 'rho_g_100_epsr1_20_mf',
        'label': r'COMSOL ($\rho_e=100, \epsilon_r=20$)',
        'marker': 'o', 's': 15, 'facecolors': 'black', 'edgecolors': 'none', 'zorder': 10
    },
    {
        'base_key': 'rho_g_500_epsr1_1_mf',
        'label': r'COMSOL ($\rho_e=500, \epsilon_r=1$)',
        'marker': 's', 's': 15, 'color': 'black', 'zorder': 10
    }
]

PLOT_TPL_REAL = {
    'component': 'real',
    'scale': 1e0,
    'xscale': 'log',
    'yscale': 'linear',
}

PLOT_TPL_IMAG = {
    'component': 'imag',
    'scale': 1e0,
    'xscale': 'log',
    'yscale': 'linear',
}

PLOT_TPL_RESISTANCE = {
    'component': 'real',
    'scale': 1e3,
    'xscale': 'log',
    'yscale': 'log',
}

PLOT_TPL_INDUCTANCE = {
    'component': 'imag_div_w',
    'scale': 1e6,
    'xscale': 'log',
    'yscale': 'linear',
}

PLOT_TPL_CONDUCTANCE = {
    'component': 'real',
    'scale': 1e3,
    'xscale': 'log',
    'yscale': 'linear',
}

PLOT_TPL_CAPACITANCE = {
    'component': 'imag_div_w',
    'scale': 1e9,
    'xscale': 'log',
    'yscale': 'linear',
}

PLOT_CONFIG = {
    'internal_impedance_matrix': {
        'suptitle': r'Single core cable installed in HDPE tube [Lafaia, 2015] and $J_s$ Method [Yin, 1990]',
        'series_to_plot': [
            {
                'key': '1',
                'type': {
                    'r_cc': {'p': 0, 'q': 0, 'color': 'black', 'linestyle': '-',  'linewidth': 1.5, 'label': r'$R_{cc}$: core self-resistance'},
                    'r_cs': {'p': 0, 'q': 1, 'color': 'black', 'linestyle': ':',  'linewidth': 1.5, 'label': r'$R_{cs}$: core-sheath mutual resistance'},
                    'r_ss': {'p': 1, 'q': 1, 'color': 'black', 'linestyle': '--', 'linewidth': 1.5, 'label': r'$R_{ss}$: sheath self-resistance'},
                }
            },
        ],
        'data_path': ['scenarios', '{key}', 'internal_matrices', 'impedance_matrix'],
        # 'comsol_series_to_plot': comsol_layout_template,
        # 'comsol_matrix_key': 'series_impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-3, 1E1),
            'legend': True,
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0, 0.35),
            'legend': False,
        }
    },

    'internal_impedance_elements': {
        'suptitle': r'P.u.l. series impedance of a single core cable with sheath path return [Patel, 2014; AMETANI, 2015]',
        'series_to_plot': [
            {
                'key': '1',
                'type': {
                    'Zcs': {'color': 'black',     'linestyle': '-', 'linewidth': 1.5, 'label': r'$R_{cs}=R_{11}+R_{12}+R_{2i}$'},
                    'z11': {'color': 'red',       'linestyle': ':', 'linewidth': 1.5, 'label': r'$R_{11}$: internal resistance of core outer surface'},
                    'z12': {'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.5, 'label': r'$R_{12}$: core outer insulator resistance'},
                    'z2i': {'color': 'blue',      'linestyle': ':', 'linewidth': 1.5, 'label': r'$R_{2i}$: internal resistance of sheath inner surface'},
                }
            },
        ],
        'data_path': ['scenarios', '{key}', 'internal_parameters', 'zcs'],
        # 'comsol_series_to_plot': comsol_layout_template,
        # 'comsol_matrix_key': 'series_impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-3, 1E1),
            'legend': True,
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0, 0.35),
            'legend': False,
        }
    },
}