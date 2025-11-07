import numpy as np

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
    'partial_internal_impedance': {
        'suptitle': r'Series Impedance Matrix of a single core cable',
        'series_to_plot': [
            {
                'key': '1',
                'type': {
                    'cc': {'p': 0, 'q': 0, 'color': 'black', 'linestyle': '-',  'linewidth': 1.5, 'label': r'$Z_{cc}$: core self-resistance'},
                    'cs': {'p': 0, 'q': 1, 'color': 'black', 'linestyle': ':',  'linewidth': 1.5, 'label': r'$Z_{cs}$: core-sheath mutual resistance'},
                    'ss': {'p': 1, 'q': 1, 'color': 'black', 'linestyle': '--', 'linewidth': 1.5, 'label': r'$Z_{ss}$: sheath self-resistance'},
                }
            },
        ],
        'path': ['internal_matrices', 'impedance_matrix'],
        'mom_so_series_to_plot': [
            {
                'key': '1',
                'type': {
                    'cc': {'p': 0, 'q': 0, 'marker': 'o', 's': 15, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': 'MoM-SO'},
                    'cs': {'p': 0, 'q': 1, 'marker': 'o', 's': 15, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': ''},
                    'ss': {'p': 1, 'q': 1, 'marker': 'o', 's': 15, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': ''},
                }
            },
        ],
        'mom_so_matrix_key': 'partial_internal_impedance',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-3, 1E1),
            'legend': True,
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0, 1.0),
            'legend': False,
        }
    },

    'coaxial_cable_impedance': {
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
        'path': ['internal_parameters', 'zcs'],
        'comsol_series_to_plot': [
            {
                'key': '1',
                'type': {
                    'Zcs': {'marker': 'o', 's': 12, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': r'COMSOL (mf)'},
                    'z11': {'marker': 'o', 's': 12, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': ''},
                    'z12': {'marker': 'o', 's': 12, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': ''},
                    'z2i': {'marker': 'o', 's': 12, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': ''},
                }
            },
        ],
        'comsol_matrix_key': 'coaxial_cable_impedance',
        'mom_so_series_to_plot': [
            {
                'key': '1',
                'type': {
                    'Zcs': {'marker': 'o', 's': 50, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 3, 'label': 'MoM-SO'},
                }
            },
        ],
        'mom_so_matrix_key': 'coaxial_cable_impedance',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-2, 1E1),
            'legend': True,
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0, 0.2),
            'legend': False,
        }
    },

    'internal_impedance_elements': {
        'suptitle': r'Series Impedance Matrix of a single core cable',
        'series_to_plot': [],
        'path': [],
        'comsol_series_to_plot': [
            {
                'key': '1',
                'type': {
                    'core_js':    {'label': r'$Z_{cc}$ ($J_s$ Method)', 'marker': 's', 's': 45, 'facecolors': 'none',  'zorder': 2, 'edgecolors': 'black'},
                    'sheath_js':  {'label': r'$Z_{ss}$ ($J_s$ Method)', 'marker': 's', 's': 45, 'facecolors': 'none',  'zorder': 2, 'edgecolors': 'blue'},
                    'mutual1_js': {'label': r'$Z_{cs}$ ($J_s$ Method)', 'marker': 's', 's': 45, 'facecolors': 'none',  'zorder': 2, 'edgecolors': 'red'},
                    
                    'core_energy':   {'label': r'$Z_{cc}$ (Energy Method)', 'marker': 'x', 's': 10, 'facecolors': 'black', 'zorder': 3},
                    'sheath_energy': {'label': r'$Z_{ss}$ (Energy Method)', 'marker': 'x', 's': 10, 'facecolors': 'blue', 'zorder': 3},
                    'mutual_energy': {'label': r'$Z_{cs}$ (Energy Method)', 'marker': 'x', 's': 10, 'facecolors': 'red', 'zorder': 3},

                }
            },
        ],
        'comsol_matrix_key': 'internal_impedance_elements',
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