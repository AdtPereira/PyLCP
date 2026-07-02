import numpy as np

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

PLOT_CONFIG = {

    'internal_impedance_matrix': {
        'suptitle': r'Series Impedance Matrix of a single core cable installed in HDPE tube [Lafaia, 2015; Yin, 1990]',
        'data_series': [
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['internal_impedance_matrix', 'js_method'],
                'plot_style': 'scatter',
                'series': {
                    'cc': {'p': 0, 'q': 0, 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': r'$Z_{cc}$ (COMSOL $J_s$)'},
                    'cs': {'p': 0, 'q': 1, 'marker': 's', 's': 20, 'facecolors': 'none', 'edgecolors': 'red',   'zorder': 2, 'label': r'$Z_{cs}$ (COMSOL $J_s$)'},
                    'ss': {'p': 1, 'q': 1, 'marker': '^', 's': 20, 'facecolors': 'none', 'edgecolors': 'blue',  'zorder': 2, 'label': r'$Z_{ss}$ (COMSOL $J_s$)'},
                },
            },
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['internal_impedance_matrix', 'energy_method'],
                'plot_style': 'line',
                'series': {
                    'cc': {'p': 0, 'q': 0, 'color': 'black', 'linestyle': '-',  'linewidth': 1.5, 'label': r'$Z_{cc}$ (Energy method)'},
                    'cs': {'p': 0, 'q': 1, 'color': 'red',   'linestyle': '--', 'linewidth': 1.5, 'label': r'$Z_{cs}$ (Energy method)'},
                    'ss': {'p': 1, 'q': 1, 'color': 'blue',  'linestyle': ':',  'linewidth': 1.5, 'label': r'$Z_{ss}$ (Energy method)'},
                },
            },
        ],
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1e-3, 1e1),
            'legend': True,
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0, 0.35),
            'legend': False,
        },
    },

    'coaxial_cable': {
        'suptitle': r'P.u.l. series impedance of a single core cable with sheath path return [Patel, 2014; Ametani, 2015]',
        'data_series': [
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['coaxial_cable_impedance'],
                'plot_style': 'scatter',
                'series': {
                    'Zcs': {'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black',    'zorder': 2, 'label': r'$Z_{cs}=Z_{11}+Z_{12}+Z_{2i}$ (COMSOL)'},
                    'z11': {'marker': 's', 's': 20, 'facecolors': 'none', 'edgecolors': 'red',      'zorder': 2, 'label': r'$Z_{11}$: core outer surface'},
                    'z12': {'marker': '^', 's': 20, 'facecolors': 'none', 'edgecolors': 'darkgreen','zorder': 2, 'label': r'$Z_{12}$: core insulator'},
                    'z2i': {'marker': 'v', 's': 20, 'facecolors': 'none', 'edgecolors': 'blue',     'zorder': 2, 'label': r'$Z_{2i}$: sheath inner surface'},
                },
            },
        ],
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1e-3, 1e1),
            'legend': True,
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0, 0.35),
            'legend': False,
        },
    },

    'internal_impedance_elements': {
        'suptitle': r'Series Impedance Matrix of a single core cable installed in HDPE tube [Lafaia, 2015; Yin, 1990]',
        'data_series': [
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['internal_impedance_elements'],
                'plot_style': 'scatter',
                'series': {
                    'self_core_js':    {'marker': 'o', 's': 45, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': r'$Z_{cc}$ ($J_s$ Method)'},
                    'self_sheath_js':  {'marker': 'o', 's': 45, 'facecolors': 'none', 'edgecolors': 'blue',  'zorder': 2, 'label': r'$Z_{ss}$ ($J_s$ Method)'},
                    'mutual_js':       {'marker': 'o', 's': 45, 'facecolors': 'none', 'edgecolors': 'red',   'zorder': 2, 'label': r'$Z_{cs}$ ($J_s$ Method)'},
                    'self_core_energy':   {'marker': 'x', 's': 20, 'color': 'black', 'zorder': 3, 'label': r'$Z_{cc}$ (Energy Method)'},
                    'self_sheath_energy': {'marker': 'x', 's': 20, 'color': 'blue',  'zorder': 3, 'label': r'$Z_{ss}$ (Energy Method)'},
                    'mutual_energy':      {'marker': 'x', 's': 20, 'color': 'red',   'zorder': 3, 'label': r'$Z_{cs}$ (Energy Method)'},
                },
            },
        ],
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1e-3, 1e1),
            'legend': True,
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0, 0.35),
            'legend': False,
        },
    },
}
