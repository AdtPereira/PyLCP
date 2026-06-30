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

_MATRIX_ELEMENTS = {
    'cc': {'p': 0, 'q': 0},
    'cs': {'p': 0, 'q': 1},
    'ss': {'p': 1, 'q': 1},
}

PLOT_CONFIG = {

    'internal_impedance_matrix': {
        'suptitle': r'Series Impedance Matrix — GMD case 3.1 ($r_0 = r_5$) vs. COMSOL ($J_s$) [Lafaia, 2015; Yin, 1990]',
        'data_series': [
            # --- Analytical: Scenario 1 (underground, ignoring HDPE) ---
            # {
            #     'source': 'analytical',
            #     'scenario_key': '1',
            #     'data_path': ['internal_matrices', 'impedance_matrix'],
            #     'plot_style': 'line',
            #     'series': {
            #         'cc': {'p': 0, 'q': 0, 'color': 'black', 'linestyle': '-',  'linewidth': 1.5, 'label': r'$Z_{cc}$ Sc.1 (ignoring HDPE)'},
            #         'cs': {'p': 0, 'q': 1, 'color': 'black', 'linestyle': ':',  'linewidth': 1.5, 'label': r'$Z_{cs}$ Sc.1'},
            #         'ss': {'p': 1, 'q': 1, 'color': 'black', 'linestyle': '--', 'linewidth': 1.5, 'label': r'$Z_{ss}$ Sc.1'},
            #     },
            # },
            # --- Analytical: Scenario 2 (area-weighted equiv. permittivity) ---
            # {
            #     'source': 'analytical',
            #     'scenario_key': '2',
            #     'data_path': ['internal_matrices', 'impedance_matrix'],
            #     'plot_style': 'line',
            #     'series': {
            #         'cc': {'p': 0, 'q': 0, 'color': 'blue', 'linestyle': '-',  'linewidth': 1.5, 'label': r'$Z_{cc}$ Sc.2 (area-weighted $\varepsilon_r$)'},
            #         'cs': {'p': 0, 'q': 1, 'color': 'blue', 'linestyle': ':',  'linewidth': 1.5, 'label': r'$Z_{cs}$ Sc.2'},
            #         'ss': {'p': 1, 'q': 1, 'color': 'blue', 'linestyle': '--', 'linewidth': 1.5, 'label': r'$Z_{ss}$ Sc.2'},
            #     },
            # },
            # --- Analytical: Scenario 3 (GMD case 3.1) ---
            {
                'source': 'analytical',
                'scenario_key': '3',
                'data_path': ['internal_matrices', 'impedance_matrix'],
                'plot_style': 'line',
                'series': {
                    'cc': {'p': 0, 'q': 0, 'color': 'red', 'linestyle': '-',  'linewidth': 1.5, 'label': r'$Z_{cc}$ GMD case 3.1'},
                    'cs': {'p': 0, 'q': 1, 'color': 'red', 'linestyle': ':',  'linewidth': 1.5, 'label': r'$Z_{cs}$ GMD case 3.1'},
                    'ss': {'p': 1, 'q': 1, 'color': 'red', 'linestyle': '--', 'linewidth': 1.5, 'label': r'$Z_{ss}$ GMD case 3.1'},
                },
            },
            # --- COMSOL: Js method matrix ---
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['internal_impedance_matrix', 'js_method'],
                'plot_style': 'scatter',
                'series': {
                    'cc': {'p': 0, 'q': 0, 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': r'COMSOL ($J_s$)'},
                    'cs': {'p': 0, 'q': 1, 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': ''},
                    'ss': {'p': 1, 'q': 1, 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 2, 'label': ''},
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
                    'Zcs': {'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black',     'zorder': 2, 'label': r'$Z_{cs}=Z_{11}+Z_{12}+Z_{2i}$ (COMSOL)'},
                    'z11': {'marker': 's', 's': 20, 'facecolors': 'none', 'edgecolors': 'red',       'zorder': 2, 'label': r'$Z_{11}$: core outer surface'},
                    'z12': {'marker': '^', 's': 20, 'facecolors': 'none', 'edgecolors': 'darkgreen', 'zorder': 2, 'label': r'$Z_{12}$: core insulator'},
                    'z2i': {'marker': 'v', 's': 20, 'facecolors': 'none', 'edgecolors': 'blue',      'zorder': 2, 'label': r'$Z_{2i}$: sheath inner surface'},
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
        'suptitle': r'Series Impedance elements of a single core cable installed in HDPE tube [Lafaia, 2015; Yin, 1990]',
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
