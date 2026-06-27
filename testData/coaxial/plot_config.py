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

PLOT_TPL_REACTANCE = {
    'component': 'imag',
    'scale': 1e3,
    'xscale': 'log',
    'yscale': 'log',
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

SCATTER_STYLES = {
    'analytical_1': {
        'color': 'blue',
        'linestyle': '--',
        'linewidth': 1.0,
        'zorder': 1,
    },
    'mom_so_1': {
        'marker': 'o',
        's': 40,  
        'edgecolors': 'black',
        'facecolors': 'none',
        'zorder': 2,
    },
    'comsol_1': {
        'marker': 'o',
        's': 6,
        'edgecolors': 'black',
        'facecolors': 'black',
        'zorder': 3,
    },
}

PLOT_CONFIG = {
    'internal_impedance_matrix_js_method': {
        'plot_function': '_plot_matricial_parameters',
        
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-3, 1E1),
            'legend': True,
        },

        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0.0, 0.2),
            'legend': False,
        },
        
        'suptitle': r'P.u.l. Internal impedance matrix of a coaxial (SC) cable with $J_s$ Method [Yin, 1990]',
        
        'data_series': [
            # --- COMSOL Formulation ---
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['internal_impedance_matrix', 'js_method'], 
                'plot_style': 'scatter',
                'series': { 
                    'cc': {'p': 0, 'q': 0, **SCATTER_STYLES['comsol_1'], 'label': 'Comsol, mf'},
                    'cs': {'p': 0, 'q': 1, **SCATTER_STYLES['comsol_1']},
                    'ss': {'p': 1, 'q': 1, **SCATTER_STYLES['comsol_1']},
                }
            },
            
            # --- Analytical Formulation ---
            {
                'source': 'analytical',
                'scenario_key': '1',
                'data_path': ['internal_impedance_matrix'], 
                'plot_style': 'line',
                'series': { 
                    'cc': {'p': 0, 'q': 0, 'color': 'black', 'linestyle': '-',  'linewidth': 1.0, 'label': r'$R_{cc}$: core self-resistance'},
                    'cs': {'p': 0, 'q': 1, 'color': 'black', 'linestyle': ':',  'linewidth': 1.0, 'label': r'$R_{cs}$: core-sheath mutual resistance'},
                    'ss': {'p': 1, 'q': 1, 'color': 'black', 'linestyle': '--', 'linewidth': 1.0, 'label': r'$R_{ss}$: sheath self-resistance'},
                }
            },
        ],
    },
    
    'internal_impedance_matrix_energy_method': {
        'plot_function': '_plot_matricial_parameters',
        
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-3, 1E1),
            'legend': True,
        },

        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0.0, 0.2),
            'legend': False,
        },
        
        'suptitle': r'P.u.l. Internal impedance matrix of a coaxial (SC) cable with Energy Method [Yin, 1990]',
        
        'data_series': [
            # --- COMSOL Formulation ---
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['internal_impedance_matrix', 'energy_method'], 
                'plot_style': 'scatter',
                'series': { 
                    'cc': {'p': 0, 'q': 0, **SCATTER_STYLES['comsol_1'], 'label': 'Comsol, mf'},
                    'cs': {'p': 0, 'q': 1, **SCATTER_STYLES['comsol_1']},
                    'ss': {'p': 1, 'q': 1, **SCATTER_STYLES['comsol_1']},
                }
            },
            
            # --- Analytical Formulation ---
            {
                'source': 'analytical',
                'scenario_key': '1',
                'data_path': ['internal_impedance_matrix'], 
                'plot_style': 'line',
                'series': { 
                    'cc': {'p': 0, 'q': 0, 'color': 'black', 'linestyle': '-',  'linewidth': 1.0, 'label': r'$R_{cc}$: core self-resistance'},
                    'cs': {'p': 0, 'q': 1, 'color': 'black', 'linestyle': ':',  'linewidth': 1.0, 'label': r'$R_{cs}$: core-sheath mutual resistance'},
                    'ss': {'p': 1, 'q': 1, 'color': 'black', 'linestyle': '--', 'linewidth': 1.0, 'label': r'$R_{ss}$: sheath self-resistance'},
                }
            },
        ],
    },

    'internal_impedance_matrix_comparison': {
        'plot_function': '_plot_non_matricial_parameters',
        
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-3, 1E1),
            'legend': True,
        },

        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0.0, 0.2),
            'legend': False,
        },
        
        'suptitle': r'P.u.l. Internal impedance matrix of a coaxial (SC) cable [Yin, 1990]',
        
        'data_series': [
            # --- COMSOL Formulation ---
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['internal_impedance_matrix'], 
                'plot_style': 'scatter',
                'series': { 
                    'self_core_js':       {**SCATTER_STYLES['comsol_1'], 'label': r'$J_s$ Method'},
                    'self_sheath_js':     {**SCATTER_STYLES['comsol_1']},
                    'mutual_js':          {**SCATTER_STYLES['comsol_1']},
                    'self_core_energy':   {**SCATTER_STYLES['mom_so_1'], 'label': r'Energy Method'},
                    'self_sheath_energy': {**SCATTER_STYLES['mom_so_1']},
                    'mutual_energy':      {**SCATTER_STYLES['mom_so_1']},
                }
            },            
        ],
    },

    'series_impedance_matrix': {
        'plot_function': '_plot_matricial_parameters',
        
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance, $R_{11}$ $(\Omega/km)$',
            'y_lim': (1E-2, 1E1),
            'legend': True,
        },

        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance, $L_{11}$ $(mH/km)$',
            'y_lim': (0.11, 0.19),
            'legend': False,
        },
        
        'suptitle': r'P.u.l. series impedance of a single core cable with sheath path return [Patel, 2014; Ametani, 2015]',
        
        'data_series': [
            # --- MoM-SO Formulation ---
            {
                'source': 'mom_so',
                'scenario_key': '1',
                'data_path': ['series_impedance_matrix'], 
                'plot_style': 'scatter',
                'series': { 
                    'Zs': {'p': 0, 'q': 0, **SCATTER_STYLES['mom_so_1'], 'label': 'MoM-SO'},
                }
            },
            
            # --- COMSOL Formulation ---
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['series_impedance_matrix'], 
                'plot_style': 'scatter',
                'series': { 
                    'Zs': {'p': 0, 'q': 0, **SCATTER_STYLES['comsol_1'], 'label': 'Comsol, mf'},
                }
            },
            
            # --- Analytical Formulation ---
            {
                'source': 'analytical',
                'scenario_key': '1',
                'data_path': ['series_impedance_matrix'], 
                'plot_style': 'line',
                'series': { 
                    'Zs': {'p': 0, 'q': 0, **SCATTER_STYLES['analytical_1'], 'label': 'Analytical'},
                }
            },
        ],
    },
    
    'coaxial_cable_parameters': {
        'plot_function': '_plot_non_matricial_parameters',
        
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-2, 1E1),
            'legend': True,
        },

        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            'y_lim': (0.0, 0.2),
            'legend': False,
        },
        
        'suptitle': r'P.u.l. impedance terms of a coaxial cable with sheath path return [Patel, 2014; Ametani, 2015]',
        
        'data_series': [
            # --- MoM-SO Formulation ---
            {
                'source': 'mom_so',
                'scenario_key': '1',
                'data_path': ['coaxial_cable_parameters'], 
                'plot_style': 'scatter',
                'series': { 
                    'Zcs': {**SCATTER_STYLES['mom_so_1'], 'label': 'MoM-SO'},
                }
            },
            
            # --- COMSOL Formulation ---
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['coaxial_cable_parameters'], 
                'plot_style': 'scatter',
                'series': { 
                    'Zcs': {**SCATTER_STYLES['comsol_1'], 'label': 'Comsol, mf'},
                    'z11': {**SCATTER_STYLES['comsol_1']},
                    'z12': {**SCATTER_STYLES['comsol_1']},
                    'z2i': {**SCATTER_STYLES['comsol_1']},
                }
            },
            
            # --- Analytical Formulation ---
            {
                'source': 'analytical',
                'scenario_key': '1',
                'data_path': ['coaxial_cable_parameters'], 
                'plot_style': 'line',
                'series': { 
                    'Zcs': {'color': 'black',     'linestyle': '-', 'linewidth': 1.0, 'label': r'$R_{cs}=R_{11}+R_{12}+R_{2i}$'},
                    'z11': {'color': 'red',       'linestyle': ':', 'linewidth': 1.5, 'label': r'$R_{11}$: internal resistance of core outer surface'},
                    'z12': {'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.5, 'label': r'$R_{12}$: core outer insulator resistance'},
                    'z2i': {'color': 'blue',      'linestyle': ':', 'linewidth': 1.5, 'label': r'$R_{2i}$: internal resistance of sheath inner surface'},
                }
            },
        ],
    },
}