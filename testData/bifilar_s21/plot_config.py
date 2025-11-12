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
        'linewidth': 1.5,
        'zorder': 1,
    },
    'analytical_2': {
        'color': 'red',
        'linestyle': ':',
        'linewidth': 1.5,
        'zorder': 1,
    },
    'mom_so_1': {
        'marker': 'o',
        's': 8,  
        'color': 'black',
        'linestyle': '-',
        'facecolors': 'black',
        'zorder': 2,
    },
    'mom_so_2': {
        'marker': 's',
        's': 8,  
        'color': 'red',
        'linestyle': '-',
        'facecolors': 'red',
        'zorder': 2,
    },
    'comsol_1': {
        'marker': 'o',
        's': 40,
        'color': 'black',
        'linestyle': '-',
        'facecolors': 'none',
        'zorder': 3,
    },
    'comsol_2': {
        'marker': 's',
        's': 40,
        'color': 'red',
        'linestyle': '-',
        'facecolors': 'none',
        'zorder': 3,
    },
}

PLOT_CONFIG = {
    'partial_impedance_matrix': {
        'plot_function': '_plot_matricial_upper_triangular',

        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance $(\Omega/km)$',
            'y_lim': (1E-2, 1E1),
            'legend': True,
        },

        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance $(mH/km)$',
            # 'y_lim': (0.1, 0.15),
            'legend': False,
        },
        
        # --- Configurações Específicas ---
        'suptitle': r'P.u.l. Partial Impedance Matrix of a bifilar transmission line',
        'data_series': [
            # --- Série MoM-SO ---
            {
                'source': 'mom_so',
                'scenario_key': '1',
                'data_path': ['partial_impedance_matrix'], 
                'plot_style': 'scatter',
                'series': { 
                    'cc': {'p': 0, 'q': 0, **SCATTER_STYLES['mom_so_1'], 'label': r'$R_{11}=R_{22}$ (MoM-SO)'},
                    'cs': {'p': 0, 'q': 1, **SCATTER_STYLES['mom_so_2'], 'label': r'$R_{12}=R_{21}$ (MoM-SO)'},
                }
            },
            
            # --- Série COMSOL ---
            {
                'source': 'comsol',
                'scenario_key': '1',
                'data_path': ['partial_impedance_matrix'], 
                'plot_style': 'scatter',
                'series': { 
                    'cc': {'p': 0, 'q': 0, **SCATTER_STYLES['comsol_1'], 'label': r'$R_{11}=R_{22}$ (Comsol, mf)'},
                    'cs': {'p': 0, 'q': 1, **SCATTER_STYLES['comsol_2'], 'label': r'$R_{12}=R_{21}$ (Comsol, mf)'},
                }
            },
        ],
    },
    
    'series_impedance_matrix': {
        'plot_function': '_plot_matricial_upper_triangular',
        
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'Resistance, $R_{11}$ $(\Omega/km)$',
            'y_lim': (1E-2, 1E2),
            'legend': True,
        },

        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'Inductance, $L_{11}$ $(mH/km)$',
            'y_lim': (0.1, 0.45),
            'legend': True,
        },
        
        # --- Configurações Específicas ---
        'suptitle': r'P.u.l. Series Impedance Matrix of a bifilar transmission line',
        
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
            {
                'source': 'analytical',
                'scenario_key': '1',
                'data_path': ['high_frequency_limit'], 
                'plot_style': 'line',
                'series': { 
                    'Zs': {'p': 0, 'q': 0, **SCATTER_STYLES['analytical_2'], 'label': 'High-Frequency Limit'},
                }
            },
        ],
    },
}

