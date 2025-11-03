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

XUE_SERIES = [
    {
        'key': 'p100_er1',
        'type': {
            'main': {'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
        }
    },
    {
        'key': 'p100_er20',
        'type': {
            'main': {'label': r'$\rho_e=100 \;\Omega m, \epsilon_r=20$', 'color': 'black', 'linestyle': '--', 'linewidth': 1.5},
        }
    },
    {
        'key': 'p500_er1',
        'type': {
            'main': {'label': r'$\rho_e=500 \;\Omega m, \epsilon_r=1$', 'color': 'black', 'linestyle': '-.', 'linewidth': 1.5},
        }
    },
] 

VANCE_SERIES = [
    {
        'key': 'p100_er1_vance',
        'type': {
            'main': {'label': 'Vance Form. (1978)', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0},
        }
    },
    {
        'key': 'p100_er20_vance',
        'type': {
            'main': {'label': '', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0},
        }
    },
    {
        'key': 'p500_er1_vance',
        'type': {
            'main': {'label': '', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.0},
        }
    },
] 

DECONTI_SERIES = [
    {
        'key': 'p100_er1_deconti',
        'type': {
            'main': {'label': 'De Conti et al. Form. (2023)', 'color': 'red', 'linestyle': ':', 'linewidth': 1.0},
        }
    },
    {
        'key': 'p100_er20_deconti',
        'type': {
            'main': {'label': '', 'color': 'red', 'linestyle': ':', 'linewidth': 1.0},
        }
    },
    {
        'key': 'p500_er1_deconti',
        'type': {
            'main': {'label': '', 'color': 'red', 'linestyle': ':', 'linewidth': 1.0},
        }
    },
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
    'fig419': {
        'suptitle': 'Figure 4.19: P.u.l. Self-impedance of phase - a sheath [Xue, 2018]',
        'p': 1, 'q': 1,
        'series_to_plot': XUE_SERIES + DECONTI_SERIES,
        'data_path': ['scenarios', '{key}', 'quasi_tem_matrices', 'series_impedance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'series_impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rs_{22} \, (\Omega/km)$',
            'y_lim': (1E1, 1E5),
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Ls_{22} \, (mH/km)$',
            'y_lim': (0.5, 2.0),
        }
    },
    
    'fig421a': {
        'suptitle': 'Figure 4.21a: P.u.l. Mutual-impedance between phase - a and phase - b sheaths [Xue, 2018]',
        'p': 1, 'q': 3,
        'series_to_plot': XUE_SERIES + DECONTI_SERIES,
        'data_path': ['scenarios', '{key}', 'quasi_tem_matrices', 'series_impedance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'series_impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rs_{24} \, (\Omega/km)$',
            'y_lim': (1E1, 1E5),
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Ls_{24} \, (mH/km)$',
            'y_lim': (0.0, 1.5),
        }
    },
    
    'fig421b': {
        'suptitle': 'Figure 4.21b: P.u.l. Mutual-impedance between phase - a and phase - c sheaths [Xue, 2018]',
        'p': 1, 'q': 5,
        'series_to_plot': XUE_SERIES + DECONTI_SERIES,
        'data_path': ['scenarios', '{key}', 'quasi_tem_matrices', 'series_impedance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'series_impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rs_{26} \, (\Omega/km)$',
            'y_lim': (1E1, 1E5),
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Ls_{26} \, (mH/km)$',
            'y_lim': (0.0, 1.5),
        }
    },
    
    'fig423': {
        'suptitle': 'Figure 4.23: P.u.l. Self-admittance of phase - a sheath [Xue, 2018]',
        'p': 1, 'q': 1,
        'series_to_plot': XUE_SERIES + VANCE_SERIES + DECONTI_SERIES,
        'data_path': ['scenarios', '{key}', 'quasi_tem_matrices', 'shunt_admittance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'shunt_admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G_{22} \, (S/km)$',
            'y_lim': (0.0, 20),
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C_{22} \, (nF/km)$',
            'y_lim': (0.0, 3.0),
        }
    },

    'fig425a': {
        'suptitle': 'Figure 4.25: P.u.l. Mutual-admittance between phase - a and phase - b sheaths [Xue, 2018]',
        'p': 1, 'q': 3,
        'series_to_plot': XUE_SERIES + VANCE_SERIES + DECONTI_SERIES,
        'data_path': ['scenarios', '{key}', 'quasi_tem_matrices', 'shunt_admittance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'shunt_admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G_{24} \, (S/km)$',
            'y_lim': (-8.0, 2.0),
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C_{24} \, (nF/km)$',
            'y_lim': (-0.8, 0.2),
        }
    },

    'fig425b': {
        'suptitle': 'Figure 4.25a: P.u.l. Mutual-admittance between phase - a and phase - c sheaths [Xue, 2018]',
        'p': 1, 'q': 5,
        'series_to_plot': XUE_SERIES + VANCE_SERIES + DECONTI_SERIES,
        'data_path': ['scenarios', '{key}', 'quasi_tem_matrices', 'shunt_admittance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'shunt_admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G_{26} \, (S/km)$',
            'y_lim': (-4.0, 2.0),
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C_{26} \, (nF/km)$',
            'y_lim': (-0.6, 0.2),
        }
    },

    'earth_propagation_constant': {
        'suptitle': r'Earth propagation constant, $\gamma_1$',
        'series_to_plot': XUE_SERIES,
        'data_path': ['scenarios', '{key}', 'earth_return_parameters', 'gamma_earth'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'gamma_earth',
        'left_plot': {
            **PLOT_TPL_REAL,
            'label': r'Attenuation Const., $\alpha \, (Np/m)$',
        },
        'right_plot': {
            **PLOT_TPL_IMAG,
            'label': r'Phase Const., $\beta \, (rad/m)$',
        }
    },

    'earth_return_impedance_self': {
        'suptitle': 'P.u.l. Self earth-return impedance of phase - a [Xue, 2018]',
        'p': 0, 'q': 0,
        'series_to_plot': XUE_SERIES,
        'data_path': ['scenarios', '{key}', 'earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rg_{11} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lg_{11} \, (mH/km)$',
        }
    },

    'earth_return_impedance_mutual_ab': {
        'suptitle': 'P.u.l. Mutual earth-return impedance between phase - a and phase - b [Xue, 2018]',
        'p': 0, 'q': 1,
        'series_to_plot': XUE_SERIES,
        'data_path': ['scenarios', '{key}', 'earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rg_{12} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lg_{12} \, (mH/km)$',
        }
    },

    'earth_return_impedance_mutual_ac': {
        'suptitle': 'P.u.l. Mutual earth-return impedance between phase - a and phase - c [Xue, 2018]',
        'p': 0, 'q': 2,
        'series_to_plot': XUE_SERIES,
        'data_path': ['scenarios', '{key}', 'earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rg_{13} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lg_{13} \, (mH/km)$',
        }
    },

    'earth_return_impedance_mutual_bc': {
        'suptitle': 'P.u.l. Mutual earth-return impedance between phase - b and phase - c [Xue, 2018]',
        'p': 1, 'q': 2,
        'series_to_plot': XUE_SERIES,
        'data_path': ['scenarios', '{key}', 'earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rg_{23} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lg_{23} \, (mH/km)$',
        }
    },
    
    'earth_return_admittance_self': {
        'suptitle': 'P.u.l. Self earth-return admittance of phase - a [Xue, 2018]',
        'p': 0, 'q': 0,
        'series_to_plot': XUE_SERIES + VANCE_SERIES,
        'data_path': ['scenarios', '{key}', 'earth_return_parameters', 'admittance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$Gg_{11} \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$Cg_{11} \, (nF/km)$',
        }
    },

    'earth_return_admittance_mutual_ab': {
        'suptitle': 'P.u.l. Mutual earth-return admittance between phase - a and phase - b [Xue, 2018]',
        'p': 0, 'q': 1,
        'series_to_plot': XUE_SERIES + VANCE_SERIES,
        'data_path': ['scenarios', '{key}', 'earth_return_parameters', 'admittance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$Gg_{12} \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$Cg_{12} \, (nF/km)$',
        }
    },

    'earth_return_admittance_mutual_ac': {
        'suptitle': 'P.u.l. Mutual earth-return admittance between phase - a and phase - c [Xue, 2018]',
        'p': 0, 'q': 2,
        'series_to_plot': XUE_SERIES + VANCE_SERIES,
        'data_path': ['scenarios', '{key}', 'earth_return_parameters', 'admittance_matrix'],
        'comsol_series_to_plot': comsol_layout_template,
        'comsol_matrix_key': 'admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$Gg_{13} \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$Cg_{13} \, (nF/km)$',
        }
    }
}