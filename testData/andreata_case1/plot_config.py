import numpy as np

COMSOL_TEMPLATE = [
    {
        'key': 'rho_g_100_epsr1_1_mf',
        'label': r'COMSOL ($\rho_e=100, \epsilon_r=1$)',
        'marker': 'o', 's': 35, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 10
    },
    {
        'key': 'rho_g_100_epsr1_20_mf',
        'label': r'COMSOL ($\rho_e=100, \epsilon_r=20$)',
        'marker': 'o', 's': 15, 'facecolors': 'black', 'edgecolors': 'none', 'zorder': 10
    },
    {
        'key': 'rho_g_500_epsr1_1_mf',
        'label': r'COMSOL ($\rho_e=500, \epsilon_r=1$)',
        'marker': 's', 's': 15, 'color': 'black', 'zorder': 10
    }
]

XUE_TEMPLATE = [
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

VANCE_TEMPLATE = [
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

MATLAB_TEMPLATE = [
    {
        'key': 'measured',
        'label': 'MATLAB (ref.)',
        'marker': 'x', 's': 12, 'color': 'blue', 'linewidths': 1.0, 'zorder': 11
    }
]

DECONTI_TEMPLATE = [
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
    # 'z11_internal_vs_matlab': {
    #     'suptitle': r'$Z_{11}$: internal-only analytical (Zi) vs. MATLAB reference (Z)',
    #     'p': 0, 'q': 0,
    #     'internal_style': {'label': r'Analytical, internal only ($Zi_{11}$)', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
    #     'matlab_series_to_plot': MATLAB_TEMPLATE,
    #     'matlab_matrix_key': 'internal_impedance_matrix',
    #     'left_plot': {
    #         **PLOT_TPL_RESISTANCE,
    #         'label': r'$Rs_{11} \, (\Omega/km)$',
    #     },
    #     'right_plot': {
    #         **PLOT_TPL_INDUCTANCE,
    #         'label': r'$Ls_{11} \, (mH/km)$',
    #     }
    # },

    # 'z12_internal_vs_matlab': {
    #     'suptitle': r'$Z_{12}$: internal-only analytical (Zi) vs. MATLAB reference (Z)',
    #     'p': 0, 'q': 1,
    #     'internal_style': {'label': r'Analytical, internal only ($Zi_{12}$)', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
    #     'matlab_series_to_plot': MATLAB_TEMPLATE,
    #     'matlab_matrix_key': 'internal_impedance_matrix',
    #     'left_plot': {
    #         **PLOT_TPL_RESISTANCE,
    #         'label': r'$Rs_{12} \, (\Omega/km)$',
    #     },
    #     'right_plot': {
    #         **PLOT_TPL_INDUCTANCE,
    #         'label': r'$Ls_{12} \, (mH/km)$',
    #     }
    # },

    # 'z22_internal_vs_matlab': {
    #     'suptitle': r'$Z_{22}$: internal-only analytical (Zi) vs. MATLAB reference (Z)',
    #     'p': 1, 'q': 1,
    #     'internal_style': {'label': r'Analytical, internal only ($Zi_{22}$)', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
    #     'matlab_series_to_plot': MATLAB_TEMPLATE,
    #     'matlab_matrix_key': 'internal_impedance_matrix',
    #     'left_plot': {
    #         **PLOT_TPL_RESISTANCE,
    #         'label': r'$Rs_{22} \, (\Omega/km)$',
    #     },
    #     'right_plot': {
    #         **PLOT_TPL_INDUCTANCE,
    #         'label': r'$Ls_{22} \, (mH/km)$',
    #     }
    # },

    # 'z22_matlab_only': {
    #     'suptitle': 'MATLAB reference only — Z_22 (phase-a sheath self-impedance)',
    #     'p': 1, 'q': 1,
    #     'series_to_plot': [],
    #     'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
    #     'matlab_series_to_plot': MATLAB_TEMPLATE,
    #     'matlab_matrix_key': 'series_impedance_matrix',
    #     'left_plot': {
    #         **PLOT_TPL_RESISTANCE,
    #         'label': r'$Rs_{22} \, (\Omega/km)$',
    #     },
    #     'right_plot': {
    #         **PLOT_TPL_INDUCTANCE,
    #         'label': r'$Ls_{22} \, (mH/km)$',
    #     }
    # },

    # 'z22_full_vs_matlab': {
    #     'suptitle': r'$Z_{22}$: full analytical (Zi+Zg, quasi-TEM) vs. MATLAB reference (Z)',
    #     'p': 1, 'q': 1,
    #     'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
    #     'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
    #     'matlab_series_to_plot': MATLAB_TEMPLATE,
    #     'matlab_matrix_key': 'series_impedance_matrix',
    #     'left_plot': {
    #         **PLOT_TPL_RESISTANCE,
    #         'label': r'$Rs_{22} \, (\Omega/km)$',
    #     },
    #     'right_plot': {
    #         **PLOT_TPL_INDUCTANCE,
    #         'label': r'$Ls_{22} \, (mH/km)$',
    #     }
    # },

    'internal_impedance_core_sheath': {
        'suptitle': 'Figure 4.19 (internal-only): $Z_i$ core/sheath analytical vs. MATLAB reference',
        'xlim': (1E-2, 1E7),
        'components': [
            {
                'p': 0, 'q': 0,
                'internal_style': {'label': r'Analytical, internal only ($Zi_{11}$)', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB ($Z_{11}$)', 'marker': 'x', 's': 12, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                'internal_style': {'label': r'Analytical, internal only ($Zi_{12}$)', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB ($Z_{12}$)', 'marker': '+', 's': 20, 'color': 'tab:red', 'linewidths': 1.2, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                'internal_style': {'label': r'Analytical, internal only ($Zi_{22}$)', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB ($Z_{22}$)', 'marker': 'o', 's': 15, 'facecolors': 'none', 'edgecolors': 'tab:blue', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
        ],
        'matlab_matrix_key': 'internal_impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rs \, (\Omega/km)$',
            'y_lim': (1E-4, 1E2),
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Ls \, (mH/km)$',
        }
    },

    'internal_admittance_core_sheath': {
        'suptitle': 'Figure 4.23 (internal-only): $Y_i$ core/sheath analytical vs. MATLAB reference',
        'internal_matrix_key': 'shunt_admittance_matrix',
        'xlim': (1E-2, 1E7),
        'components': [
            {
                'p': 0, 'q': 0,
                'internal_style': {'label': r'Analytical, internal only ($Yi_{11}$)', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB ($Y_{11}$)', 'marker': 'x', 's': 12, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                'internal_style': {'label': r'Analytical, internal only ($Yi_{12}$)', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB ($Y_{12}$)', 'marker': '+', 's': 20, 'color': 'tab:red', 'linewidths': 1.2, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                'internal_style': {'label': r'Analytical, internal only ($Yi_{22}$)', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB ($Y_{22}$)', 'marker': 'o', 's': 15, 'facecolors': 'none', 'edgecolors': 'tab:blue', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
        ],
        'matlab_matrix_key': 'internal_admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C \, (nF/km)$',
        }
    },

    'series_impedance_all_scenarios': {
        'suptitle': 'Figure 4.19: P.u.l. Self-impedance of phase - a sheath [Xue, 2018]',
        'p': 1, 'q': 1,
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'series_impedance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'series_impedance_matrix',
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rs_{22} \, (\Omega/km)$',
            # 'y_lim': (1E1, 1E5),
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Ls_{22} \, (mH/km)$',
            # 'y_lim': (0.5, 2.0),
        }
    },
    
    'fig421a': {
        'suptitle': 'Figure 4.21a: P.u.l. Mutual-impedance between phase - a and phase - b sheaths [Xue, 2018]',
        'p': 1, 'q': 3,
        'xlim': (1E4, 1E7),
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
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
        'xlim': (1E4, 1E7),
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
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
    
    'shunt_admittance_all_scenarios': {
        'suptitle': 'Figure 4.23: P.u.l. Self-admittance of phase - a sheath [Xue, 2018]',
        'p': 1, 'q': 1,
        'xlim': (1E4, 1E7),
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'shunt_admittance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'shunt_admittance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'shunt_admittance_matrix',
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
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'shunt_admittance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
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
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'shunt_admittance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
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
        'series_to_plot': XUE_TEMPLATE,
        'path': ['earth_return_parameters', 'gamma_earth'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
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
        'suptitle': 'P.u.l. Self earth-return impedance of phase - a [Xue, 2018] vs. MATLAB reference',
        'p': 0, 'q': 0,
        'series_to_plot': XUE_TEMPLATE,
        'path': ['earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'impedance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_impedance_matrix',
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
        'suptitle': 'P.u.l. Mutual earth-return impedance between phase - a and phase - b [Xue, 2018] vs. MATLAB reference',
        'p': 0, 'q': 1,
        'series_to_plot': XUE_TEMPLATE,
        'path': ['earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'impedance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_impedance_matrix',
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
        'suptitle': 'P.u.l. Mutual earth-return impedance between phase - a and phase - c [Xue, 2018] vs. MATLAB reference',
        'p': 0, 'q': 2,
        'series_to_plot': XUE_TEMPLATE,
        'path': ['earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'impedance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_impedance_matrix',
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
        'series_to_plot': XUE_TEMPLATE,
        'path': ['earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
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
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE,
        'path': ['earth_return_parameters', 'admittance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
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
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE,
        'path': ['earth_return_parameters', 'admittance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
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
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE,
        'path': ['earth_return_parameters', 'admittance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$Gg_{13} \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$Cg_{13} \, (nF/km)$',
        }
    },

    'earth_return_potential_coefficient_self': {
        'suptitle': 'P.u.l. Self earth-return potential coefficient of phase - a vs. MATLAB reference',
        'p': 0, 'q': 0,
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE,
        'path': ['earth_return_parameters', 'potential_coefficient'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_potential_coefficient_matrix',
        'left_plot': {
            **PLOT_TPL_REAL,
            'label': r'$Re\{Pg_{11}\} \, (m/F)$',
        },
        'right_plot': {
            **PLOT_TPL_IMAG,
            'label': r'$Im\{Pg_{11}\} \, (m/F)$',
        }
    },

    'earth_return_potential_coefficient_mutual_ab': {
        'suptitle': 'P.u.l. Mutual earth-return potential coefficient between phase - a and phase - b vs. MATLAB reference',
        'p': 0, 'q': 1,
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE,
        'path': ['earth_return_parameters', 'potential_coefficient'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_potential_coefficient_matrix',
        'left_plot': {
            **PLOT_TPL_REAL,
            'label': r'$Re\{Pg_{12}\} \, (m/F)$',
        },
        'right_plot': {
            **PLOT_TPL_IMAG,
            'label': r'$Im\{Pg_{12}\} \, (m/F)$',
        }
    },

    'earth_return_potential_coefficient_mutual_ac': {
        'suptitle': 'P.u.l. Mutual earth-return potential coefficient between phase - a and phase - c vs. MATLAB reference',
        'p': 0, 'q': 2,
        'series_to_plot': XUE_TEMPLATE + VANCE_TEMPLATE,
        'path': ['earth_return_parameters', 'potential_coefficient'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_potential_coefficient_matrix',
        'left_plot': {
            **PLOT_TPL_REAL,
            'label': r'$Re\{Pg_{13}\} \, (m/F)$',
        },
        'right_plot': {
            **PLOT_TPL_IMAG,
            'label': r'$Im\{Pg_{13}\} \, (m/F)$',
        }
    },
}