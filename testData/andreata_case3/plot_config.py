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
        'label': 'MATLAB',
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
    'internal_impedance_matrix': {
        'suptitle': 'P.u.l. Internal Impedance Matrix [Ametani, 2015]',
        'components': [
            {
                'p': 0, 'q': 0,
                'internal_style': {'label': r'$Zi_{11}$', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': 'x', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                'internal_style': {'label': r'$Zi_{12}$', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': 'x', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                'internal_style': {'label': r'$Zi_{22}$', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': 'x', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 6, 'q': 6,
                'internal_style': {'label': r'$Zi_{77}$', 'color': 'tab:green', 'linestyle': ':', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB', 'marker': 'x', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
        ],
        'matlab_matrix_key': 'internal_impedance_matrix',
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rs \, (\Omega/km)$',
            'y_lim': (1E-2, 1E2),
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Ls \, (mH/km)$',
        }
    },

    'internal_admittance_matrix': {
        'suptitle': 'P.u.l. Internal Admittance Matrix [Ametani, 2015]',
        'internal_matrix_key': 'shunt_admittance_matrix',
        'components': [
            {
                'p': 0, 'q': 0,
                'internal_style': {'label': r'$Yi_{11}$', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': 'x', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                'internal_style': {'label': r'$Yi_{12}$', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': 'x', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                'internal_style': {'label': r'$Yi_{22}$', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': 'x', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 6, 'q': 6,
                'internal_style': {'label': r'$Yi_{77}$', 'color': 'tab:green', 'linestyle': ':', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB', 'marker': 'x', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
        ],
        'matlab_matrix_key': 'internal_admittance_matrix',
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C \, (nF/km)$',
        }
    },

    'mutual_impedance_phase_a_sheath_ecc': {
        'suptitle': 'P.u.l. Mutual-impedance of phase-a sheath and ECC [Xue, 2018]',
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'series_impedance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'series_impedance_matrix',
        'p': 1, 'q': 6,
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rs_{27} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Ls_{27} \, (mH/km)$',
        }
    },

    'self_impedance_ecc': {
        'suptitle': 'P.u.l. Self-impedance of ECC [Xue, 2018]',
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'series_impedance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'series_impedance_matrix',
        'p': 6, 'q': 6,
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rs_{77} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Ls_{77} \, (mH/km)$',
        }
    },

    'self_admittance_ecc': {
        'suptitle': 'P.u.l. Self-admittance of ECC [Xue, 2018]',
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'shunt_admittance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'shunt_admittance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'shunt_admittance_matrix',
        'p': 3, 'q': 3,
        'xlim': (1E4, 1E7),
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G_{77} \, (S/km)$',
            # 'y_lim': (0.0, 20),
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C_{77} \, (nF/km)$',
            # 'y_lim': (0.0, 3.0),
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

    'earth_return_impedance_ecc': {
        'suptitle': 'P.u.l. Self earth-return impedance of ECC [Xue, 2018]',
        'p': 6, 'q': 6,
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'earth_return_impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'impedance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_impedance_matrix',
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rg_{77} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lg_{77} \, (mH/km)$',
        }
    },

    'earth_return_admittance_ecc': {
        'suptitle': 'P.u.l. Self earth-return admittance of ECC [Xue, 2018]',
        'p': 3, 'q': 3,
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['earth_return_parameters', 'admittance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'admittance_matrix',
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$Gg_{77} \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$Cg_{77} \, (nF/km)$',
        }
    },

    'earth_return_potential_coeff_ecc': {
        'suptitle': 'P.u.l. Self earth-return potential coefficient of ECC [Xue, 2018]',
        'p': 6, 'q': 6,
        'series_to_plot': XUE_TEMPLATE + DECONTI_TEMPLATE,
        'path': ['quasi_tem_matrices', 'earth_return_potential_coefficient'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_potential_coefficient_matrix',
        'left_plot': {
            **PLOT_TPL_REAL,
            'label': r'$Re\{Pg_{77}\} \, (m/F)$',
        },
        'right_plot': {
            **PLOT_TPL_IMAG,
            'label': r'$Im\{Pg_{77}\} \, (m/F)$',
        }
    },
}