from ..andreata_common.plot_templates import (
    COMSOL_TEMPLATE, XUE_TEMPLATE, MATLAB_TEMPLATE, DECONTI_TEMPLATE,
    PLOT_TPL_REAL, PLOT_TPL_IMAG, PLOT_TPL_RESISTANCE, PLOT_TPL_INDUCTANCE,
    PLOT_TPL_CONDUCTANCE, PLOT_TPL_CAPACITANCE,
)

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
        'p': 6, 'q': 6,
        'xlim': (1E4, 1E7),
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G_{77} \, (S/km)$',
            # 'y_lim': (0.0, 25),
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C_{77} \, (nF/km)$',
            # 'y_lim': (0.0, 2.0),
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