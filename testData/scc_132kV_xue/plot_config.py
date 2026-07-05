# PLOT_CONFIG for plotter.scc_models.SingleCoreCableModels — Xue-style coaxial SCC.
#
# Expects pul_data (see scc_132kV_xue.py) to additionally provide:
#   pul_data['analytical'] = {'frequencies': pul_data['frequencies'],
#                             'scenarios': pul_data['scenarios']}
# and a synthetic 'internal' scenario for the internal_* graphs:
#   pul_data['analytical']['scenarios']['internal'] = {
#       'impedance_matrix': pul_data['internal_matrices']['impedance_matrix'],
#       'shunt_admittance_matrix': pul_data['internal_matrices']['shunt_admittance_matrix'],
#       'potential_coefficient_matrix': Pi_3d,  # 3D-broadcast, see plan
#   }

R_LABEL = r'$R \, (\Omega/km)$'
L_LABEL = r'$L \, (mH/km)$'
G_LABEL = r'$G \, (S/km)$'
C_LABEL = r'$C \, (\mu F/km)$'
P_NORM_LABEL = r'$|P| \times 10^9 \, (\Omega m s^{-1})$'
P_ANGLE_LABEL = 'Angle of P (Degrees)'

# Single analytical scenario used by all 'composition' graphs (matches the
# original scc_params: only 'p100_xue' was ever populated for these plots).
COMPOSITION_SCENARIO_KEY = 'p100_xue'
COMPOSITION_STYLES = {
    'series_term': {'label': 'Series Impedance', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0, 'zorder': 1},
    'shunt_term': {'label': 'Shunt Admittance', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0, 'zorder': 1},
    'potential_term': {'label': 'Equivalent Potential Coeff.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0, 'zorder': 1},
    'internal_term': {'label': 'Internal', 'color': 'blue', 'linestyle': '--', 'linewidth': 1.5, 'zorder': 2},
    'earth_return_term': {'label': 'Earth-return (Integral Form.)', 'color': 'darkgreen', 'linestyle': '--', 'linewidth': 1.5, 'zorder': 3},
    'series_composition': {'label': 'Series Composition', 'color': 'red', 'linestyle': ':', 'linewidth': 1.5, 'zorder': 4},
}

# --- Shared series-of-scenarios building blocks (styling only) ---

_IMPEDANCE_SERIES = [
    {
        'source': 'analytical', 'scenario_key': 'p100_xue',
        'data_path': ['quasi_tem_matrices', 'series_impedance_matrix'],
        'series': {
            'self_core': {'p': 0, 'q': 0, 'label': 'Self-Core (Integral Form.)', 'color': 'red', 'linestyle': '-', 'linewidth': 2.0},
            'self_sheath': {'p': 1, 'q': 1, 'label': 'Self-Sheath (Integral Form.)', 'color': 'blue', 'linestyle': '-', 'linewidth': 2.0},
            'mutual': {'p': 0, 'q': 1, 'label': 'Mutual Core-Sheath (Integral Form.)', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 2.0},
        },
    },
    {
        'source': 'analytical', 'scenario_key': 'p100_deconti',
        'data_path': ['quasi_tem_matrices', 'series_impedance_matrix'],
        'series': {
            'self_core': {'p': 0, 'q': 0, 'label': 'De Conti Approx. Form.', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
            'self_sheath': {'p': 1, 'q': 1, 'label': '', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
            'mutual': {'p': 0, 'q': 1, 'label': '', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
        },
    },
]

_ADMITTANCE_SERIES = [
    {
        'source': 'analytical', 'scenario_key': 'p100_xue',
        'data_path': ['quasi_tem_matrices', 'shunt_admittance_matrix'],
        'series': {
            'self_core': {'p': 0, 'q': 0, 'label': 'Self-Core (Integral Form.)', 'color': 'red', 'linestyle': '-', 'linewidth': 2.0},
            'self_sheath': {'p': 1, 'q': 1, 'label': 'Self-Sheath (Integral Form.)', 'color': 'blue', 'linestyle': '-', 'linewidth': 2.0},
            'mutual': {'p': 0, 'q': 1, 'label': 'Mutual Core-Sheath (Integral Form.)', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 2.0},
        },
    },
    {
        'source': 'analytical', 'scenario_key': 'p100_vance',
        'data_path': ['quasi_tem_matrices', 'shunt_admittance_matrix'],
        'series': {
            'self_core': {'p': 0, 'q': 0, 'label': 'Vance Approx. Form.', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
            'self_sheath': {'p': 1, 'q': 1, 'label': '', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
            'mutual': {'p': 0, 'q': 1, 'label': '', 'color': 'black', 'linestyle': '--', 'linewidth': 1.0},
        },
    },
]

_EARTH_RETURN_IMPEDANCE_SERIES = [
    {
        'source': 'analytical', 'scenario_key': 'p100_xue',
        'data_path': ['earth_return_parameters', 'impedance_matrix'],
        'series': {'earth_return': {'p': 0, 'q': 0, 'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0}},
    },
    {
        'source': 'analytical', 'scenario_key': 'p100_deconti',
        'data_path': ['earth_return_parameters', 'impedance_matrix'],
        'series': {'earth_return': {'p': 0, 'q': 0, 'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0}},
    },
]

_EARTH_RETURN_ADMITTANCE_SERIES = [
    {
        'source': 'analytical', 'scenario_key': 'p100_xue',
        'data_path': ['earth_return_parameters', 'admittance_matrix'],
        'series': {'shunt': {'p': 0, 'q': 0, 'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0}},
    },
    {
        'source': 'analytical', 'scenario_key': 'p100_deconti',
        'data_path': ['earth_return_parameters', 'admittance_matrix'],
        'series': {'shunt': {'p': 0, 'q': 0, 'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0}},
    },
    {
        'source': 'analytical', 'scenario_key': 'p100_vance',
        'data_path': ['earth_return_parameters', 'admittance_matrix'],
        'series': {'shunt': {'p': 0, 'q': 0, 'label': 'Vance Approx. Form.', 'color': 'blue', 'linestyle': '-.', 'linewidth': 1.0}},
    },
]

_EARTH_RETURN_POTENTIAL_SERIES = [
    {
        'source': 'analytical', 'scenario_key': 'p100_xue',
        'data_path': ['earth_return_parameters', 'potential_coefficient'],
        'series': {'potential': {'p': 0, 'q': 0, 'label': 'Magalhaes/Xue Integral Form.', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0}},
    },
    {
        'source': 'analytical', 'scenario_key': 'p100_deconti',
        'data_path': ['earth_return_parameters', 'potential_coefficient'],
        'series': {'potential': {'p': 0, 'q': 0, 'label': 'De Conti et al. Approx. Form.', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0}},
    },
]

_INTERNAL_ELEMENTS = {
    'self_core': {'p': 0, 'q': 0, 'color': 'black', 'linestyle': '-', 'label': 'Self-Core'},
    'mutual': {'p': 0, 'q': 1, 'color': 'red', 'linestyle': '--', 'label': 'Mutual Core-Sheath'},
    'self_sheath': {'p': 1, 'q': 1, 'color': 'blue', 'linestyle': ':', 'label': 'Self-Sheath'},
}

PLOT_CONFIG = {
    'series_impedance_matrix': {
        'suptitle': r'P.u.l. series impedance matrix of single core coaxial cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'yscale': 'log', 'label': R_LABEL},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e6, 'xscale': 'log', 'label': L_LABEL},
        'data_series': _IMPEDANCE_SERIES,
    },
    'shunt_admittance_matrix': {
        'suptitle': r'P.u.l. shunt admittance matrix of single core coaxial cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'label': G_LABEL},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e9, 'xscale': 'log', 'label': C_LABEL},
        'data_series': _ADMITTANCE_SERIES,
    },
    'earth_return_impedance': {
        'suptitle': r'P.u.l. earth-return impedance matrix of the single-core cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'yscale': 'log', 'label': R_LABEL, 'title': 'P.u.l. series resistance'},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e6, 'xscale': 'log', 'label': L_LABEL, 'title': 'P.u.l. series inductance'},
        'data_series': _EARTH_RETURN_IMPEDANCE_SERIES,
    },
    'earth_return_admittance': {
        'suptitle': r'P.u.l. earth-return admittance matrix of the single-core cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'label': G_LABEL, 'title': 'P.u.l. shunt conductance'},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e9, 'xscale': 'log', 'label': C_LABEL, 'title': 'P.u.l. shunt capacitance'},
        'data_series': _EARTH_RETURN_ADMITTANCE_SERIES,
    },
    'earth_return_potential': {
        'suptitle': r'P.u.l. earth-return potential coefficient of the single-core cable ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'left_plot': {'component': 'abs', 'scale': 1e-9, 'xscale': 'log', 'label': P_NORM_LABEL, 'title': 'Absolute Value'},
        'right_plot': {'component': 'angle_deg', 'xscale': 'log', 'label': P_ANGLE_LABEL, 'title': 'Angle of Potential Coefficient'},
        'data_series': _EARTH_RETURN_POTENTIAL_SERIES,
    },
    'internal_impedance': {
        'suptitle': r'P.u.l. internal impedance matrix of the single-core cable',
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'yscale': 'log', 'label': R_LABEL, 'title': 'P.u.l. series resistance'},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e6, 'xscale': 'log', 'label': L_LABEL, 'title': 'P.u.l. series inductance'},
        'data_series': [
            {'source': 'analytical', 'scenario_key': 'internal', 'data_path': ['impedance_matrix'], 'series': _INTERNAL_ELEMENTS},
            {'source': 'comsol', 'scenario_key': 'measured', 'data_path': ['impedance_matrix'], 'plot_style': 'scatter',
             'series': {
                'self_core': {'p': 0, 'q': 0, 'color': 'black', 'marker': 'o', 's': 8, 'facecolors': 'none', 'zorder': 10, 'label': 'COMSOL'},
                'mutual': {'p': 0, 'q': 1, 'color': 'red', 'marker': 'o', 's': 8, 'facecolors': 'none', 'zorder': 10},
                'self_sheath': {'p': 1, 'q': 1, 'color': 'blue', 'marker': 'o', 's': 8, 'facecolors': 'none', 'zorder': 10},
             }},
        ],
    },
    'internal_admittance': {
        'suptitle': 'P.u.l. internal admittance matrix of the single-core cable',
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'label': G_LABEL, 'title': 'P.u.l. shunt conductance'},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e9, 'xscale': 'log', 'label': C_LABEL, 'title': 'P.u.l. shunt capacitance'},
        'data_series': [
            {'source': 'analytical', 'scenario_key': 'internal', 'data_path': ['shunt_admittance_matrix'], 'series': _INTERNAL_ELEMENTS},
        ],
    },
    'internal_potential': {
        'suptitle': 'P.u.l. internal potential coefficient of the single-core cable',
        'left_plot': {'component': 'abs', 'scale': 1e-9, 'xscale': 'log', 'label': P_NORM_LABEL, 'title': 'Absolute Value'},
        'right_plot': {'component': 'angle_deg', 'xscale': 'log', 'label': P_ANGLE_LABEL, 'title': 'Angle of Potential Coefficient'},
        'data_series': [
            {'source': 'analytical', 'scenario_key': 'internal', 'data_path': ['potential_coefficient_matrix'], 'series': {
                'p00': {'p': 0, 'q': 0, 'color': 'black', 'linestyle': '-', 'label': r'$P_{00}$'},
                'p01': {'p': 0, 'q': 1, 'color': 'red', 'linestyle': '--', 'label': r'$P_{01}$'},
                'p11': {'p': 1, 'q': 1, 'color': 'blue', 'linestyle': ':', 'label': r'$P_{11}$'},
            }},
        ],
    },

    # --- Composition graphs (rendered by SingleCoreCableModels._plot_composition,
    #     not by BasePlotter._plot_matricial_parameters — see plotter/scc_models.py) ---
    'series_impedance_composition_core': {
        'suptitle': r'P.u.l. core self-impedance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 0, 'q': 0, 'kind': 'series_impedance', 'comsol': True,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'yscale': 'log', 'label': R_LABEL, 'title': 'P.u.l. series resistance'},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e6, 'xscale': 'log', 'label': L_LABEL, 'title': 'P.u.l. series inductance'},
    },
    'series_impedance_composition_sheath': {
        'suptitle': r'P.u.l. sheath self-impedance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 1, 'q': 1, 'kind': 'series_impedance', 'comsol': True,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'yscale': 'log', 'label': R_LABEL, 'title': 'P.u.l. series resistance'},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e6, 'xscale': 'log', 'label': L_LABEL, 'title': 'P.u.l. series inductance'},
    },
    'series_impedance_composition_core_sheath': {
        'suptitle': r'P.u.l. core-sheath mutual impedance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 0, 'q': 1, 'kind': 'series_impedance', 'comsol': True,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'yscale': 'log', 'label': R_LABEL, 'title': 'P.u.l. series resistance'},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e6, 'xscale': 'log', 'label': L_LABEL, 'title': 'P.u.l. series inductance'},
    },
    'shunt_admittance_composition_core': {
        'suptitle': r'P.u.l. core self-admittance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 0, 'q': 0, 'kind': 'shunt_admittance', 'comsol': False,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'label': G_LABEL, 'title': 'P.u.l. shunt conductance', 'y_lim': (0, 20)},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e9, 'xscale': 'log', 'label': C_LABEL, 'title': 'P.u.l. shunt capacitance', 'y_lim': (0, 3)},
    },
    'shunt_admittance_composition_sheath': {
        'suptitle': r'P.u.l. sheath self-admittance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 1, 'q': 1, 'kind': 'shunt_admittance', 'comsol': False,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'label': G_LABEL, 'title': 'P.u.l. shunt conductance', 'y_lim': (0, 20)},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e9, 'xscale': 'log', 'label': C_LABEL, 'title': 'P.u.l. shunt capacitance', 'y_lim': (0, 3)},
    },
    'shunt_admittance_composition_core_sheath': {
        'suptitle': r'P.u.l. core-sheath mutual admittance ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 0, 'q': 1, 'kind': 'shunt_admittance', 'comsol': False,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'real', 'scale': 1e3, 'xscale': 'log', 'label': G_LABEL, 'title': 'P.u.l. shunt conductance', 'y_lim': (0, 20)},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e9, 'xscale': 'log', 'label': C_LABEL, 'title': 'P.u.l. shunt capacitance', 'y_lim': (0, 3)},
    },
    'potential_coefficients_composition_core': {
        'suptitle': r'P.u.l. core potential coefficient ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 0, 'q': 0, 'kind': 'potential_coefficient', 'comsol': False,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'abs', 'scale': 1e-9, 'xscale': 'log', 'label': P_NORM_LABEL},
        'right_plot': {'component': 'angle_deg', 'xscale': 'log', 'label': P_ANGLE_LABEL},
    },
    'potential_coefficients_composition_sheath': {
        'suptitle': r'P.u.l. sheath potential coefficient ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 1, 'q': 1, 'kind': 'potential_coefficient', 'comsol': False,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'abs', 'scale': 1e-9, 'xscale': 'log', 'label': P_NORM_LABEL},
        'right_plot': {'component': 'angle_deg', 'xscale': 'log', 'label': P_ANGLE_LABEL},
    },
    'potential_coefficients_composition_core_sheath': {
        'suptitle': r'P.u.l. core-sheath potential coefficient ($\rho_e=100 \;\Omega m, \epsilon_r=1$)',
        'p': 0, 'q': 1, 'kind': 'potential_coefficient', 'comsol': False,
        'scenario_key': COMPOSITION_SCENARIO_KEY, 'styles': COMPOSITION_STYLES,
        'left_plot': {'component': 'abs', 'scale': 1e-9, 'xscale': 'log', 'label': P_NORM_LABEL},
        'right_plot': {'component': 'angle_deg', 'xscale': 'log', 'label': P_ANGLE_LABEL},
    },
}
