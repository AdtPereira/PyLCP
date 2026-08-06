COMSOL_TEMPLATE = [
    {
        'key': 'rho_g_100_epsr1_1_mf',
        'label': r'COMSOL ($\rho_e=100, \epsilon_r=1$)',
        'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 10
    },
]

MATLAB_TEMPLATE = [
    {
        'key': 'measured',
        'label': 'MATLAB',
        'marker': 'x', 's': 12, 'color': 'blue', 'linewidths': 1.0, 'zorder': 11
    }
]

# Referência de validação cruzada (não é dado do duto): MATLAB do
# andreata_case1, mesma geometria de 3 cabos flat sem duto -- deve concordar
# com o cenário '1' ("Underground", duto ignorado) plotado junto (linha preta
# sólida do DUCT_MODEL_TEMPLATE abaixo).
MATLAB_CASE1_NO_DUCT_TEMPLATE = [
    {
        'key': 'case1_no_duct',
        'label': 'MATLAB (andreata_case1, sem duto)',
        'marker': '+', 's': 30, 'color': 'darkorange', 'linewidths': 1.2, 'zorder': 12
    }
]

# Os três modelos analíticos de duto comparados neste caso (Configuração 2 /
# Figura 5.2): ignorar o duto por completo, substituí-lo por uma isolação
# equivalente ponderada por área (ERS), ou pelo método GMD (Lafaia, 2015).
DUCT_MODEL_TEMPLATE = [
    {
        'key': '1',
        'type': {
            'main': {'label': 'Underground (ignoring HDPE)', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
        }
    },
    {
        'key': '2',
        'type': {
            'main': {'label': 'Area-weighted ERS', 'color': 'tab:blue', 'linestyle': '--', 'linewidth': 1.5},
        }
    },
    {
        'key': '3',
        'type': {
            'main': {'label': 'GMD case 3.1', 'color': 'tab:red', 'linestyle': '-.', 'linewidth': 1.5},
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

    # --- Parâmetros internos (núcleo + blindagem, fase A): sensíveis ao
    # modelo de duto (bare/ERS/GMD) porque ERS/GMD alteram a espessura e/ou
    # permissividade da isolação externa da blindagem. ---

    'core_self_impedance': {
        'suptitle': 'P.u.l. Core Self-Impedance, phase A ($Z_{cc}$) — HDPE duct modeling comparison',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['internal_matrices', 'impedance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE + MATLAB_CASE1_NO_DUCT_TEMPLATE,
        'matlab_matrix_key': 'internal_impedance_matrix',
        'p': 0, 'q': 0,
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rcc \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lcc \, (mH/km)$',
        }
    },

    'mutual_impedance_core_sheath': {
        'suptitle': 'P.u.l. Core-Sheath Mutual Impedance, phase A ($Z_{cs}$) — HDPE duct modeling comparison',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['internal_matrices', 'impedance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE + MATLAB_CASE1_NO_DUCT_TEMPLATE,
        'matlab_matrix_key': 'internal_impedance_matrix',
        'p': 0, 'q': 1,
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rcs \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lcs \, (mH/km)$',
        }
    },

    'sheath_self_impedance': {
        'suptitle': 'P.u.l. Sheath Self-Impedance, phase A ($Z_{ss}$) — HDPE duct modeling comparison',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['internal_matrices', 'impedance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE + MATLAB_CASE1_NO_DUCT_TEMPLATE,
        'matlab_matrix_key': 'internal_impedance_matrix',
        'p': 1, 'q': 1,
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rss \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lss \, (mH/km)$',
        }
    },

    'core_self_admittance': {
        'suptitle': 'P.u.l. Core Self-Admittance, phase A ($Y_{cc}$) — HDPE duct modeling comparison',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['internal_matrices', 'shunt_admittance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE + MATLAB_CASE1_NO_DUCT_TEMPLATE,
        'matlab_matrix_key': 'internal_admittance_matrix',
        'p': 0, 'q': 0,
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$Gcc \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$Ccc \, (nF/km)$',
        }
    },

    'mutual_admittance_core_sheath': {
        'suptitle': 'P.u.l. Core-Sheath Mutual Admittance, phase A ($Y_{cs}$) — HDPE duct modeling comparison',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['internal_matrices', 'shunt_admittance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE + MATLAB_CASE1_NO_DUCT_TEMPLATE,
        'matlab_matrix_key': 'internal_admittance_matrix',
        'p': 0, 'q': 1,
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$Gcs \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$Ccs \, (nF/km)$',
        }
    },

    'sheath_self_admittance': {
        'suptitle': 'P.u.l. Sheath Self-Admittance, phase A ($Y_{ss}$) — sensitive to the duct dielectric model',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['internal_matrices', 'shunt_admittance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE + MATLAB_CASE1_NO_DUCT_TEMPLATE,
        'matlab_matrix_key': 'internal_admittance_matrix',
        'p': 1, 'q': 1,
        'xlim': (1E4, 1E7),
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$Gss \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$Css \, (nF/km)$',
        }
    },

    # --- Retorno à terra (acoplamento entre fases), fase A: aproximação em
    # todos os três cenários -- a formulação de Zg/Yg não enxerga o duto,
    # apenas a posição/raio externo de cada cabo. A comparação contra
    # COMSOL/MATLAB é a única validação de fato disponível aqui. ---

    'earth_return_impedance_phase_a': {
        'suptitle': 'P.u.l. Self earth-return impedance of phase-a [Xue, 2018] — duct ignored by the formulation',
        'p': 0, 'q': 0,
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'impedance_matrix',
        'matlab_series_to_plot': MATLAB_TEMPLATE + MATLAB_CASE1_NO_DUCT_TEMPLATE,
        'matlab_matrix_key': 'earth_return_impedance_matrix',
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rg_{11} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lg_{11} \, (mH/km)$',
        }
    },

    'earth_return_admittance_phase_a': {
        'suptitle': 'P.u.l. Self earth-return admittance of phase-a [Xue, 2018] — duct ignored by the formulation',
        'p': 0, 'q': 0,
        'series_to_plot': DUCT_MODEL_TEMPLATE,
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

    'earth_return_potential_coeff_phase_a': {
        'suptitle': 'P.u.l. Self earth-return potential coefficient of phase-a [Xue, 2018]',
        'p': 0, 'q': 0,
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['earth_return_parameters', 'potential_coefficient'],
        'matlab_series_to_plot': MATLAB_TEMPLATE + MATLAB_CASE1_NO_DUCT_TEMPLATE,
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
}
