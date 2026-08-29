from ..andreata_common.plot_templates import (
    MATLAB_TEMPLATE,
    PLOT_TPL_REAL, PLOT_TPL_IMAG, PLOT_TPL_RESISTANCE, PLOT_TPL_INDUCTANCE,
    PLOT_TPL_CONDUCTANCE, PLOT_TPL_CAPACITANCE,
)

# pyLCP FEM-hybrid path: Zi/Yi from COMSOL (exact eccentric geometry,
# air + HDPE tube) + analytical earth return (Magalhaes/Xue) with the tube
# outer radius. This is the line compared against MATLAB (Andreata's FEM).
FEM_HYBRID_TEMPLATE = [
    {
        'key': 'fem',
        'type': {
            'main': {'label': 'pyLCP (FEM-hybrid)', 'color': 'black',
                     'linestyle': '-', 'linewidth': 1.5},
        }
    }
]

COMSOL_TEMPLATE = [
    {
        'key': 'rho_g_100_epsr1_1_mf',
        'label': r'COMSOL ($\rho_e=100, \epsilon_r=1$)',
        'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 10
    },
]

# Cross-validation reference (not duct data): andreata_case1 MATLAB, same
# geometry of 3 flat cables without a duct -- it must agree with scenario '1'
# ("Underground", duct ignored) plotted alongside (the solid black line of
# the DUCT_MODEL_TEMPLATE below).
MATLAB_CASE1_NO_DUCT_TEMPLATE = [
    {
        'key': 'case1_no_duct',
        'label': 'MATLAB (andreata_case1, no duct)',
        'marker': '+', 's': 30, 'color': 'darkorange', 'linewidths': 1.2, 'zorder': 12
    }
]

# The three analytical duct models compared in this case (Configuration 2 /
# Figure 5.2): ignore the duct entirely, replace it with an area-weighted
# equivalent insulation (ERS), or with the GMD method (Lafaia, 2015).
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

# COMSOL internal impedance (legacy measurement, single cable, Js method) --
# reused from the previous single-phase study (hdpe_300mm2): the
# core+sheath+duct cross-section is identical for each phase, only the earth
# return changes with the number of phases (see andreata_case2.py). One color
# per component (cc/cs/ss) to tell them apart visually -- unlike the original
# hdpe_300mm2, which uses the same black marker for all three.
INTERNAL_COMSOL_TEMPLATE = {
    'cc': {'key': 'measured', 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black',
           'zorder': 10, 'label': r'COMSOL ($J_s$ method)'},
    'cs': {'key': 'measured', 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'tab:red',
           'zorder': 10, 'label': ''},
    'ss': {'key': 'measured', 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'tab:blue',
           'zorder': 10, 'label': ''},
}

# Same style/colors as INTERNAL_COMSOL_TEMPLATE (impedance, Js method),
# reused for the internal admittance -- but with its own legend label, since
# the admittance comes from the direct charge method
# (cmsl_internal_admittance_charge_method.txt), not the Js method.
ADMITTANCE_COMSOL_TEMPLATE = {
    key: {**style, 'label': r'COMSOL (charge method)'} if key == 'cc' else style
    for key, style in INTERNAL_COMSOL_TEMPLATE.items()
}

PLOT_CONFIG = {

    # --- Combined internal parameters (core + sheath, phase A), a single
    # duct model (GMD case 3.1, the most complete of the three) against the
    # COMSOL/MATLAB references -- same plot already present in the single-phase
    # case hdpe_300mm2 (internal_impedance_matrix.png), with the MATLAB reading
    # ('measured', prefix 'andreata_hdpe') already prepared for when the duct's
    # own data arrives (see andreata_case2.py, Item 9 of BUGS_AND_FIXES.md).
    # For the per-duct-model comparison (Underground/ERS/GMD), see
    # 'core_self_impedance' etc. below. ---

    'internal_impedance_matrix': {
        'suptitle': r'P.u.l. Internal Impedance Matrix, phase A — GMD case 3.1 [Lafaia, 2015] vs. COMSOL ($J_s$) [Yin, 1990]',
        'components': [
            {
                'p': 0, 'q': 0,
                # 'internal_style': {'label': r'$Zi_{cc}$ GMD case 3.1', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                'comsol_series_to_plot': [INTERNAL_COMSOL_TEMPLATE['cc']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB', 'marker': '.', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                # 'internal_style': {'label': r'$Zi_{cs}$ GMD case 3.1', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'comsol_series_to_plot': [INTERNAL_COMSOL_TEMPLATE['cs']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:red', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                # 'internal_style': {'label': r'$Zi_{ss}$ GMD case 3.1', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'comsol_series_to_plot': [INTERNAL_COMSOL_TEMPLATE['ss']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:blue', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
        ],
        'comsol_matrix_key': 'impedance_matrix',
        'matlab_matrix_key': 'internal_impedance_matrix',
        'xlim': (1E0, 1E7),
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

    'internal_admittance_matrix': {
        'suptitle': r'P.u.l. Internal Admittance Matrix, phase A — GMD case 3.1 [Lafaia, 2015] vs. COMSOL',
        'internal_matrix_key': 'shunt_admittance_matrix',
        'components': [
            {
                'p': 0, 'q': 0,
                # 'internal_style': {'label': r'$Yi_{cc}$ GMD case 3.1', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                'comsol_series_to_plot': [ADMITTANCE_COMSOL_TEMPLATE['cc']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB', 'marker': '.', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                # 'internal_style': {'label': r'$Yi_{cs}$ GMD case 3.1', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'comsol_series_to_plot': [ADMITTANCE_COMSOL_TEMPLATE['cs']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:red', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                # 'internal_style': {'label': r'$Yi_{ss}$ GMD case 3.1', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'comsol_series_to_plot': [ADMITTANCE_COMSOL_TEMPLATE['ss']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:blue', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
        ],
        'comsol_matrix_key': 'admittance_matrix',
        'matlab_matrix_key': 'internal_admittance_matrix',
        'xlim': (1E3, 1E6),
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C \, (nF/km)$',
        }
    },

    # --- Internal parameters (core + sheath, phase A): sensitive to the duct
    # model (bare/ERS/GMD) because ERS/GMD change the thickness and/or
    # permittivity of the sheath's outer insulation. ---

    'core_self_impedance': {
        'suptitle': 'P.u.l. Core Self-Impedance, phase A ($Z_{cc}$) — HDPE duct modeling comparison',
        # 'series_to_plot': DUCT_MODEL_TEMPLATE,
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
        # 'series_to_plot': DUCT_MODEL_TEMPLATE,
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
        # 'series_to_plot': DUCT_MODEL_TEMPLATE,
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
        # 'series_to_plot': DUCT_MODEL_TEMPLATE,
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
        # 'series_to_plot': DUCT_MODEL_TEMPLATE,
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
        # 'series_to_plot': DUCT_MODEL_TEMPLATE,
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

    # --- FEM-hybrid (pyLCP) vs. MATLAB (Andreata's FEM) --------------------
    # Now the HDPE duct IS modeled: Zi/Yi come from COMSOL (exact eccentric
    # geometry) and the earth return is analytical with the tube outer radius
    # (Part A). Equivalent to the self_*/earth_return_* plots of case1. ---

    'self_impedance_phase_a_sheath': {
        'suptitle': 'P.u.l. self-impedance of phase-a sheath — pyLCP (FEM-hybrid) vs. MATLAB',
        'series_to_plot': FEM_HYBRID_TEMPLATE,
        'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'series_impedance_matrix',
        'p': 1, 'q': 1,
        'xlim': (1E-2, 1E7),
        'left_plot': {**PLOT_TPL_RESISTANCE, 'label': r'$Rs_{22} \, (\Omega/km)$'},
        'right_plot': {**PLOT_TPL_INDUCTANCE, 'label': r'$Ls_{22} \, (mH/km)$'},
    },

    'self_admittance_phase_a_sheath': {
        'suptitle': 'P.u.l. self-admittance of phase-a sheath — pyLCP (FEM-hybrid) vs. MATLAB',
        'series_to_plot': FEM_HYBRID_TEMPLATE,
        'path': ['quasi_tem_matrices', 'shunt_admittance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'shunt_admittance_matrix',
        'p': 1, 'q': 1,
        'xlim': (1E-2, 1E7),
        'left_plot': {**PLOT_TPL_CONDUCTANCE, 'label': r'$G_{22} \, (S/km)$'},
        'right_plot': {**PLOT_TPL_CAPACITANCE, 'label': r'$C_{22} \, (nF/km)$'},
    },

    'earth_return_impedance_phase_a': {
        'suptitle': 'P.u.l. self earth-return impedance of phase-a — pyLCP (FEM-hybrid, tube radius) vs. MATLAB',
        'p': 0, 'q': 0,
        'series_to_plot': FEM_HYBRID_TEMPLATE,
        'path': ['earth_return_parameters', 'impedance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
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
        # 'series_to_plot': DUCT_MODEL_TEMPLATE,
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
        'suptitle': 'P.u.l. self earth-return potential coefficient of phase-a — pyLCP (FEM-hybrid) vs. MATLAB',
        'p': 0, 'q': 0,
        'series_to_plot': FEM_HYBRID_TEMPLATE,
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
}
