from ..andreata_common.plot_templates import (
    MATLAB_TEMPLATE,
    PLOT_TPL_REAL, PLOT_TPL_IMAG, PLOT_TPL_RESISTANCE, PLOT_TPL_INDUCTANCE,
    PLOT_TPL_CONDUCTANCE, PLOT_TPL_CAPACITANCE,
)

# There is no COMSOL earth-return data for this case -- no file
# 'cmsl_ground_return_impedance.txt' exists in Results/ (the internal
# impedance and admittance, on the other hand, do have their own data -- see
# INTERNAL_COMSOL_TEMPLATE below). Template kept on hold (same single-soil
# profile as andreata_case2), ready for when/if the data arrives; today no
# config below references it actively.
COMSOL_TEMPLATE = [
    {
        'key': 'rho_g_100_epsr1_1_mf',
        'label': r'COMSOL ($\rho_e=100, \epsilon_r=1$)',
        'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 10
    },
]

# The three analytical duct models compared in this case (Configuration 4 /
# Figure 5.4): ignore the duct entirely (same heterogeneous physics as
# andreata_case3), replace it with an area-weighted equivalent insulation
# (ERS), or with the GMD method (Lafaia, 2015) -- applied to the 3 SCC
# sheaths only; the ECC is never affected by these three scenarios (see
# andreata_case4.py / README.md, Limitations section).
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


# COMSOL internal impedance (3 conductors core/sheath/ECC, specific to
# case4, Js method -- see README.md). One color per component (cc/cs/ss/ecc)
# to tell them apart visually.
INTERNAL_COMSOL_TEMPLATE = {
    'cc': {'key': 'measured', 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black',
           'zorder': 10, 'label': r'COMSOL ($J_s$ method)'},
    'cs': {'key': 'measured', 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'tab:red',
           'zorder': 10, 'label': ''},
    'ss': {'key': 'measured', 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'tab:blue',
           'zorder': 10, 'label': ''},
    'ecc': {'key': 'measured', 'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'tab:green',
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

    # --- Combined internal parameters (phase A core + sheath + ECC),
    # a single duct model (GMD case 3.1, the most complete of the three)
    # against the MATLAB reference -- same pattern as andreata_case3, extended
    # with the 6,6 component (ECC). For the per-duct-model comparison
    # (Underground/ERS/GMD), see 'core_self_impedance' etc. below. ---

    'internal_impedance_matrix': {
        'suptitle': 'P.u.l. Internal Impedance Matrix, phase A + ECC — GMD case 3.1 [Lafaia, 2015] vs. MATLAB',
        'components': [
            {
                'p': 0, 'q': 0,
                # 'internal_style': {'label': r'$Zi_{cc}$', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                'comsol_series_to_plot': [INTERNAL_COMSOL_TEMPLATE['cc']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB', 'marker': '.', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                # 'internal_style': {'label': r'$Zi_{cs}$', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'comsol_series_to_plot': [INTERNAL_COMSOL_TEMPLATE['cs']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:red', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                # 'internal_style': {'label': r'$Zi_{ss}$', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'comsol_series_to_plot': [INTERNAL_COMSOL_TEMPLATE['ss']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:blue', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 6, 'q': 6,
                # 'internal_style': {'label': r'$Zi_{77}$ (ECC)', 'color': 'tab:green', 'linestyle': ':', 'linewidth': 1.5},
                'comsol_series_to_plot': [INTERNAL_COMSOL_TEMPLATE['ecc']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:green', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
        ],
        'comsol_matrix_key': 'impedance_matrix',
        'matlab_matrix_key': 'internal_impedance_matrix',
        'xlim': (1E1, 1E7),
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
        'suptitle': 'P.u.l. Internal Admittance Matrix, phase A + ECC — GMD case 3.1 [Lafaia, 2015] vs. MATLAB',
        'internal_matrix_key': 'shunt_admittance_matrix',
        'components': [
            {
                'p': 0, 'q': 0,
                # 'internal_style': {'label': r'$Yi_{cc}$', 'color': 'black', 'linestyle': '-', 'linewidth': 1.5},
                'comsol_series_to_plot': [ADMITTANCE_COMSOL_TEMPLATE['cc']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB', 'marker': '.', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                # 'internal_style': {'label': r'$Yi_{cs}$', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'comsol_series_to_plot': [ADMITTANCE_COMSOL_TEMPLATE['cs']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:red', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                # 'internal_style': {'label': r'$Yi_{ss}$', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'comsol_series_to_plot': [ADMITTANCE_COMSOL_TEMPLATE['ss']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:blue', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 6, 'q': 6,
                # 'internal_style': {'label': r'$Yi_{77}$ (ECC)', 'color': 'tab:green', 'linestyle': ':', 'linewidth': 1.5},
                'comsol_series_to_plot': [ADMITTANCE_COMSOL_TEMPLATE['ecc']],
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:green', 'linewidths': 1.0, 'zorder': 11},
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
    # permittivity of the sheath's outer insulation. Same comparison axis as
    # andreata_case2, here on the heterogeneous model (3 SCC + ECC) of
    # andreata_case4. ---

    'core_self_impedance': {
        'suptitle': 'P.u.l. Core Self-Impedance, phase A ($Z_{cc}$) — HDPE duct modeling comparison',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['internal_matrices', 'impedance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
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
        'matlab_series_to_plot': MATLAB_TEMPLATE,
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
        'matlab_series_to_plot': MATLAB_TEMPLATE,
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

    # --- Earth return (inter-cable coupling), phase A: an approximation in
    # all three scenarios -- the Zg/Yg formulation does not see the duct,
    # only the position / outer radius of each cable (see the Architecture
    # Decision, README.md). The comparison against MATLAB is the actual
    # validation available here. ---

    'earth_return_impedance_phase_a': {
        'suptitle': 'P.u.l. Self earth-return impedance of phase-a [Xue, 2018] — duct ignored by the formulation',
        'p': 0, 'q': 0,
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['earth_return_parameters', 'impedance_matrix'],
        'comsol_series_to_plot': COMSOL_TEMPLATE,
        'comsol_matrix_key': 'impedance_matrix',
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

    # --- ECC conductor: the three analytical curves (Underground/ERS/GMD)
    # tend to coincide here, since none of the three scenarios changes the ECC
    # itself (only the SCC sheaths) -- see the note in andreata_case4.py. The
    # actual reference for the duct's effect on the ECC is the comparison with
    # MATLAB. ---

    'self_impedance_ecc': {
        'suptitle': 'P.u.l. Self-impedance of ECC — HDPE duct modeling comparison',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['quasi_tem_matrices', 'series_impedance_matrix'],
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
        'suptitle': 'P.u.l. Self-admittance of ECC — HDPE duct modeling comparison',
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['quasi_tem_matrices', 'shunt_admittance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'shunt_admittance_matrix',
        'p': 6, 'q': 6,
        'xlim': (1E4, 1E7),
        'left_plot': {
            **PLOT_TPL_CONDUCTANCE,
            'label': r'$G_{77} \, (S/km)$',
        },
        'right_plot': {
            **PLOT_TPL_CAPACITANCE,
            'label': r'$C_{77} \, (nF/km)$',
        }
    },

    'earth_return_impedance_ecc': {
        'suptitle': 'P.u.l. Self earth-return impedance of ECC — HDPE duct modeling comparison',
        'p': 6, 'q': 6,
        'series_to_plot': DUCT_MODEL_TEMPLATE,
        'path': ['quasi_tem_matrices', 'earth_return_impedance_matrix'],
        'matlab_series_to_plot': MATLAB_TEMPLATE,
        'matlab_matrix_key': 'earth_return_impedance_matrix',
        'xlim': (1E-2, 1E7),
        'left_plot': {
            **PLOT_TPL_RESISTANCE,
            'label': r'$Rg_{77} \, (\Omega/km)$',
        },
        'right_plot': {
            **PLOT_TPL_INDUCTANCE,
            'label': r'$Lg_{77} \, (mH/km)$',
        }
    },

    'earth_return_potential_coeff_ecc': {
        'suptitle': 'P.u.l. Self earth-return potential coefficient of ECC — HDPE duct modeling comparison',
        'p': 6, 'q': 6,
        'series_to_plot': DUCT_MODEL_TEMPLATE,
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
