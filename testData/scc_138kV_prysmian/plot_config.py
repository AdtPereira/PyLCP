# PLOT_CONFIG for plotter.scc_models.PrysmianCableModels — Prysmian 138kV
# core+sheath cable, validated against the exact Bessel decomposition
# [Ametani, 2015] and COMSOL (Js method, core/sheath excitation).
#
# 'core'/'sheath'/'core_sheath'/'internal_parameters' are rendered by
# PrysmianCableModels._internal_impedance_subplots (named z-component dict,
# not a (p,q)-indexed matrix — doesn't fit BasePlotter's generic engine).
# 'ground_return_impedance' uses the generic BasePlotter._plot_matricial_parameters
# engine and expects pul_data['analytical'] = {'frequencies': ..., 'scenarios': pul_data['scenarios']}.

_NAKAGAWA_STYLES = {
    'magalhaes_xue': {'label': 'Magalhães/Xue', 'color': 'black', 'linestyle': '-', 'linewidth': 2},
    'sunde': {'label': 'Sunde', 'color': 'gray', 'linestyle': ':', 'linewidth': 2},
    'pollaczek': {'label': 'Pollaczek', 'color': 'darkblue', 'linestyle': '-', 'linewidth': 2},
    'ametani': {'label': 'Ametani', 'color': 'cyan', 'linestyle': '-', 'linewidth': 1},
    'deconti': {'label': 'De Conti', 'color': 'red', 'linestyle': '--', 'linewidth': 1},
    'saad': {'label': 'Saad', 'color': 'red', 'linestyle': '--', 'linewidth': 1},
}

PLOT_CONFIG = {
    'core': {
        'suptitle': r'Prysmian 138 kV SCC P.u.l. parameters of core conductor, $z_{cs}$ [Ametani, 2015]',
        'resistance_title': 'P.u.l. resistance',
        'inductance_title': 'P.u.l. inductance',
        'series_to_plot': [
            {'key': 'zcs', 'label': r'$z_{cs} = z_{11} + z_{12} + z_{2i}$', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
            {'key': 'z11', 'label': r'$z_{11}$: Internal impedance of core outer surface', 'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.0},
            {'key': 'z12', 'label': r'$z_{12}$: Core outer insulator impedance', 'color': 'darkblue', 'linestyle': '--', 'linewidth': 1.0},
            {'key': 'z2i', 'label': r'$z_{2i}$: Internal impedance of sheath inner surface', 'color': 'red', 'linestyle': '-.', 'linewidth': 1.0},
        ],
    },
    'sheath': {
        'suptitle': r'Prysmian 138 kV SCC P.u.l. parameters of sheath conductor, $z_{s3}$ [Ametani, 2015]',
        'resistance_title': 'P.u.l. resistance',
        'inductance_title': 'P.u.l. inductance',
        'series_to_plot': [
            {'key': 'zs3', 'label': r'$z_{s3} = z_{20} + z_{23}$', 'color': 'black', 'linestyle': '-', 'linewidth': 2.0},
            {'key': 'z20', 'label': r'$z_{20}$: Internal impedance of sheath outer surface', 'color': 'red', 'linestyle': '--', 'linewidth': 1.0},
            {'key': 'z2m', 'label': r'$z_{2m}$: Sheath mutual impedance', 'color': 'darkblue', 'linestyle': '-.', 'linewidth': 1.0},
            {'key': 'z23', 'label': r'$z_{23}$: Sheath outer insulator impedance', 'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.0},
        ],
    },
    'core_sheath': {
        'suptitle': r'Prysmian 138 kV SCC P.u.l. Internal Impedance Matrix, $[z_i]$ [Ametani, 2015]',
        'resistance_title': 'P.u.l. resistance',
        'inductance_title': 'P.u.l. inductance',
        'resistance_ylim_bottom': 1e-4,
        'series_to_plot': [
            {'key': 'Zcc', 'label': r'$z_{cc} = z_{cs} + z_{s3} - 2z_{2m}$: Core self-impedance', 'color': 'black', 'linestyle': '-', 'linewidth': 1.0},
            {'key': 'Zcs', 'label': r'$z_{cs} = z_{20} + z_{23} - z_{2m}$: Mutual impedance between the core and sheath', 'color': 'darkblue', 'linestyle': '--', 'linewidth': 1.0},
            {'key': 'Zss', 'label': r'$z_{ss} = z_{20} + z_{23}$: Sheath self-impedance', 'color': 'darkgreen', 'linestyle': ':', 'linewidth': 1.0},
        ],
    },
    'internal_parameters': {
        'suptitle': r'Prysmian 138 kV SCC P.u.l. parameters of SCC [Ametani, 2015]',
        'resistance_title': 'P.u.l. resistance',
        'inductance_title': 'P.u.l. inductance',
        'series_to_plot': [
            {'key': 'z11', 'label': r'$z_{11}$: Internal impedance of core outer surface', 'color': 'black', 'linestyle': '-', 'linewidth': 1.0},
            {'key': 'z2m', 'label': r'$z_{2m}$: Sheath mutual impedance', 'color': 'darkblue', 'linestyle': '-', 'linewidth': 1.0},
            {'key': 'z2i', 'label': r'$z_{2i}$: Internal impedance of sheath inner surface', 'color': 'darkgreen', 'linestyle': '-', 'linewidth': 1.0},
            {'key': 'z20', 'label': r'$z_{20}$: Internal impedance of sheath outer surface', 'color': 'darkgray', 'linestyle': '-', 'linewidth': 1.0},
        ],
    },
    'ground_return_impedance': {
        'suptitle': r'Prysmian 138 kV SCC P.u.l. ground-return impedance for $\rho_1=1000 \;\Omega$ m and $\epsilon_{r1}=1$',
        'left_plot': {'component': 'real', 'xscale': 'log', 'label': r'$R_g \, (\Omega/m)$'},
        'right_plot': {'component': 'imag_div_w', 'scale': 1e6, 'xscale': 'log', 'label': r'$L_g \, (\mu H/m)$'},
        'data_series': [
            {
                'source': 'analytical', 'scenario_key': key,
                'data_path': ['earth_return_parameters', 'impedance_matrix'],
                'series': {'main': {'p': 0, 'q': 0, **style}},
            }
            for key, style in _NAKAGAWA_STYLES.items()
        ],
    },
}
