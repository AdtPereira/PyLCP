# PLOT_CONFIG for plotter.ohtl_plotter.OHTLPlotter — Xue (2018) Sec. 4.3 overhead line.
#
# Expects pul_data in FLAT format (see ohtl_xue_sec43.py):
#   pul_data[scenario_key][matrix_key] -> ndarray (n_freq, N, N)
#   pul_data['frequencies']            -> ndarray (n_freq,)
#   pul_data['comsol']                 -> dict[str, pandas.DataFrame] | None

R_LABEL = r'$R_s\;(\Omega/\mathrm{km})$'
L_LABEL = r'$L_s\;(\mathrm{mH/km})$'
G_LABEL = r'$G\;(\mathrm{S/km})$'
C_LABEL = r'$C\;(\mathrm{nF/km})$'
ATTEN_LABEL = r'Attenuation Constant, $\alpha_\nu$ (Np/km)'
PHASE_VEL_LABEL = r'Phase Velocity, $c_\nu/c_0$'

# Nakagawa-only scenarios (fig42, fig45, fig47)
OVERHEAD_LINE_SERIES = [
    {'key': 'p100',      'type': {'main': {'label': r'$\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=1$',  'color': 'black', 'linestyle': '-'}}},
    {'key': 'p100_er20', 'type': {'main': {'label': r'$\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=20$', 'color': 'black', 'linestyle': '--'}}},
    {'key': 'p2000',     'type': {'main': {'label': r'$\rho_e=2000\;\Omega\mathrm{m},\;\epsilon_r=1$', 'color': 'black', 'linestyle': '-.'}}},
]

# Nakagawa vs. Carson comparison at rho_e=100 Ohm.m (fig43, fig46)
NAKAGAWA_CARSON_SERIES = [
    {'key': 'p100',        'type': {'main': {'label': r'$\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=1$ (Nakagawa)',  'color': 'black', 'linestyle': '-'}}},
    {'key': 'p100_er20',   'type': {'main': {'label': r'$\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=20$ (Nakagawa)', 'color': 'black', 'linestyle': '-.'}}},
    {'key': 'p100_carson', 'type': {'main': {'label': r'$\rho_e=100\;\Omega\mathrm{m},\;\epsilon_r=1$ (Carson)',    'color': 'red',   'linestyle': '--'}}},
]

# Nakagawa vs. Carson comparison at rho_e=100 and 2000 Ohm.m (fig48)
ATTENUATION_CONSTANT_SERIES = [
    {'key': 'p100',         'type': {'main': {'label': r'$\rho_e=100\;\Omega\mathrm{m}$ (Nakagawa)',  'color': 'black', 'linestyle': '-'}}},
    {'key': 'p100_carson',  'type': {'main': {'label': r'$\rho_e=100\;\Omega\mathrm{m}$ (Carson)',    'color': 'red',   'linestyle': '-'}}},
    {'key': 'p2000',        'type': {'main': {'label': r'$\rho_e=2000\;\Omega\mathrm{m}$ (Nakagawa)', 'color': 'black', 'linestyle': '--'}}},
    {'key': 'p2000_carson', 'type': {'main': {'label': r'$\rho_e=2000\;\Omega\mathrm{m}$ (Carson)',   'color': 'red',   'linestyle': '--'}}},
]

PLOT_CONFIG = {
    'fig42': {
        'suptitle': 'Figure 4.2: P.u.l. series impedance with Nakagawa formulation [Xue, 2018]',
        'matrix_key': 'series_impedance_matrix', 'p': 0, 'q': 0,
        'series_to_plot': OVERHEAD_LINE_SERIES,
        'comsol_key': 'cmsl_ground_return_impedance_h5',
        'left_plot':  {'component': 'real',        'scale': 1e3, 'xscale': 'log', 'yscale': 'log', 'label': R_LABEL, 'title': 'P.u.l. series resistance',  'y_lim': (1E0, 1E5)},
        'right_plot': {'component': 'imag_div_w',   'scale': 1e6, 'xscale': 'log',                  'label': L_LABEL, 'title': 'P.u.l. series inductance', 'y_lim': (1, 2.5)},
    },
    'fig43': {
        'suptitle': 'Figure 4.3: P.u.l. series impedance comparison [Xue, 2018]',
        'matrix_key': 'series_impedance_matrix', 'p': 0, 'q': 0,
        'series_to_plot': NAKAGAWA_CARSON_SERIES,
        'comsol_key': 'cmsl_ground_return_impedance_h5',
        'left_plot':  {'component': 'real',        'scale': 1e3, 'xscale': 'log', 'yscale': 'log', 'label': R_LABEL, 'title': 'P.u.l. series resistance',  'y_lim': (1E0, 1E4)},
        'right_plot': {'component': 'imag_div_w',   'scale': 1e6, 'xscale': 'log',                  'label': L_LABEL, 'title': 'P.u.l. series inductance', 'y_lim': (1.4, 2.2)},
    },
    'fig45': {
        'suptitle': 'Figure 4.5: P.u.l. shunt admittance with Nakagawa formulation [Xue, 2018]',
        'matrix_key': 'shunt_admittance_matrix', 'p': 0, 'q': 0,
        'series_to_plot': OVERHEAD_LINE_SERIES,
        'left_plot':  {'component': 'real',        'scale': 1e3,  'xscale': 'log', 'label': G_LABEL, 'title': 'P.u.l. shunt conductance',  'y_lim': (-0.3, 0.1)},
        'right_plot': {'component': 'imag_div_w',   'scale': 1e12, 'xscale': 'log', 'label': C_LABEL, 'title': 'P.u.l. shunt capacitance', 'y_lim': (6.6, 7.4)},
    },
    'fig46': {
        'suptitle': 'Figure 4.6: P.u.l. shunt admittance comparison [Xue, 2018]',
        'matrix_key': 'shunt_admittance_matrix', 'p': 0, 'q': 0,
        'series_to_plot': NAKAGAWA_CARSON_SERIES,
        'left_plot':  {'component': 'real',        'scale': 1e3,  'xscale': 'log', 'label': G_LABEL, 'title': 'P.u.l. shunt conductance',  'y_lim': (-0.06, 0.02)},
        'right_plot': {'component': 'imag_div_w',   'scale': 1e12, 'xscale': 'log', 'label': C_LABEL, 'title': 'P.u.l. shunt capacitance', 'y_lim': (7.22, 7.32)},
    },
    'fig47': {
        'suptitle': 'Figure 4.7: Propagation constant with Nakagawa formulation [Xue, 2018]',
        'matrix_key': 'propagation_voltage_matrix', 'p': 0, 'q': 0,
        'series_to_plot': OVERHEAD_LINE_SERIES,
        'y_lim': {'attenuation': (1E-3, 1E1), 'phase_velocity': (0.8, 1.1)},
        'y_scale': {'attenuation': 'log'},
        'attenuation_title': 'Attenuation constant', 'phase_velocity_title': 'Normalized phase velocity',
    },
    'fig48': {
        'suptitle': r'Figure 4.8: Propagation constant comparison for $\varepsilon_r = 1$ [Xue, 2018]',
        'matrix_key': 'propagation_voltage_matrix', 'p': 0, 'q': 0,
        'series_to_plot': ATTENUATION_CONSTANT_SERIES,
        'y_lim': {'attenuation': (1E-3, 1E2), 'phase_velocity': (0.8, 1.1)},
        'y_scale': {'attenuation': 'log'},
        'attenuation_title': 'Attenuation constant', 'phase_velocity_title': 'Normalized phase velocity',
    },
}
