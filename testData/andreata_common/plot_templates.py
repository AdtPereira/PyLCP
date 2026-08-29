""" PLOT_CONFIG series/axis templates shared by the andreata_case1/2/3 cases.

Base layer (MATLAB_TEMPLATE, PLOT_TPL_*): identical across the three plot_config.py.
Soil-formulation comparison layer (COMSOL_TEMPLATE, XUE_TEMPLATE,
DECONTI_TEMPLATE): used by case1 and case3, which compare earth-return
formulations (magalhaes_xue / deconti) against 3 soil profiles. case2 does not
import this second layer -- its comparison axis is the duct model
(bare/ERS/GMD), not the soil formulation, and it defines its own local
templates (a 1-profile COMSOL_TEMPLATE, DUCT_MODEL_TEMPLATE, etc.) in
testData/andreata_case2/plot_config.py.
"""

# --- Base layer: used by the three cases ---

MATLAB_TEMPLATE = [
    {
        'key': 'measured',
        'label': 'MATLAB',
        'marker': '.', 's': 18, 'color': 'blue', 'linewidths': 1.0, 'zorder': 11
    }
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

# --- Soil-formulation comparison layer: used by case1 and case3 ---

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
