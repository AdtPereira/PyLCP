""" Templates de série/eixo de PLOT_CONFIG compartilhados pelos casos andreata_case1/2/3.

Camada base (MATLAB_TEMPLATE, PLOT_TPL_*): idêntica nos três plot_config.py.
Camada de comparação de formulação de solo (COMSOL_TEMPLATE, XUE_TEMPLATE,
DECONTI_TEMPLATE): usada por case1 e case3, que comparam formulações de
retorno à terra (magalhaes_xue / deconti) contra 3 perfis de solo. case2 não
importa essa segunda camada -- seu eixo de comparação é o modelo de duto
(bare/ERS/GMD), não a formulação de solo, e define seus próprios templates
locais (COMSOL_TEMPLATE de 1 perfil, DUCT_MODEL_TEMPLATE etc.) em
testData/andreata_case2/plot_config.py.
"""

# --- Camada base: usada pelos três casos ---

MATLAB_TEMPLATE = [
    {
        'key': 'measured',
        'label': 'MATLAB',
        'marker': '.', 's': 12, 'color': 'blue', 'linewidths': 1.0, 'zorder': 11
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

# --- Camada de comparação de formulação de solo: usada por case1 e case3 ---

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
