from ..andreata_common.plot_templates import (
    MATLAB_TEMPLATE,
    PLOT_TPL_REAL, PLOT_TPL_IMAG, PLOT_TPL_RESISTANCE, PLOT_TPL_INDUCTANCE,
    PLOT_TPL_CONDUCTANCE, PLOT_TPL_CAPACITANCE,
)

# Não há dado COMSOL para este caso (nem retorno à terra, nem interno) --
# nenhum arquivo cmsl_*.txt existe em Results/. Template mantido guardado
# (mesmo perfil de solo único do andreata_case2), pronto para quando/se o
# dado chegar; hoje nenhuma config abaixo o referencia ativamente.
COMSOL_TEMPLATE = [
    {
        'key': 'rho_g_100_epsr1_1_mf',
        'label': r'COMSOL ($\rho_e=100, \epsilon_r=1$)',
        'marker': 'o', 's': 20, 'facecolors': 'none', 'edgecolors': 'black', 'zorder': 10
    },
]

# Os três modelos analíticos de duto comparados neste caso (Configuração 4 /
# Figura 5.4): ignorar o duto por completo (mesma física heterogênea do
# andreata_case3), substituí-lo por uma isolação equivalente ponderada por
# área (ERS), ou pelo método GMD (Lafaia, 2015) -- aplicados às 3 bainhas
# SCC apenas; o ECC nunca é afetado por esses três cenários (ver
# andreata_case4.py / README.md, seção Limitações).
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


# COMSOL de impedância interna (medição legada, cabo único, método Js) --
# reaproveitado do estudo monofásico anterior (hdpe_300mm2): a seção
# transversal núcleo+blindagem+duto é idêntica em cada fase, só o retorno à
# terra muda com o número de fases (ver andreata_case2.py). Uma cor por
# componente (cc/cs/ss) para diferenciá-los visualmente -- ao contrário do
# hdpe_300mm2 original, que usa o mesmo marcador preto para os três.
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


PLOT_CONFIG = {

    # --- Parâmetros internos combinados (núcleo + blindagem da fase A + ECC),
    # um só modelo de duto (GMD case 3.1, o mais completo dos três) contra a
    # referência MATLAB -- mesmo padrão de andreata_case3, estendido com o
    # componente 6,6 (ECC). Para a comparação por modelo de duto
    # (Underground/ERS/GMD), ver 'core_self_impedance' etc. abaixo. ---

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
                'matlab_series_to_plot': [
                    {'key': 'measured', 'label': 'MATLAB', 'marker': '.', 's': 8, 'color': 'black', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 0, 'q': 1,
                # 'internal_style': {'label': r'$Yi_{cs}$', 'color': 'tab:red', 'linestyle': '--', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:red', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 1, 'q': 1,
                # 'internal_style': {'label': r'$Yi_{ss}$', 'color': 'tab:blue', 'linestyle': '-.', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:blue', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
            {
                'p': 6, 'q': 6,
                # 'internal_style': {'label': r'$Yi_{77}$ (ECC)', 'color': 'tab:green', 'linestyle': ':', 'linewidth': 1.5},
                'matlab_series_to_plot': [
                    {'key': 'measured', 'marker': '.', 's': 8, 'color': 'tab:green', 'linewidths': 1.0, 'zorder': 11},
                ],
            },
        ],
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

    # --- Parâmetros internos (núcleo + blindagem, fase A): sensíveis ao
    # modelo de duto (bare/ERS/GMD) porque ERS/GMD alteram a espessura e/ou
    # permissividade da isolação externa da blindagem. Mesmo eixo de
    # comparação do andreata_case2, aqui sobre o modelo heterogêneo (3 SCC +
    # ECC) do andreata_case4. ---

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

    # --- Retorno à terra (acoplamento entre cabos), fase A: aproximação em
    # todos os três cenários -- a formulação de Zg/Yg não enxerga o duto,
    # apenas a posição/raio externo de cada cabo (ver Decisão de
    # arquitetura, README.md). A comparação contra MATLAB é a validação de
    # fato disponível aqui. ---

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

    # --- Condutor ECC: as três curvas analíticas (Underground/ERS/GMD)
    # tendem a coincidir aqui, já que nenhum dos três cenários altera o ECC
    # em si (só as bainhas SCC) -- ver nota em andreata_case4.py. A
    # referência real do efeito do duto sobre o ECC é a comparação com
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
