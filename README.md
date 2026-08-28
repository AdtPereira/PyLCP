# PyLCP — Python Library for Cable Parameters

Framework de análise eletromagnética para cálculo de parâmetros elétricos por unidade de comprimento (PUL) em **linhas de transmissão multicondutoras** (MTLs). Combina formulações analíticas clássicas, Método dos Momentos (MoM) e validação com simulações COMSOL Multiphysics.

## Domínio de Aplicação

| Categoria | Tipos Suportados |
|-----------|-----------------|
| **Linhas aéreas (OHTL)** | Mono e multifásica com retorno pelo solo (Carson, Sunde, Nakagawa, Quasi-TEM) |
| **Cabos subterrâneos (SCC)** | Monopolar core + sheath, configurações plana e trifólio |
| **Cabos em tubos (HDPE/pipe)** | SCC excêntrico em tubo HDPE, trefoil em conduto condutor |
| **Sistemas isolados** | Fios nus (bifilar), fios revestidos, coaxiais, ribbon cable |

---

## Estrutura do Projeto

```
PyLCP/
├── mtl_main/                   # Núcleo — modelagem de MTLs
│   ├── source.py               # MulticonductorTransmissionLine (classe central)
│   ├── strategy.py             # Padrão Strategy por tipo de MTL
│   ├── propagation.py          # γ_v, γ_i, Zc, Yc no domínio de fase (compartilhado)
│   └── graphics.py             # Esquemáticos de seção transversal
│
├── analytical_forms/           # Formulações analíticas
│   ├── single_core_cable.py    # Cabos monopolares (SCC)
│   ├── overhead_lines.py       # Linhas aéreas (OHTL)
│   ├── modal_analysis.py       # Decomposição modal (Cap. 5 Andreata): α_m, v_m, Z_cm
│   └── isolated_wires.py       # Fios em meio homogêneo
│
├── mom/                        # Método dos Momentos (MoM)
│   ├── bare_wire_systems.py    # Fios nus — colocação e Galerkin
│   └── coated_wire_systems.py  # Fios revestidos
│
├── mom_so/                     # MoM com funções de Green quasi-estáticas
│   ├── quasi_static_green.py   # Matriz G quasi-estática (MoM-SO)
│   └── lossless_medium.py      # Admitância em meio sem perdas
│
├── models/                     # Geradores de modelos paramétricos
│   ├── model_generator.py      # Classe base BaseModelGenerator
│   ├── single_core_cable.py    # SingleCoreCableModelGenerator
│   ├── overhead_lines.py       # Modelos de linhas aéreas
│   ├── isolated_wires.py       # Fios isolados em arranjos planos
│   ├── hdpe.py                 # Cabos em tubo HDPE
│   └── pipe_type.py            # Cabos em conduto condutor
│
├── plotter/                    # Visualização
│   ├── models_base.py          # BasePlotter — plotagem config-driven
│   ├── scc_models.py           # SingleCoreCableModels (cabos subterrâneos)
│   ├── scc_plotter.py          # SCCPlotter(BasePlotter)
│   ├── ohtl_models.py          # OverheadLineModels (linhas aéreas)
│   ├── ohtl_plotter.py         # OHTLPlotter(BasePlotter)
│   ├── lima_models.py          # LimaModels — constante de propagação OHTL
│   ├── deConti_models.py       # DeContiModels — parâmetros matriciais
│   ├── modal_plotter.py        # ModalPropagationPlotter — Figs 5.5/5.6/5.7 (Andreata)
│   └── xue_models.py           # XueModels — impedância série Xue
│
├── utils/
│   ├── case_utils.py           # Utilitários: load_json_parameters, save_figure, ...
│   ├── passivity_check.py      # Avaliação de passividade de Z'/Y'/Yc (Gustavsen 2008, eq. 3)
│   └── comsol_data.py          # Leitura e parsing de dados COMSOL
│
├── mtl_paul/                   # Integração Fortran (RIBBON.FOR)
│   └── py_fortran.py           # FortranRunner — wrapper Python/Fortran
│
├── tulip/                      # Integração com malhas Gmsh (opcional)
│   └── py_tulip.py             # PyTulip — geração de malhas para FEM
│
├── testData/                   # Casos de teste (20 cenários)
├── run_testData.py             # Runner principal — executa todos os casos
└── environment.yml             # Ambiente conda
```

---

## Instalação

### Pré-requisitos

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) ou Anaconda
- Python 3.11

### Criar Ambiente

```bash
conda env create -f environment.yml
conda activate pylcp
```

### Dependências Principais

| Pacote | Versão | Uso |
|--------|--------|-----|
| numpy | 1.24.3 | Álgebra linear, arrays |
| scipy | 1.15.3 | Funções de Bessel, integrais |
| matplotlib | 3.10.3 | Gráficos |
| pandas | 2.2.3 | Processamento de dados COMSOL |
| sympy | 1.14.0 | Derivação simbólica |
| mpmath | 1.3.0 | Aritmética de precisão arbitrária |
| scikit-rf | 1.7.0 | Parâmetros de rede (S, Z, Y) |
| gmsh | 4.13.1 | Geração de malhas FEM *(opcional)* |
| meshpy | 2022.1.3 | Interface Python para Gmsh *(opcional)* |

> **Nota:** `gmsh` é necessário apenas para o caso `ribbon_s50` com integração Gmsh/TULIP. Os demais casos funcionam sem ele.

---

## Execução

### Todos os casos de teste

```bash
cd c:\git\PyLCP
python run_testData.py
```

### Caso individual

```bash
# Módulo como pacote Python
python -m testData.scc_132kV_xue.scc_132kV_xue
python -m testData.ohtl_single_lima.ohtl_single_lima
```

Os resultados (gráficos `.png` e esquemáticos) são salvos automaticamente em `testData/<case>/Results/`.

---

## Uso Programático

### Parâmetros PUL de um cabo monopolar

```python
import numpy as np
from mtl_main.source import MulticonductorTransmissionLine
from models.single_core_cable import SingleCoreCableModelGenerator
from analytical_forms.single_core_cable import InternalPerUnitParameters, PerUnitParameters

# 1. Carregar modelo a partir de JSON
gen = SingleCoreCableModelGenerator(__file__)
model = gen.underground_model()

# 2. Instanciar MTL
mtl = MulticonductorTransmissionLine(model)

# 3. Calcular parâmetros internos (skin effect)
f = np.logspace(0, 6, 61)            # 1 Hz a 1 MHz
internal = InternalPerUnitParameters(mtl, f)
zi = internal.approximations()['Zi_approx']

# 4. Calcular matrizes PUL completas
pul = PerUnitParameters(mtl, f)
matrices = pul.pul_matrices(zi)
# matrices['series_impedance_matrix']   → (61, N, N) Ω/m
# matrices['shunt_admittance_matrix']   → (61, N, N) S/m
# matrices['propagation_voltage_matrix'] → (61, N, N) 1/m
```

### Parâmetros PUL de uma linha aérea

```python
import numpy as np
from mtl_main.source import MulticonductorTransmissionLine
from models.overhead_lines import single_phase_model
from analytical_forms.overhead_lines import InternalPerUnitParameters, PerUnitParameters
from utils.case_utils import load_json_parameters

# 1. Carregar JSON com geometria
params = load_json_parameters(__file__)
model  = single_phase_model(params)
mtl    = MulticonductorTransmissionLine(model)

# 2. Calcular para múltiplos cenários de solo
f   = np.logspace(0, 6, 100)
pul = PerUnitParameters(mtl, f)
internal = InternalPerUnitParameters(mtl, f)
zi  = internal.approximations()['Zi_approx']

for zg_form in ['nakagawa', 'carson', 'sunde', 'quasi_tem']:
    data = pul.pul_matrices(zi, zg_form=zg_form)
    # data['series_impedance_matrix'], data['earth-return_impedance_matrix'], ...
```

---

## Casos de Teste

### Fios e Cabos Isolados

| Caso | Tipo | Descrição | Referência |
|------|------|-----------|-----------|
| `bare_wire` | `bare_wires` | Linha bifilar nua, Ø10 mm, espaçamento 21 mm | Clements (1974) |
| `bifilar_s21` | `bare_wires` | Idem, análise detalhada frequência-dependente | Clements (1974) |
| `bifilar_s25` | `bare_wires` | Linha bifilar nua, espaçamento 25 mm | Clements (1974) |
| `bifilar_s100` | `bare_wires` | Linha bifilar nua, espaçamento 100 mm | Clements (1974) |
| `coated_bifilar_s40` | `coated_wires` | Fios com revestimento isolante, esp. 40 mm | Paul (2008) |
| `coaxial` | `coaxial` | Cabo coaxial isolado em ar | Patel (2014) |
| `ribbon_s50` | `coated_wires` | Ribbon cable 3 condutores, espaçamento 50 mil | Paul (2008) |

### Linhas Aéreas (OHTL)

| Caso | Tipo | Descrição | Referência |
|------|------|-----------|-----------|
| `ohtl_single_deConti` | `overhead` | Linha monofásica, formulação De Conti | De Conti et al. |
| `ohtl_single_deConti_deri` | `overhead` | Idem com método da derivada (Deri) | Deri et al. |
| `ohtl_single_lima` | `overhead` | Propagação Quasi-TEM, Carson, Sunde, Nakagawa | Lima (2015) |
| `ohtl_single_xue` | `overhead` | Impedância série, variação de ρ_solo e ε_r | Xue (2018) |

### Cabos Subterrâneos (SCC)

| Caso | Tipo | Descrição | Referência |
|------|------|-----------|-----------|
| `scc_single_deConti` | `scc` | Monopolar simples, retorno pelo solo | De Conti |
| `scc_flat_deConti` | `scc` | 3 cabos configuração plana | De Conti |
| `scc_132kV_xue` | `scc` | Monopolar 132 kV, análise completa | Xue (2018) |
| `scc_flat_xue` | `scc` | 3 cabos configuração plana | Xue (2018) |
| `scc_trefoil_xue` | `scc` | 3 cabos configuração trifólio | Xue (2018) |
| `scc_138kV_prysmian` | `scc` | Cabo real 138 kV Prysmian | Prysmian / Patel (2014) |

### Estruturas Especiais

| Caso | Tipo | Descrição | Referência |
|------|------|-----------|-----------|
| `hdpe_220kV_2000mm2` | `hdpe` | Monopolar 220 kV, 2000 mm², em tubo HDPE 10" | Lafaia et al. (2015) |
| `hdpe_225kV_2000mm2` | `hdpe` | Idem 225 kV, com validação COMSOL | Lafaia et al. (2015) |
| `pipe_trefoil_patel` | `pipe` | 3 cabos trifásicos em conduto condutor | Patel (2014) |

---

## Arquitetura

### Fluxo de Dados

```
case.json  ──►  ModelGenerator  ──►  model dict
                                          │
                                          ▼
                              MulticonductorTransmissionLine
                              (Strategy: SCC / OHTL / Cable)
                                          │
                              ┌───────────┼───────────┐
                              ▼           ▼           ▼
                        Analytical       MoM       MoM-SO
                         (Z, Y, γ)    (C, L)     (Z quasi)
                              │           │           │
                              └───────────┼───────────┘
                                          ▼
                                      pul_data
                                    {'freq': ...,
                                     scenario_key: {
                                       'series_impedance_matrix': ...,
                                       'shunt_admittance_matrix': ...,
                                       ...}}
                                          │
                                          ▼
                              Plotter (BasePlotter subclass)
                                          │
                                          ▼
                                   Results/*.png
```

### Estratégias por Tipo de MTL

| Estratégia | Tipo (`mtl_type`) | Geometrias |
|-----------|-------------------|-----------|
| `CableStrategy` | `bare_wires`, `coated_wires`, `coaxial`, `pipe` | Sistemas isolados |
| `SingleCoreCableStrategy` | `scc` | Core + sheath enterrados |
| `SingleCoreCableInHDPEStrategy` | `hdpe` | SCC dentro de tubo HDPE |
| `SingleCoreCableWithECCInHDPEStrategy` | `shared-hdpe` | SCC + ECC em HDPE |
| `OverheadLineStrategy` | `overhead` | Fios aéreos com retorno pelo solo |

### Formato do Modelo (dict)

```python
model = {
    'type': 'scc',                   # tipo de MTL
    'idx_ref_conductor': 0,          # índice do condutor de referência (solo/retorno)
    1: {                             # condutor 1 (core)
        'line_id': 1,
        'conductor_name': 'core',
        'center_point': (0.0, 1.0),  # [m]
        'radius': (0.0, 0.02785),    # [r_in, r_out] em metros
        'conductivity': 5.8e7,       # [S/m]
        'insulation': {
            'name': 'primary_insulation',
            'thickness': 0.015,      # [m]
            'relative_permittivity': 2.5,
        },
        'fourier_order': 10,
    },
    2: { ... },                      # sheath
}
```

### Chaves retornadas por `pul_matrices()`

| Chave | Descrição | Unidade |
|-------|-----------|--------|
| `series_impedance_matrix` | Impedância série Z(f) | Ω/m |
| `shunt_admittance_matrix` | Admitância shunt Y(f) | S/m |
| `earth-return_impedance_matrix` | Contribuição do solo | Ω/m |
| `propagation_voltage_matrix` | Constante de propagação γ(f) | 1/m |
| `propagation_current_matrix` | γ em corrente | 1/m |
| `characteristic_impedance_matrix` | Impedância característica Zc | Ω |
| `characteristic_admittance_matrix` | Admitância característica Yc | S |

> Para SCC, essas quatro chaves de propagação são obtidas por
> `PerUnitParameters.propagation_matrices(quasi_tem_matrices)` (anexadas ao
> dict de `quasi_tem_approx_matrices`).

---

## Análise Modal (Cap. 5 de Andreata)

`analytical_forms/modal_analysis.py::ModalDecomposition` decompõe `Y'Z'` por
frequência, rastreia os modos ao longo da frequência (*switching-back
procedure* — Gustavsen 2008 §IV-A / Wedepohl 1996 §6) e devolve, por modo:
constante de atenuação `α_m`, constante de fase `β_m`, velocidade de fase
`v_m`, impedância/admitância características modais `Z_cm`/`Y_cm`, e as matrizes
de transformação `T_I`/`T_V`. Para a Configuração 1 (3 SCC enterrados, 6 modos)
os modos são rotulados automaticamente (`ground`, `inter_sheath_1/2`,
`coaxial_1/2/3`). `plotter/modal_plotter.py::ModalPropagationPlotter` gera as
Figs 5.5 (`α_m`), 5.6 (`v_m`) e 5.7 (`|Z_cm|`).

Detalhes de projeto, validação e limitações:
[`testData/andreata_common/DESENVOLVIMENTO_MODAL_CAP5.md`](testData/andreata_common/DESENVOLVIMENTO_MODAL_CAP5.md).

---

## Formulações de Retorno pelo Solo (OHTL)

| `zg_form` | Método | Referência |
|-----------|--------|-----------|
| `'carson'` | Integral de Carson (série infinita) | Carson (1926) |
| `'sunde'` | Integral de Sunde | Sunde (1968) |
| `'nakagawa'` | Formulação Nakagawa/Wise | Nakagawa (1981) |
| `'quasi_tem'` | Quasi-TEM (integral exata) | Lima & Paulino (2009) |
| `'quasi_tem_log'` | Quasi-TEM (aproximação logarítmica) | Lima & Paulino (2009) |
| `'deri'` | Aproximação de Deri | Deri et al. (1981) |

---

## Fluxo de Trabalho Git

O desenvolvimento contínuo ocorre na branch `ipst_2027`. A `main` recebe atualizações periódicas via merge.

### Branches

| Branch | Papel |
|--------|-------|
| `main` | Versão estável — recebe merges de `ipst_2027` |
| `ipst_2027` | Desenvolvimento ativo — commits do dia a dia |

### Sincronizar `ipst_2027` com `main` (antes de iniciar trabalho novo)

```bash
git checkout ipst_2027
git merge main
git push origin ipst_2027
```

### Promover trabalho de `ipst_2027` para `main`

```bash
git checkout main
git merge ipst_2027
git push origin main
git checkout ipst_2027
```

---

## Referências Bibliográficas

1. **Paul, C. R.** (2008). *Analysis of Multiconductor Transmission Lines*, 2nd ed. Wiley-IEEE Press.
2. **Ametani, A., Nagaoka, N., Baba, Y., Ohno, T.** (2015). *Cable System Transients*. Wiley-IEEE Press.
3. **Patel, U. R., Triverio, P., Morevev, A.** (2013-2014). Surface Admittance Approach for modeling MoM-SO internal impedance of coated conductors. *IEEE TEMC*.
4. **Xue, H.** (2018). *General Formulation for Frequency-Dependent Parameters of Underground and Overhead Transmission Lines*. Ph.D. Thesis, École Polytechnique de Montréal.
5. **Lima, A. C. S., Paulino, J. O. S.** (2009). Quasi-TEM approach for the evaluation of transmission line parameters. *IEEE TPWRD*.
6. **De Conti, A., Visacro, S.** Revisiting Carson's formulas. *IEEE TPWRD*.
7. **Carson, J. R.** (1926). Wave propagation in overhead wires with ground return. *Bell System Technical Journal*, 5(4), 539–554.
8. **Lafaia, I., et al.** (2015). Eccentric cable inside HDPE pipe modelling. *IEEE TPWRD*.
9. **Nakagawa, M.** (1981). Further studies on wave propagation along overhead transmission lines. *IEEE TPAS*, 100(7).
10. **Sunde, E. D.** (1968). *Earth Conduction Effects in Transmission Systems*. Dover Publications.
