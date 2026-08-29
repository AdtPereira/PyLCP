# PyLCP — Python Library for Cable Parameters

Electromagnetic analysis framework for computing per-unit-length (PUL) electrical parameters of **multiconductor transmission lines** (MTLs). It combines classical analytical formulations, the Method of Moments (MoM) and validation against COMSOL Multiphysics simulations.

## Application Domain

| Category | Supported Types |
|-----------|-----------------|
| **Overhead lines (OHTL)** | Single- and multi-phase with ground return (Carson, Sunde, Nakagawa, Quasi-TEM) |
| **Underground cables (SCC)** | Single-core core + sheath, flat and trefoil arrangements |
| **Cables in pipes (HDPE/pipe)** | Eccentric SCC in an HDPE pipe, trefoil in a conducting conduit |
| **Isolated systems** | Bare wires (bifilar), coated wires, coaxial, ribbon cable |

---

## Project Structure

```
PyLCP/
├── mtl_main/                   # Core — MTL modeling
│   ├── source.py               # MulticonductorTransmissionLine (central class)
│   ├── strategy.py             # Strategy pattern per MTL type
│   ├── propagation.py          # gamma_v, gamma_i, Zc, Yc in the phase domain (shared)
│   └── graphics.py             # Cross-section schematics
│
├── analytical_forms/           # Analytical formulations
│   ├── single_core_cable.py    # Single-core cables (SCC); InternalParametersFromFEM, build_pul_matrices
│   ├── overhead_lines.py       # Overhead lines (OHTL)
│   ├── modal_analysis.py       # Modal decomposition (Andreata Ch. 5): alpha_m, v_m, Z_cm
│   └── isolated_wires.py       # Wires in a homogeneous medium
│
├── mom/                        # Method of Moments (MoM)
│   ├── bare_wire_systems.py    # Bare wires — collocation and Galerkin
│   └── coated_wire_systems.py  # Coated wires
│
├── mom_so/                     # MoM with quasi-static Green's functions
│   ├── quasi_static_green.py   # Quasi-static G matrix (MoM-SO)
│   └── lossless_medium.py      # Admittance in a lossless medium
│
├── models/                     # Parametric model generators
│   ├── model_generator.py      # Base class BaseModelGenerator
│   ├── single_core_cable.py    # SingleCoreCableModelGenerator
│   ├── overhead_lines.py       # Overhead-line models
│   ├── isolated_wires.py       # Isolated wires in flat arrangements
│   ├── hdpe.py                 # Cables in an HDPE pipe
│   └── pipe_type.py            # Cables in a conducting conduit
│
├── plotter/                    # Visualization
│   ├── models_base.py          # BasePlotter — config-driven plotting
│   ├── scc_models.py           # SingleCoreCableModels (underground cables)
│   ├── scc_plotter.py          # SCCPlotter(BasePlotter)
│   ├── ohtl_models.py          # OverheadLineModels (overhead lines)
│   ├── ohtl_plotter.py         # OHTLPlotter(BasePlotter)
│   ├── lima_models.py          # LimaModels — OHTL propagation constant
│   ├── deConti_models.py       # DeContiModels — matrix parameters
│   ├── modal_plotter.py        # ModalPropagationPlotter — Figs 5.5/5.6/5.7 (Andreata)
│   └── xue_models.py           # XueModels — Xue series impedance
│
├── utils/
│   ├── case_utils.py           # Utilities: load_json_parameters, save_figure, ...
│   ├── passivity_check.py      # Passivity assessment of Z'/Y'/Yc (Gustavsen 2008, eq. 3)
│   └── comsol_data.py          # Reading and parsing of COMSOL data
│
├── mtl_paul/                   # Fortran integration (RIBBON.FOR)
│   └── py_fortran.py           # FortranRunner — Python/Fortran wrapper
│
├── tulip/                      # Gmsh mesh integration (optional)
│   └── py_tulip.py             # PyTulip — mesh generation for FEM
│
├── testData/                   # Test cases (20 scenarios)
├── run_testData.py             # Main runner — runs every case
└── environment.yml             # conda environment
```

---

## Installation

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- Python 3.11

### Create Environment

```bash
conda env create -f environment.yml
conda activate pylcp
```

### Main Dependencies

| Package | Version | Use |
|--------|--------|-----|
| numpy | 1.24.3 | Linear algebra, arrays |
| scipy | 1.15.3 | Bessel functions, integrals |
| matplotlib | 3.10.3 | Plots |
| pandas | 2.2.3 | COMSOL data processing |
| sympy | 1.14.0 | Symbolic derivation |
| mpmath | 1.3.0 | Arbitrary-precision arithmetic |
| scikit-rf | 1.7.0 | Network parameters (S, Z, Y) |
| gmsh | 4.13.1 | FEM mesh generation *(optional)* |
| meshpy | 2022.1.3 | Python interface to Gmsh *(optional)* |

> **Note:** `gmsh` is required only for the `ribbon_s50` case with Gmsh/TULIP integration. The other cases work without it.

---

## Running

### All test cases

```bash
cd c:\git\PyLCP
python run_testData.py
```

### Individual case

```bash
# Module as a Python package
python -m testData.scc_132kV_xue.scc_132kV_xue
python -m testData.ohtl_single_lima.ohtl_single_lima
```

The results (`.png` plots and schematics) are saved automatically to `testData/<case>/Results/`.

---

## Programmatic Use

### PUL parameters of a single-core cable

```python
import numpy as np
from mtl_main.source import MulticonductorTransmissionLine
from models.single_core_cable import SingleCoreCableModelGenerator
from analytical_forms.single_core_cable import InternalPerUnitParameters, PerUnitParameters

# 1. Load model from JSON
gen = SingleCoreCableModelGenerator(__file__)
model = gen.underground_model()

# 2. Instantiate the MTL
mtl = MulticonductorTransmissionLine(model)

# 3. Compute the internal parameters (skin effect)
f = np.logspace(0, 6, 61)            # 1 Hz to 1 MHz
internal = InternalPerUnitParameters(mtl, f)
zi = internal.approximations()['Zi_approx']

# 4. Compute the full PUL matrices
pul = PerUnitParameters(mtl, f)
matrices = pul.pul_matrices(zi)
# matrices['series_impedance_matrix']   -> (61, N, N) Ohm/m
# matrices['shunt_admittance_matrix']   -> (61, N, N) S/m
# matrices['propagation_voltage_matrix'] -> (61, N, N) 1/m
```

### PUL parameters of an overhead line

```python
import numpy as np
from mtl_main.source import MulticonductorTransmissionLine
from models.overhead_lines import single_phase_model
from analytical_forms.overhead_lines import InternalPerUnitParameters, PerUnitParameters
from utils.case_utils import load_json_parameters

# 1. Load JSON with the geometry
params = load_json_parameters(__file__)
model  = single_phase_model(params)
mtl    = MulticonductorTransmissionLine(model)

# 2. Compute for multiple soil scenarios
f   = np.logspace(0, 6, 100)
pul = PerUnitParameters(mtl, f)
internal = InternalPerUnitParameters(mtl, f)
zi  = internal.approximations()['Zi_approx']

for zg_form in ['nakagawa', 'carson', 'sunde', 'quasi_tem']:
    data = pul.pul_matrices(zi, zg_form=zg_form)
    # data['series_impedance_matrix'], data['earth-return_impedance_matrix'], ...
```

---

## Test Cases

### Isolated Wires and Cables

| Case | Type | Description | Reference |
|------|------|-----------|-----------|
| `bare_wire` | `bare_wires` | Bare bifilar line, 10 mm dia., 21 mm spacing | Clements (1974) |
| `bifilar_s21` | `bare_wires` | Same, detailed frequency-dependent analysis | Clements (1974) |
| `bifilar_s25` | `bare_wires` | Bare bifilar line, 25 mm spacing | Clements (1974) |
| `bifilar_s100` | `bare_wires` | Bare bifilar line, 100 mm spacing | Clements (1974) |
| `coated_bifilar_s40` | `coated_wires` | Wires with an insulating coating, 40 mm spacing | Paul (2008) |
| `coaxial` | `coaxial` | Coaxial cable isolated in air | Patel (2014) |
| `ribbon_s50` | `coated_wires` | 3-conductor ribbon cable, 50 mil spacing | Paul (2008) |

### Overhead Lines (OHTL)

| Case | Type | Description | Reference |
|------|------|-----------|-----------|
| `ohtl_single_deConti` | `overhead` | Single-phase line, De Conti formulation | De Conti et al. |
| `ohtl_single_deConti_deri` | `overhead` | Same with the derivative method (Deri) | Deri et al. |
| `ohtl_single_lima` | `overhead` | Quasi-TEM propagation, Carson, Sunde, Nakagawa | Lima (2015) |
| `ohtl_single_xue` | `overhead` | Series impedance, variation of rho_soil and eps_r | Xue (2018) |

### Underground Cables (SCC)

| Case | Type | Description | Reference |
|------|------|-----------|-----------|
| `scc_single_deConti` | `scc` | Simple single-core, ground return | De Conti |
| `scc_flat_deConti` | `scc` | 3 cables, flat arrangement | De Conti |
| `scc_132kV_xue` | `scc` | Single-core 132 kV, full analysis | Xue (2018) |
| `scc_flat_xue` | `scc` | 3 cables, flat arrangement | Xue (2018) |
| `scc_trefoil_xue` | `scc` | 3 cables, trefoil arrangement | Xue (2018) |
| `scc_138kV_prysmian` | `scc` | Real 138 kV Prysmian cable | Prysmian / Patel (2014) |

### Special Structures

| Case | Type | Description | Reference |
|------|------|-----------|-----------|
| `hdpe_220kV_2000mm2` | `hdpe` | Single-core 220 kV, 2000 mm^2, in a 10" HDPE pipe | Lafaia et al. (2015) |
| `hdpe_225kV_2000mm2` | `hdpe` | Same 225 kV, with COMSOL validation | Lafaia et al. (2015) |
| `pipe_trefoil_patel` | `pipe` | 3 three-phase cables in a conducting conduit | Patel (2014) |

---

## Architecture

### Data Flow

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
                         (Z, Y, gamma) (C, L)     (Z quasi)
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

### Strategies per MTL Type

| Strategy | Type (`mtl_type`) | Geometries |
|-----------|-------------------|-----------|
| `CableStrategy` | `bare_wires`, `coated_wires`, `coaxial`, `pipe` | Isolated systems |
| `SingleCoreCableStrategy` | `scc` | Buried core + sheath |
| `SingleCoreCableInHDPEStrategy` | `hdpe` | SCC inside an HDPE pipe |
| `SingleCoreCableWithECCInHDPEStrategy` | `shared-hdpe` | SCC + ECC in HDPE |
| `OverheadLineStrategy` | `overhead` | Overhead wires with ground return |

### Model Format (dict)

```python
model = {
    'type': 'scc',                   # MTL type
    'idx_ref_conductor': 0,          # index of the reference conductor (soil/return)
    1: {                             # conductor 1 (core)
        'line_id': 1,
        'conductor_name': 'core',
        'center_point': (0.0, 1.0),  # [m]
        'radius': (0.0, 0.02785),    # [r_in, r_out] in meters
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

### Keys returned by `pul_matrices()`

| Key | Description | Unit |
|-------|-----------|--------|
| `series_impedance_matrix` | Series impedance Z(f) | Ohm/m |
| `shunt_admittance_matrix` | Shunt admittance Y(f) | S/m |
| `earth-return_impedance_matrix` | Soil contribution | Ohm/m |
| `propagation_voltage_matrix` | Propagation constant gamma(f) | 1/m |
| `propagation_current_matrix` | gamma in current | 1/m |
| `characteristic_impedance_matrix` | Characteristic impedance Zc | Ohm |
| `characteristic_admittance_matrix` | Characteristic admittance Yc | S |

> For SCC, these four propagation keys are obtained via
> `PerUnitParameters.propagation_matrices(quasi_tem_matrices)` (appended to the
> `quasi_tem_approx_matrices` dict).

---

## Modal Analysis (Andreata Ch. 5)

`analytical_forms/modal_analysis.py::ModalDecomposition` decomposes `Y'Z'` per
frequency, tracks the modes across frequency (*switching-back procedure* —
Gustavsen 2008 sec. IV-A / Wedepohl 1996 sec. 6) and returns, per mode:
attenuation constant `alpha_m`, phase constant `beta_m`, phase velocity
`v_m`, modal characteristic impedance/admittance `Z_cm`/`Y_cm`, and the
transformation matrices `T_I`/`T_V`. For the 6-conductor systems (3 SCC, with or
without an HDPE duct) the modes are labeled automatically in 2 stages
(`ground`, `inter_sheath_1/2`, `coaxial_1/2/3`).
`plotter/modal_plotter.py::ModalPropagationPlotter` generates the plots of
`alpha_m`, `v_m` and `|Z_cm|`:

| Case | Andreata Figures |
|---|---|
| `andreata_case1` — 3 directly buried SCC (Config. 1) | 5.5 / 5.6 / 5.7 |
| `andreata_case2` — 3 SCC in individual HDPE ducts (Config. 2), **FEM-hybrid** pipeline | 5.8 / 5.9 / 5.10 |

---

## Hybrid pipeline (FEM/COMSOL internal + analytical ground return)

For geometries with no internal analytical solution (SCC inside an HDPE pipe —
Config. 2), `analytical_forms/single_core_cable.py::build_pul_matrices(...,
internal_source='fem', fem_internal=InternalParametersFromFEM(...))` assembles:

- **`Zi` / `Yi`** — from COMSOL (exact eccentric geometry, air + HDPE pipe);
- **`Zg` / `Yg`** — analytical, **closed-form** expressions of De Conti/Duarte/Alipio
  2023 (`zg_form='deconti'`; eqs. 4.59/4.63 of Andreata sec. 5.4 / 6.1), with the
  self term using the **pipe outer radius**
  (`mtl_main/strategy._cable_external_geometry` recognizes the `enclosure` field).

It reproduces `Z'`/`Y'` of Andreata's FEM reference to within **< 2 %**, without
Lafaia's equivalent-permittivity GMD method. Details:
[`testData/andreata_common/HYBRID_PIPELINE_PLAN.md`](testData/andreata_common/HYBRID_PIPELINE_PLAN.md).

Design, validation and limitations:
[`testData/andreata_common/MODAL_CH5_DEVELOPMENT.md`](testData/andreata_common/MODAL_CH5_DEVELOPMENT.md).

---

## Ground-Return Formulations (OHTL)

| `zg_form` | Method | Reference |
|-----------|--------|-----------|
| `'carson'` | Carson integral (infinite series) | Carson (1926) |
| `'sunde'` | Sunde integral | Sunde (1968) |
| `'nakagawa'` | Nakagawa/Wise formulation | Nakagawa (1981) |
| `'quasi_tem'` | Quasi-TEM (exact integral) | Lima & Paulino (2009) |
| `'quasi_tem_log'` | Quasi-TEM (logarithmic approximation) | Lima & Paulino (2009) |
| `'deri'` | Deri approximation | Deri et al. (1981) |

---

## Git Workflow

Ongoing development happens on the `ipst_2027` branch. `main` receives periodic updates via merge.

### Branches

| Branch | Role |
|--------|-------|
| `main` | Stable version — receives merges from `ipst_2027` |
| `ipst_2027` | Active development — day-to-day commits |

### Sync `ipst_2027` with `main` (before starting new work)

```bash
git checkout ipst_2027
git merge main
git push origin ipst_2027
```

### Promote work from `ipst_2027` to `main`

```bash
git checkout main
git merge ipst_2027
git push origin main
git checkout ipst_2027
```

---

## Bibliographic References

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
