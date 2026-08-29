# scc_34kV_andreata

## Purpose

Detailed simulation and analysis script for a **single buried single-core coaxial cable (SCC)**, focused on decomposing the PUL parameters into their individual contributions: internal, ground return and the resulting composition. It computes series impedance, shunt admittance and per-unit-length potential coefficients over a wide frequency range (1 Hz to 10 MHz), comparing three analytical formulations.

It differs from the `scc_flat_andreata` case by analyzing **a single isolated cable**, covering a wider frequency range, and presenting decomposition plots that make the relative contribution of each part of the parameters visible.

---

## System Configuration (JSON)

Parameters loaded from `scc_34kV_andreata.json`:

| Parameter | Value |
|---|---|
| Arrangement | Single, 1 cable |
| Burial depth | 1.2 m |
| Soil — conductivity | 0.01 S/m (rho = 100 Ohm.m) |
| Soil — relative permittivity | 1.0 |
| Core — outer radius | 10.325 mm — conductivity: 38 MS/m |
| Core insulation | XLPE, 10.675 mm, eps_r = 2.2 |
| Sheath — inner/outer radii | 21.0 / 21.8 mm — conductivity: 60 MS/m |
| Sheath insulation | PVC, 2.2 mm, eps_r = 2.8 |
| Fourier order | 10 |

The generated model contains **3 conductors**: 1 soil return (`line_id=0`) + 1 core/sheath pair (`line_id=1,2`), placed at (0, -1.2) m.

---

## Structure of the `main()` Function

### 1. Construction of the MTL model (lines 24–25)

A single model is created via `SingleCoreCableModelGenerator` + `MulticonductorTransmissionLine`:

| Variable | Soil |
|---|---|
| `mtl_model` | rho = 100 Ohm.m, eps_r = 1 (reference) |

### 2. The `pul_data` dictionary (lines 28–48)

Central structure that aggregates all of the simulation data:

```
pul_data
├── 'comsol'            -> None (no COMSOL data in this case)
├── 'frequencies'       -> 121 log-spaced points in [10^0, 10^7] Hz
├── 'internal_matrices' -> internal matrices (filled in step 3)
└── 'scenarios'         -> 3 analytical scenarios (filled in step 4)
```

**3 analytical scenarios** — same soil, different ground-return formulations:

| Key | Soil | Zg formulation | Yg formulation |
|---|---|---|---|
| `p100_xue` | 100 Ohm.m, eps_r=1 | Magalhaes/Xue (integral) | Magalhaes/Xue (integral) |
| `p100_deconti` | 100 Ohm.m, eps_r=1 | De Conti (2023) | De Conti (2023) |
| `p100_vance` | 100 Ohm.m, eps_r=1 | De Conti (2023) | Vance (1978) |

### 3. Analytical internal parameters (lines 50–53)

```python
pul = InternalPerUnitParameters(mtl_model, pul_data['frequencies'])
internal_matrices = pul.matrices()
```

Computes the internal cable matrices (impedance, admittance and potential coefficients of the core and sheath) for the 121 frequencies using the standard formulation. The result is stored in `pul_data['internal_matrices']`.

> **Difference from `scc_flat_andreata`:** the hybrid formulation (`internal_form='hybrid'`) is not used here; the `matrices()` method is called with the default parameters.

### 4. PUL parameters per analytical scenario (lines 55–63)

For each of the 3 scenarios, via `PerUnitParameters`:

1. `earth_return_parameters(zg_form, yg_form)` — computes the impedance **Zg**, potential coefficient **Pg** and ground-return admittance **Yg**.
2. `quasi_tem_approx_matrices(internal_matrices, earth_return)` — assembles the full quasi-TEM matrices: series impedance **Zs = Zi + Zg**, potential coefficient **Psh = Pi + Pg** and shunt admittance **Ysh**.

Both results are stored in each scenario's dictionary.

### 5. Plot generation (lines 66–85)

`SingleCoreCableModels` generates and automatically saves 22 figures to `Results/`, organized in three groups:

#### Group 1 — Potential Coefficients

| Method | Output file | Content |
|---|---|---|
| `potential_coefficients_composition('core_sheath')` | `potential_coefficients_composition_core_sheath.png` | Core-sheath mutual potential coeff.: internal + ground return + composition |
| `potential_coefficients_composition('core')` | `potential_coefficients_composition_core.png` | Core self potential coeff.: internal + ground return + composition |
| `potential_coefficients_composition('sheath')` | `potential_coefficients_composition_sheath.png` | Sheath self potential coeff.: internal + ground return + composition |
| `potential_coefficients_earth_return()` | `earth_return_potential.png` | Ground-return potential coefficient Pg |
| `potential_coefficients_internal()` | `internal_potential.png` | Internal potential coefficient Pi (P00, P01, P11) |

#### Group 2 — Series Impedance

| Method | Output file | Content |
|---|---|---|
| `series_impedance_composition('core_sheath')` | `series_impedance_composition_core_sheath.png` | Core-sheath mutual impedance: Zi + Zg + Zs + composition |
| `series_impedance_composition('core')` | `series_impedance_composition_core.png` | Core self impedance: Zi + Zg + Zs + composition |
| `series_impedance_composition('sheath')` | `series_impedance_composition_sheath.png` | Sheath self impedance: Zi + Zg + Zs + composition |
| `series_impedance_earth_return()` | `earth_return_impedance.png` | Ground-return impedance Zg (element [0,0]) |
| `series_impedance_internal()` | `internal_impedance.png` | Internal impedance Zi (core/sheath self and mutual) |
| `series_impedance_matrix()` | `series_impedance_matrix.png` | Full Zs matrix (core self, sheath self, mutual) |

#### Group 3 — Shunt Admittance

| Method | Output file | Content |
|---|---|---|
| `shunt_admittance_composition('core_sheath')` | `shunt_admittance_composition_core_sheath.png` | Core-sheath mutual admittance: Yi + Yg + Ysh + composition |
| `shunt_admittance_composition('core')` | `shunt_admittance_composition_core.png` | Core self admittance: Yi + Yg + Ysh + composition |
| `shunt_admittance_composition('sheath')` | `shunt_admittance_composition_sheath.png` | Sheath self admittance: Yi + Yg + Ysh + composition |
| `shunt_admittance_earth_return()` | `earth_return_admittance.png` | Ground-return admittance Yg (element [0,0]) |
| `shunt_admittance_internal()` | `internal_admittance.png` | Internal admittance Yi (core/sheath self and mutual) |
| `shunt_admittance_matrix()` | `shunt_admittance_matrix.png` | Full Ysh matrix (core self, sheath self, mutual) |

#### System schematic

| Method | Output file | Content |
|---|---|---|
| `GroundReturnMTLRepresentation(...)` | `system_schematic.png` | Circuit schematic of the MTL system (units in cm) |

---

## Decomposition Logic

The composition plots overlay four curves for each matrix element:

| Curve | Color | Meaning |
|---|---|---|
| **Internal** | Blue, dashed | Internal contribution of the cable (Zi or Yi or Pi) |
| **Earth-return** | Dark green, dashed | Ground-return contribution (Zg or Yg or Pg) |
| **Series** (quasi-TEM) | Black, solid | Full quasi-TEM result (Zs or Ysh or Psh) |
| **Composition** | Red, dotted | Direct sum internal + ground return (Zi+Zg or 1/(1/Yi+1/Yg)) |

This overlay lets you visually check in which frequency range each contribution dominates.

---

## Differences from `scc_flat_andreata`

| Aspect | `scc_34kV_andreata` | `scc_flat_andreata` |
|---|---|---|
| Arrangement | 1 single cable | 3 cables in a flat arrangement |
| MTL conductors | 3 (soil + core + sheath) | 7 (soil + 3 cores + 3 sheaths) |
| Frequency range | 1 Hz to 10 MHz | 10 kHz to 10 MHz |
| Number of scenarios | 3 | 9 |
| COMSOL data | Not used | Optional (when available) |
| Internal formulation | Standard | Hybrid (`hybrid`) |
| Plot focus | Decomposition per part | Comparison between soil scenarios |
| Plotter | `SingleCoreCableModels` | `SCCPlotter` |

---

## Dependencies

| Module | Role |
|---|---|
| `models.single_core_cable.SingleCoreCableModelGenerator` | Builds the MTL model from the JSON |
| `mtl_main.source.MulticonductorTransmissionLine` | Represents the multiconductor MTL system |
| `analytical_forms.single_core_cable.InternalPerUnitParameters` | Internal cable matrices (core + sheath) |
| `analytical_forms.single_core_cable.PerUnitParameters` | Full PUL matrices with ground return |
| `plotter.scc_models.SingleCoreCableModels` | Generation and saving of the decomposition plots |
| `mtl_main.graphics.GroundReturnMTLRepresentation` | MTL circuit schematic |

---

## Running

```bash
python -m testData.scc_34kV_andreata.scc_34kV_andreata
```

**Expected output:** 22 `.png` files saved to `testData/scc_34kV_andreata/Results/`.

---

## References

- **Xue, H.** (2018). *General Formulation and Accurate Evaluation of Earth-Return Parameters for Overhead/Underground Cables*. PhD thesis, École Polytechnique de Montréal.
- **Vance, E. F.** (1978). *Coupling to Shielded Cables*. Wiley-Interscience.
- **De Conti, A. et al.** (2023). Formulation for ground-return parameters of underground cables.
- **Ametani, A., Ohno, T. & Nagaoka, N.** (2015). *Cable System Transients: Theory, Modeling and Simulation*. Wiley-IEEE Press.
