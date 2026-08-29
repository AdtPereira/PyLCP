# scc_flat_andreata

## Purpose

Simulation and validation script for a system of **three single-core coaxial cables (SCC)** buried in a flat arrangement. It computes the per-unit-length (PUL) parameters — series impedance and shunt admittance — for multiple soil scenarios and analytical formulations, and generates the reference plots comparing against COMSOL results (when available).

The case reproduces **Figures 4.19, 4.21a/b and 4.23, 4.25a/b** of Xue's thesis (2018).

---

## System Configuration (JSON)

Parameters loaded from `scc_flat_andreata.json`:

| Parameter | Value |
|---|---|
| Arrangement | Flat, 3 cables |
| Burial depth | 1.2 m |
| Spacing between cables | 0.2 m |
| Soil — conductivity | 0.01 S/m (rho = 100 Ohm.m) |
| Soil — relative permittivity | 1.0 |
| Core — outer radius | 10.325 mm — conductivity: 38 MS/m |
| Core insulation | XLPE, 10.675 mm, eps_r = 2.2 |
| Sheath — inner/outer radii | 21.0 / 21.8 mm — conductivity: 60 MS/m |
| Sheath insulation | PVC, 2.2 mm, eps_r = 2.8 |
| Fourier order | 10 |

The generated model contains **7 conductors**: 1 soil return (`line_id=0`) + 3 core/sheath pairs (`line_id=1..6`) placed at (0,0), (0.2, 0) and (0.4, 0) m, all at depth -1.2 m.

---

## Structure of the `main()` Function

### 1. Construction of the MTL models (lines 27–36)

Three model variants are created via `SingleCoreCableModelGenerator` + `MulticonductorTransmissionLine`, representing the soil scenarios:

| Variable | Soil |
|---|---|
| `mtl_model_a` | rho = 100 Ohm.m, eps_r = 1 (reference) |
| `mtl_model_b` | rho = 100 Ohm.m, eps_r = 20 |
| `mtl_model_c` | rho = 500 Ohm.m, eps_r = 1 |

### 2. The `pul_data` dictionary (lines 38–94)

Central structure that aggregates all of the simulation data:

```
pul_data
├── 'frequencies'       -> 121 log-spaced points in [10^4, 10^7] Hz
├── 'internal_matrices' -> internal matrices (filled in step 4)
├── 'comsol'            -> COMSOL data (filled if the file is available)
│   └── 'scenarios'     -> 3 COMSOL scenarios (rho_g_100/500, epsr1_1/20)
└── 'scenarios'         -> 9 analytical scenarios (filled in step 5)
```

**9 analytical scenarios** — combination of 3 soils x 3 formulations:

| Key | Soil | Zg formulation | Yg formulation |
|---|---|---|---|
| `p100_er1` | 100 Ohm.m, eps_r=1 | Magalhaes/Xue | Magalhaes/Xue |
| `p100_er20` | 100 Ohm.m, eps_r=20 | Magalhaes/Xue | Magalhaes/Xue |
| `p500_er1` | 500 Ohm.m, eps_r=1 | Magalhaes/Xue | Magalhaes/Xue |
| `p100_er1_vance` | 100 Ohm.m, eps_r=1 | De Conti | Vance (1978) |
| `p100_er20_vance` | 100 Ohm.m, eps_r=20 | De Conti | Vance (1978) |
| `p500_er1_vance` | 500 Ohm.m, eps_r=1 | De Conti | Vance (1978) |
| `p100_er1_deconti` | 100 Ohm.m, eps_r=1 | De Conti | De Conti (2023) |
| `p100_er20_deconti` | 100 Ohm.m, eps_r=20 | De Conti | De Conti (2023) |
| `p500_er1_deconti` | 500 Ohm.m, eps_r=1 | De Conti | De Conti (2023) |

### 3. COMSOL processing — optional (lines 96–114)

Attempts to load `Results/cmsl_ground_return_impedance.txt` via `ComsolPostProcessor`. If the file exists, it builds the internal matrices at the COMSOL frequencies and computes the ground-return parameters for each of the 3 COMSOL scenarios. If it does not exist, the block is skipped with a warning and `pul_data['comsol']` is emptied.

### 4. Analytical internal parameters (lines 116–119)

```python
pul = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
internal_matrices = pul.matrices(internal_form='hybrid')
```

Computes, for the 121 frequencies, the internal cable matrices (core impedance, sheath impedance and couplings) using the hybrid formulation. The result is stored in `pul_data['internal_matrices']`.

### 5. PUL parameters per analytical scenario (lines 121–128)

For each of the 9 scenarios, via `PerUnitParameters`:

1. `earth_return_parameters(zg_form, yg_form)` — computes the ground-return impedance and admittance **Zg** and **Yg** using the specified formulation.
2. `quasi_tem_approx_matrices(internal_matrices, earth_return)` — assembles the full quasi-TEM matrices: series impedance **Zs = Zi + Zg** and shunt admittance **Ysh**.

Both results are stored in the scenario's dictionary.

### 6. Plot generation (lines 131–138)

`SCCPlotter` generates and automatically saves the following plots to `Results/`:

| Method | Figures | Content |
|---|---|---|
| `scc_earth_propagation_constant()` | `earth_propagation_constant.png` | Soil propagation constant gamma_1 (alpha and beta) |
| `scc_earth_return_impedance_matrix()` | `earth_return_impedance_self/mutual_ab/ac.png` | Zg — self-impedance and mutual terms AB, AC |
| `scc_earth_return_admittance_matrix()` | `earth_return_admittance_self/mutual_ab/ac.png` | Yg — self-admittance and mutual terms AB, AC |
| `scc_series_impedance_matrix(...)` | `fig419`, `fig421a`, `fig421b` | Zs — sheath self-impedance and mutual terms |
| `scc_shunt_admittance_matrix(...)` | `fig423`, `fig425a`, `fig425b` | Ysh — sheath self-admittance and mutual terms |
| `GroundReturnMTLRepresentation(...)` | `system_schematic.png` | Circuit schematic of the MTL system |

---

## Dependencies

| Module | Role |
|---|---|
| `models.single_core_cable.SingleCoreCableModelGenerator` | Builds the MTL model from the JSON |
| `mtl_main.source.MulticonductorTransmissionLine` | Represents the multiconductor MTL system |
| `analytical_forms.single_core_cable.InternalPerUnitParameters` | Internal cable matrices (core + sheath) |
| `analytical_forms.single_core_cable.PerUnitParameters` | Full PUL matrices with ground return |
| `utils.comsol_data.ComsolPostProcessor` | Reading and processing of the COMSOL data |
| `plotter.scc_plotter.SCCPlotter` | Generation and saving of the plots |
| `mtl_main.graphics.GroundReturnMTLRepresentation` | MTL circuit schematic |
| `plot_config.PLOT_CONFIG` | Layout, series and limits of every plot |

---

## Running

```bash
python -m testData.scc_flat_andreata.scc_flat_andreata
```

**Expected output:** 11 `.png` files saved to `testData/scc_flat_andreata/Results/`. The COMSOL data is optional — the analytical simulation runs normally without it.

---

## References

- **Xue, H.** (2018). *General Formulation and Accurate Evaluation of Earth-Return Parameters for Overhead/Underground Cables*. PhD thesis, École Polytechnique de Montréal.
- **Vance, E. F.** (1978). *Coupling to Shielded Cables*. Wiley-Interscience.
- **De Conti, A. et al.** (2023). Formulation for ground-return parameters of underground cables.
