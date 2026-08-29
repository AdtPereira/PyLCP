# andreata_case2

## Purpose

Simulation and validation script for **Configuration 2** (Figure 5.2) of the reference:
three single-core coaxial cables (SCC — core + sheath) in a flat arrangement,
**each installed inside its own HDPE duct**, buried in the soil. It computes
the per-unit-length (PUL) parameters — internal (core/sheath) and ground return —
comparing three distinct analytical ways of representing the dielectric effect
of the duct, and cross-validates against the `andreata_case1` MATLAB data.

Unlike `andreata_case1` (Configuration 1, directly buried cables) and
`andreata_case3` (Configuration 3, SCC + ECC without a duct), this case deals with a
structural limitation of pyLCP: **there is no analytical ground-return solution that
models the duct** — the `Zg`/`Yg` formulation sees only the position and outer radius of
each cable, never the presence of an HDPE duct around it. The actual validation of the
duct's effect on the ground return depends entirely on duct-specific COMSOL/MATLAB data
that **does not yet exist** for this case (see [Limitations and Pending
Data](#limitations-and-pending-data)).

Reference: Andreata, Luis Eduardo Batista. *Análise das Características de Propagação
e de Transitórios Eletromagnéticos em Cabos Subterrâneos Instalados em Tubos Não
Metálicos no Contexto de Parques Eólicos.* Graduate Program in Electrical
Engineering, UFMG, 2025. https://hdl.handle.net/1843/2086

---

## System Configuration (JSON)

Parameters loaded from `andreata_case2.json`:

| Parameter | Value |
|---|---|
| Arrangement | Flat, 3 cables |
| Burial depth | 1.2 m (cable center) |
| Spacing between cables | 0.2 m (center to center) |
| Soil — conductivity | 0.01 S/m (rho = 100 Ohm.m) |
| Soil — relative permittivity | 1.0 |
| Core — outer radius | 10.325 mm — conductivity: 38 MS/m |
| Core insulation | XLPE, 10.675 mm, eps_r = 2.2 |
| Sheath — inner/outer radii | 21.0 / 21.8 mm — conductivity: 60 MS/m |
| Sheath insulation | PVC, 2.2 mm, eps_r = 2.8 |
| Duct (HDPE) — inner/outer radii | 50.8 / 57.8 mm — eps_r = 2.35, eccentric (cable resting on the bottom) |
| Fourier order | 10 |

The duct is defined in `cable_definition.sheath.enclosure` — each of the 3 phases gets
an identical duct, generated eccentrically (`vertical_offset = enclosure_inner_radius
- cable_outer_radius`) so that the cable rests on the lower inner wall of the duct,
with the cable center (not the duct center) at depth `burial_depth`.

---

## Script Structure

### Helper function `_override_sheath_insulation` (lines 28–42)

Replaces the physical duct (`enclosure`) of each sheath with an equivalent
homogeneous insulation layer on the `sheath` conductor itself — the ERS/GMD technique
used to fold the HDPE duct's dielectric effect into the standard core+sheath internal-
parameter formulas, which have no notion of a duct layer on their own. It scans
`model.values()` looking for `conductor_name == 'sheath'` and applies the replacement to
every sheath found — generalizing to any number of phases (at the fixed indices
`2, 4, 6` of the 3-cable model).

### Function `compare_capacitance` (lines 44–70)

Prints a table comparing `C_11` (core) and `C_22` (sheath) of phase A, in uF/km,
for the three duct models (`'1'` Underground, `'2'` ERS, `'3'` GMD case 3.1) and,
when available, the COMSOL reference (energy method, `cmsl_shunt_params.txt`).
The `[0,0]`/`[1,1]` indices of the capacitance matrix always correspond to
phase-A core/sheath, regardless of the number of phases (per-cable ordering:
`np.kron(identity(N), Zij)`).

### Function `validate_against_case1_reference` (lines 72–102)

Checks whether scenario `'1'` ("Underground", duct ignored) agrees with the
`andreata_case1` MATLAB reference — same physical geometry of 3 SCC cables in a flat
arrangement, without a duct. It does not validate the duct effect itself (no counterpart
in case1); it serves as a code cross-check for the part of the pipeline shared between
the two cases (`underground_flat_model`, `InternalPerUnitParameters`, `PerUnitParameters`).
It computes the maximum relative error (inf norm) between the full quasi-TEM `Zs`/`Ysh`
of scenario `'1'` and the case1 MATLAB matrix.

### Function `main()` (lines 104–248)

#### 1. Construction of the models (lines 107–140)

Four MTL models are built from `SingleCoreCableModelGenerator`:

| Model | Geometry | Role |
|---|---|---|
| `model_0` / `mtl_0` | `flat_hdpe_enclosed_model()` — actual physical geometry, with duct | Basis of the schematic and of the ERS/GMD (provides `r4`..`r7` via `EquivalentRadiiSystems`). Not used directly as a scenario. |
| `model_1` / `mtl_1` | `underground_flat_model()` — ignores the duct entirely | Scenario `'1'` — "Underground" baseline, direct counterpart of `andreata_case1` |
| `model_2` / `mtl_2` | area-weighted equivalent insulation (ERS) via `_override_sheath_insulation` | Scenario `'2'` — "Area-weighted ERS" |
| `model_3` / `mtl_3` | GMD equivalent insulation, case 3.1 (Lafaia, 2015), keeping `r0 = r5` | Scenario `'3'` — "GMD case 3.1" |

`EquivalentRadiiSystems(mtl_0)` computes the equivalent radii (`equiv_rel_permittivity_epsr_area_weighted`,
`equivalent_parameters_from_gmd`) used to parametrize `model_2` and `model_3`.

#### 2. The `pul_data` dictionary (lines 142–154)

```
pul_data
├── 'frequencies'  -> 90 log-spaced points in [10^-2, 10^7] Hz
├── 'comsol'       -> {'scenarios': {'rho_g_100_epsr1_1_mf': {}}}  (filled if the file is available)
└── 'scenarios'    -> '1' (Underground), '2' (ERS), '3' (GMD case 3.1)
    (same soil formulation for all three: zg_form = yg_form = 'magalhaes_xue' —
    the comparison axis here is the duct model, not the soil formulation, which is
    the subject of andreata_case1)
```

#### 3. MATLAB data (lines 156–181)

Two distinct `MatlabDataReader`:

- **Duct-specific data** (`prefix='andreata_case2'`) — points at
  `Results/*.mat` of `andreata_case2` itself (files named
  `andreata_case2_*.mat`); stored as scenario `'measured'`.
- **`andreata_case1` data** (`prefix='andreata_case1'`) — points explicitly
  at the `testData/andreata_case1/andreata_case1.py` script, reading
  `andreata_case1/Results/andreata_case1_*.mat` (already existing and validated).
  Stored as scenario `'case1_no_duct'`, alongside `'measured'`, for a cross-check
  with scenario `'1'`.

Each case (`andreata_case1`, `andreata_case2`, `andreata_case3`) uses its own
file prefix (`andreata_case{N}_*.mat`) to avoid ambiguity between the
`Results/` directories of the three cases — previously they all used the generic prefix
`andreata_*`, which caused name collisions when copying/generating data from one case to
another.

Both `MatlabDataReader` use the same conductor permutation
`conductor_order=[0, 3, 1, 4, 2, 5]` (type grouping in MATLAB -> per-cable grouping in
pyLCP), inherited from `andreata_case1`.

#### 4. COMSOL processing — optional (lines 183–203)

Attempts to load `Results/cmsl_ground_return_impedance.txt` (ground-return family,
`get_earth_return_parameters`) via `ComsolPostProcessor`. It **does not yet exist** for
this case (see [Limitations](#limitations-and-pending-data)) — the block is skipped with
a warning and `pul_data['comsol']` is emptied. When present, the internal parameters used
to assemble the COMSOL quasi-TEM matrix come from the GMD case 3.1 model (`mtl_3`, the
most complete of the three variants) — only the ground return is actually measured in
COMSOL.

#### 5. PUL parameters per scenario (lines 205–229)

For each of the 3 scenarios (`'1'`, `'2'`, `'3'`), via `InternalPerUnitParameters` (form
`'hybrid'`) and `PerUnitParameters`:

1. Core+sheath internal matrices — depend on the duct model, not on the soil.
2. `earth_return_parameters(zg_form, yg_form)` + `quasi_tem_approx_matrices(...)` —
   **important note:** the `Zg`/`Yg` formulation does not model the HDPE duct in any
   of the three scenarios; the difference between them is only in the effective outer
   radius used in the self image term (which ERS/GMD extend out to `r7`). There is no
   "analytical solution with a duct" here — the three are approximations, and the actual
   validation of the duct's effect on the ground return depends on the pending
   comparison with COMSOL/MATLAB.

#### 6. Tables and validation printed to the console (lines 231–235)

- `compare_capacitance(pul_data, c_comsol)` — `C_11`/`C_22` table for the three duct
  models (+ COMSOL, if `cmsl_shunt_params.txt` is present).
- `validate_against_case1_reference(pul_data)` — maximum relative Zs/Ysh error of
  scenario `'1'` vs. `andreata_case1` MATLAB.

#### 7. Plot generation (lines 238–245)

`SCCPlotter.compare_complete_matrices` generates, for the `key_list` currently active in
the script:

| Key (`plot_config.py`) | Content |
|---|---|
| `core_self_impedance` | `Z_cc` — core self-impedance, phase A |
| `mutual_impedance_core_sheath` | `Z_cs` — core-sheath mutual, phase A |
| `sheath_self_impedance` | `Z_ss` — sheath self-impedance, phase A |
| `earth_return_impedance_phase_a` | `Zg` — ground-return self-impedance, phase A |

Each plot overlays the three analytical curves (`DUCT_MODEL_TEMPLATE`: Underground/ERS/GMD)
and, when available, the duct-specific COMSOL/MATLAB markers
(`COMSOL_TEMPLATE`/`MATLAB_TEMPLATE`) and the `andreata_case1` cross-reference
(`MATLAB_CASE1_NO_DUCT_TEMPLATE`, orange `+` marker).

`plot_config.py` also defines 5 other keys (`core_self_admittance`,
`mutual_admittance_core_sheath`, `sheath_self_admittance`,
`earth_return_admittance_phase_a`, `earth_return_potential_coeff_phase_a`), not included
in the current `key_list` of `main()` — they can be added to the
`compare_complete_matrices(key_list=[...])` call as needed.

`GroundReturnMTLRepresentation(__file__, mtl_0, units='centimeter').system_schematic()`
generates `Results/system_schematic.png` from the actual physical geometry (`model_0`, with
duct).

---

## Dependencies

| Module | Role |
|---|---|
| `models.single_core_cable.SingleCoreCableModelGenerator` | Builds the MTL models (`flat_hdpe_enclosed_model`, `underground_flat_model`) from the JSON |
| `mtl_main.source.MulticonductorTransmissionLine` | Represents the multiconductor MTL system |
| `analytical_forms.single_core_cable.InternalPerUnitParameters` | Internal cable matrices (core + sheath) |
| `analytical_forms.single_core_cable.PerUnitParameters` | Full PUL matrices with ground return |
| `analytical_forms.single_core_cable.EquivalentRadiiSystems` | ERS (area) and GMD equivalent radii/permittivities for the duct |
| `analytical_forms.single_core_cable.apply_semiconducting_layer_correction` | Semiconducting-layer correction (core and sheath) |
| `utils.comsol_data.ComsolPostProcessor` | Reading and processing of the COMSOL data (ground-return family) |
| `utils.matlab_data.MatlabDataReader` | Reading of the MATLAB reference matrices (own and `andreata_case1`) |
| `plotter.scc_plotter.SCCPlotter` | Plot generation (`compare_complete_matrices`) |
| `mtl_main.graphics.GroundReturnMTLRepresentation` | MTL circuit schematic |
| `plot_config.PLOT_CONFIG` | Layout, series and limits of every plot |

---

## Limitations and Pending Data

- **Duct COMSOL data does not yet exist.** `get_general_parameters('cmsl_ground_return_impedance')`
  returns `None` — the block is skipped gracefully. When added, it must follow the same
  column format as `andreata_case1` (`cmsl_ground_return_impedance.txt`, suffixes
  `_rho_g_100_epsr1_1_mf_vcoil_1/2/3`); the scenario key `'rho_g_100_epsr1_1_mf'` is
  already fixed in the code, matching the soil of the JSON (rho=100 Ohm.m, eps_r=1).
- **Duct MATLAB data is already available.** `MatlabDataReader(__file__, ...)` reads `.mat` in
  `Results/` with the prefix `'andreata_case2'` (e.g. `andreata_case2_frequency_range.mat`,
  `andreata_case2_series_impedance_matrix.mat` — the same set of variables that
  `andreata_case1` uses, each case with its own `andreata_case{N}_*` prefix to
  avoid name collisions between the `Results/` directories).
- **`Results/cmsl_shunt_params.txt` is from the old single-cable study** — reused
  only for the "COMSOL (energy method)" line of the capacitance table; it is not a
  reference for the three-phase system and should be replaced when the new data arrives
  (same file name expected).
- **`Results/cmsl_coaxial_cable_impedance.txt`, `cmsl_series_impedance_core_excitation.txt`,
  `cmsl_series_impedance_sheath_excitation.txt`** are leftovers from the old lineage
  (isolated-cable API) and are no longer read by the current script — dead data in
  `Results/`, not blocking execution.
- **The ground return never models the duct**, in any of the 3 scenarios — the only way
  to validate the real effect of the duct on the ground return is the COMSOL/MATLAB
  comparison above, still pending.

Full chronological details of the fixes that led to the current state of the script
(generalization to 3 phases, concentric vs. eccentric duct bug, `center_point` aliasing
bug, etc.) are in `BUGS_AND_FIXES.md`.

---

## Running

```bash
python -m testData.andreata_case2.andreata_case2
```

**Expected output:** capacitance table and cross-validation printed to the console, 4
plots (per the current `key_list`) plus `system_schematic.png` saved/displayed from
`Results/`. The duct-specific COMSOL/MATLAB data is optional — the analytical simulation
and the cross-validation against `andreata_case1` run normally without them.

---

## References

- **Andreata, L. E. B.** (2025). *Análise das Características de Propagação e de
  Transitórios Eletromagnéticos em Cabos Subterrâneos Instalados em Tubos Não Metálicos
  no Contexto de Parques Eólicos.* UFMG. https://hdl.handle.net/1843/2086
- **Xue, H.** (2018). *General Formulation and Accurate Evaluation of Earth-Return
  Parameters for Overhead/Underground Cables*. PhD thesis, École Polytechnique de
  Montréal. (`magalhaes_xue` ground-return formulations)
- **Lafaia, I.** (2015). GMD method for equivalent duct permittivity (case 3.1).
