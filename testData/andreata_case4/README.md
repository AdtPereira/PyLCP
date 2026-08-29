# andreata_case4

## Purpose

Simulation and validation script for **Configuration 4** (Figure 5.4) of the
reference: three single-core coaxial cables (SCC — core + sheath) in a flat
arrangement, each installed inside its own HDPE duct (identical to
`andreata_case2`), **plus one isolated earth conductor (ECC) sharing the duct of
the third cable** — the same position relative to cable C that the ECC already
has in `andreata_case3` (Configuration 3), only now inside a duct. It computes
the per-unit-length (PUL) parameters — internal (core/sheath/ECC) and ground
return — comparing three distinct analytical ways of representing the
dielectric effect of the duct, and cross-validates against the `andreata_case3`
MATLAB data.

Reference: Andreata, Luis Eduardo Batista. *Análise das Características de
Propagação e de Transitórios Eletromagnéticos em Cabos Subterrâneos
Instalados em Tubos Não Metálicos no Contexto de Parques Eólicos.* Graduate
Program in Electrical Engineering, UFMG, 2025.
https://hdl.handle.net/1843/2086

---

## System Configuration (JSON)

Parameters loaded from `andreata_case4.json`:

| Parameter | Value |
|---|---|
| Arrangement | Flat, 3 SCC cables + 1 ECC |
| Burial depth | 1.2 m (center of each SCC cable) |
| Spacing between SCC cables | 0.2 m (center to center) |
| Soil — conductivity | 0.01 S/m (rho = 100 Ohm.m) |
| Soil — relative permittivity | 1.0 |
| Core — outer radius | 10.325 mm — conductivity: 38 MS/m |
| Core insulation | XLPE, 10.675 mm, eps_r = 2.2 |
| Sheath — inner/outer radii | 21.0 / 21.8 mm — conductivity: 60 MS/m |
| Sheath insulation | PVC, 2.2 mm, eps_r = 2.8 |
| Duct (HDPE) — inner/outer radii | 50.8 / 57.8 mm — eps_r = 2.35, eccentric (cable resting on the bottom), **identical to `andreata_case2`** |
| ECC — outer radius | 3.3 mm — conductivity: 38 MS/m, XLPE insulation 2.0 mm |
| ECC position | `ecc_horizontal_gap = 0.02802` m / `ecc_vertical_gap = 0.00886` m (center-to-center relative to cable C) — **the same values as `andreata_case3.json`**, with no recomputation: the position of the ECC relative to cable C does not change between Configuration 3 and 4, only the presence of the duct |
| Fourier order | 10 |

The heterogeneous model (used in the 3 analytical scenarios) contains **7
conductors**: 1 soil return (`line_id=0`) + 3 core/sheath pairs
(`line_id=1..6`) + 1 ECC (`line_id=7`) — same ID convention as
`andreata_case3`.

---

## Architecture decision: two separate physical models

`EquivalentRadiiSystems` (used to compute the ERS/GMD equivalent
radii/permittivities) requires `self.model.scc` as a **flat** dict (a single
cross-section) — incompatible with the **per-cable grouped** dict that the
heterogeneous strategy (`SingleCoreCableWithECCStrategy`, needed because of the
ECC) produces. For that reason this case uses **two distinct physical models**,
each only for what it needs:

- **`flat_hdpe_enclosed_model()`** (the same method as `andreata_case2`,
  homogeneous, without ECC) — used **only** to feed `EquivalentRadiiSystems`
  and obtain `r4..r7`/equivalent permittivities. The duct's dielectric effect
  does not depend on the presence of the ECC inside it, so this homogeneous
  model is sufficient for that purpose.
- **`flat_hdpe_enclosed_with_shared_ecc_model()`** (new method,
  `models/single_core_cable.py`) — full physical model (real duct on the 3
  phases + ECC at the center-to-center position inherited from `andreata_case3`)
  used **only for the schematic** (`system_schematic.png`). It **never** goes
  through `InternalPerUnitParameters`/`PerUnitParameters` — see the limitation
  below.

The 3 analytical scenarios (Underground/ERS/GMD) start instead from
`flat_scc_with_ecc_cable_model()` (the same method as `andreata_case3`,
heterogeneous, without a duct), with the duct effect approximated via an
equivalent-insulation substitution on the 3 SCC sheaths
(`_override_sheath_insulation`, the same technique as `andreata_case2`).

---

## Limitations and Pending Data

- **There is no pyLCP analytical curve for the SCC<->ECC coupling inside the
  shared duct.** The geometry of two eccentric conductors (SCC and ECC)
  sharing an HDPE duct is **non-canonical** — there is no closed-form
  analytical formulation for that specific configuration in pyLCP's
  quasi-TEM/GMD method. It is not an implementation gap to be fixed: it is a
  physical limit of the method. For that reason the sibling single-cable cases
  `hdpe_ecc_300mm2`/`hdpe_ecc_2000mm2` (which model exactly this shared
  geometry, via `hdpe_shared_enclosed_model()`) also never call
  `InternalPerUnitParameters`/`PerUnitParameters` on that model — they only use
  it for the schematic and for direct comparison with finite-element data
  (COMSOL). This case follows the same pattern: the only real reference for
  the duct's effect on the ECC is the direct comparison with MATLAB (see the
  `*_ecc` plots).
- **The 3 scenarios (Underground/ERS/GMD) never change the ECC itself** — the
  equivalent-insulation substitution (`_override_sheath_insulation`) only acts
  on the 3 SCC sheaths. For that reason, in the plots `self_impedance_ecc`,
  `self_admittance_ecc`, `earth_return_impedance_ecc` and
  `earth_return_potential_coeff_ecc`, the three analytical curves tend to
  coincide almost completely — expected behavior, not a bug.
- **The ground return never models the duct**, in any of the 3 scenarios — the
  same limitation already documented in `andreata_case2`: the `Zg`/`Yg`
  formulation sees only the position / outer radius of each cable, never the
  presence of the duct. The difference between the 3 scenarios in the SCC
  phases' ground return is only in the effective outer radius used in the self
  image term (extended out to `r7` by ERS/GMD).
- **COMSOL internal impedance data already available**
  (`Results/cmsl_internal_impedance_matrix.txt`) — case4's own simulation,
  3 conductors (core/sheath/ECC), including the ECC self term
  (`p=6, q=6` in the `internal_impedance_matrix.png` plot). Read by
  `ComsolPostProcessor.get_scc_internal_impedance_matrix_combined()`, which
  detects this 3-conductor format (vs. the legacy 2-conductor format of
  `andreata_case2`/`hdpe_300mm2`) by the presence of the `data1_2` column and
  assembles the matrix already in the global 7-conductor space (core=0,
  sheath=1, ECC=6) — see Item 3 of `BUGS_AND_FIXES.md`.
- **There is still no COMSOL ground-return data or internal admittance data**
  (`cmsl_ground_return_impedance.txt`, `cmsl_internal_admittance_matrix.txt`)
  — the `ComsolPostProcessor.load_scc_earth_return_and_internal_scenarios(...)`
  block is called anyway ("guarded" pattern from `andreata_case3`), it only
  emits missing-data warnings for those two parts, without blocking execution.
- **Duct MATLAB data is already available.** `Results/andreata_case4_*.mat` — same
  type-grouped `(7,7,90)` format as the other cases, `conductor_order=[0, 3,
  1, 4, 2, 5, 6]` (the ECC has no pair to swap, it stays last in both
  conventions — same order as `andreata_case3`).

---

## Running

```bash
python -m testData.andreata_case4.andreata_case4
```

**Expected output:** cross-validation against `andreata_case3` printed to the
console (max relative error ~0.05%/0.12% in `Zs`/`Ysh`, Underground scenario),
`internal_impedance_matrix.png`, `internal_admittance_matrix.png` and
`system_schematic.png` saved to `testData/andreata_case4/Results/`. The scope is
deliberately reduced for now (`main()` keeps the rest of the plots — per-duct-
model comparison and ground return — commented out in
`plotter.compare_complete_matrices(...)`, ready to re-enable); the COMSOL
internal impedance data already feeds the two active plots, the rest will be
switched on as more COMSOL data arrives (ground return, internal admittance).

---

## References

- **Andreata, L. E. B.** (2025). *Análise das Características de Propagação e de
  Transitórios Eletromagnéticos em Cabos Subterrâneos Instalados em Tubos Não Metálicos
  no Contexto de Parques Eólicos.* UFMG. https://hdl.handle.net/1843/2086
- **Xue, H.** (2018). *General Formulation and Accurate Evaluation of Earth-Return
  Parameters for Overhead/Underground Cables*. PhD thesis, École Polytechnique de
  Montréal. (`magalhaes_xue` ground-return formulation)
- **Lafaia, I.** (2015). GMD method for equivalent duct permittivity (case 3.1).
