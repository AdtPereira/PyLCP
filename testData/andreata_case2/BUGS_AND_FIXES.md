# Diagnostics and fixes — `andreata_case2` (Configuration 2, Figure 5.2)

**Context:** `andreata_case2` models Configuration 2 of the reference (Figure 5.2): three
power SCC cables (core + sheath) in a flat arrangement, each installed inside its
own HDPE duct, buried in the soil. Unlike `andreata_case1`
(Configuration 1, directly buried cables) and `andreata_case3` (Configuration 3, SCC +
ECC without a duct), here **there is no analytical ground-return solution that models the
duct** — the pyLCP Zg/Yg formulation sees only the position / outer radius of each cable,
never the presence of a duct. Validating the duct effect depends entirely on
COMSOL/MATLAB.

The script found at the start of this round of work came from an **older development
lineage**: it modeled a single cable inside a duct (no three phases, no coupling between
cables, no ground return), using a COMSOL API specific to an isolated cable and a plotter
(`HDPEPlotter`) and a `plot_config.py` in a format different from the one adopted by
`andreata_case1`/`andreata_case3`. That same old lineage still serves the cases
`hdpe_2000mm2`, `hdpe_300mm2`, `hdpe_ecc_2000mm2`, `hdpe_ecc_300mm2` and was not changed.

This document records, in chronological order, what was done:
1. generalization of the model generator to 3 phases with an individual duct (Figure 5.2);
2. positioning bug — concentric instead of eccentric duct;
3. aliasing bug between the ducts of the 3 phases;
4. update of the JSON to the three-phase arrangement;
5. rewrite of `andreata_case2.py` in the mold of `andreata_case1`/`andreata_case3`;
6. generalization of the fixed sheath indices (`model_2[2]`/`model_3[2]`) to N phases;
7. rewrite of `plot_config.py` in the `SCCPlotter` config-driven format;
8. cross-validation against the `andreata_case1` MATLAB (no-duct scenario);
9. known limitations and pending external data.

---

## Item 1 — New model generator: `flat_hdpe_enclosed_model`

**Where:** `models/single_core_cable.py`

**What:** there was no method that generated N cables in a flat arrangement, each
inside its own HDPE duct — `eccentric_hdpe_enclosed_model`/
`concentric_hdpe_enclosed_model` were hardcoded for a single cable (the origin of the
old script). The MTL layer (`SingleCoreCableInHDPEStrategy` in `mtl_main/strategy.py`), by
contrast, was already generic for N cables — it is essentially a copy of
`SingleCoreCableStrategy` (used by case1's three-phase `type: "scc"`), grouping
conductors by `center_point`. In other words, the limitation was only in the geometry
generation, not in the computation infrastructure.

**Fix:** new method `flat_hdpe_enclosed_model(host_conductor='sheath')`,
generalizing `eccentric_hdpe_enclosed_model` (1 cable) to N cables, the same way
`underground_flat_model` generalizes `conventional_single_phase`: same
`burial_depth`/`spacing`/`num_conductors` from the `arrangement`, an identical duct
injected into each cable.

---

## Item 2 — Duct generated concentric; it should be eccentric (cable on the bottom of the duct)

**Where:** `models/single_core_cable.py` (`flat_hdpe_enclosed_model`)

**What:** the first version of the method used **concentric** geometry (cable and duct
sharing the same center) — by direct analogy with `concentric_hdpe_enclosed_model`
and an initial reading of Figure 5.2. The user pointed out that this is not the convention
used in the rest of the repository: `schematic_original_hdpe.png` (case `hdpe_300mm2`,
already existing) shows the cable **resting on the lower inner surface of the duct** — the
same convention as `eccentric_hdpe_enclosed_model` — and not centered in it.

**Impact:** with `burial_depth` (`h = 1.20 m`) fixed at the center of the **cable** (not
the duct), the concentric geometry placed the duct center at exactly the same level as the
cable; the correct one requires shifting the duct center **upward**, preserving `h` at the
cable center.

**Fix:** rewritten for eccentricity, replicating the logic of
`eccentric_hdpe_enclosed_model`: `vertical_offset = enclosure_inner_radius -
cable_outer_radius`; each phase's duct center is `(cable_center[0], cable_center[1] +
vertical_offset)`. Applied identically to the 3 phases (same duct, same offset).

**Validated:** `system_schematic.png` regenerated — the 3 cables appear resting on the
bottom of their respective ducts, matching `schematic_original_hdpe.png` visually and
Figure 5.2 (cables at x = 0 / 0.2 / 0.4 m, cable-center depth h = 1.20 m).

---

## Item 3 — Aliasing bug: the 3 ducts would share the same `center_point`

**Where:** `models/single_core_cable.py` (`flat_hdpe_enclosed_model`)

**What:** the pattern used by the 1-cable methods (`eccentric_/concentric_hdpe_enclosed_model`)
injects `enclosure_data` **by direct reference** into the host conductor found (there is
only one). Naively generalized to N cables, the same `enclosure_data` dict (obtained a
single time from `self.sheath['enclosure']`) would be shared by the 3 sheaths — each
assignment of `v['enclosure']['center_point']` in the loop would overwrite the same
instance, so that **all 3 phases would end up with the `center_point` of the duct of the
last processed phase**.

**Fix:** each conductor gets its own copy via `copy.deepcopy(enclosure_data)`
before setting its individual `center_point`. (`import copy` added at the top of the
file.)

---

## Item 4 — `andreata_case2.json`: three-phase arrangement

**Where:** `testData/andreata_case2/andreata_case2.json`

**What:** `arrangement` described 1 cable (`"type": "single"`, `"num_conductors": 1`,
`"spacing": null`, `"fourier_order": 4`).

**Fix:** updated to `"type": "flat"`, `"num_conductors": 3`, `"spacing": 0.2`
(m, center-to-center, per Figure 5.2), `"fourier_order": 10` (aligned to
`andreata_case1`, previously divergent with no apparent reason). `burial_depth` (1.2 m) and
the `cable_definition` definition (cable + duct) already matched Figure 5.2 and were not
changed. `name`/`note` updated to describe the three-phase Configuration 2.

---

## Item 5 — Rewrite of `andreata_case2.py` in the mold of case1/case3

**Where:** `testData/andreata_case2/andreata_case2.py`

**What:** the original script reflected the old development lineage (see the Context
section) and could not simply be "adapted" — several whole pieces were missing because,
with a single cable, they did not make sense:
- it never called `apply_semiconducting_layer_correction`, even though the JSON already
  defined `semiconducting_layer` for the core and sheath (correction never applied — a
  latent inconsistency, independent of the move to 3 phases);
- it did not compute ground return (`PerUnitParameters`) or quasi-TEM matrices — without
  coupling between phases, this had never been necessary;
- it did not use `MatlabDataReader` (case1/case3 do);
- it read COMSOL data through the isolated-cable API (`cmsl_coaxial_cable_impedance`/
  `get_coaxial_cable_parameters`/`get_scc_internal_impedance_elements`), incompatible
  with the three-phase ground-return family (`cmsl_ground_return_impedance`/
  `get_earth_return_parameters`) used by case1/case3;
- it plotted with `HDPEPlotter` (its own "data_series" API, incompatible with the
  config-driven `plot_config.py` format used by `SCCPlotter`).

**Fix:** rewritten following the skeleton of `andreata_case1.py`/`andreata_case3.py`:
- `model_0 = flat_hdpe_enclosed_model()` -> `apply_semiconducting_layer_correction(...)` ->
  `MulticonductorTransmissionLine(model_0)` — actual physical geometry, used for the
  schematic and as the basis of the ERS/GMD (via `EquivalentRadiiSystems(mtl_0)`, which
  provides `r4`..`r7`);
- `model_1 = underground_flat_model()` (+ same semiconducting correction) — "Underground"
  scenario, ignores the duct entirely; serves as a baseline and a direct counterpart of
  `andreata_case1` (see Item 8);
- `model_2`/`model_3` — the two equivalent duct models already existing in the old
  script (area-weighted ERS; GMD case 3.1) were preserved as the physical distinctive of
  the case, now extracted into a shared function `_override_sheath_insulation`
  (see Item 6) and feeding the same 3 scenarios (`'1'`, `'2'`, `'3'`) used in the rest of
  the pipeline;
- ground-return block added to the scenario loop: `PerUnitParameters(...).earth_return_parameters(...)`
  + `quasi_tem_approx_matrices(...)`, analogous to case1, but with a single soil
  formulation (`magalhaes_xue`) — the relevant comparison axis here is bare/ERS/GMD, not
  soil formulations (that is the subject of case1);
- `MatlabDataReader` added (prefix `andreata_hdpe`, see Item 9 about the file-name
  assumption);
- COMSOL switched to the ground-return family (`cmsl_ground_return_impedance` /
  `get_earth_return_parameters('rho_g_100_epsr1_1_mf')` — key chosen because it already
  matches the soil of the JSON, rho=100 Ohm.m / eps_r=1); the internal parameters used to
  assemble the COMSOL quasi-TEM matrix come from the GMD case 3.1 model (`mtl_3`), as it
  is the most complete of the three variants — only the ground return is actually
  measured in COMSOL;
- `HDPEPlotter` replaced by `SCCPlotter` (`compare_complete_matrices`), with
  `plot_config.py` rewritten (see Item 7);
- `compare_capacitance()` (printed C_11/C_22 table for phase A) preserved almost
  unchanged — it already indexed by scenario (`'1'`/`'2'`/`'3'`), so it generalizes with
  no changes to N=3 phases (the `[0,0]`/`[1,1]` indices of the capacitance matrix remain
  phase-A core/sheath, since the per-cable ordering does not change);
- the `get_shunt_capacitance_elements()` (COMSOL) read is now explicitly gated by
  `'cmsl_shunt_params' in cmsl_processor.cmsl_reader.data` — in the old script, that
  call (which accesses `self.cmsl_reader.data['cmsl_shunt_params']` directly, without
  absence handling) was gated by the presence of a *different* COMSOL file
  (`cmsl_coaxial_cable_impedance`), which was already imprecise and would have become an
  unhandled `KeyError` as soon as the COMSOL family was switched.

---

## Item 6 — Fixed sheath indices (`model_2[2]`) did not generalize to N phases

**Where:** `testData/andreata_case2/andreata_case2.py`

**What:** the old (1-cable) script located the sheath to apply the equivalent
permittivity by a fixed position: `model_2[2]['enclosure'] = None`,
`model_2[2]['insulation']['thickness'] = ...` (index `2` = the only existing sheath,
since the insertion order is `0=soil, 1=core, 2=sheath`). With 3 phases, the sheaths
are at indices `2, 4, 6` — the old code, if simply reused, would have changed only the
phase-A sheath and left B/C with the physical duct intact (a silent inconsistency, no
error).

**Fix:** extracted into the function `_override_sheath_insulation(model, thickness,
relative_permittivity)`, which scans `model.values()` for `conductor_name == 'sheath'` and
applies the replacement to every sheath found — correct for any N.
`model_3` keeps the original insulation thickness (`model_0[2]['insulation']['thickness']`,
read before the substitution — Case 3.1 keeps r0 = r5), changing only the permittivity;
`model_2` uses the thickness extended out to r7.

---

## Item 7 — Rewrite of `plot_config.py` in the `SCCPlotter` format

**Where:** `testData/andreata_case2/plot_config.py`

**What:** the old file used a `data_series` scheme (list of dicts with
`source`/`scenario_key`/`data_path`/`plot_style`/`series`) specific to `HDPEPlotter`,
incompatible with the `series_to_plot`/`path`/`p`/`q`/`components` format that `SCCPlotter`
expects (also used by `andreata_case1`/`andreata_case3`).

**Fix:** rewritten in that second format. Series templates:
`DUCT_MODEL_TEMPLATE` (scenarios `'1'`/`'2'`/`'3'` — Underground/ERS/GMD, one line each),
`COMSOL_TEMPLATE` (`rho_g_100_epsr1_1_mf`), `MATLAB_TEMPLATE` (`'measured'`, case-specific
data, not yet existing). 9 plots via `compare_complete_matrices`:
6 internal (core/sheath/mutual x impedance/admittance, all sensitive to the duct model
because ERS/GMD change the thickness and/or permittivity of the sheath's outer insulation)
+ 3 ground return per phase (impedance, admittance, potential coefficient — see the note in
Item 5 about the formulation not modeling the duct).

---

## Item 8 — Cross-validation with the `andreata_case1` MATLAB

**Where:** `testData/andreata_case2/andreata_case2.py` (`validate_against_case1_reference`),
`testData/andreata_case2/plot_config.py` (`MATLAB_CASE1_NO_DUCT_TEMPLATE`)

**Motivation:** scenario `'1'` ("Underground", duct ignored, generated by
`underground_flat_model()`) describes exactly the same physical geometry as
`andreata_case1` (3 SCC cables in a flat arrangement, without a duct). Since
`andreata_case1` already has validated MATLAB reference data (see
`andreata_case1/BUGS_AND_FIXES.md`), it serves as a code cross-check for the part of
the pipeline the two cases share (`InternalPerUnitParameters`, `PerUnitParameters`),
without depending on the duct-specific COMSOL/MATLAB data (which does not yet exist).

**Implementation:** a second `MatlabDataReader`, pointed at the script path of
`andreata_case1` (not case2's own `__file__`) — so it reads
`testData/andreata_case1/Results/*.mat` instead of `testData/andreata_case2/Results/*.mat`
— loads the same prefix (`'andreata'`) and `conductor_order` (`[0,3,1,4,2,5]`) that
`andreata_case1.py` uses for itself. The result is stored as the scenario
`'case1_no_duct'` inside `pul_data['matlab']['scenarios']`, alongside the scenario
`'measured'` (case2's own data, with duct, still absent) — the two can
coexist and be plotted together without conflict, since `SCCPlotter` indexes by scenario
key.

`validate_against_case1_reference()` prints the maximum relative error (inf norm) between
the full quasi-TEM matrix of scenario `'1'` (`series_impedance_matrix`,
`shunt_admittance_matrix`) and the case1 MATLAB reference.

**Result (with the `andreata_case1` MATLAB data already existing in its `Results/`):**
```
Zs  (full series impedance):  max relative error (inf norm) = 0.06%
Ysh (full shunt admittance):  max relative error (inf norm) = 0.15%
```
Excellent agreement — it confirms that the generation of the bare three-phase model and
the internal/ground-return parameter computation are correct in this part of the pipeline.
It does **not** validate the duct effect itself (`'2'`/`'3'`), which has no counterpart in
`andreata_case1`.

`plot_config.py` gained `MATLAB_CASE1_NO_DUCT_TEMPLATE` (orange `+` marker, key
`'case1_no_duct'`), added to `MATLAB_TEMPLATE` in the internal and ground-return plots,
to overlay this reference visually onto the black line of scenario `'1'`.

---

## Item 9 — Known limitations and pending external data

- **Duct COMSOL data does not yet exist.** `cmsl_processor.get_general_parameters('cmsl_ground_return_impedance')`
  returns `None` (file missing in `Results/`) — the whole block is skipped
  gracefully (`pul_data['comsol'] = {}`), same as `andreata_case1`. When added,
  it must follow the same column format as `andreata_case1`
  (`cmsl_ground_return_impedance.txt`, suffixes `_rho_g_100_epsr1_1_mf_vcoil_1/2/3`) — the
  scenario key `'rho_g_100_epsr1_1_mf'` is already fixed in the code, matching the
  soil of the JSON (rho=100 Ohm.m, eps_r=1); if the real file uses a different
  `rho_g`/`epsr1`, adjust the key in `andreata_case2.py` and `COMSOL_TEMPLATE`
  (`plot_config.py`).
- **Duct MATLAB data does not yet exist.** `MatlabDataReader(__file__, ...)` looks for
  `.mat` in `testData/andreata_case2/Results/` with the prefix `'andreata_hdpe'` (e.g.
  `andreata_hdpe_frequency_range.mat`, `andreata_hdpe_series_impedance_matrix.mat`, etc.
  — the same set of variables that `andreata_case1` uses, only with a different prefix).
  This prefix was chosen so as not to collide with the plain `'andreata'` reserved for
  case1 (see Item 8) — if the real prefix comes different, adjust the string in
  `andreata_case2.py`.
- **The `cmsl_shunt_params.txt` present in `Results/` is from the old single-cable
  study**, reused only for the "COMSOL (energy method)" line of the capacitance table —
  it is not a reference for the three-phase system and should be replaced when the new
  data arrives (same file name expected, `get_shunt_capacitance_elements()` does not
  change). The other legacy files in `Results/` (`cmsl_coaxial_cable_impedance.txt`,
  `cmsl_series_impedance_core_excitation.txt`, `cmsl_series_impedance_sheath_excitation.txt`)
  are no longer read by the current script (isolated-cable API, abandoned — see Item 5);
  they remain as dead data in `Results/` until a future cleanup, not blocking execution.
- **The ground return never models the duct**, in any of the 3 scenarios — see the
  extensive note in the code (`andreata_case2.py`, inside the scenario loop) and in the
  Context item above. The only way to validate the real effect of the duct on the ground
  return is the pending COMSOL/MATLAB comparison.

---

## Legacy modules not changed

`models/hdpe.py`, `models/model_generator.py`, `models/isolated_conductors.py` and the
`HDPEPlotter` class (`plotter/scc_plotter.py`) still serve the cases `hdpe_2000mm2`,
`hdpe_300mm2`, `hdpe_ecc_2000mm2`, `hdpe_ecc_300mm2` — they are the old lineage in fact, and
`andreata_case2` no longer depends on them (it uses exclusively
`models.single_core_cable.SingleCoreCableModelGenerator`, `analytical_forms.single_core_cable`
and `plotter.scc_plotter.SCCPlotter`, the same base as `andreata_case1`/`andreata_case3`).
