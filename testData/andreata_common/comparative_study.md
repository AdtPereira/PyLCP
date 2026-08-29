# Comparative study — andreata_case1 / case2 / case3

## 1. Common architecture (the "skeleton" of the three `andreata_caseN.py`)

The three scripts follow exactly the same pipeline mold:

```
imports (case_utils, ComsolPostProcessor, MatlabDataReader, SCCPlotter,
         SingleCoreCableModelGenerator, MulticonductorTransmissionLine,
         InternalPerUnitParameters/PerUnitParameters) + PLOT_CONFIG
  ↓
cable_generator = SingleCoreCableModelGenerator(__file__)
model = cable_generator.<case-specific method>()
model = apply_semiconducting_layer_correction(model, core, sheath)
mtl_model_a = MulticonductorTransmissionLine(model)
  ↓
pul_data = {'frequencies': logspace(-2,7,90), 'comsol': {...}, 'scenarios': {...}}
  ↓
[optional] ComsolPostProcessor(__file__) -> earth_return + quasi_tem per scenario
MatlabDataReader(__file__).get_scc_scenario_data(prefix, conductor_order)
  ↓
InternalPerUnitParameters + PerUnitParameters per scenario (loop)
  ↓
SCCPlotter(__file__, pul_data, PLOT_CONFIG).compare_internal_matrices(...)
GroundReturnMTLRepresentation(...).system_schematic(); plt.show()
```

This uniformity is a direct result of the rewrite work already recorded in the `BUGS_AND_FIXES.md` files — `andreata_case2.py` was deliberately rewritten "in the mold of case1/case3" (Item 5 of its changelog), and that is what is confirmed by reading the code today.

## 2. Comparative table

| Aspect | case1 | case2 | case3 |
|---|---|---|---|
| Geometry | 3 SCC directly buried (`underground_flat_model`) | 3 SCC, each in an individual HDPE duct (`flat_hdpe_enclosed_model`) | 3 SCC + 1 ECC, no duct (`flat_scc_with_ecc_cable_model`) |
| No. conductors | 6 (homogeneous) | 6 (homogeneous) | 7 (heterogeneous: 3x2 + 1x1) |
| Comparison axis | soil formulation (`magalhaes_xue` x `deconti`) x 3 profiles (rho, eps_r) = 6 scenarios | duct model (Underground / ERS / GMD) = 3 scenarios, a single soil formulation | same pattern as case1: 6 scenarios (formulation x soil profile) |
| `ComsolPostProcessor` called? | **yes** (ground return + combined internal) | **yes** (identical to case1) | **no** — the whole block is absent |
| `MatlabDataReader` | 1 instance, `conductor_order=[0,3,1,4,2,5]` | **2 instances** — its own + one pointing at case1's `Results` (cross-validation) | 1 instance, `conductor_order=[0,3,1,4,2,5,6]` |
| Local helper functions | none | `_override_sheath_insulation`, `compare_capacitance`, `validate_against_case1_reference` | none |
| `plot_config.py` | 328 lines | 333 lines | 352 lines |

## 3. Verified redundancies (with exact location)

**3.1 — `plot_config.py`: 131 lines byte-identical between case1 and case3**
Confirmed by `diff`: lines 1–131 of `testData/andreata_case1/plot_config.py` and `testData/andreata_case3/plot_config.py` are **exactly equal** — `COMSOL_TEMPLATE`, `XUE_TEMPLATE`, `VANCE_TEMPLATE`, `MATLAB_TEMPLATE`, `DECONTI_TEMPLATE` and the six `PLOT_TPL_*`. case2 duplicates a subset of this (`MATLAB_TEMPLATE` + the 6 `PLOT_TPL_*`, ~53 lines) verbatim.

**3.2 — `VANCE_TEMPLATE` is dead code in two files**
Defined in `andreata_case1/plot_config.py:40` and `andreata_case3/plot_config.py:40`, but never referenced in any `PLOT_CONFIG` of the three cases (confirmed via grep). A copy-paste leftover from another case that uses Vance's formulation.

**3.3 — `utils/case_utils.py`: a duplicated function + a broken function, both unused by the three cases**
- `matrix_to_string` is defined **twice** (`case_utils.py:35` and `case_utils.py:108`) — the second silently overrides the first (dead code, ~20 lines).
- `print_real_matrix` (`case_utils.py:116`) calls `_matrix_to_string(...)` — a name that **does not exist** in the module (the real function is called `matrix_to_string`, without the underscore). If anything ever calls `print_real_matrix`, it raises `NameError`. Since `from utils.case_utils import *` is used in the three `andreata_caseN.py`, the error would silently enter the namespace of each without being noticed — it would only fire if someone tried to use the function.
- Confirmed by grep: none of the three `andreata_*` cases uses `print_real_matrix`, `verify_kelvin_functions`, `matrix_viewer`, `format_complex_number` or `format_scientific_notation` — they are utilities from other test lineages (`mom_models.py`, `deConti_models.py`, the `_s40`/`_s50`/`_s100` cases) imported for free via `import *`.

**3.4 — Almost-identical comment blocks (copy-paste) in `andreata_case1.py` and `andreata_case2.py`**
The "COMSOL internal impedance (legacy measurement...)" block in `andreata_case1.py:106-114` and the equivalent in `andreata_case2.py:205-216` are the same text with small adaptations — including repeating the same risk note ("shares `pul_data['comsol']['frequencies']`... harmless today... if the two come to coexist"). It signals that this logic (loading the combined `cmsl_internal_impedance_matrix.txt`) is a natural candidate to become a single method in `ComsolPostProcessor` instead of a block replicated in each case's script.

## 4. Most relevant structural divergence: case3 does not load COMSOL

`andreata_case3.py` **never instantiates `ComsolPostProcessor`** — the whole block that case1/case2 have (`"Building COMSOL matrices..."` + `"Loading COMSOL internal impedance data..."`, ~35 lines) simply does not exist in case3. This is consistent with the fact that there is no `Results/cmsl_ground_return_impedance.txt` for that case (documented in Item 9, "pending", of case3's `BUGS_AND_FIXES.md`).

The problem is that **`andreata_case3/plot_config.py`** already comes with `comsol_series_to_plot: COMSOL_TEMPLATE` configured in 5 plots (`mutual_impedance_phase_a_sheath_ecc`, `self_impedance_ecc`, `self_admittance_ecc`, `earth_return_impedance_ecc`, `earth_propagation_constant`) — dead configuration today (never populated, `self.cmsl.get('frequencies')` always `None` -> the guard in `scc_plotter.py` skips silently), but which will **break exactly the way Items 4/5 already broke** the day someone adds the COMSOL file: Item 9 already documents that `get_quasi_tem_approx_matrices` (`utils/comsol_data.py:731`, `M = 2` hardcoded on line 738) does not have the `block_sizes` generalization that the analytical side got in Item 5 — it will break in `np.kron(4x4, ones((2,2)))` = 8x8 against a 7x7 `Zi`.

## 5. Cleanup/optimization opportunities, prioritized

| # | Action | Gain |
|---|---|---|
| 1 | Extract `COMSOL_TEMPLATE`/`XUE_TEMPLATE`/`VANCE_TEMPLATE`/`MATLAB_TEMPLATE`/`DECONTI_TEMPLATE`/`PLOT_TPL_*` into a shared module (`testData/andreata_common/plot_templates.py` or similar) imported by the three `plot_config.py` | Eliminates ~130 duplicated lines x 2 files; any style adjustment (color, marker, scale) applies to all three at once, with no risk of one case falling out of sync |
| 2 | Remove `VANCE_TEMPLATE` from the two files where it is dead | Reduces noise; if it is reintroduced later, add it next to the real use |
| 3 | Fix/remove `case_utils.py:108` (second definition of `matrix_to_string`) and `print_real_matrix` (broken reference to `_matrix_to_string`) | Eliminates dead code + a latent bug that only hasn't fired because nothing calls the function |
| 4 | Resolve the pending Item 9: pass an explicit `block_sizes` to `ComsolPostProcessor.get_quasi_tem_approx_matrices`, reusing `_expand_by_block_sizes` (already exists, created for Item 5) | Closes the gap before someone adds COMSOL data to case3 and hits the same bug already fixed on the MATLAB/analytical side |
| 5 | Extract the "combined internal impedance COMSOL" block (identical in case1.py/case2.py) into a convenience method in `ComsolPostProcessor`, e.g. `load_internal_and_earth_return(cmsl_processor, mtl_model, scenarios)` | Reduces ~35 duplicated lines per script to one call; the risk comment is documented once, in the code, not in two places |
| 6 | If/when case3 gets COMSOL ground-return data, copy the loading block from case1.py — today the `plot_config.py` is already ready, only the loader is missing | Prevents the "pending" from becoming another debug cycle like Items 4/5/6/7/8 |

## 6. What should **not** be touched without care

- The real architecture differences (case2 with `_override_sheath_insulation`/`compare_capacitance`/`validate_against_case1_reference`, the second `MatlabDataReader` instance pointing at case1) reflect a genuine physical need — case2 is the only one with a duct model and cross-validation against the duct-free case1 — they are not redundancy to eliminate, but also not worth generalizing to the other two cases without a concrete need.
- `MergedComsolDataReader` (in `utils/comsol_data.py`) seems redundant with `ComsolDataReader._parse_single_file`, but it actually serves a different lineage (`testData/bare_wire`) — it is not used by any of the three andreata cases, so it is not redundancy *within* the requested scope, just two parsing strategies coexisting in the same module for different cases.

## 7. Optimization plan — execution and results

The six actions of section 5 were implemented in six sequential phases, each tested (import + end-to-end run of the three `andreata_caseN.py`, plus regression on neighboring cases) before moving to the next. Two design decisions were made before starting:

- **Location of the shared templates module:** `testData/andreata_common/plot_templates.py` — a new, neutral folder, instead of `utils/` (which mixes data infrastructure with presentation config) or inside one of the cases (which would create a strange dependency on a "sibling" case).
- **Broken `print_real_matrix`:** fixed (not removed) — keeps the function available for future use.

### Phase 0 — Dead-code cleanup

- `utils/case_utils.py`: removed the second duplicated definition of `matrix_to_string` (line 108, it overrode the first); fixed the broken call in `print_real_matrix` (`_matrix_to_string` -> `matrix_to_string`).
- `andreata_case1/plot_config.py` and `andreata_case3/plot_config.py`: removed `VANCE_TEMPLATE` (dead code).

**Validation:** `print_real_matrix` tested in isolation (previously it raised `NameError`); the three `andreata_caseN.py` run end to end with no new errors.

### Phase 1 — `plot_config.py` templates extracted into a shared module

- Created `testData/andreata_common/plot_templates.py` with two layers: a base common to all three (`MATLAB_TEMPLATE`, `PLOT_TPL_REAL/IMAG/RESISTANCE/INDUCTANCE/CONDUCTANCE/CAPACITANCE`) and a case1/case3-specific layer (`COMSOL_TEMPLATE` with 3 profiles, `XUE_TEMPLATE`, `DECONTI_TEMPLATE`).
- `andreata_case1/plot_config.py` and `andreata_case3/plot_config.py`: the 131 duplicated lines became a 5-line import.
- `andreata_case2/plot_config.py`: imports only the base layer; keeps local what is specific to it (`COMSOL_TEMPLATE` with 1 profile, `MATLAB_CASE1_NO_DUCT_TEMPLATE`, `DUCT_MODEL_TEMPLATE`, `INTERNAL_COMSOL_TEMPLATE`).

**Validation:** structural check (`PLOT_CONFIG` of the three cases with the same number of keys as before, imported templates resolved by identity/equality within the dicts) + end-to-end run of the three cases, same console output as before.

### Phase 2 — Duplicated COMSOL loading block extracted

- New method `ComsolPostProcessor.load_scc_earth_return_and_internal_scenarios()` (`utils/comsol_data.py`), encapsulating the block previously duplicated in `andreata_case1.py:87-123` and `andreata_case2.py:183-224` (ground return per scenario + combined internal impedance). The risk comment about the shared `pul_data['comsol']['frequencies']` now lives once, in the method's docstring.
- `andreata_case1.py` and `andreata_case2.py` rewritten to call the new method instead of the inline block.

**Validation:** the three `andreata_caseN.py` + `scc_flat_xue.py` (another consumer of `ComsolPostProcessor`, not touched) run end to end; same output numbers for case2 (capacitance table, cross-validation errors 0.06%/0.15%). Separate finding, out of scope: `testData/scc_flat_xue/scc_flat_xue.py:129` calls `plotter.scc_earth_return_impedance_matrix()`, a method that does not exist in `SCCPlotter` — a pre-existing bug, not introduced by this work, not fixed.

### Phase 3 — Item 9 (pending) resolved: `get_quasi_tem_approx_matrices` generalized by `block_sizes`

- `utils/comsol_data.py`: hardcoded `M = 2` + per-frequency `np.kron` loop replaced by `block_sizes = internal_matrices['block_sizes']` + `_expand_by_block_sizes` (reused from `analytical_forms/single_core_cable.py`, the same function created for Item 5 on the analytical side).

**Validation:**
- With real data (`scc_flat_xue`, homogeneous case `block_sizes=[2,2,2]`): result **bit-identical** (`np.array_equal`) to the old `np.kron`, for `Zg` and `Pg`.
- Synthetically (`block_sizes=[2,2,2,1]`, ECC case format): correctly assembles a 7x7 matrix with the blocks in the right place — the old code would break with `np.kron(4x4, ones((2,2)))` -> 8x8, incompatible with the 7x7 `Zi`.
- Regression on the three cases + `scc_flat_xue`: behavior identical to before.

### Phase 4 — COMSOL enabled (guarded) in `andreata_case3.py`

- Added the `ComsolPostProcessor(__file__)` + `load_scc_earth_return_and_internal_scenarios(...)` block to `andreata_case3.py`, absent until then. Since there is no `Results/cmsl_ground_return_impedance.txt` or `cmsl_internal_impedance_matrix.txt` for this case, the block stays inert today (it only emits the two missing-data warnings), but is already protected by the Phase 3 fix and ready for the 5 plots that `plot_config.py` already has configured with `comsol_series_to_plot`.

**Validation:** `andreata_case3.py` run end to end — exit 0, two expected warnings, the same 6 scenarios computed, no new warnings.

### Phase 5 — Wide regression

Ran `scc_34kV_andreata`, `hdpe_2000mm2`, `hdpe_300mm2`, `hdpe_ecc_2000mm2`, `scc_138kV_prysmian` (in addition to the three andreata and `scc_flat_xue` already covered in the previous phases) — all with exit 0, no tracebacks or hidden errors in the log.

### Final state

All six actions of section 5 were completed and validated. No numerical behavior change was introduced in the homogeneous cases (case1, case2, `scc_flat_xue`); the heterogeneous case (case3/ECC) gained protection against a bug that had not yet been exercised, but that would break as soon as COMSOL ground-return data was added to that case.
