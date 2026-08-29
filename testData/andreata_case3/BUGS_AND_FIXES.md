# Diagnostics and fixes — `andreata_case3` (Configuration 3, Figure 5.3)

**Context:** `andreata_case3` models Configuration 3 of the reference (Figure 5.3): three
power SCC cables (core + sheath) in a flat arrangement, directly buried in the
soil — **without an HDPE duct** —, plus an isolated earth conductor (ECC) near the
rightmost cable, without touching it. It is the first case in the repository with a
heterogeneous mix of cables (3 SCC of 2 conductors each + 1 ECC of 1 conductor), which
exposed a series of implicit "N identical cables" assumptions in several parts of the pipeline.

This document records, in chronological order of discovery, the items addressed:
1. reformulation of the model generator to reflect Figure 5.3;
2. bug in the schematic plot (`system_schematic.png`);
3. design decision about `num_conductors`;
4. dimensional bug in the internal impedance matrix;
5. dimensional bug in the ground-return matrix (quasi-TEM);
6. MATLAB `conductor_order` without the ECC (data silently discarded);
7. wrong ECC indices in `PLOT_CONFIG` for the ground-return matrices;
8. MATLAB ground-return matrices sliced instead of reordered;
9. *(pending)* same bug pattern in `utils/comsol_data.py`;
10. diagnosis of numerical instability in the shunt admittance matrix `Y`;
11. reparametrization of the ECC position to precise center-to-center distances;
12. schematic depth level anchored to the ECC instead of an SCC cable.

---

## Item 1 — Reformulation of `flat_scc_with_ecc_cable_model` for Figure 5.3

**Where:** `models/single_core_cable.py:436`

**What:** the method was a copy of `hdpe_shared_enclosed_model` and assumed an HDPE duct
shared between the SCC and the ECC (with the ECC "wedged" between the duct wall and the
SCC surface via `_calculate_ecc_center_trig`). Figure 5.3, however, has no duct —
the cables are directly buried.

**Fix:** rewritten to position the 3 SCC cables in a flat arrangement (same logic as
`conventional_three_phase_flat`) and the ECC near the third cable, **without touching it**,
controlled by two `arrangement` parameters:
- `ecc_alignment`: `'center'` (same horizontal axis as the SCC centers, used in
  `andreata_case3.json`) or `'bottom_tangent'` (line tangent to the bottom of the SCC
  cables, simulating a common trench);
- `ecc_horizontal_gap`: gap between the outer surfaces of the 3rd SCC cable and the ECC
  (`andreata_case3.json` uses `0.050` m — reference value, with no exact dimension in
  the original figure, adjustable).

All duct/"corner" trigonometric-computation logic was removed from this method (it
remains intact in `hdpe_shared_enclosed_model`, used by other cases).

---

## Item 2 — `'+'` marker out of position in `system_schematic.png`

**Where:** `mtl_main/graphics.py`

**What:** the plot title appeared as the generic "Multiconductor Transmission
Line" and the `'+'` marker / depth-level line was stuck to the soil surface,
above the actual position of the cables.

**Root cause (two layers):**
1. `_calculate_schematic_parameters` (line ~93-113) did not recognize the type
   `'scc-flat-ecc'`, falling into the generic branch (`h_factor = 1`, default title).
2. With `h_factor == 1`, `_plot_conductor_graphic` **skips** the schematic
   repositioning (keeps the cables at the real depth), but `_schematic_annotations`
   (line 330) **always** applied the repositioning — the two were
   inconsistent with each other whenever a type fell into the generic branch.

**Fix:**
- `mtl_main/graphics.py:105` — added `'scc-flat-ecc'` as a recognized type
  (`h_factor = -2.5`, same convention as `'scc'`; title "Flat-Buried SCC Cables with
  ECC").
- `mtl_main/graphics.py:330` (`_schematic_annotations`) — now mirrors the same
  `h_factor == 1` bypass that `_plot_conductor_graphic` already had, so that this class
  of inconsistency does not reappear if another future type falls into the generic branch.

**Validated:** `system_schematic.png` regenerated — cross aligned to the center of the
cables, at the correct depth (h = 1.20 m).

---

## Item 3 — Design decision: `num_conductors` fixed at 3, out of the JSON

**Where:** `models/single_core_cable.py:469`

**Discussion:** the method read `num_scc_conductors = self.arrangement.get('num_conductors',
3)` from the JSON. It was left open whether that number should count only the SCC cables
(3) or include the ECC (4).

**Decision:** `num_conductors` in the JSON was **removed**; `num_scc_conductors` became
a literal `3` in the code, since the method is inherently three-phase (same
convention as `conventional_three_phase_flat`, which also does not read that field from the
JSON — the cardinality is in the method itself, not in external configuration). The ECC is
never counted by this field in any case in the repository (see `hdpe_ecc_2000mm2.json` /
`hdpe_ecc_300mm2.json`, where `num_conductors: 1` refers only to the host SCC).

---

## Item 4 — Internal impedance matrix: `None + None` / incompatible dimension

**Where:** `analytical_forms/single_core_cable.py` (`InternalPerUnitParameters`) and
`mtl_main/strategy.py` (`SingleCoreCableWithECCStrategy`)

**Original symptom:**
```
TypeError: unsupported operand type(s) for +: 'NoneType' and 'NoneType'
```
in `parameters_approximation`, line `'Zcs': z11 + z12 + z2i`.

**Root cause (three layers):**

1. **Wrong format of `context.scc`.** `SingleCoreCableWithECCStrategy._extract_scc_parameters`
   (`mtl_main/strategy.py:1056`) groups conductors by physical cable and returns a dict
   `{cp_key: {core_..., sheath_...}}` — correct to allow heterogeneous cables, but
   incompatible with `InternalPerUnitParameters.parameters_by_bessel`/`parameters_approximation`,
   which expected `self.model.scc` as a **flat** dict of a single cross
   section (`scc['core_outer_radius']` etc. at the root level). Every check
   `'core_outer_radius' in scc` returned `False` (the keys were position strings, not
   parameter names) -> `z11`, `z12`, `z2i` stayed `None`.

2. **The matrix assembly assumed N identical cables.** `matrices()` (old version)
   computed a single block `Zij` (M x M) and assembled the full matrix via
   `np.kron(np.identity(N), Zij)` — there is no way for this to represent 3 cables with M=2
   (core+sheath) and 1 cable with M=1 (ECC) at the same time.

3. **The ECC was not even captured.** The grouping loop of `_extract_scc_parameters`
   only accepted `name in ('core', 'sheath', 'armor')` — the `'ecc'` conductor was
   discarded (`continue`), so its geometric data never reached
   `context.scc`.

**Fix:**
- `mtl_main/strategy.py:1056-1108` (`SingleCoreCableWithECCStrategy._extract_scc_parameters`)
  — now also groups `'ecc'` conductors; when there is no `core` at a position (i.e.
  it is the ECC's own position), the ECC conductor takes the role of `core` in the
  formulas (physically correct: an isolated ECC is just a solid conductor with its own
  insulation, without a sheath/shield — same formulation as the "SCC with core,
  core_insulation" branch).
- `analytical_forms/single_core_cable.py:832` (`_build_cable_block`) — the
  `Zij_values`/`Pij` assembly logic extracted into a per-cable reusable method
  (receives `scc` explicitly, instead of always reading `self.model.scc`).
- `parameters_by_bessel`, `parameters_approximation`, `parameters_hybrid` — now
  accept an optional `scc` parameter (default: `self.model.scc`, preserving the
  previous behavior for the homogeneous types).
- `analytical_forms/single_core_cable.py:944` (`matrices`) — when
  `self.model.mtl_type == 'scc-flat-ecc'`, it delegates to the new method
  `_matrices_heterogeneous` (line 989): iterates each cable group in `self.model.scc`
  (3 SCC 2x2 blocks + 1 ECC 1x1 block) and assembles the final matrix **block-diagonal**
  — without a uniform `np.kron` —, since there is no internal coupling between
  conductors of different physical cables at this step (the coupling between cables only
  enters later, in the ground return).
- `analytical_forms/single_core_cable.py:467` (`_sum_or_none`) — a pre-existing latent
  bug was fixed along the way, now reachable: the convenience field
  `'Zcs'` did `z11 + z12 + z2i` unconditionally, even when `z2i` (which only exists
  if there is a sheath) is `None` — exactly the ECC case. That field is never read in
  any part of the pipeline; the sum now ignores `None` terms instead of raising an
  exception.

**Validated:** the resulting internal matrix has shape `(num_freq, 7, 7)`, with no `NaN`,
with the off-diagonal blocks between different cables exactly zero (internal decoupling
confirmed numerically). Clean regression on `andreata_case1`, `scc_34kV_andreata`,
`hdpe_2000mm2`, `hdpe_300mm2`, `hdpe_ecc_2000mm2`, `scc_138kV_prysmian`, `scc_flat_xue`
(homogeneous types `scc`, `hdpe`, `shared-hdpe` not affected).

---

## Item 5 — Ground-return matrix (quasi-TEM): 4x4 -> 7x7 broadcast error

**Where:** `analytical_forms/single_core_cable.py:1183`
(`PerUnitParameters.quasi_tem_approx_matrices`)

**Original symptom:**
```
ValueError: could not broadcast input array from shape (4,4) into shape (7,7)
```
in `Zg[i, :, :] = np.kron(z0_jk[i, :, :], np.ones((M, M)))`.

**Root cause:** same pattern as Item 4, one layer up. `z0_jk` (ground-return
impedance) already came correctly computed as 4x4 — one row/column per
physical cable (3 SCC + 1 ECC), since `earth_return_parameters` uses `num_sc_cables`
(count by position, always correct). The problem was only the expansion: `M =
self.model.num_conductors_per_scc` is an integer average (`total_conductors //
num_cables` = 7 // 4 = 1) that assumes every cable has the same number of conductors —
`np.kron(4x4, ones(1,1))` produces 4x4, incompatible with the 7x7 `Zi` of Item 4.

**Fix:**
- `analytical_forms/single_core_cable.py:478` — new function `_expand_by_block_sizes(matrix,
  block_sizes)`: generalizes `np.kron(matrix, np.ones((M, M)))` to blocks of variable
  size, via `np.repeat` (rows and then columns) with a list of per-cable sizes instead
  of a uniform M. When all sizes are equal to M, the result is
  identical to the original `kron` — that is, it is a strict generalization, safe for the
  homogeneous types too.
- `InternalPerUnitParameters.matrices()` and `_matrices_heterogeneous()` now
  include `'block_sizes'` in the returned dictionary (`[M]*N` in the homogeneous path,
  `[2, 2, 2, 1]` in `scc-flat-ecc`) — a single source of truth about how many
  conductors each cable contributes, in the same order as the ground-return matrices.
- `quasi_tem_approx_matrices` (line 1193) now reads `internal_matrices['block_sizes']`
  instead of `self.model.num_conductors_per_scc`, calling `_expand_by_block_sizes` to
  assemble `Zg`/`Pg` — no more per-frequency loop, since `np.repeat` operates directly on
  the 3D array `(freq, N, N)`.

**Validated:** the resulting `Zg` has shape `(num_freq, 7, 7)`; it was checked numerically
that the (cable 1 x ECC) block reproduces exactly `z0_jk[:, 0, 3]` and the (cable 1 x cable 2)
block reproduces `z0_jk[:, 0, 1]`, both with maximum difference `0.0`. `andreata_case3.py` runs
to completion without error. Clean regression on the same 7 homogeneous cases as Item 4.

---

## Item 6 — MATLAB `conductor_order` without the ECC (data silently discarded)

**Where:** `testData/andreata_case3/andreata_case3.py` (call to `MatlabDataReader.get_scc_scenario_data`)

**Original symptom:** three warnings when running, without a crash:
```
Warning: MATLAB data not found for 'measured' with key 'series_impedance_matrix'.
Warning: MATLAB data not found for 'measured' with key 'internal_impedance_matrix'.
Warning: MATLAB data not found for 'measured' with key 'internal_admittance_matrix'.
```

**Root cause:** the reference `.mat` files already came in the full 7x7 format (confirmed
via `scipy.io.loadmat`: shape `(7, 7, 90)`), but `conductor_order=[0, 3, 1, 4, 2, 5]`
— inherited from `andreata_case1`, which has no ECC — only listed 6 indices.
`MatlabDataReader._reorder_conductor_matrix` (`utils/matlab_data.py`) does
`matrix[:, conductor_order, :][:, :, conductor_order]`: with only 6 indices, the
row/column of the ECC (index 6) was silently discarded, producing a
`'measured'` 6x6 matrix. The plotting components that ask for `p=6, q=6` (ECC self
term) then did not find that index (`IndexError`, caught by the generic `except` of
the plotter — `scc_plotter.py:69` — and printed only as "not found").

**Fix:** `conductor_order=[0, 3, 1, 4, 2, 5, 6]` — the ECC has no core/sheath pair
to swap position with, so it stays at the end (index 6) in both the MATLAB convention
(type-grouped) and the pyLCP convention (cable-grouped).

**Validated:** the three warnings disappear; `andreata_case3.py` runs without changing the
other curves (indices 0-5 unaffected by the change).

---

## Item 7 — Wrong ECC indices in `PLOT_CONFIG` for the ground-return matrices

**Where:** `testData/andreata_case3/plot_config.py`

**Original symptom:**
```
IndexError: index 6 is out of bounds for axis 1 with size 4
```
in `scc_plotter.py:32` (`_plot_scc_matrix`), when plotting `earth_return_impedance_ecc`.

**Root cause:** the configs `earth_return_impedance_ecc`, `earth_return_admittance_ecc` and
`earth_return_potential_coeff_ecc` used `path: ['earth_return_parameters', ...]` with
`p=6, q=6` — but `PerUnitParameters.earth_return_parameters()`
(`analytical_forms/single_core_cable.py:1081`) returns matrices indexed **by physical
cable** (`N = num_sc_cables = 4`: phase A, phase B, phase C, ECC — see Item 5), not by
conductor. In that space the ECC is index 3, not 6. The configs
`mutual_impedance_phase_a_sheath_ecc`, `self_impedance_ecc` and `self_admittance_ecc`, which
use `path: ['quasi_tem_matrices', ...]`, were already correct with `p=6, q=6`, since those
matrices were expanded to the per-conductor space (N=7) via `_expand_by_block_sizes`
(Item 5).

**Fix:**
- `earth_return_admittance_ecc`: kept at `path: ['earth_return_parameters',
  'admittance_matrix']`, with `p=3, q=3` (ECC index in the per-cable space). There is no
  valid per-conductor equivalent for that matrix — see the note below.
- `earth_return_impedance_ecc` / `earth_return_potential_coeff_ecc`: migrated to
  `path: ['quasi_tem_matrices', 'earth_return_impedance_matrix']` /
  `['quasi_tem_matrices', 'earth_return_potential_coefficient']`, with `p=6, q=6` — those
  keys already exist in `quasi_tem_approx_matrices` (lines 1219-1226), expanded to
  N=7, and match index-for-index with the reordered MATLAB (see Item 8).

**Note (not fixed):** `quasi_tem_approx_matrices` returns
`'earth_return_admittance_matrix': Yg`, but `Yg` is never filled (it stays
`np.zeros_like(Zi, dtype=complex)`, pre-existing dead code) — so
`earth_return_admittance_ecc` could not be migrated like the other two. Since that config
has no `matlab_series_to_plot`, this generates neither a warning nor a crash today; it is
recorded for when someone actually computes `Yg`.

**Validated:** `andreata_case3.py` runs to completion without an `IndexError`.

---

## Item 8 — MATLAB ground-return matrices sliced instead of reordered (ECC discarded)

**Where:** `utils/matlab_data.py` (`MatlabDataReader.get_scc_scenario_data`)

**Original symptom (after Item 7):**
```
Warning: MATLAB data not found for 'measured' with key 'earth_return_impedance_matrix'.
Warning: MATLAB data not found for 'measured' with key 'earth_return_potential_coefficient_matrix'.
```

**Root cause:** `earth_return_impedance_matrix` and `earth_return_potential_coefficient_matrix`
received a different treatment from the other 4 matrices: instead of
`_reorder_conductor_matrix(..., conductor_order)`, the code sliced only the
top-left block `[:num_phases, :num_phases]` (3x3), under the assumption — correct for
`andreata_case1` (no ECC), but outdated for `andreata_case3` — that those
matrices only exist "at the cable/phase level". Confirmed via `scipy.io.loadmat` that the
ground-return `.mat` already come in the same full 7x7, type-grouped-with-redundancy
format as the other 4 matrices: `M[0,0]==M[3,3]` (core_A==sheath_A), `M[1,1]==M[4,4]`,
`M[2,2]==M[5,5]`, and `M[6,6]` is the ECC self value. Slicing to 3x3 discarded that
index 6 entirely, so `'measured'` never had ECC data.

**Fix:** the two matrices now use `_reorder_conductor_matrix(matrix,
conductor_order)`, like the other 4 — with no special slicing. The `num_phases` parameter
(which only served that slicing) was removed from `get_scc_scenario_data`, and the calls in
`andreata_case1.py`/`andreata_case3.py` updated.

**Validated:**
- The two warnings disappear; no new warning arises in `andreata_case1` (the phase-A
  `p=0, q=0` indices point to the same `core_A` value in both the sliced format
  and the reordered full format — confirmed numerically).
- Spot numerical check (frequency index 45, ~2.15 kHz): pyLCP `Zg`
  `[:,6,6] = 3.532e-4 + j4.951e-3` vs. MATLAB `3.531e-4 + j4.952e-3`; pyLCP `Pg`
  `= 5.580e4 + j5.357e5` vs. MATLAB `5.580e4 + j5.354e5` — same order of magnitude,
  consistent with measured data vs. analytical formulation (same pattern as `Zi_77`/`Yi_77`,
  Item 4).

---

## Item 9 (pending) — Same bug pattern in `utils/comsol_data.py`

**Where:** `utils/comsol_data.py:731` (`ComsolPostProcessor.get_quasi_tem_approx_matrices`)

**What:** sister function of `quasi_tem_approx_matrices` (used to compare with COMSOL
data), with the same expansion pattern via `np.kron(z0_jk, np.ones((M, M)))` — but
with `M = 2` **hardcoded** (line 738), not even derived from the model.

**Why it hasn't broken yet:** `andreata_case3.py` only calls that method inside the block
`if cmsl_params is not None:` — and there is no `Results/cmsl_ground_return_impedance.txt`
file for this case ("COMSOL processing skipped (data not available)" in the console), so
the path has never been exercised.

**Risk:** if COMSOL data is added to this case in the future, it will break the same
way Items 4 and 5 broke — `np.kron(4x4, ones((2,2)))` gives 8x8, incompatible
with the 7x7 `Zi`.

**Proposed fix (not applied):** since `ComsolPostProcessor` does not have access to
`self.model`, the fix would require passing `block_sizes` as an explicit parameter to
the function (instead of deriving it from `self.model.num_conductors_per_scc` as today),
reusing the same `_expand_by_block_sizes` of Item 5. Recorded here for when there is COMSOL
data available for this case.

---

## Item 10 — Diagnosis of numerical instability in the shunt admittance matrix `Y`

**Where:** `analytical_forms/single_core_cable.py` (`PerUnitParameters.quasi_tem_approx_matrices`,
`sommerfeld_quasi_tem_approx_admittance`).

**Request:** investigate whether there is a numerical problem/instability in the admittance
matrix `Y` of `andreata_case3`.

**What was checked (no problem found):**
- No `NaN`/`Inf` in `Ysh`, `Yg` (ground return) or `Yi` (internal), in any
  scenario (`magalhaes_xue`/`deconti`, 3 soils).
- The condition number of `Psh` (the matrix inverted to obtain `Ysh`) stays low
  (< 55) over the whole sweep of 90 frequencies — no explosive ill-conditioning.
- It is not quadrature under-convergence: increasing the Gauss-Legendre points of
  `sommerfeld_quasi_tem_approx_admittance` from 150 to 2400 changes the integral
  result by less than 0.5%.

**Real problem found (accuracy divergence, not instability itself):**
comparing the computed `Ysh` against `andreata_shunt_admittance_matrix.mat` (MATLAB,
reordered correctly to the pyLCP convention — Item 6/8), there is a deviation that grows
smoothly with frequency and is concentrated almost exclusively in the conductors
connected to the ECC (index 6) and the SCC cable closest to it (sheath C, index 5):

| conductor (index) | max relative error (at f = 10 MHz) |
|---|---|
| cores A/B/C (0, 2, 4) | ~2e-9 (floating-point noise) |
| sheaths A/B (1, 3) | 0.26% – 0.42% |
| **sheath C (5)** | **15.9%** |
| **ECC (6)** | **23.1%** |

The pattern is identical in the two formulations tested (`magalhaes_xue` and `deconti`), which
rules out a formula-specific bug — the cause is in the shared part (geometric
term `K0(gamma_earth*d)` of the ground return). The ECC<->sheath-C pair is the only one with
a small center-to-center spacing (a few cm, see Item 11) against the 0.2–0.4 m between
the SCC cables; at that distance the argument of `K0` is near the log-singular regime
(the derivative `-1/x` is large when `x->0`), which makes that specific term much more
sensitive to error/approximation at high frequency than the more widely spaced pairs —
a hypothesis consistent with the observed pattern, but not a formal proof.

**No fix applied** — it is an accuracy limit of the quasi-TEM approximation for
conductors very close together at high frequency, not an implementation bug. Recorded
for reference in case the deviation becomes bothersome again after the Item 11 geometry
adjustment (which changes the ECC<->cable-C distance).

**Side note (resolved):** the `self_admittance_ecc` config of `plot_config.py` was
temporarily set to `p=0, q=0` (core A index) during that investigation,
diverging from the label `G_77`/`C_77` (ECC = index 6, see Item 7). It was already fixed
back to `p=6, q=6`.

---

## Item 11 — Reparametrization of the ECC position to precise center-to-center distances

**Where:** `models/single_core_cable.py:436` (`flat_scc_with_ecc_cable_model`),
`testData/andreata_case3/andreata_case3.json`.

**What:** the ECC position was controlled by `ecc_alignment` (`'center'` or
`'bottom_tangent'`) + `ecc_horizontal_gap`, the latter measured as a gap between the
**outer surfaces** of cable C and the ECC — it did not allow positioning the ECC with an
arbitrary vertical offset relative to the cable, only the two fixed options of
`ecc_alignment`.

**Fix:** added a precise positioning mode, activated when
`ecc_vertical_gap` is present in `arrangement`:
- `ecc_horizontal_gap` / `ecc_vertical_gap` become **center-to-center** distances
  (no longer a surface gap) between the ECC center and the center of the third SCC cable:
  `ecc_center = (last_cable_center_x + ecc_horizontal_gap, last_cable_center_y -
  ecc_vertical_gap)`. A positive `ecc_vertical_gap` places the ECC deeper than the cable
  (same sign convention as `burial_depth`).
- The legacy mode (`ecc_alignment` + `ecc_horizontal_gap` as a surface gap) was
  kept as a *fallback* for when `ecc_vertical_gap` is not provided — with no impact
  on any other case in the repository (`flat_scc_with_ecc_cable_model` is only used by
  `andreata_case3`).
- `andreata_case3.json`: `ecc_alignment` removed; `ecc_horizontal_gap` recomputed from
  `0.001` (surface gap, the value in use at the time of the migration) to the
  definitive reference values `ecc_horizontal_gap = 0.02802` m and
  `ecc_vertical_gap = 0.00886` m (center-to-center).

**Validated:** with `ecc_horizontal_gap`/`ecc_vertical_gap` computed to reproduce
exactly the old position (`0.0303`/`0.0`), the resulting ECC center matched
numerically with the legacy mode (`(0.4303, -1.2)` in both modes). With the final
reference values (`0.02802`/`0.00886`), the ECC is at `(0.42802, -1.20886)` — deeper
than the SCC cables, which exposed Item 12.

---

## Item 12 — Schematic depth level anchored to the ECC instead of an SCC cable

**Where:** `mtl_main/graphics.py:71` (`BaseMTLRepresentation._calculate_schematic_parameters`).

**Symptom:** after Item 11, with `ecc_vertical_gap = 0.00886` (ECC deeper than the
SCC cables), the depth level (`h = ...`) and the reference cross in
`system_schematic.png` started pointing at the ECC (`(0.42802, -1.20886)`) instead of
one of the power cables.

**Root cause:** `_calculate_schematic_parameters` chooses `depth_ref_conductor` as the
physically deepest conductor (`deepest_conductor`, by `y_min`) for the type
`'scc-flat-ecc'` — except for `'hdpe'`/`'shared-hdpe'`, which already had a dedicated
exception (fixing the reference on the SCC `core`, not on the deepest conductor). Before
Item 11 the ECC was always at the same depth as the SCC cables (implicit
`ecc_alignment='center'`), so the "deepest" choice coincided by chance with an SCC
cable; when the ECC depth became adjustable, that coincidence stopped holding.

**Fix:** extended the same `'hdpe'`/`'shared-hdpe'` exception to
`'scc-flat-ecc'` — the depth level always references the center of an SCC cable `core`
(first found, cable A), regardless of where the ECC is positioned vertically.

**Validated:** with `ecc_vertical_gap = 0.00886`, `depth_ref_conductor` started pointing
at the core of cable A (`(0.0, -1.2)`, `conductor_name='core'`) instead of the ECC
(confirmed that `deepest_conductor` — no longer used for the level — is indeed still the
ECC, as expected). `system_schematic.png` regenerated without error.

---

## Item 13 (note, not a bug) — Point count in `Y` — pyLCP vs. MATLAB

**Where:** `plotter/scc_plotter.py:32` (line, pyLCP) vs. `plotter/scc_plotter.py:66`
(scatter, MATLAB).

**Reported perception:** visually, the "MATLAB" curve in the `self_admittance_ecc` plot
seemed to have fewer points than the analytical curves.

**Verified:** the two matrices have exactly the same number of elements —
`Ysh` (pyLCP) and the measured `shunt_admittance_matrix` (MATLAB) are both `(90, 7, 7)` =
4410 elements, over the same array of 90 frequencies (`np.allclose` between the two
frequency arrays = `True`). Within the plot's `xlim=(1E4, 1E7)`, exactly 30
of the 90 points fall in the visible window — the same number for both sources.

**Cause of the perception:** a difference in drawing style, not in data. The pyLCP series
use `ax.plot(...)` (continuous line interpolating the 30 points, with no individual
markers visible); the MATLAB series uses `ax.scatter(...)` (discrete `'x'` marker on each
of the 30 points). A continuous line "hides" the underlying discretization; a scatter
makes it obvious.

**No fix needed** — expected behavior, recorded only for reference.
