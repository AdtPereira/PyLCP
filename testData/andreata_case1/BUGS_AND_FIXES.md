# Bug summary — `.mat` data import (`andreata_case1` validation)

**Context:** `testData/andreata_case1/andreata_case1.py` compares pyLCP analytical results (`InternalPerUnitParameters`) against a reference exported from MATLAB (`Results/andreata_internal_impedance_matrix.mat`, variable `Z`, shape `(6,6,90)` — complex internal impedance of a 3-cable single-core system, core+sheath, flat arrangement). **4 distinct bugs** were found, all in the MATLAB data import/alignment layer — the raw `.mat` data itself is correct.

---

## Bug 1 — Non-existent `.mat` file keys (`matlab_reader.data.get`)

**Where:** `andreata_case1.py` (before the fix)

**What:** The code looked up `matlab_reader.data.get('andreata_series_impedance_matrix')` and `matlab_reader.data.get('andreata_shunt_admittance_matrix')`. `MatlabDataReader` (`utils/matlab_data.py:63-75`) indexes the data by the **file name without extension** (glob of `Results/*.mat`). Only `andreata_internal_impedance_matrix.mat` exists in the directory — those two keys never existed (leftover from a previous naming scheme). Since access is via `.get()`, the error is silent: it returns `None` without an exception.

**Impact:** `pul_data['matlab']['scenarios']['measured']` had both fields set to `None`. In the `z11/z12/z22_internal_vs_matlab` plots, `_plot_scc_internal_vs_matlab` (`plotter/scc_plotter.py:154-159`) tried to index `None[:, p, q]`, raised a `TypeError` that was caught silently — result: the MATLAB reference points simply did not appear on the plot, only a console warning.

**Fix:** use the correct key (`andreata_internal_impedance_matrix`) and name the field as the plotter expects by default (`internal_impedance_matrix`).

---

## Bug 2 — Frequency axis assumed incorrectly

**Where:** `andreata_case1.py:40` (before the fix)

**What:** The code assumed `np.logspace(0, 7, num=90)` (1 Hz–10 MHz) both for the Python analytical computation and for labeling the 90 samples of the MATLAB `Z`. The actual sweep used in MATLAB — confirmed by the `Results/andreata_frequency_range.mat` file (variable `freq1`) provided later — is `np.logspace(-2, 7, num=90)` (**0.01 Hz**–10 MHz): same number of points and same ceiling, but starting 2 decades lower. Numerical check: `freq1` matches `np.logspace(-2,7,90)` with a maximum relative difference of 4e-15.

**Impact:** since `L = Im(Z)/(2*pi*f)`, dividing each sample by the wrong `f` systematically distorts the inductance curve. With the wrong axis, `L(f)` appeared to be **monotonically increasing ~73x** over the sweep (2.2 nH -> 162 nH) — physically impossible for an internal self-impedance (it should be a low-frequency plateau followed by a skin-effect drop). With the correct axis, the curve showed exactly the expected behavior: a constant plateau (~220 nH) while R is still in the DC regime, then a smooth decay down to ~162 nH at high frequency. The real part (resistance) is not sensitive to this error (it does not depend on division by `omega`), so the problem was only visible in the inductance — and did not appear when plotting the raw `Z` directly in MATLAB (which uses its own correct native axis).

**Additional note:** `num=90` was the only case in the entire repository with that point count — every sibling case (`scc_34kV_andreata`, `scc_132kV_xue`, `scc_138kV_prysmian`, etc.) uses `num=121`. Same class of copy-paste bug already identified in `ohtl_xue_sec43.py`.

**Fix:** load `andreata_frequency_range.mat` and use that vector as `pul_data['matlab']['frequencies']`, instead of reusing the Python analytical-side vector. (Note: the Python analytical-side vector was also updated to `np.logspace(-2, 7, num=90)`, so today they coincide — but the code keeps the explicit loading of the frequency `.mat` as the source of truth, with a fallback.)

---

## Bug 3 — Inconsistent `matlab_matrix_key` in `plot_config.py`

**Where:** `testData/andreata_case1/plot_config.py:156` and `:172` (before the fix)

**What:** The three `z11/z12/z22_internal_vs_matlab` plots have the suptitle "internal-only analytical (Zi) vs. MATLAB reference (Z)", but only `z11` used `'matlab_matrix_key': 'internal_impedance_matrix'`. `z12` and `z22` still pointed to `'series_impedance_matrix'` — a key that does not exist in `pul_data['matlab']['scenarios']['measured']` (a direct consequence of Bug 1).

**Impact:** even after fixing Bug 1, the `z12` and `z22` plots would still not show the MATLAB reference (same symptom as Bug 1: `KeyError`/`TypeError` caught silently).

**Fix:** aligned the three entries to `'matlab_matrix_key': 'internal_impedance_matrix'`.

---

## Bug 4 — Conductor ordering convention incompatible between pyLCP and MATLAB

**Where:** structural — discovered while reviewing `z12_internal_vs_matlab` and `z22_internal_vs_matlab` after fixing Bugs 1–3.

**What:** MATLAB exports the 6 conductors **grouped by type**: `[core_A, core_B, core_C, sheath_A, sheath_B, sheath_C]` (indices 0–2 = cores, 3–5 = sheaths). pyLCP assembles its internal matrix **grouped by cable**: `[core_A, sheath_A, core_B, sheath_B, core_C, sheath_C]`, via `np.kron(np.identity(N), Zij)` in `InternalPerUnitParameters.matrices()` (`analytical_forms/single_core_cable.py:848`), reflecting the `conductor_id` assignment order in `models/single_core_cable.py:127-141` (core, then sheath, per cable).

**Evidence:** the diagonal elements of the MATLAB `Z` confirm this — `Z[0,0]==Z[1,1]==Z[2,2]` (~7.8536e-5 Ohm at DC) and `Z[3,3]==Z[4,4]==Z[5,5]` (~1.549406e-4 Ohm at DC), matching exactly the theoretical DC resistance (`1/(sigma*A)`) of the core and the sheath, respectively, computed from the case JSON.

**Impact:** since the code used the **same `(p,q)` pair** to index both the Python internal matrix and the MATLAB matrix:
- `z11` (p=0,q=0) coincided by chance — index 0 is `core_A` in both conventions.
- `z22` (p=1,q=1): pyLCP took the **sheath A** self-impedance; MATLAB (unfixed) delivered the **core B** self-impedance (identical to core A) — different physical quantities being compared side by side on the same plot.
- `z12` (p=0,q=1): pyLCP took the **core A <-> sheath A** mutual coupling (same cable, strong coaxial coupling); MATLAB delivered **core A <-> core B** (different cables, weak external coupling) — incompatible order of magnitude for L (as much as ~11x larger in the wrong pair).

It was validated that the correct indices (`Z[3,3]` for sheath A, `Z[0,3]` for core A <-> sheath A) produce physically consistent curves, including `L(core<->sheath) -> L(sheath, self)` converging at high frequency — the expected unit-coupling limit when the skin effect concentrates the current on the facing surfaces.

**Fix:** in `andreata_case1.py:127-130`, the MATLAB matrix is reordered with the permutation `[0,3,1,4,2,5]` (applied to both conductor dimensions) right at load time, before being stored in `pul_data['matlab']`, aligning it to the pyLCP convention once and for all — with no need for different `(p,q)` indices per plot.

---

## Recommendations for the MATLAB -> pyLCP export workflow

1. **Always export the frequency vector together with the data** (as was done for this case with `andreata_frequency_range.mat`) — this avoids recurrence of Bug 2 in new validation cases.
2. **Document/standardize the conductor ordering convention** used in the MATLAB exports (by type vs. by cable), or export directly in the same order that pyLCP uses internally (grouped by cable), eliminating the need for the reindexing permutation.
3. Bugs 1 and 3 (key/file names) suggest it may be worth an automatic sanity check in `MatlabDataReader` or in `andreata_case1.py` that warns (or fails loudly, instead of silently returning `None`) when an expected key is not found — the use of `.get()` without a check masked the first 3 bugs for a long time.
