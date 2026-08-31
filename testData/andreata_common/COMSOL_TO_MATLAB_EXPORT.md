# COMSOL (FEM) → MATLAB reference export — internal impedance / admittance

**Context:** transversal to `andreata_case1` (and reusable by `case2/3/4`). The
external developer who produces the MATLAB reference `.mat` files
(`<prefix>_internal_impedance_matrix.mat` / `..._admittance_matrix.mat`,
variables `Z` / `Y`) needs the **COMSOL/FEM** internal matrices delivered back
**in that same file format**, so a methodology comparison (analytical MATLAB ×
FEM) can be run by swapping one array for another, with zero impact on notation
or on the ordering of the internal elements.

This is the inverse direction of `MatlabDataReader` (which reads the developer's
`.mat` into pyLCP). The writer lives in the same module:
`utils/matlab_data.py` → `export_comsol_internal_matrices_to_mat()`.

---

## 1. What the method does

```python
from utils.matlab_data import export_comsol_internal_matrices_to_mat

export_comsol_internal_matrices_to_mat(
    r'testData/andreata_case1/andreata_case1.py',   # __file__ of the case (locates Results/)
    prefix='andreata_case1',
)
```

Pipeline (all steps are **structural** — no value is recomputed):

1. **Read** the COMSOL `.txt` results through `ComsolPostProcessor`
   (`utils/comsol_data.py`), which already resolves every column-naming /
   excitation-block subtlety:
   - `cmsl_internal_impedance_matrix.txt` → `get_scc_internal_impedance_matrix_combined()`
     → per-cable **2×2** `Zi` (index 0 = core, 1 = sheath), Js method;
   - `cmsl_internal_admittance_charge_method.txt` → `get_scc_internal_capacitance_matrix_combined()`
     → per-cable **2×2** nodal `C` (direct charge method).
2. **Expand** the 2×2 per-cable block into the full **6×6** system matrix
   (`_cable_block_to_type_ordered_system`): block-diagonal replication for the 3
   identical cables (`kron(I₃, block)` → pyLCP "by cable" order), then the
   inverse of `conductor_order` to land in the **MATLAB "by type" order**.
3. **Compatibility check** of the frequency grids (see §3). No interpolation.
4. **Admittance:** `Y = jω·C` with `ω` taken from the COMSOL grid itself; real
   part is exactly `+0`.
5. **Write** three `.mat` files with `scipy.io.savemat`, byte-structurally
   identical to the reference (see §2).

---

## 2. Output files (in `testData/<case>/Results/`, unless `output_dir=` is given)

| File | Variable | Structure — identical to the reference |
|---|---|---|
| `<prefix>_internal_impedance_matrix_fem.mat` | `Z` | `(6, 6, Nf)`, `complex128` |
| `<prefix>_internal_admittance_matrix_fem.mat` | `Y` | `(6, 6, Nf)`, `complex128` |
| `<prefix>_frequency_range_fem.mat` | `freq1` | `(1, Nf)`, `float64` (row vector) |

- The variable name of the frequency file is copied from the reference
  `<prefix>_frequency_range.mat` (fallback `'freq1'`).
- Only the MAT-file header metadata differs (`Platform: nt`, creation date) —
  irrelevant to loading in MATLAB.

### Conductor ordering — "by type" (the MATLAB convention)

```
index:   0      1      2      3        4        5
       core_A core_B core_C sheath_A sheath_B sheath_C
```

i.e. `Z[0:3, 0:3]` = cores, `Z[3:6, 3:6]` = sheaths, and the core↔sheath
coupling of cable *k* sits at `Z[k, k+3] = Z[k+3, k]`. This is the **inverse**
of the pyLCP "by cable" layout `[cA, sA, cB, sB, cC, sC]` produced by
`np.kron(np.identity(N), Zij)` in `InternalPerUnitParameters.matrices()`.

The permutation is the same `conductor_order` argument used by
`MatlabDataReader.get_scc_scenario_data` — default `(0, 3, 1, 4, 2, 5)` — see
`utils/matlab_data.py::_reorder_conductor_matrix` and
`andreata_case1/BUGS_AND_FIXES.md` (Bug 4). The developer can therefore drop the
`_fem` arrays straight in, with the same `(p, q)` indices as their own data.

---

## 3. Frequency grid — the "verbatim / no interpolation" contract

For a **methodology comparison** the FEM data must be preserved unrestrictedly:
no resampling, no smoothing, no recomputation. An earlier version resampled the
COMSOL sweep onto the reference grid with a PCHIP interpolation in log-frequency
— **this was removed on purpose**. Consequently the frequency grids must line up
on their own.

Two checks (`_assert_frequency_compatible`, `rtol = 1e-6`):

| # | Check | Enforcement |
|---|---|---|
| 1 | COMSOL impedance grid **==** COMSOL admittance grid | **always hard** (`ValueError`) — the two are assembled into one `(Z, Y)` pair |
| 2 | COMSOL grid **==** reference grid (`target_frequencies` or `<prefix>_frequency_range.mat`) | gated by `strict_reference` |

### `strict_reference` (default `True`)

- `True` — a mismatch in check 2 raises `ValueError`.
- `False` — **emergency mode**: the mismatch is downgraded to a warning and the
  export proceeds **on the COMSOL grid**. `<prefix>_frequency_range_fem.mat`
  then documents the real grid.

### The 90 vs 91 standardization issue (andreata_case1, 2026-08)

The standard log sweep **1 × 10⁻² … 1 × 10⁷ Hz at 10 points/decade** has

```
9 decades × 10 + 1 = 91 points   ->   np.logspace(-2, 7, 91)
```

The developer's current reference `.mat` files carry **90** points
(`np.logspace(-2, 7, 90)`, 89 intervals) — **off-standard**. The revised COMSOL
`.txt` files are correct at **91** points. So:

- `strict_reference=True` → aborts (`91 vs 90`), by design;
- `strict_reference=False` → exports `Z`/`Y` as `(6, 6, 91)` + a 91-point
  `frequency_range_fem.mat`, verbatim.

**Emergency command in use until the reference `.mat` is regenerated:**

```python
export_comsol_internal_matrices_to_mat(
    r'testData/andreata_case1/andreata_case1.py', prefix='andreata_case1',
    strict_reference=False)
```

**When the external developer regenerates the reference on the 91-point grid:**
drop `strict_reference=False`. The check goes back to strict, shapes match 1:1,
and the DC-diagonal numeric check (below) runs again.

---

## 4. Pragmatic details

### Signature

```python
export_comsol_internal_matrices_to_mat(
    script_file_path,                       # __file__ of the case script
    prefix,                                 # 'andreata_case1'
    conductor_order=(0, 3, 1, 4, 2, 5),     # by-type -> by-cable permutation
    target_frequencies=None,                # override the reference grid (array of Hz)
    output_dir=None,                        # default: the case Results/ dir
    reference_check=True,                   # print FEM-vs-reference DC-diagonal diff
    strict_reference=True,                  # False = emergency export on the COMSOL grid
) -> dict
```

Returns
`{'impedance': {'path', 'Z', 'frequencies'},
  'admittance': {'path', 'Y', 'frequencies'},
  'frequency_range': {'path', 'frequencies'}}`.

### `reference_check` output

When the reference `.mat` files exist **and have the same shape** as the output,
prints the relative difference of the **DC diagonal** (index 0), per conductor —
a sanity check, not a pass/fail. Expected order of magnitude (andreata_case1,
FEM × analytical reference):

| | max rel. diff (DC diagonal) |
|---|---|
| `Z` | ~0.5 % (cores ~0.05 %, sheaths ~0.5 %) |
| `Y` | ~0.005 % |

If the shapes differ (the 90/91 situation) the check prints a `Note:` and is
skipped.

### Errors

| Exception | Cause |
|---|---|
| `FileNotFoundError` | COMSOL `.txt` inputs missing for the case |
| `ValueError` | the two COMSOL inputs are on different grids; **or** (unless `strict_reference=False`) the COMSOL grid ≠ the reference grid |

Nothing is written when an exception is raised — the checks precede every
`savemat`.

### What to send / tell the external developer

1. The three `*_fem.mat` files.
2. That `Z` / `Y` follow **their** structure exactly (shape, dtype, variable
   name, by-type conductor order) — no re-indexing needed.
3. That the values are the **raw COMSOL FEM output** (no interpolation), on the
   standard **91-point** grid described by `<prefix>_frequency_range_fem.mat`,
   and that their own reference should move to that grid.

---

## 5. Status — andreata_case1 (2026-08-31)

- COMSOL `.txt` (impedance + admittance): **91 points**, standard grid ✔
- Reference `.mat`: **90 points** — awaiting regeneration by the external developer
- Export produced with `strict_reference=False` → `(6, 6, 91)` `Z`/`Y` +
  91-point `frequency_range_fem.mat`, values verbatim, by-type order,
  symmetric, same sparsity pattern as the reference.

---

## 6. Related — the only remaining numerical interpolation in pyLCP

`analytical_forms/single_core_cable.py::InternalParametersFromFEM._interp_zi`
(PCHIP in log-frequency, real/imag separately, linear `np.interp` fallback
outside the FEM range) — used **only** by the FEM-hybrid pipeline of
`andreata_case2` (Config. 2), to bring a coarse FEM `Zi` onto the analytical
frequency grid. It is unrelated to this export path and untouched by it. See
`HYBRID_PIPELINE_PLAN.md`.
