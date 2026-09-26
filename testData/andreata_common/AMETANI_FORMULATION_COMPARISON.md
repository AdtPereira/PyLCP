# pyLCP vs. Ametani (2015) — SC and PT cable formulation

**Date:** 2026-09.
**Reference:** Ametani, Ohno & Nagaoka, *Cable System Transients* (Wiley,
2015), Ch. 2 — sec. 2.1 (single-core coaxial, SC) and sec. 2.2
(pipe-enclosed, PT).
**Outcome:** the SC formulation is implemented as in the book. The PT
earth-return structure (eqs. 2.32/2.40) settles the ECC question left open in
[`GROUND_RETURN_ALLOCATION_STUDY.md`](GROUND_RETURN_ALLOCATION_STUDY.md), and
is now implemented. It brings Config. 4 `Z'` from 50% to 0.27% vs. MATLAB.

---

## 1. Sec. 2.1 — SC cable

| Ametani | pyLCP | Status |
|---|---|---|
| (2.5) `[Z] = [Zi] + [Z0]` | `Zs = Zi + Zg` — `PerUnitParameters.quasi_tem_approx_matrices` | ✓ |
| (2.6) `[Zi]` block-diagonal | `kron(I_N, Zij)` / `_matrices_heterogeneous` | ✓ |
| (2.8)/(2.9) core + sheath + armor | `InternalPerUnitParameters._build_cable_block`, branch 1 | ✓ |
| (2.11)/(2.12) core + sheath | branch 3 | ✓ |
| (2.13)/(2.24) core only (insulated ECC) | branch 5: `z11 + z12`, `pcj` | ✓ |
| `z11, z12, z2i, z2m, z20, z23, z3i, z3m, z30, z34`; `D1..D3` | `parameters_by_bessel` | ✓ — the book's `z3i` shows `μ2` (typo); pyLCP correctly uses `μ3` |
| (2.21)/(2.23) `[Pij]` | branches 1/3 | ✓ |
| (2.14)–(2.16) `[Z0jk]` tiled over every conductor pair of cables j, k | `_expand_by_block_sizes` | ✓ |
| `Z0jk` — Pollaczek | closed-form De Conti/Duarte/Alipio 2023 (`'deconti'`) or Sommerfeld integrals | intentional (Andreata secs. 5.4/6.1) |
| (2.18) underground `[P] = [Pi]` | `[P] = [Pi] + [Pg]` (De Conti eq. 4.63) | intentional extension — Ametani's assumption 1 neglects displacement current in the soil |

**Fixed (latent defects, prerequisite for Config. 5 — bare copper ECC):**

- Branch 6 (bare core / bare ECC) never set `Pij` → `UnboundLocalError`.
  Now `Pij = [[0]]` (eq. 2.24 with `pcj = 0`).
- Branch 4 (sheath without jacket) returned a 1×1 `Pij` with a 2×2 `Zij`.
  Now `[[pcj, 0], [0, 0]]` (eq. 2.23 with `psj = 0`).
- A zero row makes `Pi` singular: `_internal_capacitance` returns NaN for the
  inspection-only `Ci`/`Ye`; the composed `Ysh = jw (Pi + Pg)^-1` stays finite.

## 2. Sec. 2.2 — PT cable

- **`Zp`, `Pp`, `Q_jk` (eqs. 2.34, 2.35, 2.49)** assume a *conducting* pipe
  (`ρp`, `μp`; assumption 3: wall thicker than the penetration depth). The
  HDPE duct is an insulator, so they do not apply. The interior of the duct
  (eccentric cable + air + HDPE) comes from the FEM (`'fem'` scenarios) or the
  Lafaia GMD (analytical scenarios). This is not a gap in pyLCP.
- **`[Z0]` (eqs. 2.32/2.40)** — *"Z0 … is the self earth-return impedance of
  the pipe"*, repeated over **every** entry of the block of conductors inside
  it (self and mutual alike). The soil only sees the pipe.
  - Config. 2: one object per duct, self radius `D2/2`, depth of the duct
    axis (`_cable_external_geometry`). ✓
  - Config. 4 (before): the ECC inside duct C was a **4th object** with its
    own radius (0.0053 m) and position. ✗
  - Config. 4 (now): 3 objects, `ground_return_block_sizes = [2, 2, 3]`. ✓

## 3. Measured effect (Config. 4, `'fem'` scenario, vs. MATLAB)

| variant | `Z'` ECC-ECC | `Z'` C-ECC | `Z'` full | `Y'` full |
|---|---|---|---|---|
| (a) previous: 4 objects | 50.0% | 7.5% | 50.0% | 2.6% |
| (b) "correction 1": ECC inherits duct C radius, still its own object | 2.7% | 18.9% | 14.1% | 2.5% |
| **(c) Ametani 2.32/2.40: 3 objects** (implemented) | **2.7%** | **0.2%** | **0.27%**\* | **1.84%**\* |

\* inf-norm of the case's own `_validate_fem_hybrid_vs_matlab`, after the
implementation (the scratch measurement gave 2.7% / 2.3% on a max-entry
metric).

Variant (b) fixes the ECC self term but breaks the C–ECC mutual, because the
ECC keeps its own horizontal position. It also makes the ECC mode worse
(`v_m` 266% max error). Only the grouping of (c) reproduces the reference.

Modal domain (1 Hz–1 MHz max error, α / v / |Zc|): ECC mode from
117% / 28% / 346% (a) to 23% / 7% / 93% (c); spot values at 1e2/1e4/1e6 Hz are
now ≤ 2% in `v_m` for all 7 modes (see `MODAL_CH5_DEVELOPMENT.md` sec. 7.4).

## 4. Implementation

| Where | What |
|---|---|
| `mtl_main/strategy.py::_merge_groups_inside_enclosures` | Merges every duct-less conductor group whose full cross-section lies inside another group's duct (inner radius) into that group. It is geometric, so the model needs no extra field (Config. 4 ECC: 0.0441 m from the duct axis + 0.0053 m radius ≤ 0.0508 m). |
| `SingleCoreCableWithECCStrategy._cable_distance_matrices` | Uses it; the duct host is the object's representative; it checks that each object's conductors are contiguous in conductor order and exposes `ground_return_block_sizes`. |
| `PerUnitParameters.earth_return_parameters` | `N` = number of ground-return objects (`d_matrix.shape[0]`), not `num_sc_cables`. |
| `PerUnitParameters.quasi_tem_approx_matrices` | `model.ground_return_block_sizes` takes precedence over the internal `block_sizes`; checked against the conductor count. |
| `andreata_case4.py` (`'fem'`) | The explicit `block_sizes=[2, 2, 2, 1]` is removed. |

Unaffected: Configs 1/2 (other strategies), Config. 3 and the case4
analytical scenarios `1/2/3/3_deconti` (no duct in those models → nothing to
merge).
