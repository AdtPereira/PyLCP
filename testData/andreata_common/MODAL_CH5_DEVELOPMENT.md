# Development — Modal-Domain Propagation Characteristics (Andreata Ch. 5)

**Date:** 2026-08 (core); updated 2026-09 (MATLAB comparison, Configs 3/4 — see sec. 7).
**Scope:** all four Andreata cases. Core work (secs. 1–6) covers
Configurations **1 and 2** — three single-core cables (core + sheath), flat
arrangement, 34.5 kV; **directly buried** (Config. 1) or **in individual
HDPE ducts** (Config. 2), 6 conductors => 6 modes in both. Sec. 7 extends
this to the MATLAB numerical comparison (Configs 1/2) and to preparing
Configurations **3 and 4** (3 SCC + ECC, 7 conductors => 7 modes;
directly buried / in individual HDPE ducts respectively). Goal: from
`Z'(f)` and `Y'(f)` (already produced by `andreata_case1` / `andreata_case2`
/ `andreata_case3` / `andreata_case4`), compute and plot **`alpha_m`, `v_m`,
`|Z_cm|`** per mode, reproducing **Figures 5.5–5.7** (Config. 1) and
**5.8–5.10** (Config. 2).

The full plan (incl. future phases: Configs 2–5, frequency-dependent
soil) is in [`CH5_PROPAGATION_PLAN.md`](CH5_PROPAGATION_PLAN.md).

---

## 1. Files

### Created

| File | Content |
|---|---|
| `mtl_main/propagation.py` | `phase_domain_propagation(Zs, Ysh)` -> `gamma_v`, `gamma_i`, `Z_c`, `Y_c` in the phase domain (via `scipy.linalg.sqrtm`). Factored out of the loop that existed only in `overhead_lines.py`. |
| `analytical_forms/modal_analysis.py` | `ModalDecomposition`: `decompose()` -> `track()` -> `modal_parameters()` -> `_classify()`. Helpers `default_scc_roles`, `canonical_permutation`, `REFERENCE_PATTERNS_6C`. |
| `plotter/modal_plotter.py` | `ModalPropagationPlotter` — Figs 5.5 / 5.6 / 5.7. `MODE_STYLE` (fixed color/dash per mode label). |
| `utils/passivity_check.py` | `check_pul_passivity`, `print_passivity_report`, `min_real_part_eigenvalue` — passivity **assessment** (not enforcement). |

### Modified

| File | Change |
|---|---|
| `analytical_forms/overhead_lines.py` | `PerUnitParameters.pul_matrices` now calls `phase_domain_propagation` (identical result; regression checked in `ohtl_deConti_ex51`). |
| `analytical_forms/single_core_cable.py` | New `PerUnitParameters.propagation_matrices(quasi_tem_matrices)` — appends `propagation_voltage_matrix`, `propagation_current_matrix`, `characteristic_impedance_matrix`, `characteristic_admittance_matrix` to the `quasi_tem_approx_matrices` dict (backward compatible). |
| `testData/andreata_case1/andreata_case1.py` | Calls `propagation_matrices`, runs `check_pul_passivity` + report, `ModalDecomposition` on the **`p100_er1_deconti`** scenario (ground return via closed-form De Conti/Duarte/Alipio 2023 — eqs. 4.59/4.63, sec. 5.4 / 6.1), and `ModalPropagationPlotter.plot_all()` — Figs 5.5 / 5.6 / 5.7. |
| `testData/andreata_case2/andreata_case2.py` | Same on the **`fem`** scenario (FEM-hybrid pipeline — see [`HYBRID_PIPELINE_PLAN.md`](HYBRID_PIPELINE_PLAN.md); `Zi`/`Yi` from COMSOL + analytical ground return with the pipe radius, **without** Lafaia's GMD). Figs 5.8 / 5.9 / 5.10 + 4 comparison plots vs. MATLAB (`self_impedance_phase_a_sheath`, `self_admittance_phase_a_sheath`, `earth_return_impedance_phase_a`, `earth_return_potential_coeff_phase_a`). |
| *(2026-09, sec. 7)* `utils/matlab_data.py` | `get_modal_scenario_data` (MATLAB modal reference reader, mode count auto-detected). |
| *(2026-09, sec. 7)* `plotter/modal_plotter.py` | `matlab_modal` overlay on `ModalPropagationPlotter`; `print_modal_comparison_report`. |
| *(2026-09, sec. 7)* `analytical_forms/modal_analysis.py` | `default_scc_roles(..., num_ecc=0)`. |
| *(2026-09, sec. 7)* `analytical_forms/single_core_cable.py` | `InternalParametersFromFEM.from_component_blocks` (heterogeneous FEM assembly). |
| *(2026-09, sec. 7)* `testData/andreata_case3/andreata_case3.py`, `testData/andreata_case4/andreata_case4.py` | Modal decomposition wired (7 conductors); case4 additionally gets the `'fem'` FEM-hybrid scenario. |
| *(2026-09, sec. 7)* `testData/andreata_common/GROUND_RETURN_ALLOCATION_STUDY.md` | New — ground-return radius allocation study (Case 4 ECC gap). |

---

## 2. Basis (Ch. 5 equations, N = 6)

```
input:  Z'(f), Y'(f)  in C^{Nf x 6 x 6}     (canonical order [c1 c2 c3 s1 s2 s3])

1. eigendecomposition    lambda, T_I  =  eig(Y' Z')                          (eq. 5.19)
2. voltage matrix        T_V  =  (T_I^T)^-1   =>  T_V^-1 = T_I^T  (exact)    (eq. 5.20)
                         =>  Z' Y' = T_V . diag(lambda) . T_V^-1             (eq. 5.18)
3. tracking              reorders/re-phases columns of T_I(f_k) to match
                         T_I(f_{k-1})  (switching-back — see sec. 3.5)
4. modal series imped.   Z_m,j  =  [T_I^T Z' T_I]_jj                         (eq. 5.21, diagonal)
5. propagation const.    Gamma_m  =  sqrt(lambda) ,  physical branch  alpha_m = |Re|, beta_m = |Im|   (eqs. 5.24, 5.30–5.31)
6. phase velocity        v_m  =  omega / beta_m                             (eq. 5.32)
7. characteristic imped. Z_cm  =  Z_m / Gamma_m     Y_cm = 1 / Z_cm         (eqs. 5.23, 5.27)
8. mode label            projects columns of T_I(f_min) onto the patterns   (sec. 5.3, eqs. 5.34–5.36)
```

**The 6 modes of Configuration 1** (labels used in the code):

| # | label | description (Andreata sec. 5.3) | vs f |
|---|---|---|---|
| 1 | `ground` | 3 sheaths in phase, ground return | ~invariant |
| 2 | `inter_sheath_1` | 1 sheath x another sheath | ~invariant |
| 3 | `inter_sheath_2` | 1 sheath x the other two | ~invariant |
| 4 | `coaxial_1` | 3 cores x earth (LF) -> core x sheath (HF) | varies |
| 5 | `coaxial_2` | core x core (LF) -> core x sheath (HF) | varies |
| 6 | `coaxial_3` | core x 2 cores (LF) -> core x sheath (HF) | varies |

---

## 3. Design decisions (the non-obvious ones)

### 3.1 Factoring out the phase propagation (`mtl_main/propagation.py`)

The loop `sqrtm(Z'Y')` / `sqrtm(Y'Z')` / `Z_c = Y'^-1*gamma_i` / `Y_c = Z'^-1*gamma_v` existed
**only** in `analytical_forms/overhead_lines.py`; the SCC path stopped at
`Z'/Y'`. Extracted into a single module, reused by OHTL and SCC. `overhead_lines.py`
now calls it with no change of result (same LU + `sqrtm` sequence).

### 3.2 `T_V = (T_I^T)^-1` instead of a second `eig`

`Z'Y' = (Y'Z')^T` (since `Z'`, `Y'` are symmetric). Therefore, if `Y'Z' = T_I Lambda T_I^-1`,
then `Z'Y' = (T_I^-1)^T Lambda T_I^T`, i.e. `T_V = (T_I^T)^-1`. This:
- makes eq. 5.20 (`T_V^-1 = T_I^T`) **exact by construction** — with no tolerance;
- avoids a second `eig(Z'Y')` and the problem of **matching** the eigenvalues/eigenvectors
  of the two decompositions.

### 3.3 `Z_m` via `diag(T_I^T Z' T_I)`, **not** `Y_m`

Theorem (Wedepohl): `Z' t_j` is parallel to the j-th column of `T_V`, therefore
`Z_m = T_V^-1 Z' T_I = T_I^T Z' T_I` is **diagonal** in exact arithmetic.
`Y_m = T_I^-1 Y' T_V` would also be, but it depends on the **left eigenvectors**
(`T_I^-1`) and is ill-conditioned near degeneracy. Measured in
`andreata_case1`: off-diagonal of `Z_m` ~ 1e-10 (~ 0), off-diagonal of `Y_m`
reached ~38 %. For that reason `Z_cm = Z_m / Gamma_m` (and `Y_cm = 1/Z_cm`, the modal
domain is a decoupled scalar) — `Y_m` is returned only for inspection.

### 3.4 `Gamma_m = sqrt(lambda)` with physical branch `alpha = |Re|`, `beta = |Im|`

For a wave that propagates and attenuates, `Gamma = alpha + j*beta` with `alpha, beta > 0` =>
`Im(lambda) = 2*alpha*beta > 0` => the principal `sqrt(lambda)` already gives `Re > 0` and `Im > 0`. But near DC
the numerical noise can give `Im(lambda) < 0`, and `sqrt` of the plane near the negative real
axis is hypersensitive to the sign of `Im` (`sqrt(-1-0j) = -j`), generating `v_m < 0`.
Solution: `alpha_m = |Re sqrt(lambda)|`, `beta_m = |Im sqrt(lambda)|` (alpha and beta are magnitudes). Documented in
the code.

### 3.5 Tracking = *switching-back procedure*

`ModalDecomposition.track()` implements the same sub-routine that Andreata (sec. 5.2)
reused: **Gustavsen (2008) sec. IV-A**, which is the scalar-product correlation of
**Wedepohl-Nguyen-Irwin (1996) sec. 6** (the `intercheig` function of the *Matrix Fitting
Toolbox*). At each step it matches the columns of `T_I(f_k)` to those of `T_I(f_{k-1})`
by maximizing `Sum_j |<t_j(f_k), t_j(f_{k-1})>|`.

Differences from the original `intercheig` (improvements):
- assignment resolved **optimally** (`scipy.optimize.linear_sum_assignment`,
  Hungarian) instead of the greedy row-by-row maximum — more robust when several
  modes swap position at the same frequency (the HF coaxial case);
- **anchor at the middle of the band**, sweeping both ways, instead of sweeping
  from an endpoint — in Config. 1 the coaxial modes are more separated in the
  central decades.

The Wedepohl (1996) Newton-Raphson method was considered and **not adopted**:
it produces a *smooth* `T_I(omega)` fittable by rational functions, which only matters
for line models **in the modal domain** (EMTP). pyLCP + NLT (Andreata Ch. 6) is
effectively **phase domain**; to plot `alpha_m`, `v_m`, `|Z_cm|` a
consistent mode ordering is enough, which the correlation tracking
delivers. Section 6 of Wedepohl's own paper uses that correlation method
as validation of NR ("agree perfectly").

### 3.6 Refinement of degenerate eigenvalue clusters

`decompose(cluster_rtol=3e-2)`: within a cluster of eigenvalues that
agree to that relative tolerance, `np.linalg.eig` returns an **arbitrary**
basis of the eigenspace. The basis is re-chosen by diagonalizing the restricted
modal series impedance `B = U^T Z' U` (a complex symmetric matrix => its
eigenvectors satisfy `v_i^T v_j = 0`, the same property as eq. 7 of Gustavsen
2008). Without this, the 3 nearly-identical coaxial modes of Config. 1 generated
oscillation in `|Z_cm|` in the last half-decade (5e5–1e7 Hz).

### 3.7 Canonical conductor order

pyLCP assembles `Z_i = kron(I_3, Z_ij)` => rows `[c1 s1 c2 s2 c3 s3]`; Andreata
uses `[c1 c2 c3 s1 s2 s3]`. `canonical_permutation(default_scc_roles(...))`
converts (permutation `[0 2 4 1 3 5]` for the 3xSCC case). Passivity and
eigenvalues are invariant under that symmetric permutation.

### 3.8 2-stage mode classification (`_classify`)

The first version projected the columns of `T_I` **at the minimum frequency** onto
the reference patterns (eqs. 5.34–5.36). It worked for Config. 1
(similarity 0.99) but **inverted** the labels in Config. 2: with the duct, at
LF the "sheath" modes carry appreciable core current and the modes
only separate above ~100 Hz.

New labeling, robust for both configs (similarity = 1.000):

1. **separates** the 3 sheath modes from the 3 coaxial ones by the **fraction of
   `T_I` energy in the core rows**, averaged over a 100 Hz – 100 kHz band (where
   every config is decoupled). The ground/inter-sheath modes are
   frequency-invariant and have no core participation (Andreata
   sec. 5.3); the 3 with the smallest fraction => sheath group, the other 3 => coaxial group;
2. **sub-labels** each group of 3 by the corresponding row pattern
   (`[1,1,1]` = `ground`/`coaxial_1`; `[1,0,-1]` = `inter_sheath_1`/`coaxial_2`;
   `[-1,2,-1]` = `inter_sheath_2`/`coaxial_3`), via `linear_sum_assignment`,
   at an average frequency of the band.

Extra diagnostics: `classification_ref_freq_hz`, `classification_core_fraction`.

### 3.9 Configuration 2 — HDPE duct (FEM-hybrid pipeline)

The 6-mode structure is identical (Andreata eqs. 5.37–5.39 = 5.34–5.36).
The modal decomposition runs on the **`fem` scenario** of `andreata_case2`
(see [`HYBRID_PIPELINE_PLAN.md`](HYBRID_PIPELINE_PLAN.md)):

- **internal** (`Zi`, `Yi`): from COMSOL, exact eccentric geometry with air + HDPE
  pipe (`InternalParametersFromFEM`). Without Lafaia's equivalent-permittivity GMD
  trick.
- **ground return** (`Zg`, `Yg`): **closed-form expressions of De Conti/Duarte/
  Alipio 2023** (`zg_form = yg_form = 'deconti'`, eqs. 4.59 and 4.63 of Andreata —
  sec. 5.4 / 6.1), **not** the Sommerfeld integrals (`magalhaes_xue`), with the self
  term using the **pipe outer radius** (`D2/2`), not the SCC's
  (`_cable_external_geometry` in `mtl_main/strategy.py`).

Physical effect: the return current of the sheaths sees air + HDPE pipe
(large dielectric, low `eps_r`) => `Z_cm` and `v_m` of the ground/inter-sheath modes
are **much larger** than in Config. 1 (Andreata sec. 5.4.1, Figs 5.8–5.10). The
FEM-hybrid path reproduces `Z'/Y'` of Andreata's FEM reference (`.mat`) to **< 2 %**
(vs. ~10–20 % of the GMD), and the HF oscillation of the coaxial modes disappears.

### 3.10 Passivity diagnosis (`utils/passivity_check.py`)

A passive line => `eig(Re{Z'})>=0` (R' PSD), `eig(Re{Y'})>=0` (G' PSD),
`eig(Re{Y_c})>=0` (Gustavsen 2008, eq. 3) over the whole range. A negative eigenvalue
beyond numerical noise points to a **formulation problem** (ground return
with negative conductance, wrong semiconducting-layer `eps_r`, a sign flipped
in an approximation). The module **only reports** — fixing it would be the
FRP/FMP + Hamiltonian pipeline of Gustavsen (2008), deliberately out of scope:
here `Z'/Y'` are closed-form samples, not a fitted rational macromodel.
The tolerance is relative to the per-frequency scale of `|Re{M}|` (it does not flag `G' ~ 0` of a
lossless dielectric).

---

## 4. Usage

```python
from analytical_forms.modal_analysis import ModalDecomposition, default_scc_roles
from plotter.modal_plotter import ModalPropagationPlotter
from utils.passivity_check import check_pul_passivity, print_passivity_report

# qt = PerUnitParameters(...).quasi_tem_approx_matrices(...)
# qt = PerUnitParameters(...).propagation_matrices(qt)     # optional: gamma_v, Z_c, ...

print_passivity_report(check_pul_passivity(freq, qt))

md = ModalDecomposition(
    freq, qt['series_impedance_matrix'], qt['shunt_admittance_matrix'],
    conductor_roles=default_scc_roles(mtl.num_sc_cables, mtl.num_conductors_per_scc))
res = md.modal_parameters()           # dict: alpha, beta, vphase, Zcm, Ycm, TI, TV,
                                      #       mode_labels, diagnostics, ...
ModalPropagationPlotter(__file__, res, config_name='Configuration 1').plot_all()
```

`res['diagnostics']`: `max_offdiag_ratio_Zm`, `max_scalar_relation_error`,
`reciprocity`, `passivity`, `classification_similarity`.

---

## 5. Validation (constant soil rho = 100 Ohm.m; ground return = closed-form De Conti)

Diagnostics of the two configs:

```
                             Config. 1 (p100_er1_deconti)  Config. 2 (scenario 'fem', FEM-hybrid)
mode labels                  (ground, inter_sheath_1/2, coaxial_1/2/3)   (same)
max off-diag ratio Zm        ~1e-10                        ~1e-10           (~ 0 -> Z_m diagonal)
classification similarity    1.000                         1.000            (min over the 6 modes)
reconstruction T_I Lambda T_I^-1     < 1e-16               < 1e-16
passivity (Z', Y', Yc)       PASSIVE                       PASSIVE
```

*Config. 2 — `Z'`/`Y'` of the `'fem'` scenario vs. MATLAB (Andreata's FEM):*
`Z'` 0.28 % · `Y'` 1.78 % · `Zg` 0.13 % · `Pg` 0.34 %.

**Comparison with the PDF** (values at 1e6 Hz):

*Config. 1 — Figs 5.5–5.7:*

| quantity | ground mode | inter-sheath | coaxial | Andreata |
|---|---|---|---|---|
| `alpha_m` [Np/m] | 0.13 | 0.02–0.04 | 1.3e-4 | ground >> sheaths >> coaxial; crossover ~1e2 Hz OK |
| `v_m` [m/s] | 3.1e7 | 3.5–4.1e7 | 1.8e8 | OK |
| `\|Z_cm\|` [Ohm] | ~50–60 | ~16–20 | ~13 | ground ~40–60, sheaths ~15–25, coaxial ~13 OK |

*Config. 2 — Figs 5.8–5.10:*

| quantity | ground mode | inter-sheath | coaxial | Andreata |
|---|---|---|---|---|
| `alpha_m` [Np/m] | ~0.1 | ~0.02–0.04 | ~2e-4 | ground >> sheaths >> coaxial OK |
| `v_m` [m/s] | ~5e7 | **1.2–1.45e8** | ~1.8e8 | inter-sheath **much** faster than in Config. 1 OK |
| `\|Z_cm\|` [Ohm] | ~130–200 | **~55–70** | ~13 | ground ~130–150, sheaths ~60–70, coaxial ~13 OK |

The Config. 2 signature (Andreata sec. 5.4.1) — `Z_cm` and `v_m` of the
ground/inter-sheath modes **much larger** than in Config. 1 because of the air
+ HDPE pipe insulation — is reproduced.

**Validation `Z'`/`Y'`/`Zg`/`Pg` vs. MATLAB** (Andreata's FEM reference, `andreata_case2`
`fem` scenario — 4 comparison plots in `Results/`):

| quantity | max rel. error |
|---|---|
| `Z'` full series | 0.40 % |
| `Y'` full shunt | 1.85 % |
| `Zg` ground return | 0.07 % |
| `Pg` potential coeff. | 0.81 % |

---

## 6. Known limitations

1. **Constant soil** (rho = 100 Ohm.m) instead of the FD model of Salvador (2020) —
   affects the ground/inter-sheath modes at high f. Phase 3 (the
   `Z'/Y'` vs. MATLAB comparison above, < 2 %, suggests the effect is smaller than
   previously estimated).
2. **Cosmetic oscillation** of the 3 coaxial modes in Config. 1 (grid up to 10 MHz):
   genuine degeneracy — the cluster refinement (sec. 3.6) mitigates it. **In Config. 2
   with the FEM `Zi` the oscillation disappears** (the FEM `Zi` has the real external
   inductance, smooth). Next step for Config. 1: anchor the tracking on the
   asymptotic limit matrices of Wedepohl (1996) sec. 7.
3. **~2 % offset in `C_22`** of Config. 2: difference between the COMSOL FEM
   (used by pyLCP) and Andreata's own FEM — it is not a pipeline error.
4. **Configs 3/4 wired but not classified**: `ModalDecomposition` runs on both
   (7 conductors / 7 modes) and produces `alpha_m`/`v_m`/`|Z_cm|`, but
   `_classify` still bails out for `n != 6` -- modes are unlabelled
   (`mode_0`..`mode_6`). No thesis equations for a 7th "ECC mode" are in the
   codebase yet. See sec. 7.2. **Config 5 remains unimplemented.**
5. ~~`ModalPropagationPlotter` is standalone...~~ **Done (sec. 7.1)**: it now
   accepts a `matlab_modal` overlay (Configs 1/2, validated); `print_modal_
   comparison_report` gives the numeric per-mode comparison. Still standalone
   from `PLOT_CONFIG` by design (its own axis/legend logic).
6. **FEM-hybrid Case 4 `Z'` ECC self-term ~50% off** at high frequency (`Y'`,
   `Zg`, `Pg` all agree to <3%): the ground-return radius assigned to the ECC
   object does not account for it sharing phase C's duct. Root-caused, not
   yet fixed pending team confirmation -- see sec. 7.3 and
   [`GROUND_RETURN_ALLOCATION_STUDY.md`](GROUND_RETURN_ALLOCATION_STUDY.md).

---

## 7. Update (2026-09) — MATLAB comparison, Configs 3/4 prepared, FEM-hybrid Case 4

### 7.1 Numerical MATLAB overlay for the modal parameters (Configs 1/2)

The external MATLAB developer sent `alpha_m`, `beta_m`, `v_m`, `Z_cm`, `Y_cm`
per mode (6 modes, standard 91-point grid) for Configs 1 and 2 -- previously
only `Z'`/`Y'` had been cross-checked, never the modal-domain output itself.

| Delivered | File |
|---|---|
| `MatlabDataReader.get_modal_scenario_data(prefix, config_index)` -- reads `<prefix>_{alpham,betam,velocm,Zcm,Ycm}.mat`, auto-detects the mode count per file (not hardcoded to 6), assigns semantic labels (`MODAL_MODE_ORDER_6C`) only when that count is 6 | `utils/matlab_data.py` |
| `ModalPropagationPlotter(..., matlab_modal=...)` -- overlays the matching mode (joined **by label**, never by column index -- the two eigendecompositions do not share a frequency grid or eigenvector order) as open-circle markers on Figs 5.5–5.7/5.8–5.10 | `plotter/modal_plotter.py` |
| `print_modal_comparison_report(modal, matlab_modal)` -- per-mode relative error at 1e2/1e4/1e6 Hz for `alpha_m`, `v_m`, `|Z_cm|` | `plotter/modal_plotter.py` |
| Wired into `andreata_case1.py` / `andreata_case2.py` | `testData/andreata_case1/`, `testData/andreata_case2/` |

**Mode-order correspondence was verified, not assumed**: matched the
per-mode `|Z_cm|`/`alpha_m` magnitudes at 1 MHz against this document's own
sec. 5 validation table, confirming the MATLAB developer's `mode1..mode6`
follows the same thesis order used by `MODE_LABELS_6C` (ground,
inter_sheath_1/2, coaxial_1/2/3).

**Results** (pyLCP vs. MATLAB, both configs): excellent agreement for
`alpha_m`/`v_m` across the full 1e-2–1e7 Hz sweep (mostly <5%, often <1%
above 1 kHz). `|Z_cm|` agrees well up to ~1 MHz; above that the 3 coaxial
modes diverge sharply in MATLAB's own Config. 1 curve (up to ~380% at 1 MHz)
-- **the same near-degenerate-eigenvalue instability already documented in
sec. 6, limitation 2, now confirmed independently on the MATLAB side.**
Interestingly, in **Config. 2** the pyLCP coaxial `|Z_cm|` curves stay smooth
through 1e7 Hz (as sec. 6 already noted -- the FEM `Zi` removes the
oscillation) while **MATLAB's own Config. 2 reference still shows it**,
suggesting the instability tracks the *ground-truth internal impedance
model* (GMD-analytical vs. FEM) rather than being purely a tracking-algorithm
artifact -- worth raising with Alberto alongside the original question.

### 7.2 Configurations 3 and 4 prepared for modal calculation

Neither case had the modal-decomposition wiring at all before this update;
both are now 7-conductor (3 SCC + ECC) systems reusing the same
`ModalDecomposition`/`ModalPropagationPlotter` machinery as Configs 1/2.

| Delivered | File |
|---|---|
| `default_scc_roles(num_cables, conductors_per_cable, num_ecc=0)` -- appends `("ecc", i)` roles after the SCC block | `analytical_forms/modal_analysis.py` |
| `andreata_case3.py`: `propagation_matrices` + passivity check + `ModalDecomposition` on `p100_er1_deconti` (same De Conti ground return as Config 1) + `ModalPropagationPlotter`/`print_modal_comparison_report` wired (MATLAB reference not sent yet -- both no-op gracefully) | `testData/andreata_case3/andreata_case3.py` |
| `andreata_case4.py`: new `'3_deconti'` analytical scenario (De Conti ground return on the GMD case-3.1 duct model -- case4 had no `deconti` scenario at all before) + same modal wiring as case3 | `testData/andreata_case4/andreata_case4.py` |
| `MatlabDataReader.get_modal_scenario_data` generalized to auto-detect the per-file mode count (needed since Configs 3/4 will send 7 modes, not 6) | `utils/matlab_data.py` |

**Deliberately left undone**: `ModalDecomposition._classify` still only
handles `n == 6` -- there is no confirmed reference-pattern scheme for a 7th
"ECC mode," so Configs 3/4 modes come back as `mode_0`..`mode_6` (no
semantic label). This was a conscious choice over guessing: labelling by
raw index order would pair up two independently-computed, arbitrarily
ordered eigenmode sets with no verified correspondence -- exactly the
mistake sec. 7.1's mode-order verification was there to avoid for Configs
1/2. The reader/plotter/report machinery is otherwise fully ready: once
Config 3/4 MATLAB data and a 7-mode classification scheme both exist, the
overlay activates with no further code changes.

Both cases run end-to-end (passive, 7 real modes with physically sensible
`alpha_m`/`v_m`/`|Z_cm`| shapes -- see `Results/modal_*.png` in each case
folder).

### 7.3 FEM-hybrid pipeline extended to Configuration 4

Per `HYBRID_PIPELINE_PLAN.md` sec. 9 ("Config. 4: internal = FEM, 3x3 block
of the shared pipe + 2x2 for the others"), previously unimplemented.

| Delivered | File |
|---|---|
| `InternalParametersFromFEM.from_component_blocks(frequencies, component_blocks, block_sizes=None)` -- heterogeneous constructor: assembles a dense `(Nf, Ntot, Ntot)` `Zi`/`Pi` from arbitrarily-sized FEM blocks placed block-diagonally (no coupling *between* blocks); `block_sizes` (for the ground-return expansion) is accepted separately since it can group differently than the internal blocks (Case 4: 3 internal blocks vs. 4 ground-return objects). Homogeneous constructor (Config 2) untouched -- verified byte-identical regression. | `analytical_forms/single_core_cable.py` |
| `'fem'` scenario in `andreata_case4.py`: phases A/B reuse `andreata_case2`'s own FEM `Zi`/`C` (bare duct, no ECC); phase C + ECC use `andreata_case4`'s own 3-conductor COMSOL file (dense core/sheath/ECC coupling, sliced from its global 7x7 layout at indices `[0, 1, 6]`) | `testData/andreata_case4/andreata_case4.py` |
| `_validate_fem_hybrid_vs_matlab` (case4) -- same pattern as case2's function of the same name | `testData/andreata_case4/andreata_case4.py` |
| Modal decomposition switched to run on `'fem'` (falls back to `'3_deconti'` if the COMSOL data is missing) | `testData/andreata_case4/andreata_case4.py` |

**Validation**: `Y'` 2.56%, `Zg`/`Pg` (full matrix) 0.13%/0.34% -- all in the
same range as Config 2's own FEM-hybrid validation. `Z'` shows a ~50% outlier
isolated entirely to the `(ECC, ECC)` self-term above ~1 MHz (reactive part
only); root-caused to the ground-return radius assigned to the ECC object
not accounting for it sharing phase C's duct -- **not a defect in this FEM
work** (confirmed scenario-independent: identical in the purely-analytical
`'3_deconti'` path; Case 3, which has no duct, matches its own MATLAB
reference to ~1e-10). Full analysis, code citations and two candidate fixes
in [`GROUND_RETURN_ALLOCATION_STUDY.md`](GROUND_RETURN_ALLOCATION_STUDY.md)
-- intentionally left unfixed pending a team decision.

---

## 8. References

- **Andreata, L. E. B.** (2025). *Análise das Características de Propagação e de
  Transitórios Eletromagnéticos em Cabos Subterrâneos Instalados em Tubos Não
  Metálicos no Contexto de Parques Eólicos.* Thesis, PPGEE/UFMG. Ch. 5.
- **Wedepohl, L. M.; Nguyen, H. V.; Irwin, G. D.** (1996). Frequency-Dependent
  Transformation Matrices for Untransposed Transmission Lines using
  Newton-Raphson Method. *IEEE Trans. Power Delivery*, 11(3), 1538–1546.
  (sec. 6 — eigenvector correlation; sec. 7 — asymptotics.)
- **Gustavsen, B.** (2008). Fast Passivity Enforcement for Pole-Residue Models
  by Perturbation of Residue Matrix Eigenvalues. *IEEE Trans. Power Delivery*,
  23(4), 2278–2285. (sec. IV-A — switching-back; eq. 3 — passivity criterion;
  eq. 7 — eigenvalue perturbation / `v^T v` orthogonality.)
- **Paul, C. R.** (2008). *Analysis of Multiconductor Transmission Lines*, 2nd ed.
  (modal theory — Ch. 7.)
