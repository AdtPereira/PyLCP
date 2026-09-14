# Plan — Modal-Domain Propagation Characteristics (Andreata Ch. 5) — **Configurations 1 and 2**

---

## STATUS (updated)

**Core implemented; Config. 1 (Figs 5.5–5.7) and Config. 2 (Figs 5.8–5.10) validated.**

| Delivered | File |
|---|---|
| `phase_domain_propagation(Zs, Ysh)` (gamma_v, gamma_i, Zc, Yc) — factored out of the OHTL path | `mtl_main/propagation.py` |
| `overhead_lines.py` migrated to the shared module (regression OK: `ohtl_deConti_ex51`) | `analytical_forms/overhead_lines.py` |
| `PerUnitParameters.propagation_matrices()` for SCC | `analytical_forms/single_core_cable.py` |
| `ModalDecomposition` — decompose (degenerate cluster refinement) -> track (*switching-back*, anchor at the middle of the band) -> modal_parameters -> **2-stage classify** (separates sheath vs coaxial modes by the fraction of energy in the cores; sub-labels by pattern — eqs 5.34–5.39) | `analytical_forms/modal_analysis.py` |
| `ModalPropagationPlotter` — `alpha_m`, `v_m`, `\|Z_cm\|` | `plotter/modal_plotter.py` |
| `andreata_case1.py` -> Figs 5.5 / 5.6 / 5.7 (Config. 1, `p100_er1_deconti` scenario — closed-form De Conti ground return, sec. 5.4) | `testData/andreata_case1/andreata_case1.py` |
| `andreata_case2.py` -> Figs 5.8 / 5.9 / 5.10 (Config. 2, `fem` scenario — FEM-hybrid pipeline, see [`HYBRID_PIPELINE_PLAN.md`](HYBRID_PIPELINE_PLAN.md)) | `testData/andreata_case2/andreata_case2.py` |
| **Passivity diagnosis** (assessment only — Gustavsen 2008 eq. 3): `eig(Re{Z'})`, `eig(Re{Y'})`, `eig(Re{Yc})` >= 0 | `utils/passivity_check.py` |
| Tracking cited as the *switching-back procedure* (Gustavsen 2008 sec. IV-A; Wedepohl 1996 sec. 6 — the `intercheig` sub-routine that Andreata sec. 5.2 reused); optimal assignment via Hungarian | docstrings of `modal_analysis.py` |

**Results** (constant soil rho = 100 Ohm.m; ground return = **closed-form De Conti/Duarte/
Alipio 2023**, `zg_form = yg_form = 'deconti'`, eqs. 4.59/4.63 — like
Andreata sec. 5.4 / 6.1): clean labels
`(ground, inter_sheath_1, inter_sheath_2, coaxial_1, coaxial_2, coaxial_3)` in
both configurations; classification similarity **= 1.000**;
`max_offdiag_ratio_Zm ~ 0`; reconstruction `T_I Lambda T_I^-1 ~ Y'Z'` < 1e-16;
`overall: PASSIVE` in both.

| | ground mode | inter-sheath | coaxial | vs Andreata |
|---|---|---|---|---|
| **Config. 1** `\|Z_cm\|` @ 1 MHz | ~60 Ohm | ~16–20 Ohm | ~12 Ohm | Figs 5.5–5.7 OK |
| **Config. 2** `\|Z_cm\|` @ 1 MHz | ~140 Ohm | ~55–70 Ohm | ~13 Ohm | Figs 5.8–5.10 OK |
| **Config. 2** `v_m` @ 1 MHz | ~0.5e8 | ~1.2–1.45e8 | ~1.8e8 | inter-sheath **much** faster than in Config. 1 OK |

The Config. 2 signature — `Z_cm` and `v_m` of the ground/inter-sheath modes
much larger than in Config. 1 — comes from the air + HDPE pipe inside the internal
term (`Zi`/`Yi` from COMSOL) + ground return with the pipe outer radius.
`Z'`/`Y'` of the `fem` scenario match Andreata's FEM reference to **< 2 %**.

**Update (2026-09)** — see `MODAL_CH5_DEVELOPMENT.md` sec. 7 for full detail:
- ~~modal MATLAB reference overlay~~ **Done** for Configs 1/2: `MatlabDataReader.
  get_modal_scenario_data` + `ModalPropagationPlotter(matlab_modal=...)` +
  `print_modal_comparison_report` — mode order verified (not assumed) against
  this doc's own validation table. Agreement excellent for `alpha_m`/`v_m`
  everywhere; `|Z_cm|` diverges above ~1 MHz on the 3 coaxial modes **in
  MATLAB's own curves too** (Config. 1), confirming the "cosmetic oscillation"
  below is a shared, not tracking-specific, phenomenon.
- ~~Hybrid pipeline~~ **Config. 2 done** (as below); **Config. 4 done**
  (`InternalParametersFromFEM.from_component_blocks` — heterogeneous
  block-diagonal FEM assembly: phase A/B reuse Config. 2's own FEM data,
  phase C+ECC use Config. 4's own 3-conductor COMSOL file). `Y'`/`Zg`/`Pg`
  agree to <3%; `Z'` has an isolated ~50% outlier on the ECC self-term,
  root-caused to the ground-return radius not accounting for the ECC sharing
  phase C's duct — see `GROUND_RETURN_ALLOCATION_STUDY.md` (fix pending a
  team decision, deliberately not applied yet).
- Configs 3/4: modal decomposition **wired and running** (7 conductors, 7
  modes) — `alpha_m`/`v_m`/`|Z_cm|` computed, passivity OK. Mode
  *classification* (the 7th "ECC mode") is **not** implemented — no
  reference-pattern equations for it exist yet, so modes stay unlabelled
  (`mode_0`..`mode_6`) rather than guessed. **Config 5 remains untouched.**
- `f_c` (eq. 5.33) annotation: still pending, low priority.
- Frequency-dependent soil (Phase 3, ~10–20% gap at high f): still pending.

---


> Current scope: **Configurations 1 and 2** (Figs 5.1 and 5.2) — three single-core
> cables (core + sheath), flat arrangement, 34.5 kV; **directly buried** (Config. 1)
> or **in individual HDPE ducts** (Config. 2). They are **6 conductors -> 6 modes** in both.
> Goal: from `Z'(f)` and `Y'(f)` (already produced by `andreata_case1` /
> `andreata_case2`), compute and plot **`alpha_m`, `v_m`, `|Z_Cm|`** per mode,
> reproducing **Figures 5.5–5.7** (Config. 1) and **5.8–5.10** (Config. 2).
>
> Out of scope (deferred): Configurations 3–5, the earth conductor (ECC) and its
> 7th mode, shared duct / external parameters by FEM, the sec. 5.4.2
> "per mode across configurations" family (Figs 5.20–5.31).

---

## 0. Configuration 2 — what changes relative to Config. 1

Same 6-mode structure (`T_I` of Andreata eqs 5.37–5.39 = 5.34–5.36). The
HDPE duct does **not** enter the ground-return formulation (`Zg`/`Yg` only
see the position and outer radius); it enters as an **equivalent outer insulation**
on the sheath, folded via `EquivalentRadiiSystems` (the same mechanism already used
by `andreata_case2`/`andreata_case4` for the shunt parameters). `andreata_case2`
already produces 3 scenarios — `1` (duct ignored), `2` (area-weighted ERS),
`3` (GMD case 3.1) — and the modal decomposition runs on the
**`3` scenario**, the most complete (same choice as the canonical internal matrix).

**Required adjustment to the classifier** (done): the old labeling projected
`T_I` at the minimum frequency; in Config. 2 the modes only separate above
~100 Hz (with the duct, the sheath modes carry appreciable core current
at LF). The new labeling is **2-stage**: (i) separates the 3 sheath modes
from the 3 coaxial ones by the fraction of `T_I` energy in the **core**
rows, averaged over a mid band (100 Hz – 100 kHz) where every config is
decoupled; (ii) sub-labels each group by the pattern (`[1,1,1]` / `[1,0,-1]` /
`[-1,2,-1]`). Robust for both configs (similarity = 1.000).

---

## 1. Diagnosis — what already exists vs. what is missing (Config. 1)

### 1.1 Already implemented (reuse)

| Item | Where |
|---|---|
| Full `Z'(f)`, `Y'(f)` (Zi+Ze+Zg / Yi \|\| Yg) for 3 buried SCC in a flat arrangement | `analytical_forms/single_core_cable.py::PerUnitParameters.quasi_tem_approx_matrices` -> keys `series_impedance_matrix`, `shunt_admittance_matrix` |
| Internal impedance/admittance (Bessel/approx/hybrid) + semiconducting-layer correction | `InternalPerUnitParameters`, `apply_semiconducting_layer_correction` |
| Ground return: closed-form `deconti` (De Conti/Duarte/Alipio 2023, eqs. 4.59/4.63 — used in Ch. 5) and integral `magalhaes_xue` | `PerUnitParameters.earth_return_parameters` |
| Case assembled up to `Z'/Y'` with the soil `rho = 100 Ohm.m` (scenario `p100_er1`) | `testData/andreata_case1/andreata_case1.py` |
| MATLAB reference for `Z'/Y'` + conductor-order remapping (`conductor_order=[0,3,1,4,2,5]`) | `utils/matlab_data.py::MatlabDataReader.get_scc_scenario_data` |
| **PHASE propagation** `gamma_v = sqrtm(Z'Y')`, `gamma_i`, `Z_c`, `Y_c` (per-frequency loop) | **only** `analytical_forms/overhead_lines.py::PerUnitParameters.pul_matrices` (lines ~516–558) |
| Attenuation `Re gamma` x velocity `omega/Im gamma` plot for 1 element `[p,q]` | `plotter/ohtl_models.py::_propagation_subplots`, `plotter/lima_models.py` |
| Cross-section schematic | `mtl_main/graphics.py::GroundReturnMTLRepresentation` |

### 1.2 Gaps to fill

1. **Modal transformation missing.** There are no eigenvalues/eigenvectors of `Y'Z'` / `Z'Y'`, nor `T_I`, `T_V`.
2. **Phase propagation does not exist for SCC.** `single_core_cable.py::PerUnitParameters` stops at `Z'/Y'`; only the OHTL path computes `gamma`, `Z_c`.
3. **Mode tracking.** Eigenvalues swap position between neighboring frequencies (Andreata sec. 5.2). The eigenvector continuity is missing (frequency-to-frequency scalar product — Gustavsen 2008).
4. **Modal parameters.** `Gamma_m = sqrt(lambda)` (`Re >= 0` branch), `alpha_m = Re Gamma_m`, `beta_m = Im Gamma_m`, `v_m = omega/Im Gamma_m`, `Z_Cm = Gamma_m/Y_m` (Andreata eqs. 5.23–5.32).
5. **Identification/labeling of the 6 modes** (ground; inter-sheath 1 and 2; coaxial 1/2/3) by inspecting the columns of `T_I` (Andreata eqs. 5.34–5.36).
6. **Modal-domain plots.** 3 figures (`alpha_m`, `v_m`, `|Z_Cm|`) with the 6 modes overlaid -> Figs 5.5, 5.6, 5.7.
7. **Soil with frequency-dependent parameters** (Salvador et al. 2020, `rho_LF = 100 Ohm.m`). Today `sigma_g`, `eps_g` are fixed scalars. *(optional — Phase 3; the initial version uses constant `rho = 100 Ohm.m` as in `andreata_case1`)*
8. **Decoupling cutoff frequency `f_c`** (eq. 5.33) — an informative utility to annotate plots. *(optional)*

---

## 2. Fundamentals to implement (Ch. 5 equations, N = 6)

```
input:  Z'(f) [Nf,6,6] , Y'(f) [Nf,6,6]

1. products             A_i = Y'Z'   ,  A_v = Z'Y'                        (eqs. 5.14–5.15)
2. eigendecomposition   A_i = T_I . diag(lambda) . T_I^{-1}              (eq. 5.19)
                        A_v = T_V . diag(lambda) . T_V^{-1}   (same lambda)   (eq. 5.18)
                        validate  T_V^{-1} ~ T_I^{t}                     (eq. 5.20)
3. tracking             reorder/re-phase columns of T_I(f_k) to match T_I(f_{k-1})
                        (scalar products; seed from high frequency)
4. modal constant       Gamma_m(f) = sqrt(lambda(f)) ,  branch  Re{Gamma_m} >= 0   (eq. 5.24)
5. modal matrices       Z_m = T_V^{-1} Z' T_I  (diagonal)                (eq. 5.21)
                        Y_m = T_I^{-1} Y' T_V  (diagonal)                (eq. 5.22)
6. modal char. imped.   Z_Cm = Gamma_m / Y_m                            (eqs. 5.23, 5.27)
7. output               alpha_m = Re Gamma_m   beta_m = Im Gamma_m   v_m = omega / Im Gamma_m   (eqs. 5.30–5.32)
8. mode label           classify(T_I, core/sheath roles)                (sec. 5.3, eqs. 5.34–5.36)
```

**The 6 modes of Configuration 1 (Andreata sec. 5.3):**

| # | Label | Description | Behavior vs f |
|---|---|---|---|
| 1 | `ground` | current through the 3 sheaths, ground return | ~invariant |
| 2 | `inter_sheath_1` | 1 sheath x another sheath | ~invariant |
| 3 | `inter_sheath_2` | 1 sheath x the other two | ~invariant |
| 4 | `coaxial_1` | 3 cores x earth (LF) -> core x sheath (HF) | varies with f |
| 5 | `coaxial_2` | core x core (LF) -> core x sheath (HF) | varies with f |
| 6 | `coaxial_3` | core x 2 cores (LF) -> core x sheath (HF) | varies with f |

`f_c = rho_sh / [ pi mu_sh (r_sh,out - r_sh,in)^2 ]`  (eq. 5.33) — optional, annotates the start of decoupling.

---

## 3. Proposed architecture

### 3.1 New artifacts (lean)

```
mtl_main/propagation.py              # phase-domain: gamma_v, gamma_i, Zc, Yc (extracts the loop from overhead_lines.py)
analytical_forms/modal_analysis.py   # ModalDecomposition: eig, tracking, modal params + classify (6 modes)
plotter/modal_models.py              # declarative config of the 3 modal plots (ohtl_models.py style)
plotter/modal_plotter.py             # ModalPropagationPlotter(BasePlotter)
models/frequency_dependent_soil.py   # sigma_g(f), eps_g(f) — Phase 3, optional
```

Changes:
```
analytical_forms/overhead_lines.py       # use mtl_main/propagation.py (without changing the result)
analytical_forms/single_core_cable.py    # PerUnitParameters.propagation_matrices();
                                         # (Phase 3) accept vector sigma_g/eps_g
testData/andreata_case1/andreata_case1.py    # call the modal decomposition + new plotter
testData/andreata_case1/plot_config.py       # 3 entries: modal_attenuation / modal_phase_velocity / modal_char_impedance
```

### 3.2 `mtl_main/propagation.py` (phase domain)

```python
def phase_domain_propagation(Zs, Ysh):
    """Zs, Ysh: (Nf, N, N). Returns:
       propagation_voltage_matrix   gamma_v = sqrtm(Zs*Ysh)
       propagation_current_matrix   gamma_i = sqrtm(Ysh*Zs)
       characteristic_impedance_matrix  Zc = Ysh^{-1}*gamma_i
       characteristic_admittance_matrix Yc = Zs^{-1}*gamma_v
    """
```
- Replace the duplicated block in `overhead_lines.py::pul_matrices` with a call to this function;
  run `run_ohtl.py` and check that the Lima figures do not change (regression).
- In `single_core_cable.py::PerUnitParameters`, add
  `propagation_matrices(quasi_tem_matrices)` -> returns the 4 keys above appended to the dict
  (backward compatible).

### 3.3 `analytical_forms/modal_analysis.py`

```python
CANONICAL_ORDER = 'cores_then_sheaths'   # [c1,c2,c3,s1,s2,s3]

class ModalDecomposition:
    def __init__(self, freq, Zs, Ysh, conductor_roles):
        # conductor_roles: list aligned to the rows of Zs/Ysh, e.g.
        #   [('core',0),('sheath',0),('core',1),('sheath',1),('core',2),('sheath',2)]
        # derives the permutation P to the canonical order [cores | sheaths]
        ...

    def decompose(self):
        # for each f (in canonical order):
        #   lambda, TI = eig(Ysh @ Zs)
        #   _,  TV = eig(Zs @ Ysh)
        #   match the eigenvalues of TV to those of TI (nearest-lambda) and reorder TV
        #   normalize eigenvectors: inf-norm = 1; phase of the 1st non-zero component = 0
        # returns raw (without tracking)

    def track(self):
        # continuity by scalar product; seed = HIGHEST frequency (modes ~ decoupled)
        # sweep downward: at f_k, find the permutation+signs that maximize
        #   Sum_j |< TI_j(f_k) , conj(TI_j(f_{k+1})) >|   via scipy.optimize.linear_sum_assignment
        #   over the similarity matrix S = |TI(f_k)^H . TI(f_{k+1})|
        # apply the same permutation/signs to lambda, TI, TV

    def modal_parameters(self):
        Gamma = np.sqrt(lam);  Gamma *= np.where(Gamma.real < 0, -1, 1)   # Re >= 0
        Zm = einsum('fij,fjk,fkl->fil', inv(TV), Zs, TI).diagonal(axis1=1, axis2=2)
        Ym = einsum('fij,fjk,fkl->fil', inv(TI), Ysh, TV).diagonal(axis1=1, axis2=2)
        Zcm = Gamma / Ym
        return dict(gamma=Gamma, alpha=Gamma.real, beta=Gamma.imag,
                    vphase=self.w[:,None]/Gamma.imag, Zcm=Zcm, Ycm=Gamma/Zm,
                    TI=TI, TV=TV, lam=lam, mode_labels=self.classify(),
                    decoupling_freq=self._fc())   # arrays (Nf, 6)

    def classify(self, ref='high'):
        # projects the columns of TI (at high f, and reconciles with mid/low f) onto
        # REFERENCE_PATTERNS_6C (below); label = largest |projection|
```

`REFERENCE_PATTERNS_6C` — normalized `T_I` columns in the canonical order `[c1,c2,c3,s1,s2,s3]`,
transcribed from Andreata eqs. 5.34–5.36:

```python
REFERENCE_PATTERNS_6C = {
    # "sheath" modes (columns 1-3 of T_I, ~invariant with f)
    'ground'        : [ 0,   0,   0,   1/3, 1/3, 1/3],   # sheaths co-phase, ground return
    'inter_sheath_1': [ 0,   0,   0,   1/2,   0, -1/2],
    'inter_sheath_2': [ 0,   0,   0,  -1/3, 2/3, -1/3],
    # coaxial modes (columns 4-6, here in the HIGH-frequency limit: core x its own sheath)
    'coaxial_1'     : [ 1,   0,   0,  -1,   0,   0 ],
    'coaxial_2'     : [ 0,   1,   0,   0,  -1,   0 ],
    'coaxial_3'     : [ 0,   0,   1,   0,   0,  -1 ],
}
```
> Adjust the exact signs/values per pyLCP's `Z'/Y'` assembly convention
> (validate against the matrices printed by Andreata in eqs. 5.34–5.36). The classification
> uses a **normalized projection** (robust to scale/sign), not equality.

Notes:
- **Reciprocity**: check `||T_V^{-1} - T_I^{T}|| / ||T_I||` (warning, not error).
- **Degeneracy** of the coaxial modes at LF: the `linear_sum_assignment` over the
  similarity matrix resolves the label swap; if it persists, match in the subspace (block projection).
- Reuse Gustavsen (2008) eigenvalue-perturbation routine if something already exists
  in `mom_so/` — check before reimplementing.

### 3.4 Canonical conductor order

pyLCP assembles `Z_i = kron(I_3, Z_ij)` -> rows `[c1,s1,c2,s2,c3,s3]`. Andreata uses
`[c1,c2,c3,s1,s2,s3]`. Define **one** permutation matrix `P` (reuse the logic of
`conductor_order=[0,3,1,4,2,5]` already used for MATLAB), apply it to `Z'` and `Y'` **before** the
decomposition, and document it in a single place (`modal_analysis.py`).

### 3.5 `models/frequency_dependent_soil.py` *(Phase 3 — optional)*

```python
def soil_parameters(freq, rho_lf, model='alipio'):
    """Returns sigma_g(f) [S/m], eps_r_g(f). Models: 'alipio' (Alipio-Visacro 2014),
       'salvador' (Salvador et al. 2020, used by Andreata sec. 5.4), 'constant'."""
```
Integration: `PerUnitParameters.__init__` already builds `k_earth2`, `gamma_earth`, `jw_2pi_sg`
as `(Nf,)` vectors — it is enough to allow `self.sigma_1`, `self.e1` as arrays and adjust the broadcasting.
Validate in isolation (the ground and inter-sheath modes at HF are the most sensitive).

### 3.6 Plotters

`plotter/modal_models.py` — declarative config (`ohtl_models.py` style), 3 plots, series = 1 mode:
```python
MODE_STYLE = {                       # fixed colors/dashes per label (consistent across the 3 figs)
  'ground'        : {'color':'black','linestyle':'-'},
  'inter_sheath_1': {'color':'tab:blue','linestyle':'--'},
  'inter_sheath_2': {'color':'tab:blue','linestyle':'-.'},
  'coaxial_1'     : {'color':'tab:red','linestyle':'-'},
  'coaxial_2'     : {'color':'tab:green','linestyle':'-'},
  'coaxial_3'     : {'color':'tab:orange','linestyle':':'},
}
```

`plotter/modal_plotter.py::ModalPropagationPlotter(BasePlotter)`:
```python
def modal_attenuation(self):      # Fig. 5.5  — alpha_m (Np/m) x f, log-log, 6 curves
def modal_phase_velocity(self):   # Fig. 5.6  — v_m (m/s, absolute) x f, semilog-x, 6 curves
def modal_char_impedance(self):   # Fig. 5.7  — |Z_Cm| (Ohm) x f, semilog-x, 6 curves
```
- units as Andreata: `alpha_m` in **Np/m** (parametrizable Np/km), `v_m` in **absolute m/s**
  (not normalized), `|Z_Cm|` in **Ohm**; x-axis `1e-1 ... 1e6 Hz`.
- reuse `save_figure`, `_format_axis`; optional overlay of the modal MATLAB reference (Phase 4).

---

## 4. Implementation phases

### Phase 0 — Preparation (~0.5 day)
- [ ] Factor out `phase_domain_propagation` in `mtl_main/propagation.py`; migrate `overhead_lines.py`; run `run_ohtl.py` (regression: Lima figures identical).
- [ ] Define the canonical `[cores | sheaths]` permutation matrix in `modal_analysis.py`.

### Phase 1 — Phase propagation for SCC + validation (~1 day)
- [ ] `PerUnitParameters.propagation_matrices()` in `single_core_cable.py`.
- [ ] Append `gamma_v, gamma_i, Zc, Yc` to the dict of the `p100_er1` scenario in `andreata_case1.py`.
- [ ] Sanity: `Z'`, `Y'` symmetric; `||gamma_v - gamma_i^T||` small; `Re gamma_v[0,0] >= 0`.

### Phase 2 — Modal decomposition and tracking (~1.5–2 days)
- [ ] `ModalDecomposition.decompose()` + `modal_parameters()` (without tracking).
- [ ] `track()` (scalar product + `linear_sum_assignment`, seed at high frequency).
- [ ] Tests:
  - reconstruction `T_I diag(lambda) T_I^{-1} ~ Y'Z'` (rel. error < 1e-8);
  - `Gamma_m^2 ~ eigenvalues(Z'Y')` (< 1e-8);
  - `alpha_m >= 0` over the whole range; `v_m <= c` (with margin for LF noise);
  - `ground` / `inter_sheath_*` columns of `T_I` ~constant vs f.
- [ ] Reproduce **Fig. 5.5 (alpha_m)** and **Fig. 5.6 (v_m)** — 6 curves.

### Phase 3 — Mode identification + `|Z_Cm|` (~1 day)
- [ ] `classify()` with `REFERENCE_PATTERNS_6C` (eqs. 5.34–5.36); check the labels against Andreata sec. 5.3 (modes 1–3 ~invariant; 4–6 coaxial).
- [ ] Reproduce **Fig. 5.7 (|Z_Cm|)**.
- [ ] *(optional)* `frequency_dependent_soil.py` (Salvador 2020, `rho_LF = 100 Ohm.m`) and re-run the 3 figures.
- [ ] *(optional)* annotate `f_c` (eq. 5.33) on the plots.

### Phase 4 — Plots, driver and validation (~1 day)
- [ ] `modal_models.py` + `ModalPropagationPlotter`; 3 entries in `andreata_case1/plot_config.py`.
- [ ] `andreata_case1.py` calls the 3 plots; saves to `testData/andreata_case1/Results/`.
- [ ] Figure-by-figure comparison with the PDF: relative-error table per mode at `f in {1e2, 1e4, 1e6} Hz`.
- [ ] *(optional)* ingest the modal MATLAB reference (`alpha_m`, `vphase_m`, `Zcm`) as a scatter.
- [ ] Update `README.md` (Ch. 5 section, Config. 1) and `testData/andreata_case1/BUGS_AND_FIXES.md`.

**Total estimate: ~4–6 days.**

---

## 5. Validation

| Level | Reference | Criterion |
|---|---|---|
| Unit | reconstruction `T_I Lambda T_I^{-1}` vs `Y'Z'` | rel. error < 1e-8 |
| Unit | `Gamma_m^2` vs eigenvalues of `Z'Y'` | rel. error < 1e-8 |
| Physical | `alpha_m >= 0`; `v_m -> c/sqrt(eps_eff)` at high f for coaxial modes | qualitative |
| Cross | `Z'`, `Y'` (scenario `p100_er1`) vs `andreata_case1` MATLAB dump | reuse the `validate_against_case3_reference` mold |
| Cross | `coaxial_*` mode at HF vs the phase `gamma` of the isolated core (`scc_34kV_andreata`) | < 5 % |
| Figures | Figs 5.5 / 5.6 / 5.7 of the PDF | visual inspection + per-mode rel. error over 3 decades of f |
| Labels | Andreata sec. 5.3 textual description | 1:1 correspondence of the 6 mode names |

---

## 6. Risks and limitations

1. **Mode tracking at degeneracy** (nearly identical coaxial modes at LF) can swap the label;
   mitigate with high-frequency seeding + `linear_sum_assignment` on the similarity matrix.
2. **Conductor order**: an inconsistency between `[c1,s1,...]` (pyLCP) and `[c1,c2,c3,s1,...]` (Andreata)
   breaks the identification — centralize into a single permutation and test it.
3. **Constant vs FD soil**: the initial version (constant `rho = 100 Ohm.m`) deviates from Figs 5.5–5.7 at
   high f (ground and inter-sheath modes); Phase 3 closes the gap.
4. **`sqrtm` vs eigenvalues**: use the eigenvalue route as the source of truth for `Gamma_m`; `sqrtm`
   only to cross-check the phase path.
5. **Semiconducting layers**: already handled via `eps_ins_eq` in `apply_semiconducting_layer_correction`;
   confirm that the correction is applied **before** the modal decomposition.

---

## 7. File checklist

**Create**
- `mtl_main/propagation.py`
- `analytical_forms/modal_analysis.py`
- `plotter/modal_models.py`, `plotter/modal_plotter.py`
- `models/frequency_dependent_soil.py`  *(Phase 3, optional)*

**Modify**
- `analytical_forms/overhead_lines.py` — use `propagation.py` (without changing the result)
- `analytical_forms/single_core_cable.py` — `PerUnitParameters.propagation_matrices()`; *(Phase 3)* vector `sigma_g(f)`, `eps_g(f)`
- `testData/andreata_case1/andreata_case1.py` — modal decomposition + 3 plots
- `testData/andreata_case1/plot_config.py` — `modal_attenuation` / `modal_phase_velocity` / `modal_char_impedance`
- `README.md`, `testData/andreata_case1/BUGS_AND_FIXES.md`
