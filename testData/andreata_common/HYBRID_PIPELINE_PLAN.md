# Plan — Hybrid pipeline: internal parameters (FEM/COMSOL) + ground return (analytical)

---

## STATUS — implemented for Configuration 2

| Delivered | File |
|---|---|
| **Part A** — the ground return sees the outer radius of the **pipe** (`_cable_external_geometry`): when there is an `enclosure`, the self term of `Zg`/`Pg` uses `D2/2` and the depth of the pipe axis | `mtl_main/strategy.py` (`SingleCoreCableInHDPEStrategy` + `...WithECCInHDPEStrategy`) |
| **Part B** — `InternalParametersFromFEM`: `Zi` (log-log interp.) + nodal `C` -> 6-key dict, tiled `kron(I_3,.)` | `analytical_forms/single_core_cable.py` |
| Reader `get_scc_internal_capacitance_matrix_combined()` (raw nodal `C` matrix) | `utils/comsol_data.py` |
| Orchestrator `build_pul_matrices(internal_source='analytical'\|'fem', ...)` | `analytical_forms/single_core_cable.py` |
| `andreata_case2`: `'fem'` scenario (Zi/Yi from COMSOL + analytical earth return with the pipe radius), modal decomposition on it, validation vs. MATLAB, 4 comparison plots + Figs 5.8–5.10 | `testData/andreata_case2/` |

**Validation `andreata_case2` — `'fem'` scenario vs. MATLAB (Andreata's FEM):**
*(ground return = closed-form expressions of De Conti/Duarte/Alipio 2023 —
`zg_form = yg_form = 'deconti'`, eqs. 4.59/4.63, like Andreata sec. 5.4 / 6.1)*

| quantity | max rel. error (inf norm) |
|---|---|
| `Z'` (full series) | **0.28 %** |
| `Y'` (full shunt) | **1.78 %** |
| `Zg` (ground return) | **0.13 %** |
| `Pg` (potential coeff.) | **0.34 %** |

The FEM-hybrid path reproduces Andreata's FEM reference to **< 2 %** — far
better than Lafaia's GMD (kept in scenarios `1/2/3` only for comparison).
Modal: clean labels, `classif = 1.000`, `PASSIVE`; the HF oscillation of the
coaxial modes (a GMD degeneracy artifact) **disappears** with the FEM `Zi`.

Regression: `andreata_case1/3/4` unchanged (case4 uses `scc-flat-ecc`, not HDPE).

**Pending:** Configs 4–5 (the 3x3 heterogeneous block of the shared pipe — Part A
already covers the geometry; the 3x3 FEM `Zi` is missing); frequency-dependent
soil; a small ~2 % offset in `C_22` (difference between the COMSOL FEM
and Andreata's own FEM, not the pipeline).

---

## 1. Motivation and diagnosis

### The 2-term formulation (Ametani 2015) — the pyLCP one

pyLCP decomposes the PUL parameters of a buried cable into **two** terms:

```
Z'(f) = Zi(f) + Zg(f)          Y'(f) = ( Yi(f)^-1 + Yg(f)^-1 )^-1
        \_ internal _/  \_ soil _/
```

- **`Zi` / `Yi`** — everything that is **not** ground return: the conductor
  impedance (skin/proximity) **and** the external inductance/capacitance of the
  configuration's insulation up to its outer surface. For Config. 1 that
  surface is the SCC jacket (`r4`); for Config. 2 it is the **outer surface
  of the HDPE pipe** (`D2/2`), with air + HDPE inside the internal term.
- **`Zg` / `Yg`** — ground return from that surface, via integral (Magalhaes/Xue)
  or closed-form (De Conti/Duarte/Alipio 2023) expressions, which **only see
  the position and the largest external dimension** of each object
  (Andreata eq. 4.27: *"x_jk equal to the outer radius of the cables or pipes"*).

> Andreata (sec. 4.1, sec. 4.4) writes `Z' = Zi + Ze + Zg`, splitting the internal
> term into conductors (`Zi`, Bessel) + external insulation (`Ze`, FEM). It is the
> **same decomposition**: pyLCP's 2-term `Zi` = Andreata's `Zi + Ze`. The COMSOL
> FEM delivers exactly that `Zi + Ze` already summed (see sec. 5).

**Confirmed premise:** the outer boundary of the Config. 2 COMSOL model is the
outer surface of the HDPE pipe (there is no air box). Therefore `Zi_fem` and `Zg`
match exactly on the `D2/2` circle, with no overlap or gap.

### What pyLCP does today

| | Config. 1 (`andreata_case1`) | Config. 2 (`andreata_case2`) |
|---|---|---|
| **internal** | analytical (Bessel core/sheath + closed form) — correct | analytical **+ Lafaia's GMD Case 3.1 trick**: the duct is folded into an equivalent insulation of permittivity `eps_a` on the sheath (`_override_sheath_insulation`) |
| **soil** — self-term outer radius | `0.024 m` = SCC outer radius (`r4`) — **correct** | `0.024 m` (scenario `3`/GMD) or `r7` (scenario `2`/ERS) — **the HDPE pipe (`D2/2 = 0.0578 m`) is ignored in the GMD scenario** |

### Problems

1. **Ground return with the wrong radius (Part A).**
   `SingleCoreCableInHDPEStrategy._cable_distance_matrices` selects the
   largest-outer-radius component by `data['radius'][1] + insulation.thickness`
   — it **never looks at `data['enclosure']`**. The HDPE pipe disappears from the
   self term of `d_matrix`/`D_matrix`. Verified:
   `flat_hdpe_enclosed_model` -> `diag(d_matrix_ground_return) = [0.024, 0.024, 0.024]`
   while `enclosure['outer_radius'] = 0.0578`.

2. **Internal via GMD Case 3.1 (Part B).**
   The equivalent-permittivity trick is an approximation (Lafaia 2015). pyLCP
   **already reads** the FEM `Zi`/`Yi` from COMSOL for the exact eccentric
   geometry of Config. 2 — today only as a comparison *overlay*, not as a source of `Z'/Y'`.

### Goal

Restructure the pipeline to **separate the two sources**:

- **internal**: *pluggable* — `'analytical'` (canonical geometries: bare SCC,
  SCC in a metallic pipe) **or** `'fem'` (COMSOL, for non-canonical geometries:
  SCC in an HDPE pipe — Config. 2 — and later Configs 3–5);
- **soil**: **always analytical** — closed-form expressions of De Conti/Duarte/Alipio
  2023 (`'deconti'`, eqs. 4.59/4.63; what Andreata sec. 5.4 / 6.1 used) or, as an alternative,
  Sommerfeld integrals (`'magalhaes_xue'`) — using the **real largest external
  dimension** of the configuration (SCC -> `r4`; pipe -> `D2/2`).

Result: Config. 2 without Lafaia's GMD.

---

## 2. Design principle

The composition point **already exists and does not change**:
`PerUnitParameters.quasi_tem_approx_matrices(internal_matrices, earth_return)`
does `Zs = Zi + Zg` and `Ysh = jw (Pi + Pg)^-1`, reading `Zi`, `Pi` and `block_sizes`
from an **`internal_matrices` dict**. It is enough to:

- **(A)** fix the outer radius that `earth_return_parameters` uses (via
  `d_matrix_ground_return` etc.);
- **(B)** produce the `internal_matrices` dict from the FEM, with the **same
  shape** that `InternalPerUnitParameters.matrices()` returns today;
- **(C)** a thin orchestrator that chooses the internal source and calls the composition.

No change in `quasi_tem_approx_matrices`, `ModalDecomposition`,
`ModalPropagationPlotter`, `check_pul_passivity`.

---

## 3. Part A — correct outer radius in the ground return

### 3.1 Where

`mtl_main/strategy.py`:
- `SingleCoreCableInHDPEStrategy._cable_distance_matrices` (~line 586)
- `SingleCoreCableWithECCInHDPEStrategy._cable_distance_matrices` (same copy)

`SingleCoreCableStrategy` (bare SCC, ~line 360) **does not change** — there is no `enclosure` there.

### 3.2 Change

When grouping conductors by `center_point` and choosing the cable representative,
consider the `enclosure` (pipe):

```python
def get_outer_radius(conductor_tuple):
    data = conductor_tuple[1]
    encl = data.get('enclosure')
    if encl and encl.get('outer_radius'):
        return encl['outer_radius']                    # pipe dominates
    insul = (data.get('insulation') or {}).get('thickness', 0)
    return data['radius'][1] + insul
```

and in the self term (`n_tag == m_tag`):

```python
encl = n_cable.get('enclosure')
s = (encl['outer_radius'] if encl and encl.get('outer_radius')
     else n_cable['radius'][1] + (n_cable.get('insulation') or {}).get('thickness', 0))
```

**Self-term depth**: when the pipe is the representative object, use
the **pipe center** (`enclosure['center_point'][1]`, e.g. `-1.1732`) instead of
the cable center (`-1.2`) in the self `images_vertical_distance_matrix` and `D_matrix`
terms. The **mutual** terms (center-to-center, 0.2 m) **do not change** — the pipes
do not touch (`0.2 > 2 * 0.0578`).

### 3.3 Cleaner alternative (recommended)

Instead of spreading the logic through the loop, the strategy computes **one field per cable**
in the model/context:

```python
context.external_radius_ground_return = np.array([...])   # (N_cables,)
context.depth_ground_return           = np.array([...])   # (N_cables,)
```

and `_cable_distance_matrices` + `earth_return_parameters` read those vectors.
Testable in isolation; explicitly documents "which radius the soil sees".

### 3.4 Expected effect

Config. 2, self term of `Zg`/`Pg`: `x_jk` goes from `0.024` to `0.0578` m
-> smaller `K0(gamma_g * x)`, larger `ln`-term -> **smaller self `Zg`, larger self
`Pg`** (smaller return capacitance). It combines with the duct `Yi` (large
dielectric) to give the high `Z_cm`/`v_m` of the ground/inter-sheath modes
(Figs 5.9–5.10) **without** the GMD `eps_a`.

### 3.5 Regression

Config. 1 (`SingleCoreCableStrategy`, no `enclosure`) — unchanged.
`hdpe_300mm2` / `hdpe_2000mm2` — start using the pipe radius; check against
the MATLAB/COMSOL of those cases (it may **improve** the agreement).

---

## 4. Part B — internal parameters from FEM/COMSOL

### 4.1 Output contract (same as the analytical one)

`InternalPerUnitParameters.matrices()` returns today:

```python
{
  'impedance_matrix'            : Zi,   # (Nf, SumM, SumM)  complex  [Ohm/m]
  'shunt_admittance_matrix'     : Ye,   # (Nf, SumM, SumM)  complex  [S/m]
  'potential_coefficient_matrix': Pi,   # (SumM, SumM)      real, freq-indep  [m/F]
  'capacitance_matrix'          : Ci,   # (SumM, SumM)      real, freq-indep  [F/m]
  'block_sizes'                 : [M]*N # conductors per cable
}
```

Of all that, `quasi_tem_approx_matrices` only uses `impedance_matrix`,
`potential_coefficient_matrix` and `block_sizes`. The FEM builder needs to produce
at least those three.

### 4.2 New builder

`analytical_forms/single_core_cable.py` (or `utils/comsol_data.py`):

```python
class InternalParametersFromFEM:
    """Assembles the `internal_matrices` dict (same shape as
    InternalPerUnitParameters.matrices()) from per-cable matrices
    measured in FEM/COMSOL. For N identical cables, it tiles kron(I_N, block).
    For heterogeneous systems (SCC + ECC), it assembles block-diagonal, like
    _matrices_heterogeneous."""

    def __init__(self, frequencies, fem_cable_blocks):
        # fem_cable_blocks: list of dicts, one per physical cable:
        #   {'Zi': (Nf_fem, M, M) complex [Ohm/m],  'C': (M, M) real [F/m] nodal}
        # Zi is interpolated (log-log, Re/Im) on the `frequencies` grid; C is constant.
        ...

    def matrices(self):
        # Zi_full  = block-diag / kron of the per-cable Zi(interp)
        # Pi_full  = block-diag / kron of the per-cable C^-1
        # Ci_full  = block-diag / kron of the C
        # Ye_full  = jw * inv(Pi_full)         (for inspection only)
        # block_sizes per each cable's M
        return {...}   # same 6 keys
```

### 4.3 Composition (unchanged)

```python
internal = InternalParametersFromFEM(freqs, blocks).matrices()
pul      = PerUnitParameters(mtl_hdpe, freqs)           # mtl with Part A applied
earth    = pul.earth_return_parameters('deconti', 'deconti')
qt       = pul.quasi_tem_approx_matrices(internal, earth)
qt       = pul.propagation_matrices(qt)                 # gamma_v, Zc, ...
```

---

## 5. COMSOL data — the existing files already serve

**No new export is needed** and there is no COMSOL dependency in this plan.
The two files in `testData/andreata_case2/Results/` already model the **exact
eccentric geometry of Config. 2** (cable resting on the bottom of the HDPE pipe, air +
pipe in the mesh, outer boundary on the pipe surface — a confirmed premise).
Their `Zi` and `Yi` enter the FEM-hybrid path **in full**; pyLCP already
has the readers. Only a new reading method (sec. 5.3) + grid interpolation remains.

### 5.1 `cmsl_internal_impedance_matrix.txt` -> `Zi(f)` 2x2

FD magnetodynamics, current excitation. Existing reader:
`ComsolPostProcessor.get_scc_internal_impedance_matrix_combined()` ->
`{'scenarios': {'measured': {'impedance_matrix': Zi}}}` with `Zi` `(Nf, 2, 2)`
`[Ohm/m]`, indices core=0 / sheath=1:

```
Zi[:,0,0] = data1(V_core_i)     core excited, voltage on the core   (self)
Zi[:,0,1] = Zi[:,1,0] = data1(V_sheath)   core-sheath mutual
Zi[:,1,1] = data2(V_sheath_i)   sheath excited, voltage on the sheath (self)
```

This `Zi` **already contains `Ze`** — the external inductance of the air + pipe
ring up to the pipe outer surface (the model boundary). It is exactly the `Zi + Ze` of
Andreata's decomposition; the composition `Zs = Zi_fem + Zg` closes with the
analytical ground return of Part A (radius = `D2/2`).

### 5.2 `cmsl_internal_admittance_charge_method.txt` -> nodal `C` 2x2

Electrostatics, direct charge method, real `epsr` (air = 1, HDPE = 2.35, XLPE
with semiconducting correction, PVC). Existing reader:
`ComsolPostProcessor.get_scc_internal_admittance_matrix_combined()` -> `Yi = jw*C`
with **nodal** `C` `(2, 2)` `[F/m]`, frequency-independent:

```
C[0,0] = C_core_coreExc                    = Cc         (core-sheath cap.)
C[0,1] = C[1,0] = C_shIn_coreExc           = -Cc
C[1,1] = C_shIn_shExc + C_shOut_shExc      = Cc + Cs    (Cs = sheath -> pipe surface,
                                                          through PVC + air + HDPE)
```

Verified numerically: `Cs = C_shOut_shExc ~ 1.37e-10 F/m` — about 12x
smaller than the PVC-only capacitance (`~1.6e-9`), consistent with the large
air + HDPE dielectric in series. **The duct is captured.**

For the pipeline: `Pi_cable = C^-1` (2x2), and `internal_matrices['potential_
coefficient_matrix'] = kron(I_3, Pi_cable)` — the same shape as the analytical path
(where `capacitance_matrix` is also the nodal matrix, see `compare_capacitance`).

### 5.3 Actions in pyLCP (code only, no COMSOL)

- New `ComsolPostProcessor.get_scc_internal_capacitance_matrix_combined()` —
  returns the raw nodal `C` matrix `(2, 2)` (or `(7, 7)` for ECC), without the `jw`
  (today only the `*jw` version exists, for the overlay).
- `InternalParametersFromFEM` consumes `Zi` (5.1) + `C` (5.2), **interpolates `Zi`
  on the analytical grid** (log-log, Re and Im separated; `C` is constant) and tiles
  `kron(I_3, .)` -> 6-key dict.
- **Model boundary — confirmed:** the outer boundary of the two `.mph` is the
  outer surface of the HDPE pipe (no air box). `Zi_fem` = "everything up to
  `D2/2`" and closes directly with the analytical `Zg` of Part A. No pending
  verification.

### 5.4 Cross-check (optional)

`Le = Im{Zi_fem}/omega` at low frequency should match `mu0 eps0 (C0)^-1` if one day
there is a `C0` export (all `epsr = 1`) — Andreata eq. 4.48. It blocks nothing.

---

## 6. Hybrid pipeline API

`analytical_forms/single_core_cable.py` (or a new `analytical_forms/hybrid.py`):

```python
def build_pul_matrices(mtl, frequencies, *,
                       internal_source='analytical',   # | 'fem'
                       fem_blocks=None,
                       zg_form='deconti', yg_form='deconti',
                       internal_form='hybrid'):
    """Orchestrates internal (analytical OR FEM) + ground return (analytical)
    + quasi-TEM composition + phase propagation. Returns the same dict as
    quasi_tem_approx_matrices + propagation_matrices."""
    if internal_source == 'analytical':
        internal = InternalPerUnitParameters(mtl, frequencies).matrices(internal_form)
    elif internal_source == 'fem':
        internal = InternalParametersFromFEM(frequencies, fem_blocks).matrices()
    pul   = PerUnitParameters(mtl, frequencies)
    earth = pul.earth_return_parameters(zg_form, yg_form)
    qt    = pul.quasi_tem_approx_matrices(internal, earth)
    return pul.propagation_matrices(qt)
```

---

## 7. Integration — `andreata_case2` + modal

1. New `'fem'` scenario alongside `'1'/'2'/'3'`:

```python
fem_blocks = ComsolPostProcessor(__file__).get_hybrid_internal_blocks()   # 3x 2x2 block
pul_data['scenarios']['fem'] = {
    'mtl': mtl_0,                      # flat_hdpe_enclosed_model (real pipe, Part A active)
    'internal_source': 'fem', 'fem_blocks': fem_blocks,
    'zg_form': 'deconti', 'yg_form': 'deconti',
}
```

2. Scenario loop: if `internal_source == 'fem'`, use `build_pul_matrices`;
   otherwise, the current path (GMD/ERS/underground stays for comparison).

3. **The modal decomposition runs on `'fem'`** (no longer `'3'`). The
   `ModalPropagationPlotter` generates Figs 5.8/5.9/5.10 of the FEM scenario;
   optional: a dashed overlay of scenario `'3'` (GMD) for comparison.

4. `check_pul_passivity` on the `'fem'` scenario.

---

## 8. Validation

| Target | Reference | Criterion |
|---|---|---|
| Self `Zg`/`Pg` with the pipe radius | manual computation `K0(gamma_g*D2/2)` | exact |
| FEM-hybrid `Z'`, `Y'` | `andreata_case2_series_impedance_matrix.mat` / `..._shunt_admittance_matrix.mat` (Andreata's own reference, which **is** FEM) | rel. error < ~5 % (mold: `validate_against_case1_reference`) |
| FEM-hybrid `Z'`, `Y'` vs GMD | scenario `'3'` | quantify how much the GMD deviates |
| FEM-hybrid `alpha_m`, `v_m`, `\|Z_cm\|` | Figs 5.8–5.10 | shape + magnitude; **must improve** vs GMD |
| passivity `Z'`, `Y'`, `Yc` FEM-hybrid | Gustavsen eq. 3 | PSD over the whole range |

---

## 9. Generalization (Configs 3–5)

Same pattern:

| Config. | internal | outer radius for the soil |
|---|---|---|
| 3 (buried SCC + ECC) | heterogeneous analytical (already exists) | `max(r4_SCC, r_ECC)` per cable — OK today (no `enclosure`) |
| 4 (SCC in pipes + ECC sharing one pipe) | **FEM** (3x3 block of the shared pipe + 2x2 for the others) | pipe `D2/2` (Part A already covers `SingleCoreCableWithECCInHDPEStrategy`) |
| 5 (same, bare copper ECC) | **FEM** (only the ECC in the mesh changes) | pipe `D2/2` |

`InternalParametersFromFEM` is already designed for heterogeneous blocks (block-diagonal
with variable `M`) -> it serves all three.

---

## 10. Phases

### Phase A — outer radius (0.5–1 day)
- [ ] `external_radius_ground_return` / `depth_ground_return` field in the HDPE strategies.
- [ ] `_cable_distance_matrices` (HDPE + shared-HDPE) reads the field; the self term uses the pipe radius + depth.
- [ ] regression: Config. 1 unchanged; `hdpe_300mm2`/`hdpe_2000mm2` checked.
- [ ] unit test: `diag(d_matrix_ground_return)` of `flat_hdpe_enclosed_model` = `[0.0578]*3`.

### Phase B — FEM reader + builder (1 day)
- [ ] `ComsolPostProcessor.get_scc_internal_capacitance_matrix_combined()` — raw nodal `C` matrix (without `jw`).
- [ ] `ComsolPostProcessor.get_hybrid_internal_blocks()` — joins `Zi` (5.1) + `C` (5.2), returns `[{'Zi','C'}]` per cable.
- [ ] `InternalParametersFromFEM` — interpolates `Zi` on the analytical grid + tiles `kron(I_3, .)` -> 6-key dict. Anticipate the heterogeneous block (ECC) from the start.

### Phase C — orchestrator + case (1 day)
- [ ] `build_pul_matrices(internal_source='analytical'|'fem', ...)`.
- [ ] `andreata_case2`: `'fem'` scenario, modal decomposition on it, Figs 5.8–5.10 (+ dashed GMD overlay for comparison).
- [ ] `check_pul_passivity` + diagnostics on the `'fem'` scenario.

### Phase D — validation and doc (0.5–1 day)
- [ ] section 8 table filled in (`Z'`/`Y'` vs Andreata's `.mat`; GMD vs FEM).
- [ ] update `CH5_PROPAGATION_PLAN.md`, `MODAL_CH5_DEVELOPMENT.md`, README, `andreata_case2/README.md`.

**Total: ~3–4 days.** No COMSOL dependency — the two files already exist.

---

## 11. Risks

1. ~~Outer boundary of the `.mph`~~ — **confirmed**: it is the outer surface of the
   pipe, with no air box. `Zi_fem` and `Zg` match at `D2/2` (sec. 1, sec. 5.3). No risk.
2. **`Zi_fem` <-> `Zg` interface in the code.** The composition `Zs = Zi_fem + Zg` is only
   exact because the FEM boundary and the `x_jk` of `Zg` (Part A) are the same
   `D2/2` circle. Keep that contract documented in a single place
   (`InternalParametersFromFEM` + `_cable_distance_matrices`).
3. **Depth of the representative object.** Pipe centered at `-1.1732`, cable at
   `-1.2`. Define and fix it (use the pipe center when the pipe represents the cable).
4. **COMSOL frequency grid != analytical grid.** `get_general_parameters`
   already handles this for overlays; to *compose* it is necessary to interpolate the FEM
   `Zi(f)` on the analytical grid (log-log, real and imaginary parts separated).
5. **Semiconducting layers in the FEM.** Either model the real geometry, or fold the
   permittivity as on the analytical side — choose and match it with the `Zi`.
