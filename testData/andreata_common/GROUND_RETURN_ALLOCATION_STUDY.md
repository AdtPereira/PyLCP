# Study — Ground-Return (Zg/Yg) Treatment Across the 4 Andreata Cases

**Date:** 2026-09.
**Trigger:** implementing the FEM-hybrid pipeline for `andreata_case4`
(see [`HYBRID_PIPELINE_PLAN.md`](HYBRID_PIPELINE_PLAN.md)) surfaced a
50% mismatch in `Z'` against the MATLAB reference, isolated entirely to the
ECC self-term at high frequency. This document maps out how the De
Conti/Duarte/Alipio (2023) closed-form ground return is evaluated in each of
the 4 cases — in particular, **which external radius each conductor object
is assigned** — to give the team a precise basis for deciding how (or
whether) to fix it. **No code was changed as part of this study** — the
`'fem'` scenario of `andreata_case4` is implemented and running, just with
this known gap on the ECC term (see sec. 5 of `andreata_case4/andreata_case4.py`'s
FEM section and the validation printout).

---

## 1. The closed-form expressions and what they need as input

`PerUnitParameters.earth_return_parameters` (`zg_form='deconti'`),
[`analytical_forms/single_core_cable.py:1279`](../../analytical_forms/single_core_cable.py#L1279):

```python
Zg = jw_mu0_2pi * ( K0(j*kg*d_nm) + zg_term1 * exp(hnm*kg) * 2/(4 + kg**2*dnm**2) )   # eq. 4.59
Pg = jw_2pi_sg  * ( K0(j*kg*d_nm) + yg_term1 * K0(j*kg*D_nm) )                         # eq. 4.63
```

evaluated over four `(N, N)` geometry matrices — `d_matrix`, `D_matrix`
(image distance), `hnm` (image depth sum), `dnm` (horizontal separation) --
where **`N` is the number of ground-return *objects*, not the number of
conductors**. The rule used everywhere in the code:

- **self term (`n == n`)**: `dnm[n, n]` = the **external radius** of object
  `n` -- the duct's outer radius if the object sits inside one, otherwise the
  cable's own outer radius (conductor + insulation).
- **mutual term (`n != m`)**: `dnm[n, m]` = center-to-center horizontal
  distance -- no radius involved.

This radius selection is centralized in one shared helper,
[`_cable_external_geometry`](../../mtl_main/strategy.py#L7-L23) (`mtl_main/strategy.py`):

```python
def _cable_external_geometry(conductor_data):
    enclosure = conductor_data.get('enclosure')
    if enclosure and enclosure.get('outer_radius'):
        return (duct center, duct outer radius)
    return (conductor's own center, conductor radius + insulation thickness)
```

i.e. **the external radius "seen" by the ground return is decided
per-object, by looking only at that object's own `enclosure` field.** There
is no notion anywhere of "this conductor sits inside another conductor's
duct."

## 2. Who calls it, and how objects are formed

Every `_cable_distance_matrices` (there are 4 near-identical copies, one per
strategy class) runs the same algorithm:

1. Group conductors by `center_point` (same position => same physical cable);
2. Within each group, pick the **representative** conductor (largest outer
   radius -- normally the sheath, since it is the one carrying `enclosure`
   when a duct exists);
3. Call `_cable_external_geometry(representative)` to get that object's
   position + external radius;
4. Assemble `d_matrix`/`D_matrix` (`N x N`) with that radius on the diagonal.

`N` = number of groups (distinct positions) = `model.num_sc_cables`, read by
`earth_return_parameters` at
[`single_core_cable.py:1281`](../../analytical_forms/single_core_cable.py#L1281),
and later used to expand `Zg`/`Pg` back to the full per-conductor space via
[`_expand_by_block_sizes(matrix, block_sizes)`](../../analytical_forms/single_core_cable.py#L478)
inside `quasi_tem_approx_matrices`
([`single_core_cable.py:1381`](../../analytical_forms/single_core_cable.py#L1381)).

## 3. Comparison table — the 4 cases

| | **Case 1** (bare) | **Case 2** (duct) | **Case 3** (bare + ECC) | **Case 4** (duct + ECC) |
|---|---|---|---|---|
| Strategy class | `SingleCoreCableStrategy` | `SingleCoreCableInHDPEStrategy` | `SingleCoreCableWithECCStrategy` | `SingleCoreCableWithECCStrategy` (**same class as Case 3**) |
| `_cable_distance_matrices` | [strategy.py:380](../../mtl_main/strategy.py#L380) | [strategy.py:611](../../mtl_main/strategy.py#L611) | [strategy.py:1225](../../mtl_main/strategy.py#L1225) | [strategy.py:1225](../../mtl_main/strategy.py#L1225) (reused) |
| Ground-return objects (`N`) | 3 (phases A, B, C) | 3 (phases A, B, C) | 4 (A, B, C, ECC) | 4 (A, B, C, ECC) |
| `block_sizes` (conductor allocation) | `[2, 2, 2]` | `[2, 2, 2]` | `[2, 2, 2, 1]` | `[2, 2, 2, 1]` |
| Self radius -- phase A/B | cable: **0.024 m** (sheath + insulation) | **duct: 0.0578 m** (`enclosure.outer_radius`) | cable: **0.024 m** | **duct: 0.0578 m** |
| Self radius -- phase C | 0.024 m | 0.0578 m | 0.024 m | 0.0578 m |
| Self radius -- ECC | -- | -- | **0.0053 m** (ECC conductor + insulation, bare) | **0.0053 m** (identical to Case 3 -- no duct applied!) |

The matrix allocation (`_expand_by_block_sizes`) is structurally identical
across all 4 cases: each ground-return object "lends" its scalar `Zg`/`Pg`
to every conductor that belongs to it (e.g. phase A's core and sheath both
receive the same `Zg[A, A]` -- correct, since ground return does not
distinguish core from sheath, only the cable's position).

## 4. The inconsistency (Case 4)

Case 4 reuses the **exact same class** (`SingleCoreCableWithECCStrategy`)
and the same `_cable_external_geometry` helper as Case 3, applied
object-by-object. This works correctly for phases A/B/C (each has its own
`enclosure` on the sheath, so the duct radius is picked up correctly --
validated: the phase self-term `Zg`/`Pg` agree with the MATLAB reference to
0.13%/0.34%, matching andreata_case2's own FEM-hybrid validation numbers).

For the ECC, however: its conductor dict **never carries an `enclosure`**
(neither in Case 3 nor in Case 4) -- it is added as its own object, at its
own offset `center_point` relative to phase C. `_cable_external_geometry`
therefore falls through to the bare branch and returns the ECC's own radius
(0.0053 m) **identically in both cases**.

This is **physically correct in Case 3** (there is no duct at all -- the ECC
really is bare, exposed to the soil at its own radius). It is **physically
questionable in Case 4**: in the thesis's Configuration 4, the ECC sits
*inside the same duct as phase C*
(`models/single_core_cable.py::flat_hdpe_enclosed_with_shared_ecc_model`,
"ECC sharing the phase C duct"). By the stated rule -- *every ground-return
term must be evaluated per the De Conti closed-form expressions using the
external radius of the cable **or duct*** -- the ECC's self ground-return
term in Case 4 should arguably use the duct's outer radius (0.0578 m), not
the bare ECC radius (0.0053 m), because it is the duct surface that actually
separates the ECC's field from the soil.

This is the gap behind the ~50% `Z'` mismatch found while validating the
`'fem'` scenario: isolated entirely to the `(ECC, ECC)` entry, growing with
frequency, and confined to the reactive part (`Im{Zg}`; the resistive part
already agrees to <1%). Cross-checks that support this being a genuine,
scenario-independent modeling gap (not a bug in the new FEM code):

- The same discrepancy appears identically in the purely-analytical
  `'3_deconti'` scenario (no FEM involved at all).
- Case 3 (no duct) matches its own MATLAB `series_impedance_matrix.mat[ECC,
  ECC]` to ~1e-10 relative error, confirming the De Conti implementation and
  code path are internally consistent when the physical premise (no duct)
  matches the code's assumption.
- The mismatch appears *only* when Case 4's own MATLAB reference is
  compared -- consistent with Andreata's FEM reference genuinely modeling
  the ECC-inside-duct geometry, which pyLCP's `Zi (block-diagonal) + Zg
  (2-D image theory on bare radius)` decomposition currently does not.

## 5. Open question for the team

The code today has **no mechanism** to express "this conductor shares
another conductor's duct" -- the only available clue is geometric (the ECC's
position sits inside the duct circle drawn around phase C). Two candidate
fixes, not yet implemented, pending the team's decision:

1. **Case-4-specific**: when building `model_schematic`
   (`flat_hdpe_enclosed_with_shared_ecc_model`), copy phase C's `enclosure`
   (same duct center/radius) onto the ECC conductor's dict as well. A
   model-builder-only change -- `_cable_external_geometry` itself needs no
   change.
2. **General**: teach `_cable_distance_matrices`/`_cable_external_geometry`
   to detect "conductor sitting inside another conductor's enclosure"
   automatically, by geometric containment. More robust for future configs
   with other shared-duct arrangements, but a larger change surface.

Once the team confirms which radius Andreata's own formulation actually uses
for the ECC self-term in Configuration 4, resume at
`testData/andreata_case4/andreata_case4.py`'s `'fem'` scenario (already
implemented and passing everywhere else -- `Y'` 2.56%, `Zg`/`Pg` full-matrix
0.13%/0.34%) and re-validate `_validate_fem_hybrid_vs_matlab`'s `Z'` row.

---

## References

- **De Conti, A.; Duarte, A. C.; Alipio, R.** (2023). Closed-form ground-return
  impedance/admittance expressions used throughout `zg_form='deconti'`
  (eqs. 4.59, 4.63 as cited in `andreata_case1.py`/`andreata_case2.py`).
- **Andreata, L. E. B.** (2025). *Analise das Caracteristicas de Propagacao e
  de Transitorios Eletromagneticos em Cabos Subterraneos Instalados em Tubos
  Nao Metalicos no Contexto de Parques Eolicos.* Thesis, PPGEE/UFMG. Ch. 4
  (formulation), Fig. 5.4 (Configuration 4 geometry).
- [`HYBRID_PIPELINE_PLAN.md`](HYBRID_PIPELINE_PLAN.md) -- the FEM-hybrid
  pipeline this study was triggered from (sec. 9, Config. 4 generalization).
