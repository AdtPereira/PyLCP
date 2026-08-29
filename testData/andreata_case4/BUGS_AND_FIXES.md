# Diagnostics and fixes — `andreata_case4` (Configuration 4, Figure 5.4)

**Context:** `andreata_case4` combines the HDPE duct of `andreata_case2` with the
heterogeneous ECC of `andreata_case3`. Creating it required a new generator
method (`flat_hdpe_enclosed_with_shared_ecc_model`) and a new `mtl_type`
(`scc-flat-hdpe-ecc`), done together with the standardization of the `"type"` field of
the four andreata JSONs (`"scc-flat"`, `"scc-flat-hdpe"`, `"scc-flat-ecc"`,
`"scc-flat-hdpe-ecc"`). Three real problems arose in that work, all
discovered and fixed within the first implementation.

---

## Item 1 — `'HDPE'` (uppercase) not registered after making `"type"` functional

**Where:** `mtl_main/strategy.py` (`mtl_strategy_factory`), `mtl_main/graphics.py`.

**Context:** when making `underground_flat_model()`/`flat_hdpe_enclosed_model()`
sensitive to the JSON `"type"` field (previously they hardcoded `'scc'`/`'hdpe'`,
ignoring the JSON entirely), it was confirmed by grep that none of the 8 sibling
cases that share those two methods had a `"type"` field at the JSON root
level — **except** `hdpe_300mm2.json`/`hdpe_2000mm2.json`, which already
had `"type": "HDPE"` (uppercase, not `"hdpe"`) for a long time, without it
mattering (the previous hardcode ignored it). Once the field was made functional,
that value started leaking into `model['type']`, and `'HDPE'` was not
registered anywhere — `mtl_strategy_factory` raised
`ValueError: Unknown or unsupported MTL type: HDPE`.

**Root cause:** the initial safety check (grep for the literal string
`"type"` in the 8 JSONs, using a glob with braces `{a,b,c}/*.json`)
returned "no result" incorrectly — the brace glob did not expand
as expected in the search tool used, masking the presence of
`hdpe_300mm2.json`/`hdpe_2000mm2.json` in the list. A second check,
using Python's JSON parser directly file by file (slower,
but reliable), revealed the `"HDPE"` value in the two files.

**Fix:** `'HDPE'` registered as an explicit alias in
`mtl_strategy_factory` (-> `SingleCoreCableInHDPEStrategy`, same class as
`'hdpe'`) and in the same points of `mtl_main/graphics.py` where `'hdpe'` already
appeared (title/`h_factor` branch and `depth_ref_conductor` exception tuple).

**Validated:** `hdpe_300mm2.py`/`hdpe_2000mm2.py` ran to completion again
without error, with the same plots/schematics as before.

**Lesson:** when checking "no other case uses X" via textual search before
a change that depends on it, prefer a real parser (e.g. `json.load`)
over a grep with a complex glob — a false negative here almost reached
implementation without that safety net.

---

## Item 2 — Heterogeneous dispatch restricted to the exact string `'scc-flat-ecc'`

**Where:** `analytical_forms/single_core_cable.py:944`
(`InternalPerUnitParameters.matrices`).

**Original symptom:** when running `andreata_case4.py` for the first time (after
fixing Item 1), the script reached `EquivalentRadiiSystems(mtl_ers_base)`
and broke with `KeyError: 'hdpe'` in
`analytical_forms/single_core_cable.py:265` — but the real error was not on
that line, but in a prior cause:

**Root cause (two layers):**
1. **`model_ers_base`** is built via `flat_hdpe_enclosed_model()` on top of
   the `SingleCoreCableModelGenerator` of `andreata_case4.json` — whose root-level
   `"type"` is `"scc-flat-hdpe-ecc"` (describes the case as a whole, for
   the heterogeneous model). Since `flat_hdpe_enclosed_model()` now reads
   `self.input_data.get('type', 'hdpe')` (Item 0 of the standardization), and the JSON
   *has* a `"type"`, the resulting homogeneous model inherited `type =
   'scc-flat-hdpe-ecc'` instead of `'hdpe'` — mapping to
   `SingleCoreCableWithECCStrategy` (heterogeneous) instead of
   `SingleCoreCableInHDPEStrategy`. The heterogeneous `context.scc` (grouped
   by cable) has no `'hdpe'` key, hence the `KeyError`.
   **Fixed** by forcing `model_ers_base['type'] = 'hdpe'` explicitly
   in `andreata_case4.py`, right after the generation — this model is always
   homogeneous by construction (never has an ECC), so the JSON `"type"` (which
   describes the heterogeneous case as a whole) does not apply to it.
2. **More serious consequence, only noticed while investigating the first:** for the
   same reason, `model_1`/`model_2`/`model_3` (built via
   `flat_scc_with_ecc_cable_model()`, used for the 3 analytical scenarios)
   *also* inherit `type = 'scc-flat-hdpe-ecc'` from the JSON, instead of the
   `'scc-flat-ecc'` that this method uses as a default. This does not break the
   strategy registration (both types map to
   `SingleCoreCableWithECCStrategy` in `mtl_strategy_factory` — deliberate,
   see Item 2 of "Changes" in the plan), but it would silently break
   `InternalPerUnitParameters.matrices()`: the dispatch to the
   heterogeneous path (`_matrices_heterogeneous`, which knows how to assemble
   `[2,2,2,1]` blocks) checked the exact string `mtl_type == 'scc-flat-ecc'` — with
   `mtl_type == 'scc-flat-hdpe-ecc'`, it would fall into the generic/homogeneous
   branch, incompatible with the per-cable grouped `context.scc` (same class of error
   as Item 4 of `andreata_case3/BUGS_AND_FIXES.md`, now via a different path).

**Fix:** generalized the dispatch in
`analytical_forms/single_core_cable.py:944` to
`if self.model.mtl_type in ('scc-flat-ecc', 'scc-flat-hdpe-ecc'):` — the correct
check is "does this `mtl_type` map to
`SingleCoreCableWithECCStrategy`?", not a specific string; any future
`mtl_type` that reuses this same strategy inherits the correct dispatch
automatically. Also fixed in the same function that
`ComsolPostProcessor.load_scc_earth_return_and_internal_scenarios` calls
internally (same `InternalPerUnitParameters.matrices`), without needing a
separate fix in `utils/comsol_data.py` (which already reads `block_sizes` from the
result, without checking `mtl_type` — Item 5 of `andreata_case3/BUGS_AND_FIXES.md`).

**Validated:** `andreata_case4.py` runs to completion without error; cross-validation
against `andreata_case3` (scenario `'1'`, duct ignored) with a max relative
error of 0.05% (Zs) / 0.12% (Ysh) — same order of magnitude as the equivalent
`andreata_case2` validation against `andreata_case1` (0.06%/0.15%).
Clean regression on `andreata_case1`, `andreata_case2`, `andreata_case3` and
on the 8 sibling cases of `underground_flat_model`/`flat_hdpe_enclosed_model`.

**Lesson:** when registering a new `mtl_type` that reuses an existing strategy
(`mtl_strategy_factory`), always grep for exact-string checks of that `mtl_type` in
other modules (`analytical_forms/`,
`utils/comsol_data.py`, `mtl_main/graphics.py`) — registration in the factory by
itself does not guarantee the whole pipeline recognizes the new name.

---

## Item 3 — ECC self term in COMSOL silently discarded by the internal-impedance parser

**Where:** `utils/comsol_data.py` (`ComsolPostProcessor.get_scc_internal_impedance_matrix_combined`).

**Context:** when adding `Results/cmsl_internal_impedance_matrix.txt`
(real COMSOL data, specific to `andreata_case4`, 3-conductor simulation —
core, sheath, ECC — instead of the legacy 2-conductor file reused by
`andreata_case2`/`hdpe_300mm2`), the plot
`internal_impedance_matrix.png` still had no COMSOL marker for the
ECC self element (`p=6, q=6`), even though the data existed in the file.

**Symptom:** no traceback — only a single console warning, `Warning:
COMSOL data not found for 'measured' with key 'impedance_matrix'.`
(`plotter/scc_plotter.py:189`), regarding only the ECC component (the
three phase-A components — `cc`/`cs`/`ss` — plotted normally).

**Root cause (two layers):**
1. `get_scc_internal_impedance_matrix_combined()` had `N = 2` hardcoded and
   read only 3 of the 9 data columns of the new file (`data1`, `data1_1`,
   `data2_1`) — the remaining 6, including the 3 that involve the ECC (core-ECC
   mutual, sheath-ECC mutual and the ECC self term, the latter in the
   last column of the file, `data3(mf.VCoil_ecc_i0)` -> cleaned column
   `data3_2`), were silently discarded — the returned matrix
   stayed `(freq, 2, 2)` regardless of whether the file had 4 or 9
   data columns.
2. Even reading those columns, a second incompatibility: the plotter uses
   a single `(p, q)` pair to index the three sources of the same component
   (analytical/COMSOL/MATLAB — `plotter/scc_plotter.py:170-202`), and the
   case4 `plot_config.py` asks for `p=6, q=6` (**global** index of the ECC in
   the 7-conductor space, same convention as MATLAB/analytical). A "local"
   3-conductor COMSOL matrix (core=0, sheath=1, ECC=2) would have
   `Zi[:, 6, 6]` out of bounds anyway — unlike
   `cc`/`cs`/`ss`, which only work because phase A coincides by chance with
   the global indices 0/1.

**Fix:** the function now detects the file format by the presence of the
`data1_2` column (only exists when there are 3 `data1(...)` columns, i.e.
when the ECC is present) and, in that case, assembles the matrix already in the
**global** 7-conductor space (`N = 7`, core=0, sheath=1, ECC=6 — same
convention as `MatlabDataReader.conductor_order` and the `p=6, q=6` used in
all of `plot_config.py`), filling only indices 0/1/6 (the rest stays
zero) from the 9 columns: `data1_2`/`data2_2` for the core-ECC/sheath-ECC
mutuals, `data3_2` for the ECC self term. For the legacy 2-conductor
file (`data1_2` absent), the behavior stays **identical** to before
(`N = 2`, same 3 assignments).

**Validated:**
- `andreata_case4.py`: the missing-COMSOL-data warning disappears; the
  plot `internal_impedance_matrix.png` now shows a 4th COMSOL marker
  (green, ECC) at `p=6, q=6`, which closely follows the ECC MATLAB
  points (also green) in both panels (R and L) — consistent
  with the same order of magnitude and curve shape, same agreement
  pattern already observed for `cc`/`cs`/`ss`.
- Clean regression (same file/function, 2-conductor format) on
  `andreata_case1`, `andreata_case2`, `andreata_case3`, `hdpe_300mm2`,
  `scc_132kV_xue`, `scc_34kV_andreata`, `scc_138kV_prysmian` — all the
  other consumers of `load_scc_earth_return_and_internal_scenarios`/
  `get_scc_internal_impedance_matrix_combined` in the repository.

**Lesson:** a parser with a hardcoded matrix dimension (`N = 2`) does not warn
when the input file grows (more columns) — it just reads what it already
expected and discards the rest silently. It pays to detect the format by
content (presence/absence of a key column) instead of assuming a fixed
size, especially in external-data parsers (COMSOL/MATLAB)
that have more than one file layout in use in the repository.
