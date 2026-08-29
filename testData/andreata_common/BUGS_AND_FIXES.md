# Diagnostics and fixes — COMSOL internal-admittance data import (`cmsl_internal_admittance_matrix.txt`)

**Context:** this document covers work transversal to `andreata_case1` and
`andreata_case2` (shared code in `utils/comsol_data.py` and
`plotter/scc_plotter.py`): adding the COMSOL internal-admittance reference
(`Yi`, core + sheath) to the `internal_admittance_matrix.png` plots,
mirroring the internal-impedance import (`cmsl_internal_impedance_matrix.txt`)
already existing and documented in `andreata_case1/BUGS_AND_FIXES.md` (Bugs 1–4).

The most relevant item (Item 2) records a case where the **COMSOL simulation
file itself was physically incorrect** — not the Python import code — and
how that led to two incoherent fix attempts in the code before the root cause
was identified and the simulation fixed directly in COMSOL. It is left here as a
method record: facing a number that does not match the analytical reference,
first check whether the raw data (exported `.txt` file) changed between debug
sessions, before piling ever more elaborate reconstruction formulas onto the code.

Items, in chronological order:
1. new loading function (`get_scc_internal_admittance_matrix_combined`);
2. `Yi_22` (`C_ss`) diverging ~8.5x from the analytical value — root cause in the COMSOL simulation;
3. typo: invalid color `'tab:black'`;
4. COMSOL markers invisible against analytical lines of the same color;
5. analytical curves drawn behind the numerical markers (zorder);
6. `andreata_case2`: `KeyError` — `'internal_admittance_matrix'` absent from `plot_config.py`;
7. dimension error — `ValueError: x and y must be the same size` (frequency vectors conflated);
8. comparison of the pyLCP x MATLAB x COMSOL frequency vectors — no action required.

---

## Item 1 — New loading function: `get_scc_internal_admittance_matrix_combined()`

**Where:** `utils/comsol_data.py` (`ComsolPostProcessor`)

**What:** there was no reader for `cmsl_internal_admittance_matrix.txt`. The
file follows the same combined-table convention (two excitations in a single
file) already used by `cmsl_internal_impedance_matrix.txt`
(`get_scc_internal_impedance_matrix_combined()`), but with its own column
names instead of positional placeholders `data1`/`data2`: core-excitation
block (`Ccc`, `Csic`, `Csoc`, `Ccc_energy`, `Wcc`) and sheath-excitation
block (`Ccs`, `Csis`, `Csos`, `Css_energy`, `Wss`, `Wcs`).

**Fix:** new function, following the same convention of a mutual term "read
from the core-excitation side" already used for impedance:
- `C_cc` (core, self) = `Ccc` — direct reading, core excitation;
- `C_cs` (core-sheath, mutual) = `Csic` — direct reading, same excitation;
- `C_ss` (sheath, self) — see Item 2 (it was not trivial).

Result: `Yi = jw*C`, returned in the standard format
`{'frequencies': ..., 'scenarios': {'measured': {'admittance_matrix': Yi}}}`,
consumed by `SCCPlotter._plot_scc_internal_matrix()` via
`'comsol_matrix_key': 'admittance_matrix'` in `plot_config.py`.

Wired into the existing pipeline via `ComsolPostProcessor.load_scc_earth_return_and_internal_scenarios()`
(already called by `andreata_case1.py`/`andreata_case2.py`), alongside the
internal-impedance block.

---

## Item 2 — `Yi_22` (`C_ss`) diverged ~8.5x from the analytical value: root cause was the COMSOL simulation, not the code

**Where:** the simulation file `cmsl_internal_admittance_matrix.txt` (external
data source, generated in COMSOL) — not a Python function.

**What:** the analytical value of `Yi_22` (sheath self-admittance) is
`jw*(Cc+Cs)`, where `Cc` is the core<->sheath capacitance (core XLPE
insulation) and `Cs` is the sheath self-capacitance through its own
jacket (PVC) — extracted directly from the code
(`InternalPerUnitParameters._build_cable_block`, branch "sheath_insulation"):
`Pij = [[pcj+psj, psj], [psj, psj]]`, `Ye = jw*Pi^-1`, giving `Ye[0,0]=1/pcj=Cc`
and `Ye[1,1]=1/pcj+1/psj=Cc+Cs`. Numerically, for `andreata_case1`:
`Cc ~ 2.1466e-10 F/m`, `Cs ~ 1.6202e-9 F/m` (~7.5x larger than `Cc`, consistent
with the PVC jacket being much thinner than the XLPE insulation), `Cc+Cs ~ 1.8348e-9 F/m`.

In the first version of the COMSOL file, the direct reading `Csis` (sheath
charge under its own excitation) reproduced only `Cc` (~2.1467e-10 F/m) —
**not** `Cc+Cs`. `Csoc`/`Csos` (charge on the outer surface of the sheath, in the
two excitations) were ~1e-22 — numerical noise, essentially zero.

**Incorrect diagnosis #1 (incoherent operation):** it was concluded that the
COMSOL model of the PVC jacket simply was not included in that evaluation (since
no charge appeared on the outer surface in either excitation) and the
COMSOL reference was **removed** from the `Yi_22` component in
`plot_config.py`, so as not to overlay a physically incompatible quantity on the
analytical curve.

**Incorrect diagnosis #2 (incoherent operation):** revisiting the divergence,
it was noticed that the sheath-excitation block has an extra column
(`Wcs`, absent from the core-excitation block) and that
`C_ss = Css_energy + 4*Wcs` reproduced the analytical value to ~7e-4 %. That
formula (based on the energy method, with the x4 conversion factor already implicit
in `Ccc_energy = 4*Wcc`) was implemented and the COMSOL reference was
**restored** in `Yi_22` — even against the user's previous guidance not to use
the energy-method columns "for now".

**Real root cause:** between one debug session and another, the **COMSOL
file was re-exported** (timestamp changed from `Aug 10 2026, 16:44` to
`21:35`) after the user fixed the simulation itself — presumably the
boundary condition that made `Csos` not capture the charge on the outer jacket of
the sheath. In the corrected file, `Csos` started reading directly
`~1.6202e-9` (`Cs` itself) and `Css_energy` started reading directly
`~1.8349e-9` (already the full `Cc+Cs`) — the information missing in the old
file became present directly in the charge reading, without needing any
energy-based reconstruction.

This was only noticed because the formula `Css_energy + 4*Wcs` (calibrated for the
*old* file) started **adding the jacket contribution twice** over the
*new* file (which already contains it in `Css_energy`), producing
`~3.46e-9 F/m` — almost double the expected analytical value
(`~1.835e-9 F/m`) — and becoming visually obvious in the plot (COMSOL circles
of `Yi_22` off the analytical line, at twice the expected height).

**Final fix:** `C_ss = Csis + Csos` — pure direct method (sum of the charge
on the two surfaces of the sheath under its own excitation: inner, facing
the core, and outer, facing the jacket), **without** any energy
column. It matches the analytical value to ~3e-4 % (even more precise than the
energy reconstruction), and is consistent with the user's original guidance
of using only the direct method. This is the current and final implementation of
`get_scc_internal_admittance_matrix_combined()`.

**Lesson recorded:** when investigating a numerical divergence between COMSOL and the
analytical reference, first check whether the raw `.txt` file changed
(timestamp, key-column values) before assuming that the physics of the
reading method (direct vs. energy) needs adjustment — the cause may be
entirely outside the Python code, in the source simulation.

---

## Item 3 — Typo: invalid color `'tab:black'`

**Where:** `testData/andreata_case1/plot_config.py`

**What:** when trying to change the COMSOL marker color to black, `'tab:red'`
was replaced by `'tab:black'` — an invalid name in matplotlib (the prefix
`tab:` only exists for the 10 colors of the Tableau palette; pure black is just
`'black'`, with no prefix). Result: `ValueError: 'tab:black' is not a valid
color value.`, aborting the routine before saving
`internal_admittance_matrix.png`.

**Fix:** `'tab:black'` -> `'black'` in all occurrences.

---

## Item 4 — COMSOL markers invisible against analytical lines of the same color

**Where:** `testData/andreata_case1/plot_config.py` (component `Yi_11`)

**What:** the COMSOL markers use `facecolors='none'` (hollow circle) with
the edge color equal to the color of the component's own analytical line. For
`Yi_11` (solid black line, marker edge also black), a hollow circle over a solid
line of the same color is visually imperceptible — the line "passes through" the
marker. For `Yi_12`/`Yi_22` this does not happen because the lines are
dashed/dash-dot: the gaps between dashes let the colored ring show even being
the same color.

**Fix:** `facecolors='white'` (opaque fill) specifically on the COMSOL
marker of `Yi_11`, "punching" a visible gap in the solid line instead of merging
with it.

---

## Item 5 — Analytical curves drawn behind the numerical markers

**Where:** `testData/andreata_case1/plot_config.py`
(`internal_impedance_matrix` and `internal_admittance_matrix`)

**What:** the COMSOL/MATLAB `scatter` already had an explicit `zorder` (10/11);
the analytical line (`ax.plot(...)`, via `internal_style`) had no `zorder`
defined, assuming the matplotlib default (~2) — being drawn **behind**
the reference markers.

**Fix:** `'zorder': 12` added to all `internal_style` dicts
(higher than the markers' 10/11), making the analytical lines draw
on top.

---

## Item 6 — `andreata_case2`: `KeyError` — `'internal_admittance_matrix'` absent from `plot_config.py`

**Where:** `testData/andreata_case2/plot_config.py`;
`testData/andreata_case2/andreata_case2.py`

**What:** `andreata_case2.py` already called
`plotter.compare_internal_matrices(key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])`,
but `plot_config.py` never had the key `'internal_admittance_matrix'` —
a `KeyError` aborted the routine after saving only `internal_impedance_matrix.png`.
`internal_admittance_matrix.png` had never been generated successfully for
this case.

**Fix:** new entry `'internal_admittance_matrix'` mirroring the
structure/conventions already used by the sibling entry `'internal_impedance_matrix'`
in this same file (`INTERNAL_COMSOL_TEMPLATE` with a single `'measured'` scenario,
MATLAB series `'measured'`, analytical line left commented out as in the sibling
entry) — instead of literally copying the style of `andreata_case1`, which follows
its own conventions (multiple soil scenarios, active analytical line).

---

## Item 7 — Dimension error: `ValueError: x and y must be the same size`

**Where:** `utils/comsol_data.py` (`load_scc_earth_return_and_internal_scenarios`),
`plotter/scc_plotter.py` (`_plot_scc_internal_matrix`)

**What:** `pul_data['comsol']['frequencies']` was a single shared key,
written both by the impedance block and by the
admittance block (the admittance one used `setdefault`, so it silently lost
whenever the impedance one ran first). Since
`cmsl_internal_impedance_matrix.txt` and `cmsl_internal_admittance_matrix.txt`
are two independent COMSOL files, with independent parametric
sweeps, nothing guaranteed they had the same number of points — and that
premise broke as soon as the `andreata_case2` admittance file started
having 91 points while the impedance file stayed at 46. Result:
`ax1.scatter(cmsl_freq, ...)` tried to pair an x-axis of 46 points with a
y-axis of 91, raising `ValueError: x and y must be the same size` and
aborting the routine before saving `internal_admittance_matrix.png`.

**Fix:** the single key `frequencies` was replaced by
`frequencies_by_key`, a dict indexed by the matrix name
(`'impedance_matrix'`/`'admittance_matrix'`), each keeping its own
independent frequency vector. `_plot_scc_internal_matrix()` now looks up the
frequency vector corresponding to the specific `comsol_matrix_key` being
plotted, instead of assuming a single shared vector. The plots of the
ground-return family (`_plot_scc_matrix`, used by
`earth_return_*`/`self_*_phase_a_sheath`) were not changed, since there the
sharing of a single per-scenario frequency grid is legitimate —
it comes from a single COMSOL ground-return file.

**Validated:** the three `andreata_caseN.py` run end to end — exit 0,
no errors. `andreata_case1`/`andreata_case2`: impedance (46 pts) and
admittance (91 pts) now correctly independent; `andreata_case3`: without
COMSOL data (files missing), handled gracefully with warnings, without aborting
the routine.

---

## Item 8 — Comparison of the pyLCP x MATLAB x COMSOL frequency vectors — no action required

**Where:** `testData/andreata_case1/andreata_case1.py` (`pul_data['frequencies']`),
`andreata_case1_frequency_range.mat`, `cmsl_internal_admittance_matrix.txt`.

**What (informative record, no bug):**
- pyLCP (`np.logspace(-2, 7, num=90)`) and MATLAB (`.mat` file) agree to
  ~4e-15 relative difference — single-precision export rounding
  of the `.mat`, not a real difference in the sweep design. No
  adjustment needed.
- COMSOL uses its own sampling grid (46 or 91 points, depending on the
  file/case) and is plotted as an independent `scatter`, without needing
  to coincide point by point with the analytical/MATLAB line. The user decided
  (2026-08-10) to keep it that way — no interpolation or forced regeneration of
  the COMSOL grid.
- Side note: "10 points per decade" over 9 decades (1E-2 to 1E7) gives
  exactly 91 points (`9x10+1`, `logspace(-2,7,91)` verified
  numerically); "90 points" (used today by pyLCP/MATLAB) gives ~9.889
  points/decade, not exactly 9 nor exactly 10.
