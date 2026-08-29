# hdpe_225kV_2000mm2.py — Full Documentation

## Overview

Electromagnetic simulation script for a high-voltage power cable (225 kV, 2000 mm^2) buried and enclosed in an HDPE (High-Density Polyethylene) pipe. The central goal is to compute the per-unit-length (PUL) parameters — internal impedance and shunt admittance — via analytical formulations and compare them against reference results obtained in COMSOL Multiphysics.

---

## Physical Problem

A single-core cable (SCC) is buried in the soil, enclosed in an HDPE pipe. The geometry is eccentric: the cable axis does not coincide with the pipe axis. The system has:

- **Phase conductor (core):** central copper conductor (area 2000 mm^2).
- **Main insulation:** dielectric between the core and the metallic sheath.
- **Metallic sheath:** outer tubular conductor.
- **HDPE pipe:** mechanical-protection enclosure with an air gap between the sheath and the inner wall of the pipe.
- **Soil:** conducting return medium (modeled as a ground plane).

The presence of the HDPE creates a mixed region (air gap + HDPE) that must be reduced to an equivalent dielectric to allow the use of the standard analytical MTL (Multiconductor Transmission Line) formulations.

---

## Modeling Scenarios

The script defines four MTL models (`model_0` to `model_3`) representing different hypotheses about how to treat the region between the sheath and the soil:

| Model | Variable | Description |
|--------|----------|-----------|
| `model_0` | `mtl_0` | Eccentric cable with HDPE pipe — full real geometry. Used as the geometric reference for the ERS computation. |
| `model_1` | `mtl_1` | Buried cable without HDPE pipe — ignores the presence of the enclosure. Used as Scenario 1 (underground). |
| `model_2` | `mtl_2` | Cable without enclosure; the region between sheath and soil is replaced by an equivalent dielectric computed by the **area-weighting** method (area-weighted `epsr`). |
| `model_3` | `mtl_3` | Cable without enclosure; equivalent permittivity computed by the **GMD** (Geometric Mean Distance) method, case 3.1. |

Scenarios 1, 2 and 3 correspond respectively to `mtl_1`, `mtl_2` and `mtl_3` in the `pul_data['scenarios']` structure.

---

## Frequency Range

```python
frequencies = np.logspace(0, 6, num=31)  # 1 Hz to 1 MHz, 31 points on a logarithmic scale
```

---

## Dependencies

### Standard Library Modules

| Module | Use |
|--------|-----|
| `sys` | Exit on an import error |
| `os` | Terminal clear (`cls`/`clear`) |
| `copy` | Deep copy of model dictionaries |
| `time` | Execution-time measurement |
| `numpy` | Numerical operations and frequency generation |
| `matplotlib.pyplot` | Display of the plots at the end |

### PyLCP Modules

| Module | Class / Symbol | Responsibility |
|--------|------------------|-----------------|
| `utils.case_utils` | `*` (wildcard) | Formatting, matrix-printing and figure-saving utilities |
| `utils.comsol_data` | `ComsolPostProcessor` | Reading and processing of the results exported from COMSOL |
| `plotter.scc_plotter` | `HDPEPlotter` | Generation of the internal-impedance plots |
| `mtl_main.graphics` | `GroundReturnMTLRepresentation` | Generation of the cross-section schematics |
| `mtl_main.source` | `MulticonductorTransmissionLine` | Construction and validation of the MTL model |
| `models.single_core_cable` | `SingleCoreCableModelGenerator` | Generation of the model dictionaries from the case JSON |
| `analytical_forms.single_core_cable` | `InternalPerUnitParameters` | Analytical computation of the PUL matrices (impedance and admittance) |
| `analytical_forms.single_core_cable` | `EquivalentRadiiSystems` | Computation of the equivalent radii and permittivities |
| `.plot_config` | `PLOT_CONFIG` | Plot configuration (series, styles, axes) |

---

## Main Data Structure — `pul_data`

```python
pul_data = {
    'frequencies': np.ndarray,          # (31,) — frequency vector

    'comsol': {
        'frequencies': np.ndarray,       # frequencies from the COMSOL file
        'angular_frequencies': np.ndarray,
        'scenarios': {
            '1': {
                'coaxial_cable_impedance': dict,       # z11, z12, z2i, Zcs
                'internal_impedance_matrix': dict,     # js_method, energy_method
                'internal_impedance_elements': dict,   # self_core, self_sheath, mutual
            }
        }
    },

    'scenarios': {
        '1': {
            'mtl': MulticonductorTransmissionLine,     # underground model
            'internal_parameters': dict,               # output of parameters_hybrid()
            'internal_matrices': dict,                 # output of matrices()
        },
        '2': { ... },   # area-weighted ERS model
        '3': { ... },   # GMD ERS model
    },

    'analytical': {
        'frequencies': np.ndarray,
        'scenarios': { ... }            # reference to the same 'scenarios' dict
    }
}
```

---

## Execution Flow — `main()`

### 1. Generation of the Geometric Models

```python
model_generator = SingleCoreCableModelGenerator(__file__)
model_0 = model_generator.eccentric_hdpe_enclosed_model()
model_1 = model_generator.underground_model()
```

- `SingleCoreCableModelGenerator` loads the physical parameters from the `hdpe_225kV_2000mm2.json` file in the same directory.
- `eccentric_hdpe_enclosed_model()` returns a dictionary with conductors, insulations and enclosure (HDPE pipe) placed eccentrically.
- `underground_model()` returns the same cable without an enclosure, representing simple burial.

### 2. Construction of the MTL Objects

```python
mtl_0 = MulticonductorTransmissionLine(model_0)
mtl_1 = MulticonductorTransmissionLine(model_1)
```

`MulticonductorTransmissionLine` validates the model, extracts the conducting surfaces and computes the geometric properties used by the analytical formulations.

### 3. Computation of the Equivalent Radii Systems (ERS)

#### Area-Weighted Method -> `model_2`

```python
ers = EquivalentRadiiSystems(mtl_0)
epsr_area = ers.equiv_rel_permittivity_epsr_area_weighted()
r4 = epsr_area['sheath_outer_radius']
r7 = epsr_area['sheath_enclosure_outer_radius']

model_2 = copy.deepcopy(model_0)
model_2[2]['enclosure'] = None
model_2[2]['insulation']['thickness'] = r7 - r4
model_2[2]['insulation']['relative_permittivity'] = epsr_area['equivalent_relative_permittivity']
mtl_2 = MulticonductorTransmissionLine(model_2)
```

The equivalent permittivity is computed as an area-weighted average of the dielectric, air-gap and HDPE regions. The enclosure is removed and replaced by a single equivalent insulation layer of thickness `r7 - r4`.

#### GMD Method (Case 3.1) -> `model_3`

```python
ers = EquivalentRadiiSystems(mtl_0)
gmp = ers.equivalent_parameters_from_gmd()
eps_a = gmp['equivalent_relative_permittivity']['case 3.1']

model_3 = copy.deepcopy(model_0)
model_3[2]['enclosure'] = None
model_3[2]['insulation']['relative_permittivity'] = eps_a
mtl_3 = MulticonductorTransmissionLine(model_3)
```

Based on Lafaia's GMD methodology (2015). In case 3.1, the outer radius of the equivalent insulator is the cable outer radius (`r5`), without modifying the layer thickness.

### 4. Reading of the COMSOL Data

```python
cmsl_processor = ComsolPostProcessor(__file__)
cmsl_params = cmsl_processor.get_general_parameters('cmsl_coaxial_cable_impedance')
pul_data['comsol'].update(cmsl_params)

value['coaxial_cable_impedance'] = cmsl_processor.get_coaxial_cable_parameters()
scc_elements = cmsl_processor.get_scc_internal_impedance_elements()
value['internal_impedance_matrix'] = scc_elements
value['internal_impedance_elements'] = scc_elements
```

`ComsolPostProcessor` automatically locates the case's `Results/` directory and reads the exported files. The data obtained is:

| Key | Content |
|-------|----------|
| `coaxial_cable_impedance` | `z11`, `z12`, `z2i`, `Zcs` — coaxial-cable impedance elements |
| `internal_impedance_matrix` | 2x2 internal-impedance matrices (JS method and energy method) |
| `internal_impedance_elements` | Individual elements: `self_core`, `self_sheath`, `mutual` |

### 5. Analytical Computation of the PUL Parameters

```python
for key, value in pul_data['scenarios'].items():
    pul = InternalPerUnitParameters(value['mtl'], pul_data['frequencies'])
    value['internal_parameters'] = pul.parameters_hybrid()
    value['internal_matrices'] = pul.matrices()
```

For each analytical scenario (1, 2, 3):

- `parameters_hybrid()` — hybrid approach: Bessel functions below 100 kHz, approximations above. Returns the individual elements `z11`, `z12`, `z2i`, `Zcs`, etc.
- `matrices()` — assembles the full series-impedance and shunt-admittance matrices in the frequency domain.

**Return structure of `parameters_hybrid()`:**

```python
{
    'zcs': {'z11': array, 'z12': array, 'z2i': array, 'Zcs': array},
    'zs3': {'z20': array, 'z23': array},
    'z2m': array,
    'potentials': {'pcj': array, 'psj': array}
}
```

**Return structure of `matrices()`:**

```python
{
    'impedance_matrix':          np.ndarray,  # (31, N, N) — Ohm/m
    'resistance_matrix':         np.ndarray,  # (31, N, N) — Re(Z)
    'inductance_matrix':         np.ndarray,  # (31, N, N) — Im(Z)/(2*pi*f)
    'shunt_admittance_matrix':   np.ndarray,  # (31, N, N) — S/m
    'potential_coefficient_matrix': np.ndarray,  # (N, N) — frequency-independent
    'capacitance_matrix':        np.ndarray,  # (N, N) — F/m
}
```

### 6. Plot Generation

#### Internal-Impedance Plots

```python
plotter = HDPEPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=False)
plotter.hdpe_internal_impedance_matrix()
plotter.hdpe_internal_impedance_elements()
```

- `hdpe_internal_impedance_matrix()` — dual plot (resistance and inductance) comparing the three analytical scenarios with the COMSOL data for the PUL matrix elements (Z_cc, Z_cs, Z_ss).
- `hdpe_internal_impedance_elements()` — decomposes the impedance into individual contributions of each layer (z11, z12, z2i), comparing the JS method with the energy method.

#### Cross-Section Schematics

```python
schematic_configs = [
    {'mtl': mtl_0, 'filename': 'schematic_original_hdpe'},
    {'mtl': mtl_1, 'filename': 'schematic_ignored_hdpe'},
]
for config in schematic_configs:
    schematic = GroundReturnMTLRepresentation(
        __file__, config['mtl'], autoSave=True, units='millimeter'
    )
    schematic.system_schematic(base_filename=config['filename'])
```

Generates two cross-section schematics in `Results/`:
- `schematic_original_hdpe.*` — real section with the HDPE pipe.
- `schematic_ignored_hdpe.*` — section without the pipe (underground).

`autoSave=True` saves the files automatically; `autoSave=False` on the impedance plotters leaves control to `plt.show()` at the end.

---

## Generated Outputs

| File | Directory | Description |
|---------|-----------|-----------|
| `schematic_original_hdpe.png/svg` | `Results/` | Cross-section with HDPE pipe |
| `schematic_ignored_hdpe.png/svg` | `Results/` | Cross-section without pipe |
| Impedance plots | Interactive window | Displayed via `plt.show()` (not saved automatically) |

---

## Reference Geometric Parameters

The radii below identify the system's surfaces in increasing order from the cable axis:

| Symbol | Surface |
|---------|-----------|
| `r1` | Inner radius of the core |
| `r2` | Outer radius of the core |
| `r3` | Outer radius of the core insulation |
| `r4` | Outer radius of the metallic sheath |
| `r5` | Outer radius of the sheath insulation |
| `r6` | Inner radius of the HDPE pipe (start of the air gap) |
| `r7` | Outer radius of the HDPE pipe |

---

## Running

```bash
# From the project root directory (c:\git\PyLCP)
python -m testData.hdpe_225kV_2000mm2.hdpe_225kV_2000mm2
```

The terminal is cleared automatically at the start of each run (`main()`). The total time is reported at the end of the simulation.
