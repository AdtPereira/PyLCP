import sys
import os
import copy
import time
import numpy as np
import matplotlib.pyplot as plt

# --- Import custom modules ---
try:
    from utils.case_utils import *
    from utils.comsol_data import ComsolPostProcessor
    from plotter.scc_plotter import HDPEPlotter
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from models.single_core_cable import SingleCoreCableModelGenerator
    from analytical_forms.single_core_cable import InternalPerUnitParameters
    from analytical_forms.single_core_cable import EquivalentRadiiSystems
    from .plot_config import PLOT_CONFIG
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def compare_capacitance(pul_data: dict, c_comsol: dict = None) -> None:
    """
    Print a comparison table of C_11 and C_22 (capacitance matrix diagonal)
    for the four proposed models plus the COMSOL reference (energy method).
    Both values are frequency-independent, expressed in µF/km.
    """
    rows = {}
    labels = {'1': 'Underground', '2': 'Area-weighted ERS', '3': 'GMD case 3.1'}
    for key, label in labels.items():
        cap = pul_data['scenarios'][key]['internal_matrices']['capacitance_matrix']
        rows[label] = {'c11': cap[0, 0], 'c22': cap[1, 1]}

    col_w = 26
    print("\n" + "=" * 62)
    print("  Per-Unit-Length Capacitance  [uF/km]")
    print("=" * 62)
    print(f"  {'Model':<{col_w}} {'C_11':>12}  {'C_22':>12}")
    print("-" * 62)
    for label, vals in rows.items():
        print(f"  {label:<{col_w}} {vals['c11'] * 1e9:>12.4f}  {vals['c22'] * 1e9:>12.4f}")
    if c_comsol is not None:
        print("-" * 62)
        print(f"  {'COMSOL (energy method)':<{col_w}} {c_comsol['c11'] * 1e9:>12.4f}  {c_comsol['c22'] * 1e9:>12.4f}")
    print("=" * 62 + "\n")

def main():
    """
    Main function to run the simulation and plotting.
    """
    os.system('cls' if os.name == 'nt' else 'clear')
    st = time.time()    
    model_generator = SingleCoreCableModelGenerator(__file__)    
    model_0 = model_generator.eccentric_hdpe_enclosed_model()
    model_1 = model_generator.underground_model()
    mtl_0 = MulticonductorTransmissionLine(model_0)
    mtl_1 = MulticonductorTransmissionLine(model_1)

    # Equivalent Radii Systems (ERS) for shunt parameters
    print("Calculating equivalent radii systems for shunt parameters...")
    model_2 = copy.deepcopy(model_0)
    ers = EquivalentRadiiSystems(mtl_0)
    epsr_area = ers.equiv_rel_permittivity_epsr_area_weighted()
    r4, r7 = epsr_area['sheath_outer_radius'], epsr_area['sheath_enclosure_outer_radius']
    
    # Remove enclosure for Model Case 2
    model_2[2]['enclosure'] = None  
    model_2[2]['insulation']['thickness'] = r7 - r4
    model_2[2]['insulation']['relative_permittivity'] = epsr_area['equivalent_relative_permittivity']
    mtl_2 = MulticonductorTransmissionLine(model_2)

    # Case 3.1: The outer radius r0 of the equivalent insulator is taken as the cable outer radius r5 in Fig. 1 (b) 
    print("Calculating equivalent radii systems for shunt parameters...")
    model_3 = copy.deepcopy(model_0)
    ers = EquivalentRadiiSystems(mtl_0)
    gmp = ers.equivalent_parameters_from_gmd()
    eps_a = gmp['equivalent_relative_permittivity']['case 3.1']
    model_3[2]['enclosure'] = None
    model_3[2]['insulation']['relative_permittivity'] = eps_a
    mtl_3 = MulticonductorTransmissionLine(model_3)

    pul_data = {
        'frequencies': np.logspace(0, 6, num=31),
        'comsol': {
            'scenarios': {
                '1': {},
            },
        },
        'scenarios': {
            '1': {
                'mtl': mtl_1,
            },
            '2': {
                'mtl': mtl_2,
            },
            '3': {
                'mtl': mtl_3,
            },
        }
    }

    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_params = cmsl_processor.get_general_parameters('cmsl_coaxial_cable_impedance')
    pul_data['comsol'].update(cmsl_params)
    for key, value in pul_data['comsol']['scenarios'].items():
        print(f"  -> Processando COMSOL para: {key}")
        value['coaxial_cable_impedance'] = cmsl_processor.get_coaxial_cable_parameters()
        scc_elements = cmsl_processor.get_scc_internal_impedance_elements()
        value['internal_impedance_matrix'] = scc_elements
        value['internal_impedance_elements'] = scc_elements
    
    pul_data['analytical'] = {
        'frequencies': pul_data['frequencies'],
        'scenarios': pul_data['scenarios'],
    }

    print("\nCalculating per-unit-length parameters and quasi-TEM matrices for all scenarios...")
    for key, value in pul_data['scenarios'].items():
        print(f"Calculating internal parameters for Model Case {key}...")
        pul = InternalPerUnitParameters(value['mtl'], pul_data['frequencies'])
        value['internal_parameters'] = pul.parameters_hybrid()
        value['internal_matrices'] = pul.matrices()
            
    compare_capacitance(pul_data, cmsl_processor.get_shunt_capacitance_elements())
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = HDPEPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=True)
    plotter.hdpe_internal_impedance_matrix()
    plotter.hdpe_internal_impedance_elements()

    # 1. Crie uma lista de configurações para cada esquemático
    schematic_configs = [
        {'mtl': mtl_0, 'filename': 'schematic_original_hdpe'},
        {'mtl': mtl_1, 'filename': 'schematic_ignored_hdpe'},
    ]

    for config in schematic_configs:
        schematic = GroundReturnMTLRepresentation(
            __file__,
            config['mtl'], 
            autoSave=True, 
            units='millimeter'
        )
        schematic.system_schematic(base_filename=config['filename'])
    plt.show()
    
if __name__ == "__main__":
    main()