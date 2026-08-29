'''
PRYSMIAN_138kV_CORE_SHEATH Cable Data (in meters):
1. CONDUCTOR: Compact circular copper stranded conductor, per the requirements of standard NBR NM 280 (class 2). Nominal cross-section: 500 mm2
    Nominal diameter: 25.95 mm
2. CONDUCTOR TAPING: Semiconducting tape containing swelling powder plus a nylon tape, both applied helically over the conductor.
    Nominal diameter: 26.73 mm
3. CONDUCTOR SHIELD: Extruded layer of XLPE-based semiconducting compound. Nominal thickness: 1.5 mm
    Nominal diameter: 29.73 mm
4. INSULATION: Extruded layer of cross-linked polyethylene (XLPE). Nominal thickness: 13.31 mm
    Nominal diameter: 60.35 mm
    Nominal relative permittivity: 2.3
5. INSULATION SHIELD: Extruded layer of XLPE-based semiconducting compound. Nominal thickness: 1.5 mm
    Nominal diameter: 63.35 mm
6. INSULATION TAPING: Semiconducting tape containing swelling powder, applied helically over the insulation shield.
    Nominal diameter: 64.63 mm
7. METALLIC SHEATH: Extruded lead-alloy sheath. Nominal thickness: 3.00 mm
    Nominal diameter: 70.63 mm
    Nominal cross-section: 637.4 mm2
8. OUTER JACKET: Extruded layer of high-density polyethylene (HDPE) containing termite-protection additive and graphite powder.
    Nominal thickness: 4.0 mm
    Nominal diameter: 78.63 mm

ELECTRICAL PROPERTIES
1. RMS PHASE-TO-EARTH VOLTAGE (kV): 79.69
2. RMS PHASE-TO-PHASE VOLTAGE (kV): 138
3. BASIC IMPULSE LEVEL (BIL) (kV): 650
4. MAXIMUM DC CONDUCTOR RESISTANCE AT 20 deg C (ohm/km): 0.0366
5. CAPACITANCE (uF/km): 0.1805
'''

import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt

# --- Import custom modules ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    from utils.case_utils import *
    from utils.comsol_data import ComsolPostProcessor
    from plotter.scc_models import PrysmianCableModels
    from .plot_config import PLOT_CONFIG
    from models.single_core_cable import SingleCoreCableModelGenerator
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import (
        InternalPerUnitParameters, PerUnitParameters, apply_semiconducting_layer_correction,
    )
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    cable_generator = SingleCoreCableModelGenerator(__file__)
    model = cable_generator.underground_flat_model()
    model = apply_semiconducting_layer_correction(model, cable_generator.core, cable_generator.sheath)
    mtl_model = MulticonductorTransmissionLine(model)

    print("Loading COMSOL internal impedance results...")
    cmsl_processor = ComsolPostProcessor(__file__)
    scc_internal_cmsl = cmsl_processor.get_scc_internal_impedance_matrix_combined()

    pul_data = {
        'comsol': scc_internal_cmsl,
        'frequencies': np.logspace(0, 7, num=121),
        'logger_data': {
            'frequencies': [1e2, 1e4, 1e5],
            'scenario': 'ametani'
        },
        'scenarios': {
            'magalhaes_xue': {
                'mtl': mtl_model,
                'zg_form': 'magalhaes_xue'
            },
            'sunde': {
                'mtl': mtl_model,
                'zg_form': 'sunde'
            },
            'pollaczek': {
                'mtl': mtl_model,
                'zg_form': 'pollaczek'
            },
            'ametani': {
                'mtl': mtl_model,
                'zg_form': 'ametani'
            },
            'deconti': {
                'mtl': mtl_model,
                'zg_form': 'deconti'
            },
            'saad': {
                'mtl': mtl_model,
                'zg_form': 'saad'
            },
        }
    }

    print("Calculating internal parameters for all frequencies...")
    pul = InternalPerUnitParameters(mtl_model, pul_data['frequencies'])
    internal_matrices = pul.matrices()
    pul_data['internal_matrices'] = internal_matrices
    pul_data['internal_parameters'] = pul.parameters_hybrid()

    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        
        earth_return = pul.earth_return_parameters(value['zg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)

        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem

    # Alias consumed by BasePlotter's generic engine (source='analytical'),
    # used by ground_return_impedance().
    pul_data['analytical'] = {'frequencies': pul_data['frequencies'], 'scenarios': pul_data['scenarios']}

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = PrysmianCableModels(__file__, pul_data, PLOT_CONFIG)
    plotter.internal_impedance_parameters(graph_key='core')
    plotter.internal_impedance_parameters(graph_key='sheath')
    plotter.internal_impedance_parameters(graph_key='core_sheath')
    plotter.internal_impedance_parameters(graph_key='internal_parameters')
    plotter.ground_return_impedance()
    GroundReturnMTLRepresentation(__file__, mtl_model, units='centimeter').system_schematic()
    plt.show()

if __name__ == "__main__":
    main()