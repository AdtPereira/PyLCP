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
    from models.single_core_cable import SingleCoreCableModelGenerator
    from plotter.scc_plotter import CoaxialCablePlotter
    from mtl_main.graphics import IsolatedMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters
    from mom_so.quasi_static_green import QuasiStatic
    from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing
    from .plot_config import PLOT_CONFIG
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    model = SingleCoreCableModelGenerator(__file__).isolated_coaxial_cable()
    mtl = MulticonductorTransmissionLine(model)

    pul_data = {
        'frequencies': np.logspace(0, 6, num=121),
        'comsol': {
            'scenarios': {
                '1': {},
            },
        },
        'mom_so': {
            'frequencies': np.logspace(0, 6, num=31),
            'scenarios': {
                '1': {
                    'mtl': mtl,
                },
            },
        },
        'scenarios': {
            '1': {
                'mtl': mtl,
            },
        }
    }

    print("Importing COMSOL data for coaxial cable model...")
    cmsl_processor = ComsolPostProcessor(__file__)
    cmsl_params = cmsl_processor.get_general_parameters('cmsl_coaxial_cable_impedance')
    pul_data['comsol'].update(cmsl_params)
    for key, value in pul_data['comsol']['scenarios'].items():
        print(f"  -> Processando COMSOL para: {key}")
        value['coaxial_cable_impedance'] = cmsl_processor.get_coaxial_cable_parameters()

    print("\nCalculating per-unit-length parameters by Analytical Formulation (Ametani, 2015)")
    for key, value in pul_data['scenarios'].items():
        print(f"  -> Calculating internal parameters for Model Case {key}...")
        
        pul = InternalPerUnitParameters(value['mtl'], pul_data['frequencies'])
        value['internal_parameters'] = pul.parameters_hybrid()
        value['internal_matrices'] = pul.matrices()

    print("\nCalculating per-unit-length parameters by MoM-SO (Patel, 2014)")
    for key, value in pul_data['mom_so']['scenarios'].items():
        print(f"  -> Calculating internal impedance matrix for Model Case {key}...")

        green_matrix = QuasiStatic(value['mtl']).green_matrix()
        mom_so = HomogeneousLosslessMedium(value['mtl'], pul_data['mom_so']['frequencies'])
        post_processor = LosslessPostProcessing(value['mtl'])
        z_partial_stack = mom_so.z_partial(green_matrix)

        value['partial_internal_impedance'] = z_partial_stack
        value['coaxial_cable_impedance'] = post_processor.z_total(z_partial_stack)

    print(f"\nEnd of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = CoaxialCablePlotter(__file__, pul_data, PLOT_CONFIG, autoSave=False)
    plotter.coaxial_cable_internal_impedance_elements()
    IsolatedMTLRepresentation(__file__, mtl, units='millimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()