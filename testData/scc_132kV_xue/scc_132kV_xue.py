import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt

# --- Import custom modules ---
os.system('cls' if os.name == 'nt' else 'clear')
try:
    from utils.case_utils import *
    from plotter.scc_models import SingleCoreCableModels
    from models.single_core_cable import SingleCoreCableModelGenerator
    from mtl_main.graphics import GroundReturnMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters, PerUnitParameters
    print("Modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    model = SingleCoreCableModelGenerator(__file__).underground_model()
    mtl_model = MulticonductorTransmissionLine(model)

    # 1. ESTRUTURA DE DADOS CENTRALIZADA
    pul_data = {
        'comsol': None, # cmsl_reader.data,
        'frequencies': np.logspace(4, 7, num=121),
        'scenarios': {
            'p100_xue': {
                'mtl': mtl_model,
                'zg_form': 'magalhaes_xue',
                'yg_form': 'magalhaes_xue',
            },
            'p100_deconti': {
                'mtl': mtl_model,
                'zg_form': 'deconti',
                'yg_form': 'deconti',
            },
        }
    }
    
    print("Calculating internal parameters for all frequencies...")
    pul = InternalPerUnitParameters(mtl_model, pul_data['frequencies'])
    internal_matrices = pul.matrices()
    pul_data['internal_matrices'] = internal_matrices

    for key, value in pul_data['scenarios'].items():
        print(f"Calculating scenario: {key}...")
        pul = PerUnitParameters(value['mtl'], pul_data['frequencies'])
        
        earth_return = pul.earth_return_parameters(value['zg_form'], value['yg_form'])
        quasi_tem = pul.quasi_tem_approx_matrices(internal_matrices, earth_return)

        value['earth_return_parameters'] = earth_return
        value['quasi_tem_matrices'] = quasi_tem
    
    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = SingleCoreCableModels(__file__, pul_data)
    plotter.series_impedance_matrix(condutor='core')
    plotter.series_impedance_matrix(condutor='sheath')
    plotter.series_impedance_matrix(condutor='core_sheath')
    plotter.shunt_admittance_matrix(condutor='sheath')
    plotter.potential_coefficients_matrix(condutor='core')
    plotter.potential_coefficients_matrix(graph_form='real_and_imaginary', condutor='core')
    plotter.earth_return_admittance_matrix()
    plotter.earth_return_impedance_matrix()
    plotter.earth_return_potential_coefficients_matrix()
    plotter.internal_admittance_matrix()
    plotter.internal_impedance_matrix()
    plotter.internal_potential_coefficients_matrix()
    GroundReturnMTLRepresentation(__file__, mtl_model, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()