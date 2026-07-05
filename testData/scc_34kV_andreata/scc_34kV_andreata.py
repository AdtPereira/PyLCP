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
        'frequencies': np.logspace(0, 7, num=121),
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
            'p100_vance': {
                'mtl': mtl_model,
                'zg_form': 'deconti',
                'yg_form': 'vance',
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
    plotter = SingleCoreCableModels(__file__, pul_data, autoSave=True)
    plotter.potential_coefficients_composition(conductor='core_sheath')
    plotter.potential_coefficients_composition(conductor='core')
    plotter.potential_coefficients_composition(conductor='sheath')
    plotter.potential_coefficients_earth_return()
    plotter.potential_coefficients_internal()
    plotter.series_impedance_composition(conductor='core_sheath')
    plotter.series_impedance_composition(conductor='core')
    plotter.series_impedance_composition(conductor='sheath')
    plotter.series_impedance_earth_return()
    plotter.series_impedance_internal()
    plotter.series_impedance_matrix()
    plotter.shunt_admittance_composition(condutor='core_sheath')
    plotter.shunt_admittance_composition(condutor='core')
    plotter.shunt_admittance_composition(condutor='sheath')
    plotter.shunt_admittance_earth_return()
    plotter.shunt_admittance_internal()
    plotter.shunt_admittance_matrix()
    GroundReturnMTLRepresentation(__file__, mtl_model, units='centimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()