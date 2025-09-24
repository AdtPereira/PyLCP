import sys
import os
import time
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# --- Configure project root for module imports (sem alteração) ---
try:
    os.system('cls' if os.name == 'nt' else 'clear')
    project_root = Path(__file__).resolve().parents[2]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root)) 
    print(f"Project root configured at: {project_root}")
    case_name = os.path.splitext(os.path.basename(__file__))[0]
    print(f"Case name identified as: '{case_name}'")
except IndexError:
    raise RuntimeError("Could not find project root. Ensure the directory structure is correct.")

# --- Import custom modules ---
try:
    from utils.case_utils import *
    from models import single_core_cables as scc 
    from plotter.patel_models import PatelModels
    from mtl_main.graphics import IsolatedMTLRepresentation
    from mtl_main.source import MulticonductorTransmissionLine
    from analytical_forms.single_core_cable import InternalPerUnitParameters
    from mom_so.quasi_static_green import QuasiStatic
    from mom_so.lossless_medium import HomogeneousLosslessMedium, LosslessPostProcessing
    print("Core modules imported successfully.")
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

# --- Load COMSOL Data ---
COMSOL_DATA = {}
try:
    COMSOL_DATA['core_sheath_return'] = load_comsol_results(__file__, comsol_tag='')
    COMSOL_DATA['core_exc'] = load_comsol_results(__file__, comsol_tag='_core_exc')
    COMSOL_DATA['core_exc_constrains'] = load_comsol_results(__file__, comsol_tag='_core_exc_constrains')
    COMSOL_DATA['sheath_exc'] = load_comsol_results(__file__, comsol_tag='_sheath_exc')
    COMSOL_DATA['shunt_params'] = load_comsol_results(__file__, comsol_tag='_shunt_params')
    print("COMSOL data loaded successfully.")

    # Display the first few rows of the loaded data to verify
    print("--- Data Head ---")
    print(COMSOL_DATA['core_sheath_return'].head())

    # Display a concise summary of the DataFrame
    print("\n--- DataFrame Info core_sheath_return---")
    COMSOL_DATA['core_sheath_return'].info()    
    print("\n--- DataFrame Info core_exc---")    
    COMSOL_DATA['core_exc'].info()
    print("\n--- DataFrame Info core_exc_constrains---")
    COMSOL_DATA['core_exc_constrains'].info()
    print("\n--- DataFrame Info sheath_exc---")
    COMSOL_DATA['sheath_exc'].info()
    print("\n--- DataFrame Info shunt_params---")
    COMSOL_DATA['shunt_params'].info()
except FileNotFoundError as e:
    print(f"Warning: COMSOL data file not found. Skipping comparison. Details: {e}")

def verify_constrain_equation(mtl_model):
    """
    Verifies the constraint equation by comparing the expected total current
    with the sum of the calculated conduction current (from Eq. 2.20) and
    the displacement current from COMSOL.
    """
    try:
        df = COMSOL_DATA['core_exc_constrains']
        # --- Prerequisite: Ensure you have exported these columns from COMSOL ---
        freq = df['freq'].values
        Ik = df['coil_current'].values      # Expected total current in the conductor
        Vs = df['coil_voltage'].values 
        int_Az = df['core_surface_int_az'].values
        int_Jdz = df['core_surface_int_jdz'].values
    except KeyError as e:
        print(f"\n[ERROR] Missing data column in the exported file: {e}")
        print("Please ensure you have exported the surface integrals of 'mf.Az' and 'mf.Jdz' from COMSOL.")
        return

    # --- Physical Parameters ---
    jw = 1j * 2 * np.pi * freq
    sigma = 1 / mtl_model.scc['core_resistivity']
    r2 = mtl_model.scc['core_outer_radius']
    Sck = np.pi * r2**2
    
    # --- Calculation of Current Components ---
    # Source current density from the coil voltage
    model_length = 1.0 # Length of the 2D model used for Jsk calculation
    Jsk = sigma * Vs / model_length
    
    # Conduction current, calculated using the terms from the QMS paper (Eq. 2.20)
    conduction_current = (-jw * sigma * int_Az) + (Sck * Jsk)
    
    # The total current is the sum of the conduction and displacement parts
    calculated_Ik = conduction_current + int_Jdz

    print("\n--- Verifying the Full Constraint Equation (Conduction + Displacement) ---")
    print(f"{'Frequency (Hz)':<16} | {'Ik Calculated (A)':<30} | {'Ik Expected (A)':<20} | {'Error Relative (%)':<20}")
    print("-" * 95)

    # Calculate and print the verification for each frequency
    # Note: We compare against the real part of expected_Ik for error calculation
    relative_errors = np.abs((calculated_Ik - Ik) / Ik) * 100
    for i in range(len(freq)):
        print(f"{freq[i]:<16.2f} | {calculated_Ik[i]:<30.4e} | {Ik[i]:<20.2f} | {relative_errors[i]:<20.4f}")

    # Final verification using numpy.allclose for numerical precision
    if np.allclose(calculated_Ik.real, Ik):
        print("\n[SUCCESS] The constraint equation was verified: I_conduction + I_displacement = Ik")
    else:
        print("\n[FAILURE] The constraint equation was NOT verified.")
        
def main():
    """ Main function to run the simulation and plotting using vectorized calculations. """
    st = time.time()    
    input_json = load_json_parameters(__file__, show_content=True)
    model = scc.isolated_coaxial_cable(input_json, show_model=True)
    mtl_model = MulticonductorTransmissionLine(model)

    # --- VECTORIZED CALCULATION ---
    analytical_freqs = np.logspace(0, 6, num=200)
    numerical_freqs = np.logspace(0, 6, num=31)

    # Analytical Formulation (Ametani et al., 2015)
    print("Calculating internal parameters for all frequencies...")
    internal = InternalPerUnitParameters(mtl_model, analytical_freqs)
    
    # MoM-SO formulation (Patel, 2014)
    print("Vectorized numeric routine (MoM-SO)...")
    green_matrix = QuasiStatic(mtl_model).green_matrix()
    mom_so = HomogeneousLosslessMedium(mtl_model, numerical_freqs)
    post_processor = LosslessPostProcessing(mtl_model)
    z_partial_stack = mom_so.z_partial(green_matrix)    # Partial impedance matrix
    zs_stack = post_processor.z_total(z_partial_stack)  # Total series impedance matrix
    verify_constrain_equation(mtl_model)

    # Populate the pul_data dictionary 
    pul_data = {
        'analytical': {
            "frequencies": analytical_freqs,
            "internal_series_parameters": {
                "bessel": internal.parameters_by_bessel(),
                "approximation": internal.parameters_approximation(),
                "hybrid": internal.parameters_hybrid(transition_frequency=1e5)
            },
            "internal_matrices": {
                "bessel": internal.internal_matrices(internal_form='bessel'),
                "approximation": internal.internal_matrices(internal_form='approximation'),
                "hybrid": internal.internal_matrices(internal_form='hybrid')
            },
        },
        'numerical': {
            "frequencies": numerical_freqs,
            "partial_impedance_matrix": z_partial_stack,
            "series_impedance_matrix": zs_stack,
            "series_resistance_matrix": post_processor.rs_matrix(zs_stack),
            "series_inductance_matrix": post_processor.ls_matrix(zs_stack, numerical_freqs)
        },
        'comsol': COMSOL_DATA
    }

    print(f"End of the routine! Time spent on simulation: {(time.time() - st):.1f} seconds.\n")
    plotter = PatelModels(pul_data, case_name)
    plotter.internal_impedance_matrix_js_method()
    plotter.internal_impedance_matrix_energy_method()
    plotter.internal_impedance_elements()
    plotter.internal_admittance_elements()
    IsolatedMTLRepresentation(mtl_model, case_name, units='millimeter').system_schematic()
    plt.show()    

if __name__ == "__main__":
    main()