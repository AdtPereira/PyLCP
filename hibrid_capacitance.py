import numpy as np
import math
import scipy.constants as sc

def calculate_hybrid_capacitance(r4, r5, R1, R2, d, epsr1, epsr2, epsr3):
    """
    Calculates the equivalent hybrid capacitance per unit length for a
    three-layer cable with a central eccentric air gap.

    The model consists of three capacitors in series:
    1. A concentric layer from r_a to r_b.
    2. An eccentric layer from r_b to r_c.
    3. A concentric layer from r_c to r_d.

    Args:
        r_a (float): Radius of the central conductor [meters].
        r_b (float): Outer radius of the first concentric layer [meters].
        r_c (float): Outer radius of the eccentric air gap [meters].
        r_d (float): Inner radius of the outer sheath/conductor [meters].
        d (float): Eccentricity of the air gap layer [meters].
        epsilon_r1 (float): Relative permittivity of the first dielectric layer.
        epsilon_r_ar (float): Relative permittivity of the air gap (typically ~1.0).
        epsilon_r3 (float): Relative permittivity of the third dielectric layer.

    Returns:
        float: The equivalent capacitance per unit length in Farads per meter (F/m).
               Returns float('inf') if the geometry is invalid.
    """
    # --- Input Validation ---
    if not (0 < r4 < r5 < R1 < R2):
        print("Error: Radii must be positive and in increasing order (r_a < r_b < r_c < r_d).")
        return float('inf')
    if not (d >= 0 and d <= (R1 - r5)):
        print(f"Error: Eccentricity 'd' must be non-negative and less than the available gap ({R1 - r5:.4f} m).")
        return float('inf')
    if not (epsr1 >= 1 and epsr2 >= 1 and epsr3 >= 1):
        print("Error: Relative permittivities must be >= 1.")
        return float('inf')

    # --- Calculate each term of the denominator (the "elastance" of each layer) ---
    # Term 1: First concentric layer
    term1 = np.log(r5 / r4) / (epsr1 * sc.epsilon_0)

    # Term 2: Eccentric air gap layer
    # Argument for the acosh function
    acosh_arg = (r5**2 + R1**2 - d**2) / (2 * r5 * R1)
    term2 = math.acosh(acosh_arg) / (epsr2 * sc.epsilon_0)

    # Term 3: Third concentric layer
    term3 = np.log(R2 / R1) / (epsr3 * sc.epsilon_0)

    # --- Calculate the final equivalent capacitance ---
    return (2 * np.pi) / (term1 + term2 + term3)

# --- Example Usage ---
if __name__ == "__main__":
    print("Cálculo da Capacitância Equivalente Híbrida de um Cabo Exemplo")
    print("-" * 60)

    # Material properties (relative permittivity)
    er1 = 2.3
    er_ar = 1.0
    er3 = 2.3

    # --- Convert units to meters for the calculation function ---
    r4 = 0.05720 # Raio do condutor central
    r5 = 0.06220 # Raio da primeira camada de isolação
    R1 = 0.09925 # Raio interno da camada excêntrica
    R2 = 0.11250 # Raio externo da camada excêntrica
    d  = 0.03705 # Excentricidade da camada excêntrica

    # --- Print input parameters ---
    print("Parâmetros Geométricos:")
    print(f"  r4 = {r4} m, r5 = {r5} m, R1 = {R1} m, R2 = {R2} m")
    print(f"  Excentricidade d = {d} m\n")
    print("Permissividades Relativas:")
    print(f"  Camada 1: ε_r1 = {er1}")
    print(f"  Camada 2 (Ar): ε_r_ar = {er_ar}")
    print(f"  Camada 3: ε_r3 = {er3}\n")

    # --- Call the function and get the result ---
    cap_equivalent = calculate_hybrid_capacitance(r4, r5, R1, R2, d, er1, er_ar, er3)
    print("-" * 60)
    print(f"  Capacitância Equivalente = {cap_equivalent * 1e12:.3f} nF/km")

