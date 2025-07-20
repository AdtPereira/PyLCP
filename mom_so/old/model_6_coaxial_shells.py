"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

Below is a high-level overview of the script components:

Imports and Global Variables:

Required libraries and global variables are imported and defined.

BIFILAR_TL: A list containing properties of the two conductors.
Classes:

Geometry: Handles basic geometry calculations, such as distance matrices between conductor centers.
ParametersWithFrequency: Extends Geometry to include frequency-dependent parameters.
GreensMatrices: Uses the geometry to compute Green's matrices, which are essential for the method 
of moments.
MoMSuperficialOperator: Implements the method of moments, calculating matrices like U, Ys, G, and Z.
AnalyticalFormulation: Provides analytical formulations for high-frequency resistance, external 
inductance, and impedance.
Plotter: Handles plotting of series resistance and inductance against frequency.

Functions:

clear_screen: Clears the console screen.
main: The main function orchestrates the scattering calculations and plotting. It performs the 
following steps:
Clears the screen.
Initializes objects for the method of moments and analytical formulations.
Computes series resistance, external inductance, and impedance over a range of frequencies.
Plots the results using the Plotter class.
Detailed Class and Function Explanations
Geometry
__init__: Initializes the geometry of the system based on conductor properties.
distance_matrices: Calculates matrices for distances and angles between conductor centers.
ParametersWithFrequency
__init__: Extends the Geometry class to include frequency-dependent parameters such as 
conductivity, permeability, and permittivity.
ynp_operator: Calculates the surface admittance operator for a conductor.
GreensMatrices
__init__: Initializes Green's matrices using the geometry of the system.
dissertation and ieee_paper: Calculate Green's functions using different methods.
sub_matrices: Generates Green's sub-matrices.
MoMSuperficialOperator
__init__: Extends ParametersWithFrequency to initialize the method of moments parameters.
matrix_u: Constructs matrix U.
matrix_ys: Constructs matrix Ys.
matrix_g and matrix_g_ieee: Constructs matrix G using different methods.
matrix_z: Computes the impedance matrix Z.
AnalyticalFormulation
__init__: Initializes the analytical formulation based on the conductor properties and frequency.
pul_parameters: Calculates high-frequency resistance, external inductance, and impedance.
Plotter
__init__: Initializes the plotting class with frequency and impedance data.
series_resistance: Plots series resistance against frequency.
series_inductance: Plots series inductance against frequency.

REFERENCES:
[1] PATEL, Utkarsh R. A Surface Admittance Approach For Fast Calculation of the 
    Series Impedance of Cables Including Skin, Proximity, and Ground Return Effects.
    2014. University of Toronto, Graduate Department of The Edward S. Rogers Sr. 
    Department of Electrical & Computer Engineering. 

[2] U. R. Patel, B. Gustavsen and P. Triverio, "An Equivalent Surface Current Approach
    for the Computation of the Series Impedance of Power Cables with Inclusion of Skin
    and Proximity Effects," in IEEE Transactions on Power Delivery, vol. 28, no. 4, pp.
    2474-2482, Oct. 2013, doi: 10.1109/TPWRD.2013.2267098.

[3] U. R. Patel, B. Gustavsen and P. Triverio, "Application of the MoM-SO Method for 
    Accurate Impedance Calculation of Single-Core Cables Enclosed by a Conducting Pipe," 
    Proc. International Conference on Power Systems Transients (IPST 2013), Vancouver, 
    Canada July 18-20, 2013. https://www.ipstconf.org/Proc_IPST2013.php

"""

import time
import numpy as np
import mom_so_lossless_systems
from mom_so_system_data import TRANSMISSION_LINE
from mom_so_geometry import GraphicRepresentation
from mom_so_user_interface import UserInterface

# Multiconductor Transmission Line choices
# 0 - Two-Wire Transmission Line (TWT) with 100 mm spacing
# 1 - Two-Wire Transmission Line (TWT) with 25 mm spacing
# 2 - Coax Cable (Single Core Cable type)
# 3 - Coax Cable (Single Core Cable type) with asymmetrical core
# 4 - Coax Cable (Single Core Cable type) with screen wires (sw)
# 5 - Trefoil Power Cable
MTL_DICT = TRANSMISSION_LINE[6]


def main():
    """ Main function to perform the calculations and display the results."""

    # User configuration
    UserInterface().clear_screen()
    UserInterface().set_numpy_print_options()

    # Initialize the frequency points
    print("Calculations were started ... ...")
    start_time = time.time()
    frequency = np.logspace(0, 6, num=200)
    frequency_mom = np.logspace(0, 6, num=30)

    # Perform the calculations
    two_wire = mom_so_lossless_systems.TwoWire(
        MTL_DICT, frequency, frequency_mom, green_evaluation='Analytically')
    # two_wire = mom_so_results.TwoWire(
    #   MTL, frequency, frequency_mom, green_evaluation='Numerically')

    # End the timer
    elapsed_time = time.time() - start_time
    print(f"End of the routine! Time spent on simulation: {
          elapsed_time:.2f} seconds.\n")

    # Display the geometry of the transmission line
    GraphicRepresentation(MTL_DICT).wires_and_cables(line_type='bifilar')
    two_wire.plot_series_resistance()
    two_wire.plot_series_inductance()


if __name__ == "__main__":
    main()
