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
import warnings
import numpy as np
import mom_so_lossy_systems
from mom_so.mom_so_multiconductor import GraphicRepresentation as graph
from mom_so_user_interface import UserInterface
from mom_so_system_data import TRANSMISSION_LINE

# Multiconductor Transmission Line choices
# 7 - Underground Single core cable (SCC) with soil return
# 8 - Underground Double core cable (DCC) with soil return
MTL_DICT = TRANSMISSION_LINE[7]


def main():
    """ Main function to perform the calculations and display the results."""

    # User configuration
    UserInterface().clear_screen()
    UserInterface().set_numpy_print_options()

    # Warnings filter
    warnings.simplefilter("error")

    # Initialize the frequency points
    print("Calculations were started ... ...")
    start_time = time.time()
    frequency = np.logspace(0, 5, num=200)
    frequency_mom = [1E2] # np.logspace(0, 5, num=30)

    # Perform the MoM-SO calculations
    # green_type = ['Analytically', 'Numerically']
    mom_so_lossy_systems.BuriedSingleCoreCable(
        MTL_DICT, frequency, frequency_mom)

    # End the timer
    elapsed_time = time.time() - start_time
    print(f"\nEnd of the routine! Time spent on simulation: {
          elapsed_time:.2f} seconds.\n")

    # Display the geometry of the transmission line
    graph(MTL_DICT).underground_system()
    #scc.plot_series_resistance()


if __name__ == "__main__":
    main()
