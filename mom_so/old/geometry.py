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
"""

import numpy as np
from scipy.constants import epsilon_0
from mtl_data.mtl import MulticonductorTransmissionLine as MTL

# class UndergroundSystem(MTL):
#     """ This class contains the basic geometry of the system. """

#     def __init__(self, mtl):
#         super().__init__(mtl)

#         # Conductor surfaces dictionary
#         self.conductor_surfaces = []

#         # Hole surfaces dictionary
#         self.hole_surfaces = []

#         # Ground dictionary
#         self.ground = []

#         # Number of conductor and hole surfaces
#         self.surfaces_type = []

#         # Verify the items in the dictionary
#         for item in self.mtl:
#             # Check if the item is conductor
#             if item['line_type'] == 'active':
#                 self.surfaces_type.append(1)
#                 # Check if the conductor is hollow
#                 if item['radius'][0] != 0:
#                     for radius in item['radius']:
#                         self.conductor_surfaces.append({
#                             'center': item['center_point'],
#                             'fourier_order': item['surface_points'],
#                             'radius': radius
#                         })

#                 # Then, the conductor is solid
#                 else:
#                     self.conductor_surfaces.append({
#                         'center': item['center_point'],
#                         'fourier_order': item['surface_points'],
#                         'radius': item['radius'][1]
#                     })

#             # Check if the item is a hole
#             elif item['line_type'] == 'hole':
#                 self.surfaces_type.append(0)
#                 self.hole_surfaces.append({
#                     'center': item['center_point'],
#                     'fourier_order': item['surface_points'],
#                     'radius': item['radius'][1]
#                 })

#             # Check if the item is ground
#             elif item['line_type'] == 'return':
#                 self.ground.append({
#                     'conductivity': item['conductivity'],
#                     'permittivity': item['relative_permittivity'] * epsilon_0
#                 })

#         # Dimension N
#         # Equation (2.36) - PAG. 32 [1]
#         self.N = sum([2*Np['fourier_order']+1 for Np in self.conductor_surfaces])

#         # Dimension N_hat
#         # Equation (3.3) - PAG. 59 [1]
#         self.Nhat = sum([2*Nh['fourier_order']+1 for Nh in self.hole_surfaces])
