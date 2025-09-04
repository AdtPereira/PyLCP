# -*- coding: utf-8 -*-
"""
    This script provides a parametric method to generate the dictionary for a 
    pipe-type cable in a trefoil configuration, resting at the bottom of the pipe.
    The parameters are based on the reference table from Patel, 2014.

    MODEL = {
        'name': 'patel_pipe_trefoil',
        'type': 'pipe',
        'note': 'An Patel pipe-type trefoil model Fig. 2.9 [Patel 2014]',
        'idx_ref_conductor': 0,
        0: {
            'line_id': 0,
            'conductor_name': 'pipe',
            'line_type': 'return',
            'line_return': None,
            'center_point': (0.0, 0.0144636),
            'radius': [0.045, 0.050],
            'conductivity': 1.0E7,
            'subconductors': None,
            'insulation': None,
            'conductor_layers': None,
            'relative_permeability': 100,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        },
        1: {
            'line_id': 1,
            'conductor_name': 'core',
            'line_type': 'active',
            'line_return': 0,
            'center_point': (0.0, 0.01871),
            'radius': [0.0, 0.010],
            'conductivity': 5.8E7,
            'subconductors': None,
            'insulation': {'name': 'primary_insulation',
                            'type': 'XLPE',
                            'center_point': (0.0, 0.01871),
                            'thickness': 0.004,
                            'fourier_order': 0,
                            'relative_permittivity': 2.3,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        },
        2: {
            'line_id': 2,
            'conductor_name': 'sheath',
            'line_type': 'active',
            'line_return': 0,
            'center_point': (0.0, 0.01871),
            'radius': [0.0140, 0.0142],
            'conductivity': 5.8E7,
            'subconductors': None,
            'insulation': {'name': 'secondary_insulation',
                            'type': 'HDPE',
                            'center_point': (0.0, 0.01871),
                            'thickness': 0.002,
                            'fourier_order': 0,
                            'relative_permittivity': 2.3,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        },
        3: {
            'line_id': 3,
            'conductor_name': 'core',
            'line_type': 'active',
            'line_return': 0,
            'center_point': (0.0162, -0.00935),
            'radius': [0.0, 0.010],
            'conductivity': 5.8E7,
            'subconductors': None,
            'insulation': {'name': 'primary_insulation',
                            'type': 'XLPE',
                            'center_point': (0.0162, -0.00935),
                            'thickness': 0.004,
                            'fourier_order': 0,
                            'relative_permittivity': 2.3,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        },
        4: {
            'line_id': 4,
            'conductor_name': 'sheath',
            'line_type': 'active',
            'line_return': 0,
            'center_point': (0.0162, -0.00935),
            'radius': [0.0140, 0.0142],
            'conductivity': 5.8E7,
            'subconductors': None,
            'insulation': {'name': 'secondary_insulation',
                            'type': 'HDPE',
                            'center_point': (0.0162, -0.00935),
                            'thickness': 0.002,
                            'fourier_order': 0,
                            'relative_permittivity': 2.3,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        },
        5: {
            'line_id': 5,
            'conductor_name': 'core',
            'line_type': 'active',
            'line_return': 0,
            'center_point': (-0.0162, -0.00935),
            'radius': [0.0, 0.010],
            'conductivity': 5.8E7,
            'subconductors': None,
            'insulation': {'name': 'primary_insulation',
                            'type': 'XLPE',
                            'center_point': (-0.0162, -0.00935),
                            'thickness': 0.004,
                            'fourier_order': 0,
                            'relative_permittivity': 2.3,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        },
        6: {
            'line_id': 6,
            'conductor_name': 'sheath',
            'line_type': 'active',
            'line_return': 0,
            'center_point': (-0.0162, -0.00935),
            'radius': [0.0140, 0.0142],
            'conductivity': 5.8E7,
            'subconductors': None,
            'insulation': {'name': 'secondary_insulation',
                            'type': 'HDPE',
                            'center_point': (-0.0162, -0.00935),
                            'thickness': 0.002,
                            'fourier_order': 0,
                            'relative_permittivity': 2.3,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        },
    }
"""

import math
import json

def trefoil_symmetric_model(
    core_radius_mm: float,
    core_conductivity: float,
    insulation_thickness_mm: float,
    insulation_rel_permittivity: float,
    screen_thickness_mm: float,
    screen_conductivity: float,
    jacket_thickness_mm: float,
    jacket_rel_permittivity: float,
    pipe_outer_diameter_mm: float,
    pipe_thickness_mm: float,
    pipe_conductivity: float,
    pipe_rel_permeability: float,
    show_model: bool = False
) -> dict:
    """
    Generates a parametric model dictionary for a pipe-type trefoil cable system
    in a SYMMETRIC configuration where the trefoil and pipe are centered at (0,0).

    Args:
        (Arguments are the same as the asymmetric version)
        ...
        show_model (bool): If True, prints the generated MODEL dictionary
                           to the console. Defaults to False.

    Returns:
        dict: A dictionary containing the complete model definition.
    """
    # --- 1. Convert all units from mm to meters ---
    core_radius = core_radius_mm / 1000.0
    insulation_thickness = insulation_thickness_mm / 1000.0
    screen_thickness = screen_thickness_mm / 1000.0
    jacket_thickness = jacket_thickness_mm / 1000.0
    pipe_outer_radius = (pipe_outer_diameter_mm / 2.0) / 1000.0
    pipe_thickness = pipe_thickness_mm / 1000.0

    # --- 2. Calculate radii of different layers ---
    insulation_outer_radius = core_radius + insulation_thickness
    screen_inner_radius = insulation_outer_radius
    screen_outer_radius = screen_inner_radius + screen_thickness
    cable_total_radius = screen_outer_radius + jacket_thickness
    pipe_inner_radius = pipe_outer_radius - pipe_thickness
    
    # --- 3. Calculate trefoil center coordinates (relative to trefoil centroid at 0,0) ---
    # This logic remains the same as the cables are positioned relative to their own center
    triangle_side_length = 2 * cable_total_radius
    dist_centroid_to_vertex = triangle_side_length / math.sqrt(3)

    center_1 = (0.0, dist_centroid_to_vertex)
    center_2 = (triangle_side_length / 2.0, -dist_centroid_to_vertex / 2.0)
    center_3 = (-triangle_side_length / 2.0, -dist_centroid_to_vertex / 2.0)
    
    cable_centers = [center_1, center_2, center_3]

    # --- 4. Define the pipe center coordinate ---
    # For the symmetric case, the pipe center is at the origin.
    pipe_center = (0.0, 0.0)

    # --- 5. Build the MODEL dictionary ---
    MODEL = {
        'name': 'patel_pipe_trefoil_parametric_symmetric',
        'type': 'pipe',
        'note': 'An Patel pipe-type trefoil model - Symmetric Configuration',
        'idx_ref_conductor': 0,
    }

    # Add Pipe (Conductor 0)
    MODEL[0] = {
        'line_id': 0,
        'conductor_name': 'pipe',
        'line_type': 'return',
        'line_return': None,
        'center_point': pipe_center,
        'radius': [pipe_inner_radius, pipe_outer_radius],
        'conductivity': pipe_conductivity,
        'subconductors': None,
        'insulation': None,
        'conductor_layers': None,
        'relative_permeability': pipe_rel_permeability,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    }

    # Add Cores and Sheaths for the three cables
    conductor_id = 1
    for center_point in cable_centers:
        # Core
        MODEL[conductor_id] = {
            'line_id': conductor_id,
            'conductor_name': 'core',
            'line_type': 'active',
            'line_return': 0,
            'center_point': center_point,
            'radius': [0.0, core_radius],
            'conductivity': core_conductivity,
            'subconductors': None,
            'insulation': {'name': 'primary_insulation',
                            'type': 'XLPE',
                            'center_point': center_point,
                            'thickness': insulation_thickness,
                            'fourier_order': 0,
                            'relative_permittivity': insulation_rel_permittivity,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        }
        conductor_id += 1
        
        # Sheath/Screen
        MODEL[conductor_id] = {
            'line_id': conductor_id,
            'conductor_name': 'sheath',
            'line_type': 'active',
            'line_return': 0,
            'center_point': center_point,
            'radius': [screen_inner_radius, screen_outer_radius],
            'conductivity': screen_conductivity,
            'subconductors': None,
            'insulation': {'name': 'secondary_insulation',
                            'type': 'HDPE', # High-density polyethylene (Jacket material)
                            'center_point': center_point,
                            'thickness': jacket_thickness,
                            'fourier_order': 0,
                            'relative_permittivity': jacket_rel_permittivity,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        }
        conductor_id += 1

    if show_model:
        print("\n--- Generated SYMMETRIC MODEL Dictionary ---")
        print(json.dumps(MODEL, indent=2, default=str))
        print("------------------------------------------\n")

    return MODEL

def trefoil_asymmetric_model(
    core_radius_mm: float,
    core_conductivity: float,
    insulation_thickness_mm: float,
    insulation_rel_permittivity: float,
    screen_thickness_mm: float,
    screen_conductivity: float,
    jacket_thickness_mm: float,
    jacket_rel_permittivity: float,
    pipe_outer_diameter_mm: float,
    pipe_thickness_mm: float,
    pipe_conductivity: float,
    pipe_rel_permeability: float,
    show_model: bool = False
) -> dict:
    """
    Generates a parametric model dictionary for a pipe-type trefoil cable system.

    The function calculates all geometric properties based on the input parameters
    from the reference table, including the asymmetric placement of the trefoil
    at the bottom of the pipe. All dimensions are converted to meters.

    Args:
        core_radius_mm (float): Radius of the central conductor [mm].
        core_conductivity (float): Electrical conductivity of the core [S/m].
        insulation_thickness_mm (float): Thickness of the primary insulation [mm].
        insulation_rel_permittivity (float): Relative permittivity of the primary insulation.
        screen_thickness_mm (float): Thickness of the screen/sheath [mm].
        screen_conductivity (float): Electrical conductivity of the screen/sheath [S/m].
        jacket_thickness_mm (float): Thickness of the outer jacket [mm].
        jacket_rel_permittivity (float): Relative permittivity of the jacket.
        pipe_outer_diameter_mm (float): Outer diameter of the steel pipe [mm].
        pipe_thickness_mm (float): Thickness of the steel pipe wall [mm].
        pipe_conductivity (float): Electrical conductivity of the pipe [S/m].
        pipe_rel_permeability (float): Relative permeability of the pipe.
        show_model (bool): If True, prints the generated MODEL dictionary
                           to the console. Defaults to False.

    Returns:
        dict: A dictionary containing the complete model definition.
    """
    # --- 1. Convert all units from mm to meters ---
    core_radius = core_radius_mm / 1000.0
    insulation_thickness = insulation_thickness_mm / 1000.0
    screen_thickness = screen_thickness_mm / 1000.0
    jacket_thickness = jacket_thickness_mm / 1000.0
    pipe_outer_radius = (pipe_outer_diameter_mm / 2.0) / 1000.0
    pipe_thickness = pipe_thickness_mm / 1000.0

    # --- 2. Calculate radii of different layers ---
    insulation_outer_radius = core_radius + insulation_thickness
    screen_inner_radius = insulation_outer_radius
    screen_outer_radius = screen_inner_radius + screen_thickness
    cable_total_radius = screen_outer_radius + jacket_thickness
    pipe_inner_radius = pipe_outer_radius - pipe_thickness
    
    # --- 3. Calculate trefoil center coordinates (relative to trefoil centroid at 0,0) ---
    triangle_side_length = 2 * cable_total_radius
    dist_centroid_to_vertex = triangle_side_length / math.sqrt(3)

    center_1 = (0.0, dist_centroid_to_vertex)
    center_2 = (triangle_side_length / 2.0, -dist_centroid_to_vertex / 2.0)
    center_3 = (-triangle_side_length / 2.0, -dist_centroid_to_vertex / 2.0)
    
    cable_centers = [center_1, center_2, center_3]

    # --- 4. Calculate the pipe center coordinate (y_p) ---
    # The distance between the pipe center and a bottom cable center is (R_pipe - r_cable)
    dist_centers = pipe_inner_radius - cable_total_radius
    
    # Using the distance formula: D^2 = (x2-x1)^2 + (y2-y1)^2
    # D = dist_centers, (x1,y1) = (0, y_p), (x2,y2) = center_2
    # dist_centers^2 = (center_2[0])^2 + (center_2[1] - y_p)^2
    # Solving for y_p:
    y_p_squared_term = dist_centers**2 - center_2[0]**2
    # We take the positive root for the resting configuration
    y_p_term = math.sqrt(y_p_squared_term) 
    # y_p = center_2[1] + y_p_term --> This corresponds to hanging from the top
    # y_p = center_2[1] - (-y_p_term) --> This corresponds to resting on the bottom
    # After re-arranging the equation: -center_2[1] - y_p = -y_p_term
    pipe_center_y = center_2[1] + y_p_term

    pipe_center = (0.0, pipe_center_y)

    # --- 5. Build the MODEL dictionary ---
    MODEL = {
        'name': 'patel_pipe_trefoil_parametric',
        'type': 'pipe',
        'note': 'An Patel pipe-type trefoil model Fig. 2.9 [Patel 2014] - Parametric',
        'idx_ref_conductor': 0,
    }

    # Add Pipe (Conductor 0)
    MODEL[0] = {
        'line_id': 0,
        'conductor_name': 'pipe',
        'line_type': 'return',
        'line_return': None,
        'center_point': pipe_center,
        'radius': [pipe_inner_radius, pipe_outer_radius],
        'conductivity': pipe_conductivity,
        'subconductors': None,
        'insulation': None,
        'conductor_layers': None,
        'relative_permeability': pipe_rel_permeability,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    }

    # Add Cores and Sheaths for the three cables
    conductor_id = 1
    for center_point in cable_centers:
        # Core
        MODEL[conductor_id] = {
            'line_id': conductor_id,
            'conductor_name': 'core',
            'line_type': 'active',
            'line_return': 0,
            'center_point': center_point,
            'radius': [0.0, core_radius],
            'conductivity': core_conductivity,
            'subconductors': None,
            'insulation': {'name': 'primary_insulation',
                            'type': 'XLPE',
                            'center_point': center_point,
                            'thickness': insulation_thickness,
                            'fourier_order': 0,
                            'relative_permittivity': insulation_rel_permittivity,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        }
        conductor_id += 1
        
        # Sheath/Screen
        MODEL[conductor_id] = {
            'line_id': conductor_id,
            'conductor_name': 'sheath',
            'line_type': 'active',
            'line_return': 0,
            'center_point': center_point,
            'radius': [screen_inner_radius, screen_outer_radius],
            'conductivity': screen_conductivity,
            'subconductors': None,
            'insulation': {'name': 'secondary_insulation',
                            'type': 'HDPE', # High-density polyethylene (Jacket material)
                            'center_point': center_point,
                            'thickness': jacket_thickness,
                            'fourier_order': 0,
                            'relative_permittivity': jacket_rel_permittivity,
                            'relative_permeability': 1.0,
                        },
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        }
        conductor_id += 1

    if show_model:
        print("\n--- Generated MODEL Dictionary ---")
        print(json.dumps(MODEL, indent=2, default=str))
        print("--------------------------------\n")

    return MODEL
