# -*- coding: utf-8 -*-
"""
This script provides a parametric method to generate the MODEL dictionary 
for a flat arrangement of single-core cables, based on input parameters 
typically loaded from a JSON file.
This version dynamically handles the presence or absence of cable layers.
"""

import json

def _add_conductor_to_model(MODEL, conductor_id, center_point, layer_data, layer_name):
    """
    Helper function to add a single conductor layer (core, sheath, or armor) to the MODEL.
    It dynamically handles the presence of insulation.
    """
    insulation_data = layer_data.get('insulation')
    # Default to no insulation
    insulation_dict = None

    if insulation_data:
        insulation_dict = {
            'name': f"{layer_name}_insulation",
            'type': insulation_data.get('type', 'insulation'),
            'center_point': center_point,
            'thickness': insulation_data['thickness_m'],
            'relative_permittivity': insulation_data['relative_permittivity'],
            'relative_permeability': 1.0,
            'fourier_order': 0,
        }

    MODEL[conductor_id] = {
        'line_id': conductor_id,
        'conductor_name': layer_name,
        'line_type': 'active',
        'line_return': 0,
        'center_point': center_point,
        'radius': [layer_data['inner_radius_m'], layer_data['outer_radius_m']],
        'conductivity': layer_data['conductivity_S_per_m'],
        'subconductors': None,
        'insulation': insulation_dict,
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0,
        'fourier_order': 0,
    }
    return conductor_id + 1

def isolated_coaxial_cable(
    input_json: dict, 
    show_model: bool = False
) -> dict:
    """
    Generates a parametric model for a single SCC cable, dynamically
    building it based on the layers defined in the JSON.
    """
    # Safely get cable definition components
    cable_def = input_json.get('cable_definition', {})
    cable_ref = input_json.get('reference', {})
    frequency = input_json.get('frequency', {})
    core = cable_def.get('core')
    sheath = cable_def.get('sheath')
    armor = cable_def.get('armor')

    # Build the base MODEL dictionary
    MODEL = {
        'name': input_json.get('name', 'SINGLE_PHASE_SCC'),
        'type': input_json.get('type', 'coaxial'),
        'note': input_json.get('note', 'A parametric single-phase SCC model.'),
        'idx_ref_conductor': 0,
        'frequency': frequency,
        0: {
            'line_id': 0,
            'conductor_name': cable_ref.get('name', 'sheath'),
            'line_type': 'return',
            'line_return': None,
            'center_point': (0.0, 0.0),
            'radius': [cable_ref['inner_radius_m'], cable_ref['outer_radius_m']],
            'conductivity': cable_ref['conductivity_S_per_m'],
            'subconductors': None,
            'insulation': None,
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': -1.0,
            'fourier_order': 0,
        },
    }

    # Dynamically add conductors based on what's defined in the JSON
    conductor_id = 1
    for center_point in [(0.0, 0.0)]:
        if core:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, core, 'core')
        
        if sheath:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, sheath, 'sheath')
            
        if armor:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, armor, 'armor')

    if show_model:
        print("\n--- Generated Dynamic SCC MODEL Dictionary ---")
        print(json.dumps(MODEL, indent=2, default=str))
        print("---------------------------------------------------\n")

    return MODEL

def single_phase_model(
    input_json: dict, 
    show_model: bool = False
) -> dict:
    """
    Generates a parametric model for a single SCC cable, dynamically
    building it based on the layers defined in the JSON.
    """
    # Safely get cable definition components
    cable_def = input_json.get('cable_definition', {})
    core = cable_def.get('core')
    sheath = cable_def.get('sheath')
    armor = cable_def.get('armor')

    # Positions for the cable
    depth = input_json['arrangement']['burial_depth_m']
    c1 = (0.0, -depth)
    
    # Build the base MODEL dictionary
    MODEL = {
        'name': input_json.get('name', 'SINGLE_PHASE_SCC'),
        'type': 'scc',
        'note': input_json.get('note', 'A parametric single-phase SCC model.'),
        'idx_ref_conductor': 0,
        0: {
            'line_id': 0,
            'conductor_name': 'soil',
            'line_type': 'return',
            'line_return': None,
            'conductivity': input_json['soil']['conductivity_S_per_m'],
            'relative_permeability': 1.0,
            'relative_permittivity': input_json['soil']['relative_permittivity'],
            'relative_permittivity_out': 1.0,
        },
    }

    # Dynamically add conductors based on what's defined in the JSON
    conductor_id = 1
    for center_point in [c1]:
        if core:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, core, 'core')
        
        if sheath:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, sheath, 'sheath')
            
        if armor:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, armor, 'armor')

    if show_model:
        print("\n--- Generated Dynamic SCC MODEL Dictionary ---")
        print(json.dumps(MODEL, indent=2, default=str))
        print("---------------------------------------------------\n")

    return MODEL

def three_phase_flat_model(
    input_json: dict, 
    show_model: bool = False
) -> dict:
    """
    Generates a parametric model for a three-phase flat arrangement, dynamically
    building each cable based on the layers defined in the JSON.
    """
    # Safely get cable definition components
    cable_def = input_json.get('cable_definition', {})
    core = cable_def.get('core')
    sheath = cable_def.get('sheath')
    armor = cable_def.get('armor')

    # Positions for the three cables
    depth = input_json['arrangement']['burial_depth_m']
    spacing = input_json['arrangement']['spacing_m']
    c1 = (0.0, -depth)
    c2 = (spacing, -depth)
    c3 = (2 * spacing, -depth)
    
    # Build the base MODEL dictionary
    MODEL = {
        'name': input_json.get('name', 'THREE_PHASE_FLAT_SCC'),
        'type': 'scc',
        'note': input_json.get('note', 'A parametric three-phase flat SCC model.'),
        'idx_ref_conductor': 0,
        0: {
            'line_id': 0,
            'conductor_name': 'soil',
            'line_type': 'return',
            'line_return': None,
            'conductivity': input_json['soil']['conductivity_S_per_m'],
            'relative_permeability': 1.0,
            'relative_permittivity': input_json['soil']['relative_permittivity'],
            'relative_permittivity_out': 1.0,
        },
    }

    # Loop through cable positions and dynamically add conductors
    conductor_id = 1
    for center_point in [c1, c2, c3]:
        if core:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, core, 'core')
        
        if sheath:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, sheath, 'sheath')
            
        if armor:
            conductor_id = _add_conductor_to_model(MODEL, conductor_id, center_point, armor, 'armor')

    if show_model:
        print("\n--- Generated Dynamic SCC MODEL Dictionary ---")
        print(json.dumps(MODEL, indent=2, default=str))
        print("---------------------------------------------------\n")

    return MODEL