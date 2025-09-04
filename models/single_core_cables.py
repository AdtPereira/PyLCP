# -*- coding: utf-8 -*-
"""
This script provides a parametric method to generate the MODEL dictionary 
for a flat arrangement of single-core cables, based on input parameters 
typically loaded from a JSON file like 'scc_flat_xue.json'.
"""

import json

def three_phase_flat_model(
    input_json: dict, 
    show_model: bool = False
) -> dict:
    """
    Generates a parametric model dictionary for a flat arrangement of SCC cables.

    Args:
        input_json (dict): A dictionary containing the cable's physical definition 
                            (core, insulation, sheath, jacket); cable layout 
                            (burial depth, spacing); and soil properties.
        show_model (bool): If True, prints the generated MODEL dictionary.

    Returns:
        dict: The complete MODEL dictionary for the simulation.
    """
    # Cable Layers
    core = input_json['cable_definition']['core']
    core_insulation = input_json['cable_definition']['core']['insulation']
    sheath = input_json['cable_definition']['sheath']
    sheath_insulation = input_json['cable_definition']['sheath']['insulation']

    # Positions for the three cables in a flat arrangement
    depth = input_json['arrangement']['burial_depth_m']
    spacing = input_json['arrangement']['spacing_m']
    c1 = (0.0, -depth)
    c2 = (spacing, -depth)
    c3 = (2 * spacing, -depth)
    
    # --- 3. Build the MODEL dictionary ---
    MODEL = {
        'name': 'XUE_FLAT_ARRANGEMENT',
        'type': 'scc',
        'note': 'A parametric Xue flat arrangement cable model based on Fig. 4.18 (a) [Xue 2018]',
        'idx_ref_conductor': 0,
        0: {
            'line_id': 0,
            'conductor_name': 'soil',
            'line_type': 'return',
            'line_return': None,
            'conductivity': input_json['soil'].get('conductivity_S_per_m', 0.01),
            'relative_permeability': 1.0,
            'relative_permittivity': input_json['soil'].get('relative_permittivity', 1.0),
            'relative_permittivity_out': 1.0,
        },
    }

    # Loop through the cable positions to create the conductors
    conductor_id = 1
    for center_point in [c1, c2, c3]:
        # Core conductor
        MODEL[conductor_id] = {
            'line_id': conductor_id,
            'conductor_name': 'core',
            'line_type': 'active',
            'line_return': 0,
            'center_point': center_point,
            'radius': [core['inner_radius_m'], core['outer_radius_m']],
            'conductivity': core['conductivity_S_per_m'],
            'insulation': {
                'name': 'primary_insulation',
                'type': 'XLPE',
                'center_point': center_point,
                'thickness': core_insulation['thickness_m'],
                'relative_permittivity': core_insulation['relative_permittivity'],
                'fourier_order': 0,
                'relative_permeability': 1.0,
            },
            'subconductors': None,
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        }
        conductor_id += 1
        
        # Sheath conductor
        MODEL[conductor_id] = {
            'line_id': conductor_id,
            'conductor_name': 'sheath',
            'line_type': 'active',
            'line_return': 0,
            'center_point': center_point,
            'radius': [sheath['inner_radius_m'], sheath['outer_radius_m']],
            'conductivity': sheath['conductivity_S_per_m'],
            'insulation': {
                'name': 'secondary_insulation',
                'type': 'HDPE',
                'center_point': center_point,
                'thickness': sheath_insulation['thickness_m'],
                'relative_permittivity': sheath_insulation['relative_permittivity'],
                'fourier_order': 0,
                'relative_permeability': 1.0,
            },
            'subconductors': None,
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': 0,
        }
        conductor_id += 1

    if show_model:
        print("\n--- Generated SCC Flat Arrangement MODEL Dictionary ---")
        print(json.dumps(MODEL, indent=2, default=str))
        print("---------------------------------------------------\n")

    return MODEL
