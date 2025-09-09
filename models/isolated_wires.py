# -*- coding: utf-8 -*-
"""
This script provides a parametric method to generate the MODEL dictionary 
for a flat arrangement of single-core cables, based on input parameters 
typically loaded from a JSON file.
This version dynamically handles the presence or absence of cable layers.
"""

import json
from utils.case_utils import UNITS_DATA

def _add_conductor_to_model(MODEL, input_json, conductor_id):
    """
    Helper function to add a single conductor layer (core, sheath, or armor) to the MODEL.
    It dynamically handles the presence of insulation.
    """
    arrangement = input_json.get('arrangement', None)
    conductor = input_json.get('conductor', {})
    insulation_data = conductor.get('insulation', None)

    idx_ref_conductor = int(arrangement.get('idx_ref_conductor', 0))
    scale_unit = UNITS_DATA[arrangement.get('unit', 'meter')]['scale']

    insulation_dict = None
    x_coordinate = conductor_id * arrangement.get('spacing', 0.0) / scale_unit 
    conductor_outer_radius = conductor.get('outer_radius', 0.0) / scale_unit
    fourier_order = arrangement.get('fourier_order', 4)

    if insulation_data:
        insulation_dict = {
            'name': f"primary_insulation",
            'center_point': (x_coordinate, 0.0),
            'thickness': insulation_data['thickness'] / scale_unit,
            'relative_permittivity': insulation_data['relative_permittivity'],
            'fourier_order': fourier_order,
        }

    MODEL[conductor_id] = {
        'line_id': conductor_id,
        'conductor_name': "conductor_" + str(conductor_id),
        'line_type': 'active' if conductor_id != idx_ref_conductor else 'return',
        'line_return': idx_ref_conductor if conductor_id != idx_ref_conductor else None,
        'center_point': (x_coordinate, 0.0),
        'radius': [0.0, conductor_outer_radius],
        'conductivity': conductor.get('conductivity_S_per_m', 5.8E7),
        'subconductors': None,
        'insulation': insulation_dict,
        'conductor_layers': None,
        'relative_permeability': 1.0,
        'relative_permittivity': 1.0,
        'relative_permittivity_out': 1.0,
        'potential_to_infinity': 1.0 if conductor_id != idx_ref_conductor else -1.0,
        'fourier_order': fourier_order,
    }

def circular_conductor_wires(
    input_json: dict, 
    show_model: bool = False
) -> dict:
    """
    Generates a parametric model for a single SCC cable, dynamically
    building it based on the layers defined in the JSON.
    """
    # Safely get cable definition components
    arrangement = input_json.get('arrangement', None)
    mom = input_json.get('mom', {})
    num_conductors = int(arrangement.get('num_conductors', 0))

    # Build the base MODEL dictionary
    MODEL = {
        'type': input_json.get('type', None),
        'name': input_json.get('name', 'unknown'),
        'note': input_json.get('note', 'unknown'),
        'idx_ref_conductor': mom.get('idx_ref_conductor', 0),
    }

    # Dynamically add conductors based on what's defined in the JSON
    for conductor_id in range(num_conductors):
        _add_conductor_to_model(MODEL, input_json, conductor_id)

    if show_model:
        print("\n--- Generated Dynamic SCC MODEL Dictionary ---")
        print(json.dumps(MODEL, indent=2, default=str))
        print("---------------------------------------------------\n")

    return MODEL
