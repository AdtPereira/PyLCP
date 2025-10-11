# -*- coding: utf-8 -*-
"""
This script provides an object-oriented, parametric method to generate 
the MODEL dictionary for various arrangements of single-core cables.

The CableModelGenerator class takes a JSON-like dictionary as input and 
can generate models for single-phase, three-phase flat, and isolated 
coaxial arrangements, dynamically handling the presence or absence of 
cable layers.
"""

import json
from pathlib import Path
import numpy as np
from typing import Dict, Any, Tuple
from utils.case_utils import UNITS_DATA

class SingleCoreCableModelGenerator:
    """
    A class to generate parametric models for single-core cable arrangements.
    """

    def __init__(self, file_path: str, silent_mode: bool = False):
        """
        Initializes the generator with cable and environmental definitions.

        Args:
            input_json (Dict[str, Any]): A dictionary containing the definitions
                                         for the cable, soil, and arrangement.
        """
        self.script_path = Path(file_path)
        self.silent_mode = silent_mode

        self.input_data = self.load_json_parameters()
        
        self.cable_def = self.input_data.get('cable_definition', {})
        self.soil = self.input_data.get('soil', {})
        self.arrangement = self.input_data.get('arrangement', {})
        self.reference = self.input_data.get('reference', {})
        self.scale_unit = UNITS_DATA[self.arrangement.get('unit', 'meter')]['scale']
        self.fourier_order = self.arrangement.get('fourier_order', 0)

        # Cable layer definitions
        self.core = self.cable_def.get('core')
        self.sheath = self.cable_def.get('sheath')
        self.armor = self.cable_def.get('armor')

    def load_json_parameters(self) -> dict:
        """
        Dynamically loads parameters from a '.json' file.
        It assumes the '.in.json' file has the same base name as the
        calling script and is located in the same directory.

        Args:
            script_file_path (str): The __file__ attribute from the calling script.
            show_content (bool): If True, prints the content of the loaded
                                dictionary to the console. Defaults to False.

        Returns:
            dict: A dictionary with the parameters loaded from the JSON file.
        """
        
        model_file_name = f"{self.script_path.stem}.json"
        model_path = self.script_path.resolve().parent / model_file_name

        print(f"Case name identified as: '{self.script_path.stem}'")
        
        if not model_path.exists():
            raise FileNotFoundError(f"The parameter file could not be found at: {model_path}")

        with open(model_path, 'r') as f:
            parameters = json.load(f)
        
        print(f"Successfully loaded parameters from: {model_path}")
        
        if not self.silent_mode:
            print(f"\n--- Content of {model_file_name} ---")
            print(json.dumps(parameters, indent=2))
            print("--------------------------------------\n")
            
        return parameters

    def _add_single_layer(self, 
                          model: Dict[str, Any], 
                          conductor_id: int, 
                          center_point: Tuple[float, float], 
                          layer_data: Dict[str, Any], 
                          layer_name: str) -> int:
        """
        Adds a single conductor layer (e.g., core, sheath) to the model.
        """
        insulation_data = layer_data.get('insulation')
        insulation_dict = None

        if insulation_data:
            insulation_dict = {
                'name': f"{layer_name}_insulation",
                'type': insulation_data.get('type', 'insulation'),
                'center_point': center_point,
                'thickness': insulation_data['thickness'],
                'relative_permittivity': insulation_data['relative_permittivity'],
                'relative_permeability': 1.0,
                'fourier_order': self.fourier_order,
            }

        model[conductor_id] = {
            'line_id': conductor_id,
            'conductor_name': layer_name,
            'line_type': 'active',
            'line_return': 0,
            'center_point': center_point,
            'radius': [layer_data['inner_radius'], layer_data['outer_radius']],
            'conductivity': layer_data['conductivity_S_per_m'],
            'subconductors': None,
            'insulation': insulation_dict,
            'conductor_layers': None,
            'relative_permeability': 1.0,
            'relative_permittivity': 1.0,
            'relative_permittivity_out': 1.0,
            'potential_to_infinity': 1.0,
            'fourier_order': self.fourier_order,
        }
        return conductor_id + 1

    def _add_cable_conductors(self, 
                              model: Dict[str, Any], 
                              conductor_id: int, 
                              center_point: Tuple[float, float]) -> int:
        """
        Adds all defined conductive layers (core, sheath, armor) for a single cable.
        """
        if self.core:
            conductor_id = self._add_single_layer(model, conductor_id, center_point, self.core, 'core')
        if self.sheath:
            conductor_id = self._add_single_layer(model, conductor_id, center_point, self.sheath, 'sheath')
        if self.armor:
            conductor_id = self._add_single_layer(model, conductor_id, center_point, self.armor, 'armor')
        return conductor_id
    
    @staticmethod
    def _show_model(model: Dict[str, Any]):
        """Prints the generated model dictionary in a readable format."""
        print(f"\n--- Generated Model: {model.get('name', 'N/A')} ---")
        print(json.dumps(model, indent=2, default=str))
        print("---------------------------------------------------\n")

    def concentric_hdpe_enclosed_model(self, host_conductor: str = 'sheath') -> Dict[str, Any]:
        """
        Generates a model for a cable inside an HDPE enclosure.
        The enclosure is defined within a conductor in the JSON and is added
        as a property to that conductor in the final model dictionary.
        The geometry remains concentric.
        """
        
        # --- Data Retrieval ---
        # Find the conductor that defines the enclosure (typically the outermost one).
        host_conductor_data = getattr(self, host_conductor)
        host_insulation = host_conductor_data.get('insulation')
        enclosure_data = host_conductor_data.get('enclosure')

        if not enclosure_data:
            raise ValueError(f"Enclosure definition not found within conductor '{host_conductor}'.")

        # --- Concentricity Calculation ---
        # The cable's center is the reference point, defined by the burial depth.
        # The enclosure's center is the same as the cable's center for concentric arrangement.
        cable_center = (0.0, -self.arrangement['burial_depth'])        

        # --- Model Generation ---
        model = {
            'name': self.input_data.get('name', 'HDPE_Enclosed_SCC_System'),
            'type': 'hdpe',
            'note': 'A parametric model of a cable eccentrically placed inside an HDPE enclosure.',
            'idx_ref_conductor': 0,
            0: {
                'line_id': 0,
                'conductor_name': 'soil',
                'line_type': 'return',
                'line_return': None,
                'conductivity': self.soil['conductivity_S_per_m'],
                'relative_permeability': 1.0,
                'relative_permittivity': self.soil['relative_permittivity'],
                'relative_permittivity_out': 1.0,
            },
        }

        conductor_id = 1
        # Add the cable conductors (core, sheath) at their reference position
        conductor_id = self._add_cable_conductors(model, conductor_id, cable_center)

        # --- Inject Enclosure Data into Host Conductor ---
        # Find the host conductor's entry in the generated model.
        for k, v in model.items():
            if isinstance(k, int) and k > 0: 
                if v.get('conductor_name') == host_conductor:
                    # Add the enclosure dictionary, including its calculated center point.
                    v['enclosure'] = enclosure_data
                    v['enclosure']['center_point'] = cable_center
                else:
                    v['enclosure'] = None

        if not self.silent_mode:
            self._show_model(model)

        return model
    
    def eccentric_hdpe_enclosed_model(self, host_conductor: str = 'sheath') -> Dict[str, Any]:
        """
        Generates a model for a cable inside an HDPE enclosure.
        The enclosure is defined within a conductor in the JSON and is added
        as a property to that conductor in the final model dictionary.
        The geometry remains eccentric.
        """
        
        # --- Data Retrieval ---
        # Find the conductor that defines the enclosure (typically the outermost one).
        host_conductor_data = getattr(self, host_conductor)
        host_insulation = host_conductor_data.get('insulation')
        enclosure_data = host_conductor_data.get('enclosure')

        if not enclosure_data:
            raise ValueError(f"Enclosure definition not found within conductor '{host_conductor}'.")

        # --- Eccentricity Calculation ---
        # 1. The cable's center is the reference point, defined by the burial depth.
        cable_center = (0.0, -self.arrangement['burial_depth'])
        
        # 2. The enclosure's center is calculated based on the cable's position.        
        cable_outer_radius = host_conductor_data['outer_radius'] + (host_insulation['thickness'] if host_insulation else 0)        
        enclosure_inner_radius = enclosure_data['inner_radius']
        
        if cable_outer_radius > enclosure_inner_radius:
            raise ValueError("Cable does not fit inside the enclosure based on JSON dimensions.")

        vertical_offset = enclosure_inner_radius - cable_outer_radius
        
        # 3. The enclosure's center is displaced upwards by the offset.
        enclosure_center = (cable_center[0], cable_center[1] + vertical_offset)

        # --- Model Generation ---
        model = {
            'name': self.input_data.get('name', 'HDPE_Enclosed_SCC_System'),
            'type': 'hdpe',
            'note': 'A parametric model of a cable eccentrically placed inside an HDPE enclosure.',
            'idx_ref_conductor': 0,
            0: {
                'line_id': 0,
                'conductor_name': 'soil',
                'line_type': 'return',
                'line_return': None,
                'conductivity': self.soil['conductivity_S_per_m'],
                'relative_permeability': 1.0,
                'relative_permittivity': self.soil['relative_permittivity'],
                'relative_permittivity_out': 1.0,
            },
        }

        conductor_id = 1
        # Add the cable conductors (core, sheath) at their reference position
        conductor_id = self._add_cable_conductors(model, conductor_id, cable_center)

        # --- Inject Enclosure Data into Host Conductor ---
        # Find the host conductor's entry in the generated model.
        for k, v in model.items():
            if isinstance(k, int) and k > 0: 
                if v.get('conductor_name') == host_conductor:
                    # Add the enclosure dictionary, including its calculated center point.
                    v['enclosure'] = enclosure_data
                    v['enclosure']['center_point'] = enclosure_center
                else:
                    v['enclosure'] = None

        if not self.silent_mode:
            self._show_model(model)

        return model
    
    def underground_model(self) -> Dict[str, Any]:
        """
        Generates a parametric model for underground cables in a flat arrangement.
        The number of cables and their spacing is determined by the 'arrangement'
        data provided during initialization.
        """
        depth = self.arrangement['burial_depth']
        num_conductors = self.arrangement.get('num_conductors', 1)
        
        # Default spacing to 0 if not specified (for the single conductor case)
        spacing = self.arrangement.get('spacing', 0) if num_conductors > 1 else 0

        # Dynamically generate center points for any number of conductors
        center_points = [(i * spacing, -depth) for i in range(num_conductors)]
        
        model = {
            'name': self.input_data.get('name', 'Underground_SCC_System'),
            'type': 'scc',
            'note': self.input_data.get('note', 'A parametric underground SCC model.'),
            'idx_ref_conductor': 0,
            0: {
                'line_id': 0,
                'conductor_name': 'soil',
                'line_type': 'return',
                'line_return': None,
                'conductivity': self.soil['conductivity_S_per_m'],
                'relative_permeability': 1.0,
                'relative_permittivity': self.soil['relative_permittivity'],
                'relative_permittivity_out': 1.0,
            },
        }

        conductor_id = 1
        for cp in center_points:
            conductor_id = self._add_cable_conductors(model, conductor_id, cp)

        if not self.silent_mode:
            self._show_model(model)

        return model
    
    def conventional_single_phase(self) -> Dict[str, Any]:
        """
        Generates a parametric model for a single buried SCC cable.
        """
        depth = self.arrangement['burial_depth']
        center_points = [(0.0, -depth)]
        
        model = {
            'name': self.input_data.get('name', 'SINGLE_PHASE_SCC'),
            'type': 'scc',
            'note': self.input_data.get('note', 'A parametric single-phase SCC model.'),
            'idx_ref_conductor': 0,
            0: {
                'line_id': 0,
                'conductor_name': 'soil',
                'line_type': 'return',
                'line_return': None,
                'conductivity': self.soil['conductivity_S_per_m'],
                'relative_permeability': 1.0,
                'relative_permittivity': self.soil['relative_permittivity'],
                'relative_permittivity_out': 1.0,
            },
        }

        conductor_id = 1
        for cp in center_points:
            conductor_id = self._add_cable_conductors(model, conductor_id, cp)

        if not self.silent_mode:
            self._show_model(model)

        return model

    def conventional_three_phase_flat(self) -> Dict[str, Any]:
        """
        Generates a model for a three-phase flat arrangement of SCC cables.
        """
        depth = self.arrangement['burial_depth']
        spacing = self.arrangement['spacing']
        center_points = [(0.0, -depth), (spacing, -depth), (2 * spacing, -depth)]
        
        model = {
            'name': self.input_data.get('name', 'THREE_PHASE_FLAT_SCC'),
            'type': 'scc',
            'note': self.input_data.get('note', 'A parametric three-phase flat SCC model.'),
            'idx_ref_conductor': 0,
            0: {
                'line_id': 0,
                'conductor_name': 'soil',
                'line_type': 'return',
                'line_return': None,
                'conductivity': self.soil['conductivity_S_per_m'],
                'relative_permeability': 1.0,
                'relative_permittivity': self.soil['relative_permittivity'],
                'relative_permittivity_out': 1.0,
            },
        }

        conductor_id = 1
        for cp in center_points:
            conductor_id = self._add_cable_conductors(model, conductor_id, cp)

        if not self.silent_mode:
            self._show_model(model)

        return model

    def isolated_coaxial_cable(self) -> Dict[str, Any]:
        """
        Generates a model for a single isolated coaxial cable.
        """
        model = {
            'name': self.input_data.get('name', 'ISOLATED_COAXIAL_SCC'),
            'type': 'coaxial',
            'note': self.input_data.get('note', 'A parametric isolated coaxial SCC model.'),
            'idx_ref_conductor': 0,
            0: {
                'line_id': 0,
                'conductor_name': self.reference.get('name', 'sheath'),
                'line_type': 'return',
                'line_return': None,
                'center_point': (0.0, 0.0),
                'radius': [self.reference['inner_radius'], self.reference['outer_radius']],
                'conductivity': self.reference['conductivity_S_per_m'],
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

        # For a coaxial cable, there is only one center point at the origin
        self._add_cable_conductors(model, 1, (0.0, 0.0))

        if not self.silent_mode:
            self._show_model(model)

        return model

    def simple_trefoil(self) -> Dict[str, Any]:
        """
        Gera um modelo para um arranjo trifólio (trefoil) de três fases de cabos SCC.
        Assume-se que 'burial_depth' corresponde à profundidade dos centros dos
        dois condutores da base (h3 na figura de referência). O 'spacing' é a
        distância de centro a centro entre cabos adjacentes.
        """
        # Profundidade dos cabos da base (h3)
        depth_bottom = self.arrangement['burial_depth']
        spacing = self.arrangement['spacing']

        # --- Cálculo da Geometria Trifólio ---
        # A altura do triângulo equilátero formado pelos cabos.
        height = spacing * np.sqrt(3) / 2
        
        # As coordenadas dos cabos da base (b e c) são conhecidas.
        y_bottom = -depth_bottom
        x_side = spacing / 2

        # A coordenada do cabo superior (a) é calculada a partir da base.
        # Sua posição vertical é a da base mais a altura do triângulo.
        y_top = y_bottom + height
        
        center_points = [
            (0.0, y_top),         # Cabo superior (a)
            (-x_side, y_bottom),  # Cabo inferior esquerdo (b)
            (+x_side, y_bottom)   # Cabo inferior direito (c)
        ]
        
        model = {
            'name': self.input_data.get('name', 'THREE_PHASE_TREFOIL_SCC'),
            'type': 'scc',
            'note': self.input_data.get('note', 'A parametric three-phase trefoil SCC model based on bottom conductor depth.'),
            'idx_ref_conductor': self.arrangement.get('idx_ref_conductor', 0),
            0: {
                'line_id': 0,
                'conductor_name': 'soil',
                'line_type': 'return',
                'line_return': None,
                'conductivity': self.soil['conductivity_S_per_m'],
                'relative_permeability': 1.0,
                'relative_permittivity': self.soil['relative_permittivity'],
                'relative_permittivity_out': 1.0,
            },
        }

        conductor_id = 1
        for cp in center_points:
            conductor_id = self._add_cable_conductors(model, conductor_id, cp)

        if not self.silent_mode:
            self._show_model(model)

        return model