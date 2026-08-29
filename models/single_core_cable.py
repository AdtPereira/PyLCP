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
import copy
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
        self.driver_frequency = self.input_data.get('frequency_driver', {})
        
        self.scale_unit = UNITS_DATA[self.arrangement.get('unit', 'meter')]['scale']
        self.fourier_order = self.arrangement.get('fourier_order', 0)

        # Cable layer definitions
        self.core = self.cable_def.get('core', None)
        self.sheath = self.cable_def.get('sheath', None)
        self.armor = self.cable_def.get('armor', None)
        self.ecc = self.cable_def.get('ecc', None)

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
        Adds a single conductor layer (e.g., core, sheath or armor) to the model.
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
    
    def _add_ecc_conductor(self,
                           model: Dict[str, Any],
                           conductor_id: int,
                           center_point: Tuple[float, float]) -> int:
        """
        Adds the ECC (Earth Continuity Conductor) to the model
        at its specific 'center_point'.
        """
        if self.ecc:
            conductor_id = self._add_single_layer(model, conductor_id, center_point, self.ecc, 'ecc')
        return conductor_id
    
    def _calculate_ecc_center_trig(self,
                                   enclosure_center: Tuple[float, float],
                                   R_enc: float,
                                   R_scc: float,
                                   R_ecc: float,
                                   vertical_offset_c: float) -> Tuple[float, float]:
        """
        Computes the position of the ECC "wedged" between the SCC and the HDPE duct
        using the Law of Cosines.

        Args:
            enclosure_center (P_enc): Center point of the duct.
            R_enc: Inner radius of the duct.
            R_scc: Outer radius of the SCC cable.
            R_ecc: Outer radius of the ECC cable.
            vertical_offset_c (c): Distance between P_enc and P_scc.
        """

        # Geometry validation
        if R_enc <= (R_scc + R_ecc):
            raise ValueError(
                f"Impossible geometry: the ECC (R={R_ecc}) does not fit in the space "
                f"between the SCC (R={R_scc}) and the HDPE duct (R={R_enc}). "
                f"Condition (R_enc > R_scc + R_ecc) failed."
            )

        # Sides of the triangle formed by the centers (P_enc, P_scc, P_ecc)
        a = R_scc + R_ecc  # Dist (P_scc -> P_ecc)
        b = R_enc - R_ecc  # Dist (P_enc -> P_ecc)
        c = vertical_offset_c  # Dist (P_enc -> P_scc)

        # Law of Cosines to find the angle 'alpha' at vertex P_enc
        # a^2 = b^2 + c^2 - 2bc*cos(alpha)
        numerator = (b**2) + (c**2) - (a**2)
        denominator = 2 * b * c

        # Floating-point error handling
        cos_alpha = max(min(numerator / denominator, 1.0), -1.0)

        alpha = np.arccos(cos_alpha)
        sin_alpha = np.sin(alpha)

        # Compute the ECC center
        # (Assuming the ECC is in the x+, y- quadrant)
        # The ECC center is at a distance 'b' from 'enclosure_center',
        # rotated by the angle 'alpha' from the vertical line P_enc -> P_scc.

        ecc_center_x = enclosure_center[0] + b * sin_alpha
        ecc_center_y = enclosure_center[1] - b * cos_alpha # Subtract because the axis points downward

        return (ecc_center_x, ecc_center_y)
    
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

    def flat_hdpe_enclosed_model(self, host_conductor: str = 'sheath') -> Dict[str, Any]:
        """
        Generates a model for Configuration 2 (Figure 5.2): three power SCC
        cables in a flat arrangement, each installed inside its own HDPE duct
        buried in the soil. The geometry is eccentric: the cable rests on the
        bottom of its duct (same convention as
        `eccentric_hdpe_enclosed_model`, see also `schematic_original_hdpe.png`
        of the `hdpe_300mm2` case), not concentric. Since 'burial_depth' refers
        to the center of the cable (not the center of the duct), the duct center
        is displaced upward relative to the cable center.

        Generalizes `eccentric_hdpe_enclosed_model` (1 cable) to N cables, the
        same way `underground_flat_model` generalizes
        `conventional_single_phase`: same 'burial_depth'/'spacing' from the
        'arrangement', an identical duct (same vertical offset) injected into
        each cable.
        """
        # --- Cable and duct data (defined in the host conductor, e.g. 'sheath') ---
        host_conductor_data = getattr(self, host_conductor)
        host_insulation = host_conductor_data.get('insulation')
        enclosure_data = host_conductor_data.get('enclosure')

        if not enclosure_data:
            raise ValueError(f"Enclosure definition not found within conductor '{host_conductor}'.")

        # --- Vertical cable -> duct offset (same eccentric logic as
        # `eccentric_hdpe_enclosed_model`, applied identically to each phase) ---
        cable_outer_radius = host_conductor_data['outer_radius'] + (host_insulation['thickness'] if host_insulation else 0)
        enclosure_inner_radius = enclosure_data['inner_radius']

        if cable_outer_radius > enclosure_inner_radius:
            raise ValueError("Cable does not fit inside the enclosure based on JSON dimensions.")

        vertical_offset = enclosure_inner_radius - cable_outer_radius

        # --- SCC cable positions (flat arrangement); 'burial_depth' is the cable center ---
        depth = self.arrangement['burial_depth']
        spacing = self.arrangement.get('spacing', 0)
        num_conductors = self.arrangement.get('num_conductors', 3)
        cable_centers = [(i * spacing, -depth) for i in range(num_conductors)]

        # --- Model Generation ---
        model = {
            'name': self.input_data.get('name', 'FLAT_HDPE_Enclosed_SCC_System'),
            'type': self.input_data.get('type', 'hdpe'),
            'note': self.input_data.get(
                'note', 'A parametric flat arrangement of SCC cables, each individually enclosed in an HDPE duct.'),
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
        for cp in cable_centers:
            conductor_id = self._add_cable_conductors(model, conductor_id, cp)

        # --- Inject one eccentric duct (displaced upward) per phase ---
        # Each conductor gets its own copy of 'enclosure_data' -- a shared dict
        # (by reference) would make the 'center_point' of the last processed
        # cable overwrite those of the previous phases.
        for v in model.values():
            if isinstance(v, dict) and v.get('conductor_name') == host_conductor:
                cable_center = v['center_point']
                enclosure_center = (cable_center[0], cable_center[1] + vertical_offset)
                v['enclosure'] = copy.deepcopy(enclosure_data)
                v['enclosure']['center_point'] = enclosure_center

        if not self.silent_mode:
            self._show_model(model)

        return model

    def flat_hdpe_enclosed_with_shared_ecc_model(self, host_conductor: str = 'sheath') -> Dict[str, Any]:
        """
        Generates a model for Configuration 4 (Figure 5.4): three power SCC
        cables in a flat arrangement, each installed inside its own HDPE duct
        buried in the soil -- same as `flat_hdpe_enclosed_model` --, plus an
        earth continuity conductor (ECC) sharing the duct of the third cable
        (same cable near the ECC in `flat_scc_with_ecc_cable_model`,
        Configuration 3).

        The position of the ECC relative to the third SCC cable is the same as
        in Configuration 3 -- the only physical difference between the two
        configurations is the addition of the HDPE duct. For that reason this
        method reuses, without modification, the center-to-center positioning
        formula ("precise" mode) of `flat_scc_with_ecc_cable_model`:
        `ecc_center = (last_cable_center + ecc_horizontal_gap,
        last_cable_center - ecc_vertical_gap)`, reading the same
        `arrangement.ecc_horizontal_gap`/`ecc_vertical_gap` from the JSON (equal
        to those of `andreata_case3.json`). There is no geometric fit
        computation inside the duct: this geometry (SCC + ECC sharing a duct) is
        non-canonical and has no analytical formulation in pyLCP -- this model
        is only used for the schematic (see `andreata_case4/README.md`).
        """
        # --- Cable and duct data (defined in the host conductor, e.g. 'sheath') ---
        host_conductor_data = getattr(self, host_conductor)
        host_insulation = host_conductor_data.get('insulation')
        enclosure_data = host_conductor_data.get('enclosure')

        if not enclosure_data:
            raise ValueError(f"Enclosure definition not found within conductor '{host_conductor}'.")

        assert self.ecc is not None, "ECC conductor data must be provided for the shared-duct model."

        # --- Vertical cable -> duct offset (same eccentric logic as
        # `flat_hdpe_enclosed_model`, applied identically to each phase) ---
        cable_outer_radius = host_conductor_data['outer_radius'] + (host_insulation['thickness'] if host_insulation else 0)
        enclosure_inner_radius = enclosure_data['inner_radius']

        if cable_outer_radius > enclosure_inner_radius:
            raise ValueError("Cable does not fit inside the enclosure based on JSON dimensions.")

        vertical_offset = enclosure_inner_radius - cable_outer_radius

        # --- SCC cable positions (flat arrangement); 'burial_depth' is the cable center ---
        depth = self.arrangement['burial_depth']
        spacing = self.arrangement.get('spacing', 0)
        num_conductors = self.arrangement.get('num_conductors', 3)
        cable_centers = [(i * spacing, -depth) for i in range(num_conductors)]

        # --- ECC position ("precise" mode of `flat_scc_with_ecc_cable_model`,
        # center-to-center relative to the last SCC cable) ---
        last_cable_center = cable_centers[-1]
        ecc_horizontal_gap = self.arrangement.get('ecc_horizontal_gap', 0.0)
        ecc_vertical_gap = self.arrangement['ecc_vertical_gap']
        ecc_center = (
            last_cable_center[0] + ecc_horizontal_gap,
            last_cable_center[1] - ecc_vertical_gap,
        )

        # --- Model Generation ---
        model = {
            'name': self.input_data.get('name', 'FLAT_HDPE_Enclosed_SCC_System_with_shared_ECC'),
            'type': self.input_data.get('type', 'scc-flat-hdpe-ecc'),
            'note': self.input_data.get(
                'note', 'A parametric flat arrangement of SCC cables, each individually enclosed in an '
                        'HDPE duct, with an ECC sharing the third duct.'),
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

        # --- Add Conductors: 3 SCC cables (1..6), then the ECC (7) ---
        conductor_id = 1
        for cp in cable_centers:
            conductor_id = self._add_cable_conductors(model, conductor_id, cp)
        conductor_id = self._add_ecc_conductor(model, conductor_id, ecc_center)

        # --- Inject one eccentric duct (displaced upward) per phase --
        # Each conductor gets its own copy of 'enclosure_data', same as
        # `flat_hdpe_enclosed_model`. The ECC gets no duct of its own (it has
        # no 'enclosure' field), same convention as `flat_scc_with_ecc_cable_model`.
        for v in model.values():
            if isinstance(v, dict) and v.get('conductor_name') == host_conductor:
                cable_center = v['center_point']
                enclosure_center = (cable_center[0], cable_center[1] + vertical_offset)
                v['enclosure'] = copy.deepcopy(enclosure_data)
                v['enclosure']['center_point'] = enclosure_center

        if not self.silent_mode:
            self._show_model(model)

        return model

    def hdpe_shared_enclosed_model(self, host_conductor: str = 'sheath') -> Dict[str, Any]:
        """
        Generates a model for an SCC cable and an ECC sharing an HDPE duct.

        Geometry:
        1. 'burial_depth' defines the center of the HDPE duct.
        2. The SCC (Main Cable) rests on the bottom of the HDPE duct.
        3. The ECC (Earth Continuity Conductor) is positioned in the "wedge"
           space between the outer surface of the SCC and the inner wall of the HDPE duct.
        """

        # --- 1. SCC and Duct data ---
        host_conductor_data = getattr(self, host_conductor)
        host_insulation = host_conductor_data.get('insulation')
        enclosure_data = host_conductor_data.get('enclosure')

        # --- 2. ECC data ---
        ecc_insulation = self.ecc.get('insulation')
        ecc_outer_radius = self.ecc['outer_radius'] + (ecc_insulation['thickness'] if ecc_insulation else 0)

        assert enclosure_data is not None, "Enclosure definition must be provided for shared enclosure model."
        assert self.ecc is not None, "ECC conductor data must be provided for shared enclosure model."

        # --- 3. Position Computation ---

        # 3.1. 'burial_depth' defines the SCC center (standard installation convention).
        cable_outer_radius = host_conductor_data['outer_radius'] + (host_insulation['thickness'] if host_insulation else 0)
        enclosure_inner_radius = enclosure_data['inner_radius']

        # Geometry validations
        assert enclosure_inner_radius > cable_outer_radius, "Enclosure inner radius must be larger than cable outer radius."
        total_width_check = cable_outer_radius + 2 * ecc_outer_radius + cable_outer_radius
        assert total_width_check < (2 * enclosure_inner_radius), "The cable and ECC do not fit side by side within the HDPE enclosure."

        # 3.2. The SCC rests on the bottom of the duct; the duct center is vertical_offset above.
        vertical_offset = enclosure_inner_radius - cable_outer_radius
        cable_center = (0.0, -self.arrangement['burial_depth'])
        enclosure_center = (cable_center[0], cable_center[1] + vertical_offset)

        # --- 4. ECC position (Constraint 3) ---
        # Delegates the trigonometric computation to the private method
        ecc_center = self._calculate_ecc_center_trig(
            enclosure_center=enclosure_center,
            R_enc=enclosure_inner_radius,
            R_scc=cable_outer_radius,
            R_ecc=ecc_outer_radius,
            vertical_offset_c=vertical_offset
        )

        # --- 5. Model Generation ---
        model = {
            'name': self.input_data.get('name', 'generic shared HDPE enclosed SCC system'),
            'type': self.input_data.get('type', 'shared-hdpe'),
            'note': self.input_data.get('note', 'NA'),
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

        # --- 6. Add Conductors ---
        conductor_id = 1

        # Adds the SCC at 'cable_center' (computed to rest on the bottom)
        conductor_id = self._add_cable_conductors(model, conductor_id, cable_center)

        # Adds the ECC at 'ecc_center' (computed to be in the "corner")
        conductor_id = self._add_ecc_conductor(model, conductor_id, ecc_center)

        # --- 7. Inject Duct (Enclosure) Data ---
        for k, v in model.items():
            if isinstance(k, int) and k > 0:
                if v.get('conductor_name') == host_conductor:
                    # Adds the duct dictionary, including its center (primary reference).
                    v['enclosure'] = enclosure_data
                    v['enclosure']['center_point'] = enclosure_center
                else:
                    v['enclosure'] = None

        if not self.silent_mode:
            self._show_model(model)

        return model

    def flat_scc_with_ecc_cable_model(self, host_conductor: str = 'sheath') -> Dict[str, Any]:
        """
        Generates a model for Configuration 3 (Figure 5.3): three power SCC
        cables in a flat arrangement, directly buried in the soil (no HDPE
        duct), plus an earth continuity conductor (ECC) near the rightmost
        cable, without touching it.

        The ECC position is controlled by 'arrangement' parameters, in one of
        two mutually exclusive modes:

        Precise mode (used when 'ecc_vertical_gap' is present in the JSON):
          - 'ecc_horizontal_gap' / 'ecc_vertical_gap': distances (in meters,
            may be negative) between the ECC center and the center of the third
            SCC cable (the rightmost one), taken center-to-center (not gap
            between surfaces). A positive 'ecc_vertical_gap' places the ECC
            deeper than the cable (same sign convention as 'burial_depth'); a
            positive 'ecc_horizontal_gap' moves it horizontally away from the
            cable.

        Legacy mode (used when 'ecc_vertical_gap' is absent):
          - 'ecc_alignment': 'center' (default) places the ECC center on the
            same horizontal axis as the SCC cable centers (all at
            '-burial_depth'); 'bottom_tangent' places the ECC resting on the
            same horizontal line tangent to the bottom surface of the SCC
            cables (a situation where all cables rest on the bottom of a common
            trench).
          - 'ecc_horizontal_gap': horizontal gap (in meters) between the outer
            surface of the last SCC cable and the outer surface of the ECC
            (cables "near" each other, not touching).
        """

        # --- 1. SCC data ---
        host_conductor_data = getattr(self, host_conductor)
        host_insulation = host_conductor_data.get('insulation')
        cable_outer_radius = host_conductor_data['outer_radius'] + (host_insulation['thickness'] if host_insulation else 0)

        # --- 2. ECC data ---
        assert self.ecc is not None, "ECC conductor data must be provided for the flat SCC + ECC model."
        ecc_insulation = self.ecc.get('insulation')
        ecc_outer_radius = self.ecc['outer_radius'] + (ecc_insulation['thickness'] if ecc_insulation else 0)

        # --- 3. SCC cable positions (flat arrangement, directly buried) ---
        # Configuration 3 is always three-phase (3 SCC cables); the cardinality does
        # not come from the JSON, analogous to conventional_three_phase_flat.
        depth = self.arrangement['burial_depth']
        spacing = self.arrangement['spacing']
        num_scc_conductors = 3
        cable_centers = [(i * spacing, -depth) for i in range(num_scc_conductors)]

        # --- 4. ECC position (near the last SCC cable, without touching) ---
        last_cable_center = cable_centers[-1]
        ecc_vertical_gap = self.arrangement.get('ecc_vertical_gap')

        if ecc_vertical_gap is not None:
            # Precise mode: 'ecc_horizontal_gap'/'ecc_vertical_gap' are
            # center-to-center distances between the ECC and the third SCC cable
            # (not gaps between outer surfaces).
            ecc_horizontal_gap = self.arrangement.get('ecc_horizontal_gap', 0.0)
            ecc_center = (
                last_cable_center[0] + ecc_horizontal_gap,
                last_cable_center[1] - ecc_vertical_gap,
            )
        else:
            # Legacy mode: 'ecc_alignment' + 'ecc_horizontal_gap' as the gap
            # between the outer surfaces of the cable and the ECC.
            ecc_alignment = self.arrangement.get('ecc_alignment', 'center')
            ecc_horizontal_gap = self.arrangement.get('ecc_horizontal_gap', 0.0)
            horizontal_offset = cable_outer_radius + ecc_horizontal_gap + ecc_outer_radius

            if ecc_alignment == 'center':
                # ECC on the same horizontal axis as the SCC cable centers.
                ecc_center = (last_cable_center[0] + horizontal_offset, -depth)
            elif ecc_alignment == 'bottom_tangent':
                # ECC resting on the line tangent to the bottom surface of the SCC cables
                # (all cables resting on the bottom of a common trench).
                trench_floor_y = -depth - cable_outer_radius
                ecc_center = (last_cable_center[0] + horizontal_offset, trench_floor_y + ecc_outer_radius)
            else:
                raise ValueError(
                    f"Unknown 'ecc_alignment' value: '{ecc_alignment}'. Expected 'center' or 'bottom_tangent'."
                )

        # --- 5. Model Generation ---
        model = {
            'name': self.input_data.get('name', 'generic flat ECC system'),
            'type': self.input_data.get('type', 'scc-flat-ecc'),
            'note': self.input_data.get('note', 'NA'),
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

        # --- 6. Add Conductors ---
        conductor_id = 1

        # Adds the three SCC cables in the flat arrangement.
        for cp in cable_centers:
            conductor_id = self._add_cable_conductors(model, conductor_id, cp)

        # Adds the ECC at 'ecc_center' (near the last cable, without touching).
        conductor_id = self._add_ecc_conductor(model, conductor_id, ecc_center)

        if not self.silent_mode:
            self._show_model(model)

        return model

    def underground_flat_model(self) -> Dict[str, Any]:
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
            'type': self.input_data.get('type', 'scc'),
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

    def isolated_wires(self) -> Dict[str, Any]:
        """
        Generates a model for a single isolated coaxial cable.
        """
        N = self.arrangement.get('num_conductors', 2)
        idx_ref_conductor = int(self.arrangement.get('idx_ref_conductor', 0))

        model = {
            'name': self.input_data.get('name', 'unknown'),
            'type': self.input_data.get('type', 'unknown'),
            'note': self.input_data.get('note', 'unknown'),
            'idx_ref_conductor': idx_ref_conductor,
            'frequency_driver': self.driver_frequency,
        }

        for i in range(N):
            x_coordinate = i * self.arrangement.get('spacing', 0)

            if self.core is not None:
                insulation_data = self.core.get('insulation')
                insulation_dict = None
            
                if insulation_data:
                    insulation_dict = {
                        'name': f"conductor_{i}_insulation",
                        'type': insulation_data.get('type', 'insulation'),
                        'center_point': (x_coordinate, 0.0),
                        'thickness': insulation_data['thickness'],
                        'relative_permittivity': insulation_data['relative_permittivity'],
                        'relative_permeability': 1.0,
                        'fourier_order': self.fourier_order,
                    }

                model[i] = {
                    'line_id': i,
                    'conductor_name': f"conductor_{i}",
                    'line_type': 'active' if i != idx_ref_conductor else 'return',
                    'line_return':  idx_ref_conductor if i != idx_ref_conductor else None,
                    'center_point': (x_coordinate, 0.0),
                    'radius': [self.core['inner_radius'], self.core['outer_radius']],
                    'conductivity': self.core.get('conductivity_S_per_m', 0),
                    'subconductors': None,
                    'insulation': insulation_dict,
                    'conductor_layers': None,
                    'relative_permeability': 1.0,
                    'relative_permittivity': 1.0,
                    'relative_permittivity_out': 1.0,
                    'potential_to_infinity': 1.0 if i != idx_ref_conductor else -1.0,
                    'fourier_order': self.fourier_order,
                }

        if not self.silent_mode:
            self._show_model(model)

        return model
    
    def simple_trefoil(self) -> Dict[str, Any]:
        """
        Generates a model for a trefoil arrangement of three-phase SCC cables.
        It is assumed that 'burial_depth' corresponds to the depth of the centers
        of the two bottom conductors (h3 in the reference figure). 'spacing' is the
        center-to-center distance between adjacent cables.
        """
        # Depth of the bottom cables (h3)
        depth_bottom = self.arrangement['burial_depth']
        spacing = self.arrangement['spacing']

        # --- Trefoil Geometry Computation ---
        # The height of the equilateral triangle formed by the cables.
        height = spacing * np.sqrt(3) / 2

        # The coordinates of the bottom cables (b and c) are known.
        y_bottom = -depth_bottom
        x_side = spacing / 2

        # The coordinate of the top cable (a) is computed from the base.
        # Its vertical position is the base position plus the triangle height.
        y_top = y_bottom + height

        center_points = [
            (0.0, y_top),         # Top cable (a)
            (-x_side, y_bottom),  # Bottom-left cable (b)
            (+x_side, y_bottom)   # Bottom-right cable (c)
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