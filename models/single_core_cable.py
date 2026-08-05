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
        Adiciona o condutor ECC (Earth Continuity Conductor) ao modelo
        em seu 'center_point' específico.
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
        Calcula a posição do ECC "entalado" (wedge) entre o SCC e o Duto HDPE
        usando a Lei dos Cossenos.

        Args:
            enclosure_center (P_enc): Ponto central do duto.
            R_enc: Raio interno do duto.
            R_scc: Raio externo do cabo SCC.
            R_ecc: Raio externo do cabo ECC.
            vertical_offset_c (c): Distância entre P_enc e P_scc.
        """
        
        # Validação de geometria
        if R_enc <= (R_scc + R_ecc):
            raise ValueError(
                f"Geometria impossível: O ECC (R={R_ecc}) não cabe no espaço "
                f"entre o SCC (R={R_scc}) e o duto HDPE (R={R_enc}). "
                f"Condição (R_enc > R_scc + R_ecc) falhou."
            )

        # Lados do triângulo formado pelos centros (P_enc, P_scc, P_ecc)
        a = R_scc + R_ecc  # Dist (P_scc -> P_ecc)
        b = R_enc - R_ecc  # Dist (P_enc -> P_ecc)
        c = vertical_offset_c  # Dist (P_enc -> P_scc)

        # Lei dos Cossenos para encontrar o ângulo 'alpha' no vértice P_enc
        # a² = b² + c² - 2bc*cos(alpha)
        numerator = (b**2) + (c**2) - (a**2)
        denominator = 2 * b * c

        # Tratamento de erro de ponto flutuante
        cos_alpha = max(min(numerator / denominator, 1.0), -1.0)
        
        alpha = np.arccos(cos_alpha)
        sin_alpha = np.sin(alpha)

        # Calcular o centro do ECC
        # (Assumindo que o ECC está no quadrante x+, y-)
        # O centro do ECC está a uma distância 'b' do 'enclosure_center',
        # rotacionado pelo ângulo 'alpha' a partir da linha vertical P_enc -> P_scc.
        
        ecc_center_x = enclosure_center[0] + b * sin_alpha
        ecc_center_y = enclosure_center[1] - b * cos_alpha # Subtrai pois o eixo é para baixo

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

    def hdpe_shared_enclosed_model(self, host_conductor: str = 'sheath') -> Dict[str, Any]:
        """
        Gera um modelo para um cabo SCC e um ECC compartilhando um duto HDPE.
        
        Geometria:
        1. 'burial_depth' define o centro do Duto HDPE.
        2. O SCC (Cabo Principal) repousa no fundo do duto HDPE.
        3. O ECC (Condutor de Continuidade) é posicionado no espaço "entalado" (wedge)
           entre a superfície externa do SCC e a parede interna do duto HDPE.
        """
        
        # --- 1. Dados do SCC e Duto ---
        host_conductor_data = getattr(self, host_conductor)
        host_insulation = host_conductor_data.get('insulation')
        enclosure_data = host_conductor_data.get('enclosure')

        # --- 2. Dados do ECC ---
        ecc_insulation = self.ecc.get('insulation')
        ecc_outer_radius = self.ecc['outer_radius'] + (ecc_insulation['thickness'] if ecc_insulation else 0)

        assert enclosure_data is not None, "Enclosure definition must be provided for shared enclosure model."
        assert self.ecc is not None, "ECC conductor data must be provided for shared enclosure model."
        
        # --- 3. Cálculo de Posição ---

        # 3.1. 'burial_depth' define o centro do SCC (convenção padrão de instalação).
        cable_outer_radius = host_conductor_data['outer_radius'] + (host_insulation['thickness'] if host_insulation else 0)
        enclosure_inner_radius = enclosure_data['inner_radius']

        # Validações de geometria
        assert enclosure_inner_radius > cable_outer_radius, "Enclosure inner radius must be larger than cable outer radius."
        total_width_check = cable_outer_radius + 2 * ecc_outer_radius + cable_outer_radius
        assert total_width_check < (2 * enclosure_inner_radius), "The cable and ECC do not fit side by side within the HDPE enclosure."

        # 3.2. O SCC repousa no fundo do duto; o centro do duto fica vertical_offset acima.
        vertical_offset = enclosure_inner_radius - cable_outer_radius
        cable_center = (0.0, -self.arrangement['burial_depth'])
        enclosure_center = (cable_center[0], cable_center[1] + vertical_offset)

        # --- 4. Posição do ECC (Restrição 3) ---
        # Delega o cálculo trigonométrico para o método privado
        ecc_center = self._calculate_ecc_center_trig(
            enclosure_center=enclosure_center,
            R_enc=enclosure_inner_radius,
            R_scc=cable_outer_radius,
            R_ecc=ecc_outer_radius,
            vertical_offset_c=vertical_offset
        )

        # --- 5. Geração do Modelo ---
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

        # --- 6. Adicionar Condutores ---
        conductor_id = 1
        
        # Adiciona o SCC no 'cable_center' (calculado para repousar no fundo)
        conductor_id = self._add_cable_conductors(model, conductor_id, cable_center)

        # Adiciona o ECC no 'ecc_center' (calculado para estar no "canto")
        conductor_id = self._add_ecc_conductor(model, conductor_id, ecc_center)

        # --- 7. Injetar Dados do Duto (Enclosure) ---
        for k, v in model.items():
            if isinstance(k, int) and k > 0: 
                if v.get('conductor_name') == host_conductor:
                    # Adiciona o dicionário do duto, incluindo seu centro (referência primária).
                    v['enclosure'] = enclosure_data
                    v['enclosure']['center_point'] = enclosure_center
                else:
                    v['enclosure'] = None

        if not self.silent_mode:
            self._show_model(model)

        return model

    def flat_scc_with_ecc_cable_model(self, host_conductor: str = 'sheath') -> Dict[str, Any]:
        """
        Gera um modelo para a Configuração 3 (Figura 5.3): três cabos SCC de potência
        em arranjo plano (flat), diretamente enterrados no solo (sem duto HDPE), mais
        um cabo de aterramento (ECC) próximo ao cabo mais à direita, sem encostar
        nele.

        A posição do ECC é controlada por parâmetros de 'arrangement', em um de
        dois modos mutuamente exclusivos:

        Modo preciso (usado quando 'ecc_vertical_gap' está presente no JSON):
          - 'ecc_horizontal_gap' / 'ecc_vertical_gap': distâncias (em metros,
            podem ser negativas) entre o centro do ECC e o centro do terceiro
            cabo SCC (o mais à direita), tomadas centro-a-centro (não folga
            entre superfícies). 'ecc_vertical_gap' positivo posiciona o ECC
            mais fundo que o cabo (mesma convenção de sinal de
            'burial_depth'); 'ecc_horizontal_gap' positivo o afasta
            horizontalmente do cabo.

        Modo legado (usado quando 'ecc_vertical_gap' está ausente):
          - 'ecc_alignment': 'center' (padrão) posiciona o centro do ECC no mesmo
            eixo horizontal dos centros dos cabos SCC (todos a '-burial_depth');
            'bottom_tangent' posiciona o ECC apoiado na mesma linha horizontal
            tangente à superfície inferior dos cabos SCC (situação em que todos
            repousam no fundo de uma vala comum).
          - 'ecc_horizontal_gap': folga horizontal (em metros) entre a superfície externa do
            último cabo SCC e a superfície externa do ECC (cabos "próximos", não
            encostados).
        """

        # --- 1. Dados do SCC ---
        host_conductor_data = getattr(self, host_conductor)
        host_insulation = host_conductor_data.get('insulation')
        cable_outer_radius = host_conductor_data['outer_radius'] + (host_insulation['thickness'] if host_insulation else 0)

        # --- 2. Dados do ECC ---
        assert self.ecc is not None, "ECC conductor data must be provided for the flat SCC + ECC model."
        ecc_insulation = self.ecc.get('insulation')
        ecc_outer_radius = self.ecc['outer_radius'] + (ecc_insulation['thickness'] if ecc_insulation else 0)

        # --- 3. Posição dos Cabos SCC (arranjo plano, diretamente enterrados) ---
        # Configuração 3 é sempre trifásica (3 cabos SCC); a cardinalidade não
        # vem do JSON, análogo a conventional_three_phase_flat.
        depth = self.arrangement['burial_depth']
        spacing = self.arrangement['spacing']
        num_scc_conductors = 3
        cable_centers = [(i * spacing, -depth) for i in range(num_scc_conductors)]

        # --- 4. Posição do ECC (próximo ao último cabo SCC, sem encostar) ---
        last_cable_center = cable_centers[-1]
        ecc_vertical_gap = self.arrangement.get('ecc_vertical_gap')

        if ecc_vertical_gap is not None:
            # Modo preciso: 'ecc_horizontal_gap'/'ecc_vertical_gap' são distâncias
            # centro-a-centro entre o ECC e o terceiro cabo SCC (não folgas entre
            # superfícies externas).
            ecc_horizontal_gap = self.arrangement.get('ecc_horizontal_gap', 0.0)
            ecc_center = (
                last_cable_center[0] + ecc_horizontal_gap,
                last_cable_center[1] - ecc_vertical_gap,
            )
        else:
            # Modo legado: 'ecc_alignment' + 'ecc_horizontal_gap' como folga entre
            # as superfícies externas do cabo e do ECC.
            ecc_alignment = self.arrangement.get('ecc_alignment', 'center')
            ecc_horizontal_gap = self.arrangement.get('ecc_horizontal_gap', 0.0)
            horizontal_offset = cable_outer_radius + ecc_horizontal_gap + ecc_outer_radius

            if ecc_alignment == 'center':
                # ECC no mesmo eixo horizontal dos centros dos cabos SCC.
                ecc_center = (last_cable_center[0] + horizontal_offset, -depth)
            elif ecc_alignment == 'bottom_tangent':
                # ECC apoiado na linha tangente à superfície inferior dos cabos SCC
                # (todos os cabos repousando no fundo de uma vala comum).
                trench_floor_y = -depth - cable_outer_radius
                ecc_center = (last_cable_center[0] + horizontal_offset, trench_floor_y + ecc_outer_radius)
            else:
                raise ValueError(
                    f"Unknown 'ecc_alignment' value: '{ecc_alignment}'. Expected 'center' or 'bottom_tangent'."
                )

        # --- 5. Geração do Modelo ---
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

        # --- 6. Adicionar Condutores ---
        conductor_id = 1

        # Adiciona os três cabos SCC no arranjo plano.
        for cp in cable_centers:
            conductor_id = self._add_cable_conductors(model, conductor_id, cp)

        # Adiciona o ECC em 'ecc_center' (próximo ao último cabo, sem encostar).
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