from abc import ABC, abstractmethod
import numpy as np
import scipy.constants as sc

class MTLStrategy(ABC):
    """ Abstract Base Class defining the interface for MTL type-specific logic."""

    @abstractmethod
    def preprocess_mtl_data(self, mtl_data: dict) -> dict:
        """
        Filters and preprocesses the MTL dictionary based on the type.
        For example, 'overhead' lines might have different conductor indexing.
        """
        pass

    @abstractmethod
    def validate(self, mtl_data: dict, idx_ref: int):
        """
        Performs specific validations for the MTL type.
        """
        pass
    
    def apply_mtl_ref_properties(self, context, mtl_ref: dict) -> None:
        """
        Calculates and applies type-specific properties directly to the context object.
        The 'context' is an instance of MulticonductorTransmissionLine.
        """
        # Default implementation does nothing
        pass

    def apply_mtl_properties(self, context, mtl: dict) -> None:
        """
        Calculates and applies type-specific properties directly to the context object.
        The 'context' is an instance of MulticonductorTransmissionLine.
        """
        # Default implementation does nothing
        pass
    
    @staticmethod
    def _conductors_center_distance_matrix(mtl):
        """ Calculates the distance matrix between the center points of the conductors. """
        d_pq = np.zeros((len(mtl), len(mtl)))
        x_pq = np.zeros_like(d_pq)
        y_pq = np.zeros_like(d_pq)
        theta_pq = np.zeros_like(d_pq)

        for p, conductor_p in mtl.items():
            for q, conductor_q in mtl.items():
                cp = conductor_p['center_point']
                cq = conductor_q['center_point']

                # Distance between the center points of the conductors
                x_pq[p, q] = cp[0] - cq[0]
                y_pq[p, q] = cp[1] - cq[1]

                # Angle between the center points of the conductors
                theta_pq[p, q] = np.arctan2(y_pq[p, q], x_pq[p, q])

                # Distance between the center points of the conductors
                d_pq[p, q] = np.linalg.norm(np.array(cp) - np.array(cq))

        return {
            'distance_pq': d_pq,
            'x_pq': x_pq,
            'y_pq': y_pq,
            'theta_pq': theta_pq
        }

    @staticmethod
    def _surface_distance_with_ground_return(mtl: dict) -> dict:
        """
        Calculates all necessary distance matrices for overhead line analysis,
        including geometric, ground return, horizontal, and vertical separations.
        Returns a dictionary of the calculated matrices.
        """
        N = len(mtl)
        d_matrix = np.zeros((N, N))
        D_matrix = np.zeros((N, N))
        vertical_separation_matrix = np.zeros((N, N))
        horizontal_separation_matrix = np.zeros((N, N))

        # Create a sorted list of conductors to ensure consistent ordering
        conductors = sorted(mtl.items())

        for n_idx, (n_tag, n_cond_data) in enumerate(conductors):
            for m_idx, (m_tag, m_cond_data) in enumerate(conductors):
                cn = n_cond_data['center_point']
                cm = m_cond_data['center_point']

                # Horizontal and Vertical separation
                dnm = cn[0] - cm[0]
                hnm = cn[1] + cm[1]

                # Geometric Distance (d_nm)
                if n_tag == m_tag:
                    d = n_cond_data['radius'][1] # Use radius for self-distance
                else:
                    d = np.sqrt(dnm ** 2 + (cn[1] - cm[1]) ** 2)

                # Distance to Image (D_nm)
                D = np.sqrt(dnm ** 2 + hnm ** 2)

                # Populate the matrices at the correct indices
                d_matrix[n_idx, m_idx] = d
                D_matrix[n_idx, m_idx] = D
                vertical_separation_matrix[n_idx, m_idx] = hnm
                horizontal_separation_matrix[n_idx, m_idx] = dnm

        return {
            'd_matrix_ground_return': d_matrix,
            'D_matrix_ground_return': D_matrix,
            'vertical_separation_matrix': vertical_separation_matrix,
            'horizontal_separation_matrix': horizontal_separation_matrix
        }

    # @staticmethod
    # def _cable_distance_with_ground_return(mtl: dict) -> dict:
    #     """
    #     Calculates all necessary distance matrices for overhead line analysis,
    #     including geometric, ground return, horizontal, and vertical separations.
    #     Returns a dictionary of the calculated matrices.
    #     """
    #     # Create a sorted list of CORE conductors to ensure consistent ordering
    #     core_conductors = sorted({key: value for key, value in mtl.items() if value.get('conductor_name') in ['core']}.items())
        
    #     # Count number of CORE conductors 
    #     N = len(core_conductors)
        
    #     d_matrix = np.zeros((N, N))
    #     D_matrix = np.zeros((N, N))
    #     vertical_separation_matrix = np.zeros((N, N))
    #     horizontal_separation_matrix = np.zeros((N, N))

    #     for n_idx, (n_tag, n_conductor) in enumerate(core_conductors):
    #         for m_idx, (m_tag, m_conductor) in enumerate(core_conductors):
    #             cn = n_conductor['center_point']
    #             cm = m_conductor['center_point']

    #             # Horizontal and Vertical separation
    #             dnm = cn[0] - cm[0]
    #             hnm = cn[1] + cm[1]

    #             # Geometric Distance (d_nm)
    #             if n_tag == m_tag:
    #                 d = n_conductor['radius'][1] 
    #             else:
    #                 d = np.sqrt(dnm ** 2 + (cn[1] - cm[1]) ** 2)

    #             # Distance to Image (D_nm)
    #             D = np.sqrt(dnm ** 2 + hnm ** 2)

    #             # Populate the matrices at the correct indices
    #             d_matrix[n_idx, m_idx] = d
    #             D_matrix[n_idx, m_idx] = D
    #             vertical_separation_matrix[n_idx, m_idx] = hnm
    #             horizontal_separation_matrix[n_idx, m_idx] = dnm

    #     return {
    #         'd_matrix_ground_return': d_matrix,
    #         'D_matrix_ground_return': D_matrix,
    #         'vertical_separation_matrix': vertical_separation_matrix,
    #         'horizontal_separation_matrix': horizontal_separation_matrix
    #     }

    @staticmethod
    def _extract_scc_parameters(mtl: dict) -> dict:
        """
        Extracts geometric and physical parameters from a single-core cable
        data structure using descriptive names for clarity. It handles solid
        cores (inner radius = 0) and hollow layers.

        Args:
            mtl (dict): The dictionary containing conductor data for the SCC.

        Returns:
            dict: A dictionary with the calculated parameters (radii, rho, mu, epsilon).
        """
        scc = {}
        core, sheath, armor = None, None, None

        # Identify each layer by its name
        for conductor_data in mtl.values():
            name = conductor_data.get('conductor_name')
            if name == 'core':
                core = conductor_data
            elif name == 'sheath':
                sheath = conductor_data
            elif name == 'armor':
                armor = conductor_data
        
        # === CORE ===
        if core:
            scc['core_inner_radius'], scc['core_outer_radius'] = core['radius']
            if 'insulation' in core and core['insulation']:
                scc['core_insulation_outer_radius'] = scc['core_outer_radius'] + core['insulation']['thickness']
                scc['core_insulation_permittivity'] = core['insulation']['relative_permittivity'] * sc.epsilon_0
                scc['core_insulation_permeability'] = core['insulation']['relative_permeability'] * sc.mu_0
            
            # Extract physical properties for the core conductor (layer 1)
            scc['core_resistivity'] = 1 / core['conductivity']
            scc['core_permeability'] = core['relative_permeability'] * sc.mu_0
            scc['core_permittivity'] = core['relative_permittivity'] * sc.epsilon_0
        
        # === SHEATH ===
        if sheath and scc.get('core_insulation_outer_radius') is not None:
            assert np.isclose(scc['core_insulation_outer_radius'], sheath['radius'][0]), \
                (f"Geometric mismatch: Core's insulation outer radius ({scc['core_insulation_outer_radius']}) "
                     f"does not match sheath's inner radius ({sheath['radius'][0]})")

            scc['sheath_inner_radius'], scc['sheath_outer_radius'] = sheath['radius']
            if 'insulation' in sheath and sheath['insulation']:
                scc['sheath_insulation_outer_radius'] = scc['sheath_outer_radius'] + sheath['insulation']['thickness']
                scc['sheath_insulation_permittivity'] = sheath['insulation']['relative_permittivity'] * sc.epsilon_0
                scc['sheath_insulation_permeability'] = sheath['insulation']['relative_permeability'] * sc.mu_0

            # Extract physical properties for the sheath conductor (layer 2)
            scc['sheath_resistivity'] = 1 / sheath['conductivity']
            scc['sheath_permeability'] = sheath['relative_permeability'] * sc.mu_0
            scc['sheath_permittivity'] = sheath['relative_permittivity'] * sc.epsilon_0

        # === ARMOR ===
        if armor and scc.get('sheath_insulation_outer_radius') is not None:
            assert np.isclose(scc['sheath_insulation_outer_radius'], armor['radius'][0]), \
                (f"Geometric mismatch: Sheath's insulation outer radius ({scc['sheath_insulation_outer_radius']}) "
                     f"does not match armor's inner radius ({armor['radius'][0]})")
            
            scc['armor_inner_radius'], scc['armor_outer_radius'] = armor['radius']
            if 'insulation' in armor and armor['insulation']:
                scc['armor_insulation_outer_radius'] = scc['armor_outer_radius'] + armor['insulation']['thickness']
                scc['armor_insulation_permittivity'] = armor['insulation']['relative_permittivity'] * sc.epsilon_0
                scc['armor_insulation_permeability'] = armor['insulation']['relative_permeability'] * sc.mu_0

            # Extract physical properties for the armor conductor (layer 3)
            scc['armor_resistivity'] = 1 / armor['conductivity']
            scc['armor_permeability'] = armor['relative_permeability'] * sc.mu_0
            scc['armor_permittivity'] = armor['relative_permittivity'] * sc.epsilon_0

        return scc
    
class SingleCoreCableStrategy(MTLStrategy):
    """Strategy for single-core cable MTLs."""

    def preprocess_mtl_data(self, mtl_input: dict) -> dict:
        # For single-core cables (scc), we expect integer keys starting from 1.
        mtl = {key: value for key, value in mtl_input.items() if isinstance(key, int) and key > 0}
        mtl_ref = {key: value for key, value in mtl_input.items() if key == 0 and value.get('line_type') == 'return' and value.get('conductor_name') == 'soil'}
        return mtl, mtl_ref

    def validate(self, mtl: dict, mtl_ref: dict, mtl_input: dict):
        key_conductors = sorted(mtl.keys())
        key_expected = list(range(1, len(key_conductors) + 1))

        assert len(mtl) > 0, "Single Core Cable-based MTL must have at least one conductor."
        assert key_conductors == key_expected, "Cable conductor tags must be a sequence starting from 1."
        # assert mtl_data[0]['line_type'] == 'return', "Conductor with index '0' must be the return path for cables."
        # assert idx_ref in mtl_data, f"The reference conductor index {idx_ref} must be in the MTL dictionary."

    def _count_scc_and_conductors(self, mtl_input: dict) -> tuple:
        """
        Counts the number of (sc) cables (N) and conductors per cable (M)
        based on the provided data structure.
        """
        # 1. Filter to get only active conductors
        active_conductors = [
            v for k, v in mtl_input.items()
            if isinstance(k, int) and v.get('line_type') == 'active'
        ]

        if not active_conductors:
            return 0, 0

        # 2. Count the number of cables (N) by finding unique center points
        num_cables = len(set([cond['center_point'] for cond in active_conductors]))

        # 3. Count conductors per cable (M)
        total_active_conductors = len(active_conductors)
        
        if num_cables > 0:
            conductors_per_cable = total_active_conductors // num_cables
        else:
            conductors_per_cable = 0

        return num_cables, conductors_per_cable
    
    def _cable_distance_with_ground_return(self, mtl: dict) -> dict:
        """
        Calculates all necessary distance matrices for overhead line analysis,
        including geometric, ground return, horizontal, and vertical separations.
        Returns a dictionary of the calculated matrices.
        """
        # Create a sorted list of conductors to ensure consistent ordering
        cables = sorted({key: value for key, value in mtl.items() if value.get('conductor_name') in ['sheath']}.items())
        
        # Count number of cables conductors 
        N = len(cables)
        
        d_matrix = np.zeros((N, N))
        D_matrix = np.zeros((N, N))
        vertical_separation_matrix = np.zeros((N, N))
        horizontal_separation_matrix = np.zeros((N, N))

        for n_idx, (n_tag, n_conductor) in enumerate(cables):
            for m_idx, (m_tag, m_conductor) in enumerate(cables):
                cn = n_conductor['center_point']
                cm = m_conductor['center_point']

                # Horizontal and Vertical separation
                dnm = cn[0] - cm[0]
                hnm = cn[1] + cm[1]

                # Geometric Distance (d_nm)
                if n_tag == m_tag:
                    d = n_conductor['radius'][1] + n_conductor['insulation']['thickness'] # Use outer radius for self-distance
                else:
                    d = np.sqrt(dnm ** 2 + (cn[1] - cm[1]) ** 2)

                # Distance to Image (D_nm)
                D = np.sqrt(dnm ** 2 + hnm ** 2)

                # Populate the matrices at the correct indices
                d_matrix[n_idx, m_idx] = d
                D_matrix[n_idx, m_idx] = D
                vertical_separation_matrix[n_idx, m_idx] = hnm
                horizontal_separation_matrix[n_idx, m_idx] = dnm

        return {
            'd_matrix_ground_return': d_matrix,
            'D_matrix_ground_return': D_matrix,
            'vertical_separation_matrix': vertical_separation_matrix,
            'horizontal_separation_matrix': horizontal_separation_matrix
        }
    
    def apply_mtl_ref_properties(self, context, mtl_input: dict) -> None:
        # Number of single core cables (N) and conductors per cable (M)
        context.num_sc_cables, context.num_conductors_per_scc = self._count_scc_and_conductors(mtl_input)

    def apply_mtl_properties(self, context, mtl: dict) -> None:
        """
        Calculates and applies overhead-line-specific distance matrices to the MTL object
        by delegating the calculation to a static helper method.
        """
        # 1. Delegate the complex calculation to the static method
        properties = self._cable_distance_with_ground_return(mtl)
        context.d_matrix_ground_return = properties['d_matrix_ground_return']
        context.D_matrix_ground_return = properties['D_matrix_ground_return']
        context.vertical_separation_matrix = properties['vertical_separation_matrix']
        context.horizontal_separation_matrix = properties['horizontal_separation_matrix']

        # Extract and apply SCC geometric parameters
        context.scc = MTLStrategy._extract_scc_parameters(mtl)

        # Conductors Permeability [np.array]
        context.mu = np.array([sc.mu_0 * conductor['relative_permeability'] for conductor in mtl.values()]) 

        # Conductors Permittivity [np.array]
        context.epsilon = np.array([sc.epsilon_0 * conductor['relative_permittivity'] for conductor in mtl.values()]) 

        # Conductors conductivity [np.array]
        context.sigma = np.array([conductor['conductivity'] for conductor in mtl.values()])

        # Free-Space Permittivity [np.array]
        context.epsilon_out = np.array([sc.epsilon_0 * conductor['relative_permittivity_out'] for conductor in mtl.values()])

class SingleCoreCableInHDPEStrategy(MTLStrategy):
    """Strategy for single-core cable MTLs."""

    def preprocess_mtl_data(self, mtl_input: dict) -> dict:
        # For single-core cables (scc), we expect integer keys starting from 1.
        mtl = {key: value for key, value in mtl_input.items() if isinstance(key, int) and key > 0}
        mtl_ref = {key: value for key, value in mtl_input.items() if key == 0 and value.get('line_type') == 'return' and value.get('conductor_name') == 'soil'}
        return mtl, mtl_ref

    def validate(self, mtl: dict, mtl_ref: dict, mtl_input: dict):
        key_conductors = sorted(mtl.keys())
        key_expected = list(range(1, len(key_conductors) + 1))

        assert len(mtl) > 0, "Single Core Cable-based MTL must have at least one conductor."
        assert key_conductors == key_expected, "Cable conductor tags must be a sequence starting from 1."
        # assert mtl_data[0]['line_type'] == 'return', "Conductor with index '0' must be the return path for cables."
        # assert idx_ref in mtl_data, f"The reference conductor index {idx_ref} must be in the MTL dictionary."

    def _count_scc_and_conductors(self, mtl_input: dict) -> tuple:
        """
        Counts the number of (sc) cables (N) and conductors per cable (M)
        based on the provided data structure.
        """
        # 1. Filter to get only active conductors
        active_conductors = [
            v for k, v in mtl_input.items()
            if isinstance(k, int) and v.get('line_type') == 'active'
        ]

        if not active_conductors:
            return 0, 0

        # 2. Count the number of cables (N) by finding unique center points
        num_cables = len(set([cond['center_point'] for cond in active_conductors]))

        # 3. Count conductors per cable (M)
        total_active_conductors = len(active_conductors)
        
        if num_cables > 0:
            conductors_per_cable = total_active_conductors // num_cables
        else:
            conductors_per_cable = 0

        return num_cables, conductors_per_cable
    
    def _extract_hdpe_parameters(self, mtl: dict) -> dict:
        """
        Extracts geometric and physical parameters from a single-core cable
        data structure using descriptive names for clarity. It handles solid
        cores (inner radius = 0) and hollow layers.

        Args:
            mtl (dict): The dictionary containing conductor data for the SCC.

        Returns:
            dict: A dictionary with the calculated parameters (radii, rho, mu, epsilon).
        """
        hdpe = {'core': {}, 'sheath': {}, 'armor': {}}
        core_enclosure, sheath_enclosure, armor_enclosure = None, None, None

        # Identify each enclosure by its name
        for conductor_data in mtl.values():
            name = conductor_data.get('conductor_name')
            if name == 'core' and conductor_data.get('enclosure') is not None:
                core_enclosure = conductor_data['enclosure']
            elif name == 'sheath' and conductor_data.get('enclosure') is not None:
                sheath_enclosure = conductor_data['enclosure']
            elif name == 'armor' and conductor_data.get('enclosure') is not None:
                armor_enclosure = conductor_data['enclosure']
        
        # === CORE ===
        if core_enclosure:
            pass
        
        # === SHEATH ===
        if sheath_enclosure:
            hdpe['sheath']['inner_radius'] = sheath_enclosure['inner_radius']
            hdpe['sheath']['outer_radius'] = sheath_enclosure['outer_radius']
            
            if 'insulation' in sheath_enclosure and sheath_enclosure['insulation']:
                hdpe['sheath']['insulation_permittivity'] = sheath_enclosure['insulation']['relative_permittivity'] * sc.epsilon_0

            # Extract physical properties for the sheath conductor (layer 2)
            hdpe['sheath']['permittivity'] = sheath_enclosure['relative_permittivity'] * sc.epsilon_0

        # === ARMOR ===
        if armor_enclosure and hdpe.get('sheath_insulation_outer_radius') is not None:
            pass

        return hdpe
    
    def apply_mtl_ref_properties(self, context, mtl_input: dict) -> None:
        # Number of single core cables (N) and conductors per cable (M)
        context.num_sc_cables, context.num_conductors_per_scc = self._count_scc_and_conductors(mtl_input)

    def apply_mtl_properties(self, context, mtl: dict) -> None:
        """
        Calculates and applies overhead-line-specific distance matrices to the MTL object
        by delegating the calculation to a static helper method.
        """
        # 1. Delegate the complex calculation to the static method
        properties = MTLStrategy._cable_distance_with_ground_return(mtl)
        context.d_matrix_ground_return = properties['d_matrix_ground_return']
        context.D_matrix_ground_return = properties['D_matrix_ground_return']
        context.vertical_separation_matrix = properties['vertical_separation_matrix']
        context.horizontal_separation_matrix = properties['horizontal_separation_matrix']

        # Extract and apply SCC geometric parameters
        context.scc = MTLStrategy._extract_scc_parameters(mtl)
        context.scc['hdpe'] = self._extract_hdpe_parameters(mtl)

        # Conductors Permeability [np.array]
        context.mu = np.array([sc.mu_0 * conductor['relative_permeability'] for conductor in mtl.values()]) 

        # Conductors Permittivity [np.array]
        context.epsilon = np.array([sc.epsilon_0 * conductor['relative_permittivity'] for conductor in mtl.values()]) 

        # Conductors conductivity [np.array]
        context.sigma = np.array([conductor['conductivity'] for conductor in mtl.values()])

        # Free-Space Permittivity [np.array]
        context.epsilon_out = np.array([sc.epsilon_0 * conductor['relative_permittivity_out'] for conductor in mtl.values()])

class CableStrategy(MTLStrategy):
    """Strategy for cable-based MTLs like 'coaxial', 'coated_wires', etc."""

    def preprocess_mtl_data(self, mtl_input: dict) -> dict:
        # For cables, we expect integer keys starting from 0.
        mtl =  {key: value for key, value in mtl_input.items() if isinstance(key, int)}
        mtl_ref = {key: value for key, value in mtl_input.items() if key == 0 and value.get('line_type') == 'return'}
        return mtl, mtl_ref

    def validate(self, mtl: dict,  mtl_ref: dict, mtl_input: dict):
        key_conductors = sorted(mtl.keys())
        key_expected = list(range(len(key_conductors)))

        assert len(mtl) > 1, "Cable-based MTL must have at least two conductors."
        assert key_conductors == key_expected, "Cable conductor tags must be a sequence starting from 0."
        assert mtl[0]['line_type'] == 'return', "Conductor with index '0' must be the return path for cables."
        # assert idx_ref in mtl_data, f"The reference conductor index {idx_ref} must be in the MTL dictionary."

    def _count_scc_and_conductors_plus_ref(self, mtl_input: dict) -> tuple:
        """
        Counts the number of (sc) cables (N) and conductors per cable (M)
        based on the provided data structure, including the reference conductor.
        """
        # 1. Filter to get only active conductors
        all_conductors = [v for k, v in mtl_input.items() if isinstance(k, int)]

        if not all_conductors:
            return 0, 0

        # 2. Count the number of cables (N) by finding unique center points
        num_cables = len(set([cond['center_point'] for cond in all_conductors]))

        # 3. Count conductors per cable (M)
        total_active_conductors = len(all_conductors)
        
        if num_cables > 0:
            conductors_per_cable = total_active_conductors // num_cables
        else:
            conductors_per_cable = 0

        return num_cables, conductors_per_cable
    
    def apply_mtl_ref_properties(self, context, mtl_input: dict) -> None:
        
        # Number of single core cables (N) and conductors per cable (M)
        context.num_sc_cables, context.num_conductors_per_scc = self._count_scc_and_conductors_plus_ref(mtl_input)

    def apply_mtl_properties(self, context, mtl: dict) -> None:
        """Applies cable-specific distance properties to the MTL object."""
        # 1. Delegate the complex calculation to the static method
        distances = MTLStrategy._conductors_center_distance_matrix(mtl)
        context.D_pq = distances['distance_pq']
        context.x_pq = distances['x_pq']
        context.y_pq = distances['y_pq']
        context.theta_pq = distances['theta_pq']

        # Extract and apply SCC geometric parameters
        context.scc = MTLStrategy._extract_scc_parameters(mtl)

        # Conductors Permeability [np.array]
        context.mu = np.array([sc.mu_0 * conductor['relative_permeability'] for conductor in mtl.values()]) 

        # Conductors Permittivity [np.array]
        context.epsilon = np.array([sc.epsilon_0 * conductor['relative_permittivity'] for conductor in mtl.values()]) 

        # Conductors conductivity [np.array]
        context.sigma = np.array([conductor['conductivity'] for conductor in mtl.values()])

        # Free-Space Permittivity [np.array]
        context.epsilon_out = np.array([sc.epsilon_0 * conductor['relative_permittivity_out'] for conductor in mtl.values()])

class OverheadCableStrategy(MTLStrategy):
    """Strategy for cable-based MTLs like 'coaxial', 'coated_wires', etc."""

    def preprocess_mtl_data(self, mtl_input: dict) -> dict:
        # For cables, we expect integer keys starting from 0.
        mtl =  {key: value for key, value in mtl_input.items() if isinstance(key, int)}
        mtl_ref = {key: value for key, value in mtl_input.items() if key == 0 and value.get('line_type') == 'return' and value.get('conductor_name') == 'soil'}
        return mtl, mtl_ref

    def validate(self, mtl: dict,  mtl_ref: dict, mtl_input: dict):
        key_conductors = sorted(mtl.keys())
        key_expected = list(range(len(key_conductors)))

        assert len(mtl) > 1, "Cable-based MTL must have at least two conductors."
        assert key_conductors == key_expected, "Cable conductor tags must be a sequence starting from 0."
        assert mtl[0]['line_type'] == 'return', "Conductor with index '0' must be the return path for cables."
        # assert idx_ref in mtl_data, f"The reference conductor index {idx_ref} must be in the MTL dictionary."

    def apply_mtl_properties(self, context, mtl: dict) -> None:
        """Applies cable-specific distance properties to the MTL object."""
        # 1. Delegate the complex calculation to the static method
        distances = MTLStrategy._conductors_center_distance_matrix(mtl)
        context.D_pq = distances['distance_pq']
        context.x_pq = distances['x_pq']
        context.y_pq = distances['y_pq']
        context.theta_pq = distances['theta_pq']

        # Extract and apply SCC geometric parameters
        context.scc = MTLStrategy._extract_scc_parameters(mtl)

        # Conductors Permeability [np.array]
        context.mu = np.array([sc.mu_0 * conductor['relative_permeability'] for conductor in mtl.values()]) 

        # Conductors Permittivity [np.array]
        context.epsilon = np.array([sc.epsilon_0 * conductor['relative_permittivity'] for conductor in mtl.values()]) 

        # Conductors conductivity [np.array]
        context.sigma = np.array([conductor['conductivity'] for conductor in mtl.values()])

        # Free-Space Permittivity [np.array]
        context.epsilon_out = np.array([sc.epsilon_0 * conductor['relative_permittivity_out'] for conductor in mtl.values()])

class OverheadLineStrategy(MTLStrategy):
    """Strategy for overhead lines."""

    def preprocess_mtl_data(self, mtl_input: dict) -> dict:
        # For overhead lines, conductors are indexed starting from 1.
        mtl = {key: value for key, value in mtl_input.items() if isinstance(key, int) and key > 0}
        mtl_ref = {key: value for key, value in mtl_input.items() if key == 0 and value.get('line_type') == 'return' and value.get('conductor_name') == 'soil'}
        return mtl, mtl_ref

    def validate(self, mtl: dict,  mtl_ref: dict, mtl_input: dict):
        key_conductors = sorted(mtl.keys())
        key_expected = list(range(1, len(key_conductors) + 1))
        
        assert len(mtl) > 0, "Overhead MTL must have at least one conductor."
        assert key_conductors == key_expected, "Overhead conductor tags must be a sequence starting from 1."
        # assert mtl_data[0]['line_type'] == 'return', "Conductor with index '0' must be the return path for cables."
        # assert idx_ref in mtl_data, f"The reference conductor index {idx_ref} must be in the MTL dictionary."

    def _count_scc_and_conductors(self, mtl_ref: dict) -> tuple:
        """
        Counts the number of (sc) cables (N) and conductors per cable (M)
        based on the provided data structure.
        """
        # 1. Filter to get only active conductors
        active_conductors = [
            v for k, v in mtl_ref.items()
            if isinstance(k, int) and v.get('line_type') == 'active'
        ]

        if not active_conductors:
            return 0, 0

        # 2. Count the number of cables (N) by finding unique center points
        num_cables = len(set([cond['center_point'] for cond in active_conductors]))

        # 3. Count conductors per cable (M)
        total_active_conductors = len(active_conductors)
        
        if num_cables > 0:
            conductors_per_cable = total_active_conductors // num_cables
        else:
            conductors_per_cable = 0

        return num_cables, conductors_per_cable
    
    def apply_mtl_ref_properties(self, context, mtl_ref: dict) -> None:
        # Number of single core cables (N) and conductors per cable (M)
        context.num_sc_cables, context.num_conductors_per_scc = self._count_scc_and_conductors(mtl_ref)

    def apply_mtl_properties(self, context, mtl: dict) -> None:
        """
        Calculates and applies overhead-line-specific distance matrices to the MTL object
        by delegating the calculation to a static helper method.
        """
        # 1. Delegate the complex calculation to the static method
        properties = MTLStrategy._surface_distance_with_ground_return(mtl)
        context.d_matrix_ground_return = properties['d_matrix_ground_return']
        context.D_matrix_ground_return = properties['D_matrix_ground_return']
        context.vertical_separation_matrix = properties['vertical_separation_matrix']
        context.horizontal_separation_matrix = properties['horizontal_separation_matrix']

        # Conductors Permeability [np.array]
        context.mu = np.array([sc.mu_0 * conductor['relative_permeability'] for conductor in mtl.values()]) 

        # Conductors Permittivity [np.array]
        context.epsilon = np.array([sc.epsilon_0 * conductor['relative_permittivity'] for conductor in mtl.values()]) 

        # Conductors conductivity [np.array]
        context.sigma = np.array([conductor['conductivity'] for conductor in mtl.values()])

        # Free-Space Permittivity [np.array]
        context.epsilon_out = np.array([sc.epsilon_0 * conductor['relative_permittivity_out'] for conductor in mtl.values()])

def mtl_strategy_factory(mtl_type: str) -> MTLStrategy:
    """Factory function to select the appropriate MTL strategy."""
    strategies = {
        'coated_wires': CableStrategy,
        'bare_wires': CableStrategy,
        'coaxial': CableStrategy,
        'overhead': OverheadLineStrategy,
        'scc': SingleCoreCableStrategy,
        'pipe': CableStrategy,
        'hdpe': SingleCoreCableInHDPEStrategy,
    }
    
    strategy_class = strategies.get(mtl_type)
    if not strategy_class:
        raise ValueError(f"Unknown or unsupported MTL type: {mtl_type}")
        
    return strategy_class()
