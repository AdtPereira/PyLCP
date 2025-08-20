from abc import ABC, abstractmethod
import numpy as np

class MTLStrategy(ABC):
    """
    Abstract Base Class defining the interface for MTL type-specific logic.
    """

    @abstractmethod
    def preprocess_mtl_dict(self, mtl_data: dict) -> dict:
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

    # def calculate_properties(self, mtl_data: dict) -> dict:
    #     """
    #     Calculates type-specific properties. For instance, distance matrices
    #     might not apply to all configurations in the same way.
    #     Returns a dictionary of calculated properties.
    #     """
    #     # Default implementation returns no specific properties
    #     return {}

    def apply_properties(self, context, mtl_data: dict) -> None:
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
            cp = conductor_p['center_point']
            for q, conductor_q in mtl.items():
                cq = conductor_q['center_point']

                # Distance between the center points of the conductors
                x_pq[p, q] = cp[0] - cq[0]
                y_pq[p, q] = cp[1] - cq[1]

                # Angle between the center points of the conductors
                theta_pq[p, q] = np.arctan2(y_pq[p, q], x_pq[p, q])

                # Distance between the center points of the conductors
                d_pq[p, q] = np.linalg.norm(np.array(cp) - np.array(cq))

        return d_pq, x_pq, y_pq, theta_pq
    

    @staticmethod
    def _distances_with_ground_return(mtl_data: dict) -> dict:
        """
        Calculates all necessary distance matrices for overhead line analysis,
        including geometric, ground return, horizontal, and vertical separations.
        Returns a dictionary of the calculated matrices.
        """
        N = len(mtl_data)

        # Initialize all four matrices
        d_matrix = np.zeros((N, N))
        D_matrix = np.zeros((N, N))
        vertical_separation_matrix = np.zeros((N, N))
        horizontal_separation_matrix = np.zeros((N, N))

        # Create a sorted list of conductors to ensure consistent ordering
        conductors = sorted(mtl_data.items())

        for n_idx, (n_tag, n_cond_data) in enumerate(conductors):
            for m_idx, (m_tag, m_cond_data) in enumerate(conductors):
                cn = n_cond_data['center_point']
                cm = m_cond_data['center_point']

                # Horizontal and Vertical separation
                dn_dm = cn[0] - cm[0]
                hn_hm = cn[1] + cm[1]

                # Geometric Distance (d_nm)
                if n_tag == m_tag:
                    d = n_cond_data['radius'][1] # Use radius for self-distance
                else:
                    d = np.sqrt(dn_dm ** 2 + (cn[1] - cm[1]) ** 2)

                # Distance to Image (D_nm)
                D = np.sqrt(dn_dm ** 2 + hn_hm ** 2)

                # Populate the matrices at the correct indices
                d_matrix[n_idx, m_idx] = d
                D_matrix[n_idx, m_idx] = D
                vertical_separation_matrix[n_idx, m_idx] = hn_hm
                horizontal_separation_matrix[n_idx, m_idx] = dn_dm

        return {
            'd_matrix_ground_return': d_matrix,
            'D_matrix_ground_return': D_matrix,
            'vertical_separation_matrix': vertical_separation_matrix,
            'horizontal_separation_matrix': horizontal_separation_matrix
        }

class CableStrategy(MTLStrategy):
    """Strategy for cable-based MTLs like 'coaxial', 'coated_wires', etc."""

    def preprocess_mtl_dict(self, mtl_data: dict) -> dict:
        # For cables, we expect integer keys starting from 0.
        return {key: value for key, value in mtl_data.items() if isinstance(key, int)}

    def validate(self, mtl_data: dict, idx_ref: int):
        key_conductors = sorted(mtl_data.keys())
        key_expected = list(range(len(key_conductors)))
        
        assert len(mtl_data) > 1, "Cable-based MTL must have at least two conductors."
        assert key_conductors == key_expected, "Cable conductor tags must be a sequence starting from 0."
        assert mtl_data[0]['line_type'] == 'return', "Conductor with index '0' must be the return path for cables."
        assert idx_ref in mtl_data, f"The reference conductor index {idx_ref} must be in the MTL dictionary."

    def apply_properties(self, context, mtl_data: dict) -> None:
        """Applies cable-specific distance properties to the MTL object."""
        d_pq, dx_pq, dy_pq, theta_pq = MTLStrategy._conductors_center_distance_matrix(mtl_data)
        
        context.D_pq = d_pq
        context.dx_pq = dx_pq
        context.dy_pq = dy_pq
        context.theta_pq = theta_pq

class OverheadStrategy(MTLStrategy):
    """Strategy for overhead lines."""

    def preprocess_mtl_dict(self, mtl_data: dict) -> dict:
        # For overhead lines, conductors are indexed starting from 1.
        return {key: value for key, value in mtl_data.items() if isinstance(key, int) and key > 0}

    def validate(self, mtl_data: dict, idx_ref: int):
        key_conductors = sorted(mtl_data.keys())
        key_expected = list(range(1, len(key_conductors) + 1))
        
        assert len(mtl_data) >= 1, "Overhead MTL must have at least one conductor."
        assert key_conductors == key_expected, "Overhead conductor tags must be a sequence starting from 1."
        # Overhead lines might not have a specific return conductor at index 0
        assert idx_ref is not None, "The reference conductor index must be defined."

    # Matrices Distance dnm and Dnm [np.array] - Equations 2.48 and 2.49 [2]
    # def _calculate_distances_with_ground_return(self, mtl, n, m):
    #     """ 
    #     Calculates distances between conductors n and m, including ground return effect.
    #     Originalmente 'conductor_distances' em mtl.py.
    #     """
    #     cn = mtl[n]['center_point']
    #     cm = mtl[m]['center_point']
    #     rn = mtl[n]['radius'][1]

    #     # Horizontal separation
    #     if n == m:
    #         dn_dm = rn
    #     else:
    #         dn_dm = cn[0] - cm[0]

    #     # Vertical separation (for ground return image)
    #     hn_hm = cn[1] + cm[1]

    #     # Geometric Distance
    #     d = np.sqrt(dn_dm ** 2 + (cn[1] - cm[1]) ** 2)

    #     # Distance to Image
    #     D = np.sqrt(dn_dm ** 2 + hn_hm ** 2)

    #     return d, D, hn_hm, dn_dm
    
    # def calculate_properties(self, mtl_data: dict) -> dict:
    #     """
    #     Pre-calculates all necessary distance matrices for overhead line analysis,
    #     including geometric, ground return, horizontal, and vertical separations.
    #     """
    #     # The number of conductors is derived from the length of the mtl_data dictionary
    #     N = len(mtl_data)

    #     d_matrix_ground_return = np.zeros((N, N))       # Geometric Distance matrix
    #     D_matrix_ground_return = np.zeros((N, N))       # Distance to Image matrix
    #     vertical_separation_matrix = np.zeros((N, N))   # Vertical distance matrix
    #     horizontal_separation_matrix = np.zeros((N, N)) # Horizontal distance matrix

    #     # The keys in mtl_data for overhead lines start from 1, so we need a mapping
    #     # from conductor tag (e.g., 1, 2, 3) to matrix index (e.g., 0, 1, 2)
    #     # Assuming tags are sequential, e.g., 1, 2, ..., N

    #     # Create a sorted list of conductors to ensure consistent ordering
    #     conductors = sorted(mtl_data.items())

    #     # A lógica de 'conductor_distances' agora vive aqui, dentro de um laço
    #     for n_idx, (n_tag, n_cond_data) in enumerate(conductors):
    #         for m_idx, (m_tag, m_cond_data) in enumerate(conductors):
    #             cn = n_cond_data['center_point']
    #             cm = m_cond_data['center_point']

    #             # --- The core logic from the old conductor_distances method ---

    #             # Horizontal and Vertical separation
    #             dn_dm = cn[0] - cm[0]
    #             hn_hm = cn[1] + cm[1]

    #             # Geometric Distance
    #             if n_tag == m_tag:
    #                 d = n_cond_data['radius'][1]
    #             else:
    #                 d = np.sqrt(dn_dm ** 2 + (cn[1] - cm[1]) ** 2)

    #             # Populate the matrices at the correct indices
    #             d_matrix_ground_return[n_idx, m_idx] = d
    #             D_matrix_ground_return[n_idx, m_idx] = np.sqrt(dn_dm ** 2 + hn_hm ** 2)
    #             vertical_separation_matrix[n_idx, m_idx] = hn_hm
    #             horizontal_separation_matrix[n_idx, m_idx] = dn_dm

    #     return {
    #         'd_matrix_ground_return': d_matrix_ground_return,
    #         'D_matrix_ground_return': D_matrix_ground_return,
    #         'vertical_separation_matrix': vertical_separation_matrix,
    #         'horizontal_separation_matrix': horizontal_separation_matrix
    #     }

    # def apply_properties(self, context, mtl_data: dict) -> None:
    #     """Applies overhead-line-specific distance matrices to the MTL object."""

    #     # The number of conductors is derived from the length of the mtl_data dictionary
    #     N = len(mtl_data)

    #     d_matrix_ground_return = np.zeros((N, N))       # Geometric Distance matrix
    #     D_matrix_ground_return = np.zeros((N, N))       # Distance to Image matrix
    #     vertical_separation_matrix = np.zeros((N, N))   # Vertical distance matrix
    #     horizontal_separation_matrix = np.zeros((N, N)) # Horizontal distance matrix

    #     # The keys in mtl_data for overhead lines start from 1, so we need a mapping
    #     # from conductor tag (e.g., 1, 2, 3) to matrix index (e.g., 0, 1, 2)
    #     # Assuming tags are sequential, e.g., 1, 2, ..., N

    #     # Create a sorted list of conductors to ensure consistent ordering
    #     conductors = sorted(mtl_data.items())

    #     # A lógica de 'conductor_distances' agora vive aqui, dentro de um laço
    #     for n_idx, (n_tag, n_cond_data) in enumerate(conductors):
    #         for m_idx, (m_tag, m_cond_data) in enumerate(conductors):
    #             cn = n_cond_data['center_point']
    #             cm = m_cond_data['center_point']

    #             # --- The core logic from the old conductor_distances method ---

    #             # Horizontal and Vertical separation
    #             dn_dm = cn[0] - cm[0]
    #             hn_hm = cn[1] + cm[1]

    #             # Geometric Distance
    #             if n_tag == m_tag:
    #                 d = n_cond_data['radius'][1]
    #             else:
    #                 d = np.sqrt(dn_dm ** 2 + (cn[1] - cm[1]) ** 2)

    #             # Populate the matrices at the correct indices
    #             d_matrix_ground_return[n_idx, m_idx] = d
    #             D_matrix_ground_return[n_idx, m_idx] = np.sqrt(dn_dm ** 2 + hn_hm ** 2)
    #             vertical_separation_matrix[n_idx, m_idx] = hn_hm
    #             horizontal_separation_matrix[n_idx, m_idx] = dn_dm

    #     context.d_matrix_ground_return = d_matrix_ground_return
    #     context.D_matrix_ground_return = D_matrix_ground_return
    #     context.vertical_separation_matrix = vertical_separation_matrix
    #     context.horizontal_separation_matrix = horizontal_separation_matrix

    def apply_properties(self, context, mtl_data: dict) -> None:
        """
        Calculates and applies overhead-line-specific distance matrices to the MTL object
        by delegating the calculation to a static helper method.
        """
        # 1. Delegate the complex calculation to the static method
        properties = MTLStrategy._distances_with_ground_return(mtl_data)

        # 2. Apply the results to the context object
        context.d_matrix_ground_return = properties['d_matrix_ground_return']
        context.D_matrix_ground_return = properties['D_matrix_ground_return']
        context.vertical_separation_matrix = properties['vertical_separation_matrix']
        context.horizontal_separation_matrix = properties['horizontal_separation_matrix']

def mtl_strategy_factory(mtl_type: str) -> MTLStrategy:
    """Factory function to select the appropriate MTL strategy."""
    strategies = {
        'coated_wires': CableStrategy,
        'bare_wires': CableStrategy,
        'coaxial': CableStrategy,
        'scc': CableStrategy,
        'overhead': OverheadStrategy,
    }
    
    strategy_class = strategies.get(mtl_type)
    if not strategy_class:
        raise ValueError(f"Unknown or unsupported MTL type: {mtl_type}")
        
    return strategy_class()
