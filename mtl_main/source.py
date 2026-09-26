""" This module defines the MulticonductorTransmissionLine class.

REFERENCES:
[1] PATEL, Utkarsh R. A Surface Admittance Approach For Fast Calculation of the 
    Series Impedance of Cables Including Skin, Proximity, and Ground Return Effects.
    2014. University of Toronto, Graduate Department of The Edward S. Rogers Sr. 
    Department of Electrical & Computer Engineering. 

[2] U. R. Patel, B. Gustavsen and P. Triverio, "An Equivalent Surface Current Approach
    for the Computation of the Series Impedance of Power Cables with Inclusion of Skin
    and Proximity Effects," in IEEE Transactions on Power Delivery, vol. 28, no. 4, pp.
    2474-2482, Oct. 2013, doi: 10.1109/TPWRD.2013.2267098.

[3] U. R. Patel, B. Gustavsen and P. Triverio, "Application of the MoM-SO Method for 
    Accurate Impedance Calculation of Single-Core Cables Enclosed by a Conducting Pipe," 
    Proc. International Conference on Power Systems Transients (IPST 2013), Vancouver, 
    Canada July 18-20, 2013. https://www.ipstconf.org/Proc_IPST2013.php

[4] PAUL, Clayton R. Analysis of multiconductor transmission lines. 2. ed. Hoboken,
    N.J.: John Wiley & Sons, Inc., c2008.

"""

import copy
import numpy as np
from pyparsing import Dict
import scipy.constants as sc
from mtl_main.strategy import mtl_strategy_factory

class MulticonductorTransmissionLine:
    def __init__(self, model: Dict):
        """Initialize the MulticonductorTransmissionLine class.

        Key Points
        Deep Copy: When initializing the MulticonductorTransmissionLine class, 
        we use copy.deepcopy to ensure that self.conductors is an independent 
        copy of TREFOIL.

        Isolation of Modifications: Any modifications made to self.conductors 
        within the class will not affect the original TREFOIL object.
        """

        mtl_input = copy.deepcopy(model)
        self.mtl_type = mtl_input.get('type', 'unknown')
        self.mtl_frequency = mtl_input.get('frequency_driver', {})  
        self.mtl_idx_ref = mtl_input.get('idx_ref_conductor', 0)

        # MULTICONDUCTOR TRANSMISSION LINES PARAMETERS
        self.mu = None
        self.sigma = None
        self.epsilon = None    
        self.epsilon_out = None    
        self.num_sc_cables = None
        self.num_conductors_per_scc = None
        self.d_matrix_ground_return = None
        self.D_matrix_ground_return = None
        self.images_vertical_distance_matrix = None
        self.horizontal_separation_matrix = None
        self.ground_return_block_sizes = None   # conductors per ground-return object (ECC strategies)

        # Delegate preprocessing and validation to the strategy
        strategy = mtl_strategy_factory(self.mtl_type)
        self.mtl, self.mtl_ref = strategy.preprocess_mtl_data(mtl_input)
        strategy.validate(self.mtl, self.mtl_ref, mtl_input)

        # Delegate calculation of type-specific properties
        strategy.apply_mtl_ref_properties(self, mtl_input)
        strategy.apply_mtl_properties(self, self.mtl)

        # Conductor surfaces dictionary
        self._define_mtl_surfaces()

        # Number of Fourier coefficients List  
        self.NF_List = [2 * surface['fourier_order'] + 1 for surface in self.surfaces]
        
        # Dimension N - Equation (2.36) [1]
        self.N = sum(self.NF_List)

        # Number of Fourier harmonic coefficients per conductor
        self.NF = self.NF_List[self.mtl_idx_ref]

    # Method to generate the PUL data structure
    def get_pul_data_structure(self):
        """
        Generates the base data structure for the per-unit-length (PUL) parameters,
        using the frequency data from the model's 'frequency_driver'.

        The frequency resolution for the 'analytical' method is set to 4x
        the 'mom_so' resolution (steps_per_decade).

        Returns:
            dict: A structured dictionary to store PUL data.
        """

        freq = self.mtl_frequency

        # Try to generate frequencies dynamically from the model
        if freq and freq.get('spacing') == 'log':
            min_hz = freq.get('min_Hz', 1.0)
            max_hz = freq.get('max_Hz', 1e6)

            # Steps per decade for mom_so (numerical)
            steps_per_decade = freq.get('steps_per_decade', 5) # Default 5

            # UPDATED: steps per decade for analytical (4x mom_so)
            steps_per_decade_analytical = steps_per_decade * 4

            start_log = np.log10(min_hz)
            stop_log = np.log10(max_hz)
            num_decades = stop_log - start_log

            # Compute num_points for mom_so
            num_mom_so = int((num_decades * steps_per_decade) + 1)

            # Compute num_points for analytical
            num_analytical = int((num_decades * steps_per_decade_analytical) + 1)

            freq = np.logspace(start_log, stop_log, num=num_mom_so)
            freq_analytical = np.logspace(start_log, stop_log, num=num_analytical)

        else:
            # Fallback to static values (in case frequency_driver fails)
            # (Default: 5 steps/decade for mom_so -> 31 points)
            freq = np.logspace(0, 6, num=31)
            # (Default: 20 steps/decade for analytical (4*5) -> 121 points)
            freq_analytical = np.logspace(0, 6, num=121)

        pul_data = {
            'comsol': {
                'frequencies': freq,
                'scenarios': {
                    '1': {},
                },
            },
            'mom_so': {
                'frequencies': freq,
                'scenarios': {
                    '1': {'mtl': self},  # 'self' is the mtl instance
                },
            },
            'analytical': {
                'frequencies': freq_analytical,
                'scenarios': {
                    '1': {'mtl': self},  # 'self' is the mtl instance
                },
            },
        }
        return pul_data
    
    # Helper method to define MTL surfaces
    def _define_mtl_surfaces(self):
        """
        Defines the surfaces for all conductors and their insulations,
        handling geometric duplicates by prioritizing conductor surfaces.
        """
        # Step 1: Generate a temporary list of all potential surfaces, including duplicates.
        self.surfaces = []
        all_surfaces = []
        for key, conductor in self.mtl.items():
            # Check if the conductor is hollow
            if conductor['radius'][0] != 0:
                for radius in conductor['radius']:
                    all_surfaces.append({
                        'type': 'conductor',
                        'tag': key,
                        'radius': radius,
                        'center_point': conductor['center_point'],
                        'fourier_order': conductor['fourier_order'],
                    })
            
            # Then, the conductor is solid
            else:
                all_surfaces.append({
                    'type': 'conductor',
                    'tag': key,
                    'radius': conductor['radius'][1],
                    'center_point': conductor['center_point'],
                    'fourier_order': conductor['fourier_order'],
                })

            # Check if the conductor has insulation
            if conductor['insulation'] is not None:
                insulation = conductor['insulation']
                # Add the insulation surface
                all_surfaces.append({
                    'type': 'primary_insulation',
                    'tag': key,
                    'radius': conductor['radius'][1] + insulation['thickness'],
                    'center_point': insulation['center_point'],
                    'fourier_order': insulation['fourier_order'],
                    'relative_permittivity': insulation['relative_permittivity'],
                })

        # Step 2: Use a dictionary to identify unique surfaces based on geometry,
        # applying the preference rule for conductors.
        unique_surfaces_map = {}
        for surface in all_surfaces:
            # A unique key is defined by the center point and radius.
            key = (surface['center_point'], surface['radius'])
            
            # If this geometry is new, or if the new surface is a conductor
            # and the existing one is not, we add/update the map.
            if key not in unique_surfaces_map or \
               (surface['type'] == 'conductor' and unique_surfaces_map[key]['type'] != 'conductor'):
                unique_surfaces_map[key] = surface

        # Step 3: Reconstruct the final list, preserving the original order of appearance
        # of the unique geometries.
        added_geometries = set()
        for surface in all_surfaces:
            key = (surface['center_point'], surface['radius'])
            if key not in added_geometries:
                # Append the definitive version of the surface from our map
                self.surfaces.append(unique_surfaces_map[key])
                added_geometries.add(key)

    # Contour Position Vector [np.array] - Equation (2.2) [1]
    def contour_vector_position(self, surfaces_list, theta, p):
        """ The boundary cp can be traced by the position vector rp(ap, θp) """

        xp, yp = surfaces_list[p]['center']
        ap = surfaces_list[p]['radius']

        return np.array([xp + ap * np.cos(theta), yp + ap * np.sin(theta)])

    # Create a diagonal matrix by concatenating a list of matrices along the diagonal
    def create_diagonal_matrix(self, matrix_list):
        """
        Create a diagonal matrix by concatenating a list of matrices along the diagonal.

        Parameters:
        matrix_list (list): A list of numpy arrays representing matrices.

        Returns:
        numpy.ndarray: A diagonal matrix created by concatenating the matrices along the diagonal.
        """
        # Check if all elements in the list are square matrices
        for i, matrix in enumerate(matrix_list):
            if matrix.shape[0] != matrix.shape[1]:
                raise ValueError(f"Matrix at index {i} is not square: shape {matrix.shape}")

        # Determine the size of the resulting matrix Y
        size_y = sum(mat.shape[0] for mat in matrix_list)

        # Create a zero matrix of appropriate size
        matrix_y = np.zeros((size_y, size_y), dtype=complex)

        # Variables to keep track of the position in the matrix Y
        row_start = 0
        col_start = 0

        for matrix in matrix_list:
            matrix_size = matrix.shape[0]

            # Insert the matrix into the diagonal block of Y
            matrix_y[row_start:row_start + matrix_size,
                     col_start:col_start + matrix_size] = matrix

            # Update the position variables
            row_start += matrix_size
            col_start += matrix_size

        return matrix_y

    # Create a block matrix from a list of sub-matrices
    def create_block_matrix(self, matrix_list, num_rows, num_cols):
        """
        Create a block matrix from a list of sub-matrices.

        Args:
            matrix_list (list): A list of numpy arrays representing the sub-matrices.
            numRows (int): The number of rows in the block matrix.
            numColumns (int): The number of columns in the block matrix.

        Returns:
            numpy.ndarray: The block matrix created from the sub-matrices.

        Raises:
            ValueError: If the number of sub-matrices does not match numRows * numColumns.

        """

        # Verify the total number of sub-matrices matches the product of numRows and numColumns
        n = len(matrix_list)
        if num_rows * num_cols != n:
            raise ValueError(
                "The number of sub-matrices must match numRows * numColumns")

        # Determine the row and column sizes for the final matrix
        row_sizes = [matrix_list[i * num_cols].shape[0]
                     for i in range(num_rows)]
        col_sizes = [matrix_list[i].shape[1] for i in range(num_cols)]

        # Create the final matrix filled with zeros
        matrix = np.zeros((sum(row_sizes), sum(col_sizes)),
                          dtype=matrix_list[0].dtype)

        # Fill the final matrix with sub-matrices
        current_row = 0
        for i in range(num_rows):
            current_col = 0
            for j in range(num_cols):
                sub_matrix = matrix_list[i * num_cols + j]
                rows, cols = sub_matrix.shape
                matrix[current_row:current_row + rows,
                       current_col:current_col + cols] = sub_matrix
                current_col += cols
            current_row += row_sizes[i]

        return matrix
