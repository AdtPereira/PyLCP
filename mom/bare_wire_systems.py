import numpy as np
import scipy.integrate as spi
import scipy.constants as spc

from utils.case_utils import *
from mtl_main.source import MulticonductorTransmissionLine

@staticmethod
def galerkin_integrand(phi_p, p_idx, q_idx, test_func_idx, basis_func_idx, self_instance):
    """
    Static method that computes the integrand of the Galerkin Method.

    Computes f_pa(phi_p) * g_qb(phi_p), where 'f' is the test function and 'g' is the
    potential of the basis function.

    Args:
        phi_p (float): Angle on the observation conductor 'p' (integration variable).
        p_idx (int): Index of the observation surface.
        q_idx (int): Index of the source surface.
        test_func_idx (int): Index of the test function 'a' on conductor 'p'.
        basis_func_idx (int): Index of the basis function 'b' on conductor 'q'.
        self_instance (object): The class instance used to access the model data.

    Returns:
        float: The value of the integrand.
    """
    # Get surface data
    field_surface = self_instance.model.surfaces[p_idx]
    source_surface = self_instance.model.surfaces[q_idx]
    epsilon = self_instance.model.epsilon_out[source_surface['tag']]

    # --- 1. Compute the value of the test function f_pa(phi_p) ---
    is_test_cos = (test_func_idx % 2 != 0)
    k_test = (test_func_idx + 1) // 2 if is_test_cos else test_func_idx // 2

    if k_test == 0:
        f_pa = 1.0
    elif is_test_cos:
        f_pa = np.cos(k_test * phi_p)
    else:
        f_pa = np.sin(k_test * phi_p)

    # --- 2. Compute the potential g_qb(phi_p) ---
    # Coordinates of the observation point on conductor 'p'
    obs_point = np.array(field_surface['center_point']) + \
                field_surface['radius'] * np.array([np.cos(phi_p), np.sin(phi_p)])

    # Vector from the center of source 'q' to the observation point
    rho_b_vector = obs_point - np.array(source_surface['center_point'])
    rho_b = np.linalg.norm(rho_b_vector)
    theta_b = np.arctan2(rho_b_vector[1], rho_b_vector[0])

    is_basis_cos = (basis_func_idx % 2 != 0)
    k_basis = (basis_func_idx + 1) // 2 if is_basis_cos else basis_func_idx // 2

    g_qb = 0.0
    if k_basis == 0:  # Constant term of the source
        if p_idx == q_idx:
            g_qb = (-source_surface['radius'] / epsilon) * np.log(field_surface['radius'])
        else:
            g_qb = (-source_surface['radius'] / epsilon) * np.log(rho_b)
    else:  # Harmonic terms of the source
        if p_idx == q_idx:
            term = np.cos(k_basis * phi_p) if is_basis_cos else np.sin(k_basis * phi_p)
            g_qb = (source_surface['radius'] / (2 * k_basis * epsilon)) * term
        else:
            term = np.cos(k_basis * theta_b) if is_basis_cos else np.sin(k_basis * theta_b)
            g_qb = (source_surface['radius'] / (2 * k_basis * epsilon)) * ((source_surface['radius'] / rho_b)**k_basis) * term

    return f_pa * g_qb

@staticmethod
def maxwellian_capacitance(model, mom_data):
    """
    Computes the physical (Maxwellian) capacitance matrix of dimension (N-1)x(N-1)
    from the generalized capacitance matrix of dimension NxN.

    This process happens in two steps:
    1.  First, a full Maxwellian matrix (NxN) is computed using
        Equation 5.21, which is given by:
        C_full_ij = c_ij - (row_sum_i * col_sum_j) / total_sum
    2.  Then, the matrix is reduced to (N-1)x(N-1) by removing the row and the
        column corresponding to the reference conductor, whose index is
        specified by the class attribute `self.idx_ref`.
    """
    cgen = mom_data['generalized_capacitance']
    idx_ref = model.mtl_idx_ref

    # --- Validations ---
    assert isinstance(cgen, np.ndarray), "The generalized capacitance matrix must be a NumPy array."
    assert cgen.ndim == 2 and cgen.shape[0] == cgen.shape[1], "The generalized capacitance matrix must be square."
    assert cgen.shape[0] > 1, "Computing the Maxwellian capacitance requires at least 2 conductors."
    assert 0 <= idx_ref < cgen.shape[0], f"The reference index self.idx_ref ({idx_ref}) is outside the valid range [0, {cgen.shape[0]-1}]."

    # --- Step 1: Compute the full Maxwellian matrix (NxN) ---
    total_sum = np.sum(cgen)

    # Avoid division by zero
    assert np.abs(total_sum) > 1e-15, "The sum of the elements of the generalized capacitance matrix is zero, resulting in division by zero."

    correction_matrix = np.outer(np.sum(cgen, axis=1), np.sum(cgen, axis=0)) / total_sum
    C_full = cgen - correction_matrix

    # --- Step 2: Reduce the matrix to (N-1)x(N-1) ---
    # Uses np.delete to remove the row (axis=0) and the column (axis=1)
    # corresponding to the reference conductor index `self.idx_ref`.
    return np.delete(np.delete(C_full, idx_ref, axis=0), idx_ref, axis=1)

class BareWireMoMSolver:
    """
    Computes the capacitance and charge distribution for bare-wire systems
    using the Method of Moments (MoM) with a harmonic series expansion.

    This class uses an instance of MulticonductorTransmissionLine (MTL)
    to obtain the geometric and electrical parameters of the system. It then
    runs the full MoM simulation, populating its own result
    attributes.

    It runs the full MoM simulation, populating every result
    attribute. The construction of the D matrix now includes the constant,
    cosine and sine expansion terms, following expressions (20a), (20b)
    and (20c) of Clements (1975).

    In this class, the maximum harmonic order is defined by 'k',
    while NF (number of coefficients) is derived as 2*k + 1.
    """
    def __init__(self, model: MulticonductorTransmissionLine):
        # MTL Geometry Model
        self.model = model

        # Number of Fourier harmonic coefficients per conductor
        self.NF = [2*surface['fourier_order']+1 for surface in self.model.surfaces][0]

        # Result attributes
        self.mom_data = {'collocation': {}, 'galerkin': {}}

        self.DR_ratio = (model.D_pq[0, 1]) / (model.surfaces[0]['radius'])
        assert self.DR_ratio > 2, "The D/R ratio must be greater than 2 to ensure convergence of the solution."
        self.C_exact_bare_wires = np.pi * spc.epsilon_0 / np.arccosh(0.5 * self.DR_ratio)

    def _collocation_points(self):
        """
        Computes and stores the collocation points, classifying them in a dictionary
        nested by the conductor 'tag' and by the surface type ('conductor', 'sheath').
        """
        # Initialize the main dictionary that will be the class attribute.
        collocation_data = {}

        # Equation (A.4b): Rotation angle for the set of points.
        delta = np.pi / (2 * self.NF)

        # Compute the base angles, which are rotated by delta to obtain
        # the angles of the observation points (match points).
        base_angles = np.linspace(0, 2 * np.pi, self.NF, endpoint=False)
        match_angles = base_angles + delta

        # Iterate over every surface defined in the base MTL class.
        for surface in self.model.surfaces:
            if surface['tag'] not in collocation_data:
                collocation_data[surface['tag']] = {}

            # Compute the Cartesian coordinates for the source and observation points.
            match_points = np.array(surface['center_point']) + surface['radius'] * np.array([np.cos(match_angles), np.sin(match_angles)]).T

            # Populate the dictionary for the specific surface with its data.
            collocation_data[surface['tag']][surface['type']] = {
                'observation': {
                    'cartesian': match_points,
                    'angles_rad': match_angles
                }
            }

        self.mom_data['collocation']['data'] = collocation_data

    def _generalized_capacitance_clements(self):
        """
        Computes the generalized capacitance matrix C from the T matrix (D^-1).
        """
        moment_matrix = self.mom_data['collocation']['moment_matrix']
        T_matrix = np.linalg.inv(moment_matrix)
        C_matrix = np.zeros((2, 2))

        for n in range(2):      # Index of the charge conductor
            r_i = self.model.surfaces[n]['radius']
            for m in range(2):  # Index of the potential conductor
                sum_of_T_elements = np.sum(T_matrix[(n * self.NF), (m * self.NF):((m + 1) * self.NF)])
                C_matrix[n, m] = 2 * np.pi * r_i * sum_of_T_elements

        self.mom_data['collocation']['generalized_capacitance'] = C_matrix

    def _generalized_capacitance_savage(self):
        """
        Computes the generalized capacitance matrix C from the T matrix (D^-1)
        following the Savage (1993) formulation for the Galerkin Method.

        The formulation is given by: C_ij = (2*pi)^2 * r_i * T_ij[0,0], where T_ij[0,0]
        is the top-left element of the corresponding submatrix of the inverse matrix T.
        """
        # 1. Invert the D matrix to obtain the T matrix
        moment_matrix = self.mom_data['galerkin']['moment_matrix']
        T_matrix = np.linalg.inv(moment_matrix)

        # 2. Get the number of conductors (surfaces)
        num_conductors = len(self.model.surfaces)
        C_matrix = np.zeros((num_conductors, num_conductors))

        # 3. Iterate over each element of the capacitance matrix to be computed
        for i in range(num_conductors):      # Index 'i' for the charge conductor (row)
            for j in range(num_conductors):  # Index 'j' for the potential conductor (column)

                # 4. Get the radius of the charge conductor 'i'
                # This is more robust than using self.R, since it accounts for different radii.
                r_i = self.model.surfaces[i]['radius']

                # 5. Locate element (0,0) of the T_ij submatrix
                # This is the top-left element of the block relating the observation
                # on conductor 'i' with the source on conductor 'j'.
                T_ij_00 = T_matrix[i * self.NF, j * self.NF]

                # 6. Compute the capacitance element C_ij according to Equation 4.37
                C_matrix[i, j] = (2 * np.pi)**2 * r_i * T_ij_00

        self.mom_data['galerkin']['generalized_capacitance'] = C_matrix

    def run_collocation_method(self):
        """
        Runs the full MoM simulation, assembling the system of equations for
        every surface (conducting and dielectric) based on the new
        data structures.
        """
        self._collocation_points()
        collocation_data = self.mom_data['collocation']['data']
        moment_matrix = np.zeros((self.model.N, self.model.N))
        V_vector = np.zeros(self.model.N)

        # 1. Prepare the system indices and vectors
        # Pre-compute the number of coefficients (NF) for each surface
        nfs_per_surface = [2 * surface['fourier_order'] + 1 for surface in self.model.surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        # 2. Assemble the [D] Matrix and the [V] Vector
        # Loop over the OBSERVATION surfaces p (rows of the matrix)
        for p, field_surface in enumerate(self.model.surfaces):
            tag_p = field_surface['tag']
            type_p = field_surface['type']
            radius_p = field_surface['radius']

            # Get the observation points for surface p
            match_points = collocation_data[tag_p][type_p]['observation']['cartesian']

            # Fill the potential vector V for the row block of surface p
            if type_p == 'conductor':
                V_vector[offsets[p] : offsets[p] + nfs_per_surface[p]] = self.model.mtl[tag_p]['potential_to_infinity']

            # The boundary condition on the dielectric sheath results in 0 on the right-hand side of the equation
            elif type_p == 'primary_insulation':
                V_vector[offsets[p] : offsets[p] + nfs_per_surface[p]] = 0.0

            # Loop over the SOURCE surfaces q (columns of the matrix)
            for q, source_surface in enumerate(self.model.surfaces):
                radius_q = source_surface['radius']
                epsilon = self.model.epsilon_out[source_surface['tag']]

                # Loop over each observation point m on surface p
                for m in range(nfs_per_surface[p]):
                    row_idx = offsets[p] + m

                    # Angle of the observation point relative to the center of its OWN surface
                    rho_i_vector = match_points[m] - np.array(field_surface['center_point'])
                    theta_i = np.arctan2(rho_i_vector[1], rho_i_vector[0])

                    # Loop over each basis function n on surface q
                    for n in range(nfs_per_surface[q]):
                        col_idx = offsets[q] + n

                        # Local harmonic index of the source
                        harmonic_idx = n
                        is_cosine_term = (harmonic_idx % 2 != 0)
                        k = (harmonic_idx + 1) // 2 if is_cosine_term else harmonic_idx // 2

                        # Source angle and vector 'b' relative to the center of the SOURCE surface 'q'
                        rho_b_vector = match_points[m] - np.array(source_surface['center_point'])
                        rho_b = np.linalg.norm(rho_b_vector)
                        theta_b = np.arctan2(rho_b_vector[1], rho_b_vector[0])

                        # ========================================================================
                        # ==== START OF THE D MATRIX ELEMENT COMPUTATION LOGIC ==================
                        # ========================================================================

                        # === BLOCK 1: POTENTIAL COMPUTATION (phi) ==============================
                        # === Applies the boundary condition V = Vm on the conducting surfaces. =

                        # Self-interaction (Observer ON the source boundary)
                        # Constant term (k=0)
                        if harmonic_idx == 0:
                            if p == q:
                                moment_matrix[row_idx, col_idx] = (-radius_q / epsilon) * np.log(radius_p)

                            # Mutual interaction
                            else:
                                moment_matrix[row_idx, col_idx] = (-radius_q / epsilon) * np.log(rho_b)

                        # Harmonic terms (k>0)
                        else:
                            # Self-interaction
                            if p == q:
                                term = np.cos(k * theta_i) if is_cosine_term else np.sin(k * theta_i)
                                moment_matrix[row_idx, col_idx] = (radius_q / (2 * k * epsilon)) * term

                            # Mutual interaction
                            else:
                                term = np.cos(k * theta_b) if is_cosine_term else np.sin(k * theta_b)
                                moment_matrix[row_idx, col_idx] = (radius_q / (2 * k * epsilon)) * ((radius_q / rho_b)**k) * term

                        # ========================================================================
                        # ==== END OF THE D MATRIX ELEMENT COMPUTATION LOGIC ===================
                        # ========================================================================

        # 4. Store the results in the mom_data dictionary
        self.mom_data['collocation']['V_vector'] = V_vector
        self.mom_data['collocation']['moment_matrix'] = moment_matrix
        self.mom_data['collocation']['sigma_coeffs'] = np.linalg.solve(moment_matrix, V_vector)

        self._generalized_capacitance_clements()
        cap_matrix = maxwellian_capacitance(self.model, self.mom_data['collocation'])
        self.mom_data['collocation']['maxwellian_capacitance'] = cap_matrix

    def run_galerkin_method(self):
        """
        Runs the full MoM simulation using the Galerkin Method.
        """
        moment_matrix = np.zeros((self.model.N, self.model.N))
        V_vector = np.zeros(self.model.N)

        nfs_per_surface = [2 * surface['fourier_order'] + 1 for surface in self.model.surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        # Loop over the OBSERVATION surfaces p (rows of the matrix)
        for p, field_surface in enumerate(self.model.surfaces):
            # Loop over the SOURCE surfaces q (columns of the matrix)
            for q, source_surface in enumerate(self.model.surfaces):

                # Loop over the TEST FUNCTIONS 'm' on surface 'p'
                for m in range(nfs_per_surface[p]):
                    row_idx = offsets[p] + m

                    # Loop over the BASIS FUNCTIONS 'n' on surface 'q'
                    for n in range(nfs_per_surface[q]):
                        col_idx = offsets[q] + n

                        # --- Numerical Integration with scipy.integrate.quad ---
                        integral_value, _ = spi.quad(
                            galerkin_integrand, 0, 2 * np.pi,
                            args=(p, q, m, n, self),
                            limit=350
                        )
                        moment_matrix[row_idx, col_idx] = integral_value

            # Filling the V Vector according to the Galerkin formulation
            if field_surface['type'] == 'conductor':
                potential = self.model.mtl[field_surface['tag']]['potential_to_infinity']
                # Only the constant term (m=0) of the right-hand side integral is non-zero
                V_vector[offsets[p]] = 2 * np.pi * potential

        # Solve the system and get the results
        self.mom_data['galerkin']['V_vector'] = V_vector
        self.mom_data['galerkin']['moment_matrix'] = moment_matrix
        self.mom_data['galerkin']['sigma_coeffs'] = np.linalg.solve(moment_matrix, V_vector)

        self._generalized_capacitance_savage()
        cap_matrix = maxwellian_capacitance(self.model, self.mom_data['galerkin'])
        self.mom_data['galerkin']['maxwellian_capacitance'] = cap_matrix
