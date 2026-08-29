import numpy as np
import pandas as pd
import scipy.constants as sc
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from utils.case_utils import *
from mtl_main.source import MulticonductorTransmissionLine

class MulticonductorCoatedWireSystems:
    """
    Computes the capacitance and charge distribution for bare-wire systems
    using the Method of Moments (MoM) with a harmonic series expansion.

    This class inherits from MulticonductorTransmissionLine (MTL) and specializes it
    for the two-wire case, deriving its parameters from an 'mtl' configuration.

    It runs the full MoM simulation, populating every result attribute.
    The construction of the D matrix now includes the constant, cosine and sine
    expansion terms, following expressions (20a), (20b) and (20c) of Clements (1975).

    In this class, the maximum harmonic order is defined by 'k',
    while NF (number of coefficients) is derived as 2*k + 1.
    """
    def __init__(self, model: MulticonductorTransmissionLine):
        # MTL Geometry Model
        self.model = model

        # Radius of conductor 'p' (first conductor)
        self.R = self.model.surfaces[0]['radius']
        self.D = self.model.D_pq[0, 1]
        self.DR_ratio = self.D / self.R
        assert self.DR_ratio > 2, "The D/R ratio must be greater than 2 to ensure convergence of the solution."
        self.C_exact_bare_wires = np.pi * sc.epsilon_0 / np.arccosh(0.5 * self.DR_ratio)

        # Conductor and insulation surfaces
        self.conductor_surfaces = [s for s in self.model.surfaces if s['type'] == 'conductor']
        self.insulation_surfaces = [s for s in self.model.surfaces if s['type'] == 'primary_insulation']
        self.ordered_surfaces = self.conductor_surfaces + self.insulation_surfaces

        # Result attributes
        self.collocation_data = None
        self.D_matrix = None
        self.T_matrix = None
        self.V_vector = None
        self.sigma_coeffs = None
        self.C_generalized = None
        self.C_maxwellian = None

    def _calculate_collocation_points(self):
        """
        Computes and stores the collocation points, classifying them in a
        dictionary nested by the conductor 'tag' and by the surface type
        ('conductor', 'primary_insulation').
        """
        # Initialize the main dictionary that will be the class attribute.
        self.collocation_data = {}

        # Equation (A.4a): Angular separation between collocation points.
        theta = 2 * np.pi / self.model.NF

        # Equation (A.4b): Rotation angle for the set of points.
        delta = np.pi / (2 * self.model.NF)

        # Compute the base angles, which are rotated by delta to obtain
        # the angles of the observation points (match points).
        base_angles = np.linspace(0, 2 * np.pi, self.model.NF, endpoint=False)
        field_angles = base_angles + delta

        # The source points are placed halfway between the observation
        # points to ensure numerical stability.
        source_angles = field_angles - (theta / 2)

        # Iterate over every surface defined in the base MTL class.
        for surface in self.model.surfaces:
            tag = surface['tag']
            surface_type = surface['type']
            center = np.array(surface['center_point'])
            radius = surface['radius']

            # Create the dictionary for the conductor 'tag' if it does not exist yet.
            if tag not in self.collocation_data:
                self.collocation_data[tag] = {}

            # Compute the Cartesian coordinates for the source and observation points.
            source_points = center + radius * np.array([np.cos(source_angles), np.sin(source_angles)]).T
            field_points = center + radius * np.array([np.cos(field_angles), np.sin(field_angles)]).T

            # Populate the dictionary for the specific surface with its data.
            self.collocation_data[tag][surface_type] = {
                'source': {
                    'cartesian': source_points,
                    'angles_rad': source_angles,
                    'angles_deg': np.degrees(source_angles)
                },
                'observation': {
                    'cartesian': field_points,
                    'angles_rad': field_angles,
                    'angles_deg': np.degrees(field_angles)
                }
            }

    def _calculate_generalized_capacitance(self):
        """
        Computes the generalized capacitance matrix C in a robust way.

        This revised version fixes the flaw of the previous implementation, ensuring
        that the C matrix is built correctly regardless of the values or the
        order of the conductor 'tags'.
        """
        self.T_matrix = np.linalg.inv(self.D_matrix)

        num_conductors = len(self.conductor_surfaces)
        nfs_per_surface = [2 * s['fourier_order'] + 1 for s in self.ordered_surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        # 2. Create a map for easy access to the properties and offsets of each surface
        surface_map = {}
        for i, surface in enumerate(self.ordered_surfaces):
            tag = surface['tag']
            if tag not in surface_map:
                surface_map[tag] = {}
            surface_map[tag][surface['type']] = {
                'radius': surface['radius'],
                'offset': offsets[i],
                'nf': nfs_per_surface[i]
            }

        # --- START OF REVISED LOGIC ---

        # 3. Ensure a consistent ordering for the capacitance matrix
        ### Get a sorted list of the conductor tags. Essential for consistency.
        sorted_conductor_tags = sorted([s['tag'] for s in self.conductor_surfaces])

        ### Create a map from 'tag' to the matrix index (0, 1, 2...).
        tag_to_idx = {tag: i for i, tag in enumerate(sorted_conductor_tags)}

        # 4. Compute the capacitance matrix
        C_matrix = np.zeros((num_conductors, num_conductors))

        ### Loop over the conductor TAGS, not over generic indices.
        for i_tag in sorted_conductor_tags:
            for j_tag in sorted_conductor_tags:

                # Get the correct matrix indices from the tags
                row = tag_to_idx[i_tag]
                col = tag_to_idx[j_tag]

                # Information about the column block of conductor 'j_tag'
                info_cond_j = surface_map[j_tag]['conductor']
                col_start_j = info_cond_j['offset']
                col_end_j = col_start_j + info_cond_j['nf']

                # Term 1 (Eq. 5.48): Contribution of the surface of conductor 'i_tag'.
                info_cond_i = surface_map[i_tag]['conductor']
                row_idx_cond_i = info_cond_i['offset']
                radius_cond_i = info_cond_i['radius']
                sum_bij = np.sum(self.T_matrix[row_idx_cond_i, col_start_j:col_end_j])
                term1 = 2 * np.pi * radius_cond_i * sum_bij

                # Term 2 (Eq. 5.48): Contribution of the sheath surface of 'i_tag'.
                term2 = 0.0
                if 'primary_insulation' in surface_map[i_tag]:
                    surface_i = surface_map[i_tag]['primary_insulation']
                    row_idx_i = surface_i['offset']
                    sum_b_prime_ij = np.sum(self.T_matrix[row_idx_i, col_start_j:col_end_j])
                    term2 = 2 * np.pi * surface_i['radius'] * sum_b_prime_ij

                ### Assign the value to the correct position in the matrix using the mapped indices.
                C_matrix[row, col] = term1 + term2

        self.C_generalized = C_matrix

    def _calculate_maxwellian_capacitance(self):
        """
        Computes the physical (Maxwellian) capacitance matrix of dimension (N-1)x(N-1),
        faithfully replicating the logic and ordering of the RIBBON.FOR code.

        The ordering of the final matrix is based on the original sequence of the
        conductors, simply removing the row/column of the reference conductor,
        as implemented in the reference source code.
        """
        gc = self.C_generalized

        # Assuming self.model.idx_ref has already been converted to 0-base in __init__
        ref_idx = self.model.mtl_idx_ref

        # --- Validations ---
        num_conductors = gc.shape[0]
        assert 0 <= ref_idx < num_conductors, f"Invalid reference index ({ref_idx})."

        # --- Step 1: Compute the required sums, as in RIBBON.FOR ---
        total_sum = np.sum(gc)
        if np.abs(total_sum) < 1e-15:
            raise ValueError("The sum of the elements of the generalized capacitance matrix is close to zero.")

        row_sums = np.sum(gc, axis=1)
        col_sums = np.sum(gc, axis=0)

        # --- Step 2: Compute the (N-1)x(N-1) matrix with the simple RIBBON.FOR ordering ---

        # Get the original conductor indices, except the reference one.
        # The order is the natural index order (0, 1, 2, ... N-1), which corresponds to I=1,N in FORTRAN.
        final_indices = [i for i in range(num_conductors) if i != ref_idx]

        # Initialize the final (N-1)x(N-1) matrix
        C_maxwellian = np.zeros((num_conductors - 1, num_conductors - 1))

        # Fill the final matrix iterating over the indices while preserving the original order
        for i_new, i_orig in enumerate(final_indices):
            for j_new, j_orig in enumerate(final_indices):

                # Apply the RIBBON.FOR / Eq. (5.21) formula
                # Indices i_orig and j_orig map directly to the rows/columns of gc, row_sums and col_sums
                correction_term = (row_sums[i_orig] * col_sums[j_orig]) / total_sum
                C_maxwellian[i_new, j_new] = gc[i_orig, j_orig] - correction_term

        self.C_maxwellian = C_maxwellian

    def run_simulation(self):
        """
        Runs the full MoM simulation, implementing the physics for
        the conducting and dielectric boundaries.

        This revised version distinguishes between observation points internal and
        external to a source boundary, implementing the potential formulas
        of Tables II.a and II.b of Clements (1975).

        NOTE: The boundary condition of the electric displacement vector on the
        sheath surface still needs to be implemented. This version computes
        the potential on every boundary.
        """
        self._calculate_collocation_points()

        # 1. Prepare the system indices and vectors
        nfs_per_surface = [2 * surface['fourier_order'] + 1 for surface in self.ordered_surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        self.D_matrix = np.zeros((self.model.N, self.model.N))
        self.V_vector = np.zeros(self.model.N)

        # 2. Assemble the [D] Matrix and the [V] Vector
        # Loop over the OBSERVATION surfaces p (rows of the matrix)
        for p, field_surface in enumerate(self.ordered_surfaces):
            tag_p = field_surface['tag']
            type_p = field_surface['type']
            radius_p = field_surface['radius']
            center_p = np.array(field_surface['center_point'])
            nf_p = nfs_per_surface[p]
            offset_p = offsets[p]

            # Get the observation points for surface p
            observation_points = self.collocation_data[tag_p][type_p]['observation']['cartesian']

            # Fill the potential vector V for the row block of surface p
            if type_p == 'conductor':
                self.V_vector[offset_p : offset_p + nf_p] = self.model.mtl[tag_p]['potential_to_infinity']

            # The boundary condition on the dielectric sheath results in 0 on the right-hand side [cite: 222]
            elif type_p == 'primary_insulation':
                self.V_vector[offset_p : offset_p + nf_p] = 0.0

            # Loop over the SOURCE surfaces q (columns of the matrix)
            for q, source_surface in enumerate(self.ordered_surfaces):
                tag_q = source_surface['tag']
                type_q = source_surface['type']
                epsilon = self.model.epsilon_out[source_surface['tag']]
                center_q = np.array(source_surface['center_point'])
                nf_q = nfs_per_surface[q]
                offset_q = offsets[q]

                # Get the source points for surface q
                source_points = self.collocation_data[tag_q][type_q]['source']['cartesian']

                # Loop over each observation point 'm' on surface 'p'
                for m in range(nf_p):
                    row_idx = offset_p + m

                    # Vector pointing from the center of the SOURCE surface 'q' to the OBSERVATION point 'm'.
                    rho_i_vector = observation_points[m] - center_q
                    rho_i = np.linalg.norm(rho_i_vector)
                    theta_i = np.arctan2(rho_i_vector[1], rho_i_vector[0])

                    # Unit vector (un_rho_i) from the center of surface 'q' to the OBSERVATION point 'm'.
                    un_rho_i = rho_i_vector / rho_i

                    # Unit normal vector (un_p) from the center of surface 'p' to the OBSERVATION point 'm'.
                    un_p = (observation_points[m] - center_p) / radius_p

                    # Loop over each basis function 'n' on surface 'q'
                    for n in range(nf_q):
                        col_idx = offset_q + n
                        harmonic_ord = n

                        # Trigonometric Term at observation point
                        is_cosine_term = (harmonic_ord % 2 != 0)
                        k = (harmonic_ord + 1) // 2 if is_cosine_term else harmonic_ord // 2
                        harmonic_term = np.cos(k * theta_i) if is_cosine_term else np.sin(k * theta_i)
                        k2epsilon = k * 2 * epsilon

                        # Source vector 'rho_b' relative to the center of the SOURCE surface 'q'
                        rho_b = np.linalg.norm(source_points[m] - center_q)

                        # ====================================================================================
                        # ==== START OF THE D MATRIX ELEMENT COMPUTATION LOGIC ==============================
                        # ====================================================================================
                        is_observer_inside = (rho_i < rho_b) and not np.isclose(rho_i, rho_b)

                        # === BLOCK 1: POTENTIAL COMPUTATION (phi) ==========================================
                        # === Applies the boundary condition V = Vm on the conducting surfaces. =============

                        if type_p == 'conductor':
                            # --- TABLE II.b: rho_i < rho_b (Interaction for an observer INSIDE the source boundary) ---
                            if is_observer_inside:
                                if harmonic_ord == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = - rho_b * np.log(rho_b) / epsilon

                                else: # Harmonic Terms (k>0)
                                    self.D_matrix[row_idx, col_idx] = rho_i**k / k2epsilon / rho_b**(k-1) * harmonic_term

                            # --- TABLE II.a: rho_i >= rho_b (Interaction for an observer OUTSIDE or ON the source boundary) ---
                            else:
                                if harmonic_ord == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = - rho_b * np.log(rho_i) / epsilon

                                else: # Harmonic Terms (k>0)
                                    self.D_matrix[row_idx, col_idx] = rho_b**(k+1) / k2epsilon / rho_i**k * harmonic_term

                        # === BLOCK 2: DISPLACEMENT VECTOR BOUNDARY CONDITION (epsilon*E) ===================
                        # === Applies continuity of the normal component of D across the dielectric sheath ==

                        elif type_p == 'primary_insulation':
                            er = field_surface['relative_permittivity']

                            # Dot product of the unit vectors in RIBBON.FOR: COS(TH - ANG)
                            RDN = np.dot(un_p, un_rho_i)

                            # 4. RIBBON.FOR: TDN = -sin(TH-ANG) = sin(ANG-TH)
                            TDN = np.cross(un_p, un_rho_i)

                            # --- TABLE II.a: rho_i = rho_b (Interaction for an observer ON the dielectric boundary) ---
                            if type_q == 'primary_insulation' and tag_p == tag_q:
                                if harmonic_ord == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = (0 - 1) * (rho_b / rho_i) * RDN

                                else: # Harmonic Terms (k>0)
                                    self.D_matrix[row_idx, col_idx] = - 0.5 * (er + 1) * (rho_b / rho_i)**(k-1) * RDN * harmonic_term

                            # --- TABLE II.a: rho_i >= rho_b (Interaction for an observer OUTSIDE the dielectric boundary) ---
                            else:
                                if harmonic_ord == 0: # Constant Term (k=0)
                                    self.D_matrix[row_idx, col_idx] = (er - 1) * (rho_b / rho_i) * RDN

                                else: # Harmonic Terms (k>0)
                                    if is_cosine_term:
                                        harmonic_term = np.cos(k * theta_i) * RDN - np.sin(k * theta_i) * TDN
                                    else:
                                        harmonic_term = np.sin(k * theta_i) * RDN + np.cos(k * theta_i) * TDN
                                    self.D_matrix[row_idx, col_idx] = 0.5 * (er - 1) * (rho_b / rho_i)**(k+1) * harmonic_term

                        # ====================================================================================
                        # ==== END OF THE D MATRIX ELEMENT COMPUTATION LOGIC ===============================
                        # ====================================================================================

        # 3. Solve the system and get the results
        self.sigma_coeffs = np.linalg.solve(self.D_matrix, self.V_vector)
        self._calculate_generalized_capacitance()
        self._calculate_maxwellian_capacitance()

    def print_results(self):
        """Prints a summary of the simulation results."""
        if self.C_maxwellian is None:
            print("Running simulation first...")
            self.run_simulation()

        if self.model.NF < 4:
            matrix_viewer(self.D_matrix, "D Matrix")
            # matrix_viewer(self.sigma_coeffs, "Sigma Coefficients")
        else:
            print(f"\nD Matrix Shape: {self.D_matrix.shape}.")

        matrix_viewer(self.C_generalized, "MoM Generalized Capacitance Matrix (F/m)")
        matrix_viewer(self.C_maxwellian, "Maxwellian Bifilar Capacitance (MoM) (F/m)")

    def plot_collocation_points(self):
        """
        Generates an interactive plot of the collocation points using Plotly,
        reflecting the new dictionary structure of self.collocation_data.
        """
        if self.collocation_data is None:
            self._calculate_collocation_points()

        # 1. Prepare the data for Plotly from the new nested structure
        plot_data = []
        # Iterate over each conductor 'tag' in the dictionary (e.g. 0, 1)
        for tag, conductor_surfaces in self.collocation_data.items():
            # Iterate over each surface of that conductor (e.g. 'conductor', 'primary_insulation')
            for surface_type, surface_data in conductor_surfaces.items():

                # Look up the corresponding radius in the self.model.surfaces list,
                # since it is not stored in self.collocation_data.
                matching_surface = next(s for s in self.model.surfaces if s['tag'] == tag and s['type'] == surface_type)
                radius = matching_surface['radius']

                # Extract and add the source point data
                for i, pt in enumerate(surface_data['source']['cartesian']):
                    plot_data.append({
                        'x': pt[0], 'y': pt[1],
                        'type': 'Source',
                        'tag': tag,
                        'surface': surface_type,
                        'radius': radius,
                        'angle_rad': surface_data['source']['angles_rad'][i]
                    })

                # Extract and add the observation point data
                for i, pt in enumerate(surface_data['observation']['cartesian']):
                    plot_data.append({
                        'x': pt[0], 'y': pt[1],
                        'type': 'Observation',
                        'tag': tag,
                        'surface': surface_type,
                        'radius': radius,
                        'angle_rad': surface_data['observation']['angles_rad'][i]
                    })

        df = pd.DataFrame(plot_data)
        fig = go.Figure()

        # 3. Add the circle shapes (this part does not change, as it already iterates over self.model.surfaces)
        for surface in self.model.surfaces:
            fig.add_shape(type="circle",
                        xref="x", yref="y",
                        x0=surface['center_point'][0] - surface['radius'], y0=surface['center_point'][1] - surface['radius'],
                        x1=surface['center_point'][0] + surface['radius'], y1=surface['center_point'][1] + surface['radius'],
                        line_color="Black", fillcolor="LightGray", opacity=0.7)

        # 4. Add the collocation points from the DataFrame
        for pt_type, color, symbol in [('Source', 'blue', 'circle'), ('Observation', 'red', 'x-thin')]:
            df_subset = df[df['type'] == pt_type]
            fig.add_trace(go.Scatter(
                x=df_subset['x'], y=df_subset['y'],
                mode='markers',
                marker=dict(color=color, symbol=symbol, size=8, line=dict(width=1, color='DarkSlateGrey')),
                name=pt_type,
                # Update customdata and hovertemplate to display the new information
                customdata=df_subset[['tag', 'surface', 'radius', 'angle_rad']],
                hovertemplate=(
                    f"<b>{pt_type}</b><br>"
                    "Conductor (tag): %{customdata[0]}<br>"
                    "Surface: %{customdata[1]}<br>"
                    "Coord X: %{x:.4f} m<br>"
                    "Coord Y: %{y:.4f} m<br>"
                    "Angle: %{customdata[3]:.3f} rad<br>"
                    "Radius: %{customdata[2]:.4f} m"
                    "<extra></extra>"
                )
            ))

        # 5. Configure the plot layout (does not change)
        fig.update_layout(
            title='Interactive Map of Collocation Points',
            xaxis_title='X Coordinate (m)',
            yaxis_title='Y Coordinate (m)',
            yaxis_scaleanchor="x",
            yaxis_scaleratio=1,
            legend_title_text='Point Type',
            template='plotly_white'
        )
        fig.show()

    def plot_charge_density(self, tag_to_plot=1):
        """
        Plots the charge density for a specific conductor, dynamically aligning
        the exact solution with the actual geometry of the system.

        Args:
            tag_to_plot (int): The 'tag' of the conductor for which the charge
                            density will be plotted.
        """
        if self.sigma_coeffs is None:
            self.run_simulation()

        # 1. Get data for the conductor to be plotted and for its pair
        all_tags = list(self.model.mtl.keys())
        if len(all_tags) != 2:
            print("Error: plot_charge_density was designed for 2-conductor systems.")
            return
        other_tag = next(tag for tag in all_tags if tag != tag_to_plot)

        center_plot = np.array(self.model.mtl[tag_to_plot]['center_point'])
        center_other = np.array(self.model.mtl[other_tag]['center_point'])

        conductor_surface = next(s for s in self.model.surfaces if s['tag'] == tag_to_plot and s['type'] == 'conductor')
        R = conductor_surface['radius']
        D = np.linalg.norm(center_plot - center_other)
        DR_ratio = D / R

        # 2. Compute the Analytical Solution with Correct Alignment and Sign
        theta_plot = np.linspace(0, 2 * np.pi, 360)

        vec_to_other = center_other - center_plot
        angle_of_max_charge = np.arctan2(vec_to_other[1], vec_to_other[0])
        denominator = DR_ratio - 2 * np.cos(theta_plot - angle_of_max_charge)

        delta_v = self.model.mtl[tag_to_plot]['potential_to_infinity'] - self.model.mtl[other_tag]['potential_to_infinity']
        numerator = (DR_ratio**2 / 4) - 1

        # FINAL FIX: Remove the abs() to preserve the sign of the charge
        charge_density_exact = (self.C_exact_bare_wires * delta_v / R) * (numerator / denominator)

        # 3. Reconstruct the MoM Solution for the Correct Conductor
        nfs_per_surface = [2 * s['fourier_order'] + 1 for s in self.model.surfaces]
        offsets = np.cumsum([0] + nfs_per_surface)

        surface_index = next(i for i, s in enumerate(self.model.surfaces) if s['tag'] == tag_to_plot and s['type'] == 'conductor')

        offset = offsets[surface_index]
        nf = nfs_per_surface[surface_index]
        coeffs_to_plot = self.sigma_coeffs[offset : offset + nf]

        charge_density_mom = np.zeros_like(theta_plot)
        charge_density_mom += coeffs_to_plot[0]
        max_k = (nf - 1) // 2
        for k in range(1, max_k + 1):
            cos_coeff = coeffs_to_plot[2 * k - 1]
            sin_coeff = coeffs_to_plot[2 * k]
            charge_density_mom += cos_coeff * np.cos(k * theta_plot) + sin_coeff * np.sin(k * theta_plot)

        # 4. Plot Generation
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(np.rad2deg(theta_plot), charge_density_exact, 'r-', label='Exact Solution')
        ax.plot(np.rad2deg(theta_plot), charge_density_mom, 'k-.', label=f'MoM (tag={tag_to_plot})')
        ax.set_title(f'Charge Distribution (Conductor {tag_to_plot}) with D/R = {DR_ratio:.2f}')
        ax.set_xlabel('Angle (Degrees)'); ax.set_ylabel('Charge Density (C/m^2)')
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.set_xticks(np.arange(0, 361, 90)); ax.set_xlim(0, 360)
        ax.legend()
        plt.tight_layout()

    def plot_harmonic_coefficients(self):
        """
        Reproduces and extends Figure 4(c) of Clements (1975), showing the magnitude
        of every coefficient of the harmonic series with adjusted indexing.

        The plot shows the ratio between the magnitude of each harmonic coefficient
        and the magnitude of the constant coefficient. The plot follows the convention:
        - j=1: Constant Term
        - j=2, 4, 6,...: Cosine Coefficients
        - j=3, 5, 7,...: Sine Coefficients
        """
        if self.sigma_coeffs is None:
            print("Running the simulation to obtain the coefficients...")
            self.run_simulation()

        # Isolate the coefficients of the first conductor
        coeffs_conductor1 = self.sigma_coeffs[:self.model.NF]

        # The constant coefficient (alpha_n1) is at index 0 in the code
        constant_term = coeffs_conductor1[0]
        if np.abs(constant_term) < 1e-15: # Avoid division by zero
            print("The constant coefficient is close to zero. Normalization is not possible.")
            return

        # --- Index Mapping for Plotting ---
        # j=1: Constant Term. Its normalized ratio is 1.0.
        plot_j_const = [1]
        ratio_const = [1.0]

        # j=2, 4, 6,...: Cosine Coefficients (indices 1, 3, 5,... in the code)
        code_indices_cos = np.arange(1, self.model.NF, 2)
        plot_j_cos = code_indices_cos + 1 # Maps [1, 3, 5] to [2, 4, 6]
        coeffs_cos = coeffs_conductor1[code_indices_cos]
        ratio_cos = np.abs(coeffs_cos) / np.abs(constant_term)

        # j=3, 5, 7,...: Sine Coefficients (indices 2, 4, 6,... in the code)
        code_indices_sin = np.arange(2, self.model.NF, 2)
        plot_j_sin = code_indices_sin + 1 # Maps [2, 4, 6] to [3, 5, 7]
        coeffs_sin = coeffs_conductor1[code_indices_sin]
        ratio_sin = np.abs(coeffs_sin) / np.abs(constant_term)

        # --- Plot Generation ---
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(12, 8))

        # Plot each series with a distinct marker
        ax.plot(plot_j_const, ratio_const, marker='s', markersize=6, linestyle='none',
                color='blue', label='Constant Term (j=1)')

        ax.plot(plot_j_cos, ratio_cos, marker='^', markersize=6, linestyle='none',
                fillstyle='none', markeredgecolor='black', label='Cosine Coefficients (j=2, 4, ...)')

        ax.plot(plot_j_sin, ratio_sin, marker='o', markersize=6, linestyle='none',
                color='red', label='Sine Coefficients (j=3, 5, ...)')

        ax.set_title(f'Normalized Magnitude of the Harmonic Coefficients (d/a = {self.DR_ratio:.1f})', fontsize=14)
        ax.set_ylabel(r'$|\alpha_{nj} / \alpha_{n1}|$', fontsize=12)
        ax.set_xlabel('Coefficient Index (j)', fontsize=12)
        ax.legend()
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=-0.05)
        ax.set_xticks(np.arange(1, self.model.NF + 1))
        ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
