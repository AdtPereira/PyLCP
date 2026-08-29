"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

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

import numpy as np
from scipy.special import iv, kv, jv, jvp
from scipy.linalg import lu_factor, lu_solve
from mtl_main.source import MulticonductorTransmissionLine

# 4.2.2 Per-Unit-Length Inductance and Capacitance for Wire-Type Lines [4]
class WiresHomogeneousMedia:
    """
    This class contains the analytical formulation for wire-type lines in a
    homogeneous medium.

    It calculates:
    1. Frequency-dependent series loop impedance for 2-wire (bifilar) systems.
    2. Frequency-independent (static) L/C parameters for 2-wire systems.
    3. Frequency-independent (static) L/C matrices for N-wire systems (with 1 reference).
    
    This class follows a structure similar to InternalPerUnitParameters,
    using composition (taking an MTL model) and vectorizing over frequency.

    References:
    [1] Clayton R. Paul, "Introduction to Electromagnetic Compatibility", 2nd Edition, Wiley, 2007.
        4.2.2 Per-Unit-Length Inductance and Capacitance for Wire-Type Lines
        5.2.1 Wide-Separation Approximations for Wires in Homogeneous Media 
    
    """

    def __init__(self, model: MulticonductorTransmissionLine, f: np.ndarray):
        """
        Initializes the calculator with the MTL model and frequency array.

        Args:
            model (MulticonductorTransmissionLine): The MTL object containing
                geometry and material properties.
            f (np.ndarray): Array of frequencies (in Hz) for calculations.
        """
        self.model = model
        self.f = np.asarray(f)
        self.num_freq = len(self.f)
        self.num_conductors = len(model.surfaces)
        self.jw = 1j * 2 * np.pi * self.f
        
        # Store model properties for easier access
        self.mtl = self.model.mtl
        self.mu = self.model.mu
        self.epsilon_out = self.model.epsilon_out

        # Constants for Kelvin functions
        self.kelvin_exp = np.exp(1j * 3 * np.pi / 4)

    # --- Kelvin Function Helpers ---

    def ber(self, xi):
        """ Kelvin ber(xi) function """
        return np.real(jv(0, xi * self.kelvin_exp))

    def bei(self, xi):
        """ Kelvin bei(xi) function """
        return np.imag(jv(0, xi * self.kelvin_exp))

    def ber_prime(self, xi):
        """ Kelvin ber'(xi) derivative function """
        return np.real(self.kelvin_exp * jvp(0, xi * self.kelvin_exp, 1))

    def bei_prime(self, xi):
        """ Kelvin bei'(xi) derivative function """
        return np.imag(self.kelvin_exp * jvp(0, xi * self.kelvin_exp, 1))

    def bifilar_series_impedance(self):
        """
        Computes the series impedance of the bifilar loop, vectorized over frequency.
        Assumes a system with one active conductor (tag=1) and one return conductor (tag=0),
        and that both are identical (uses the parameters of conductor 1).

        The output shape was changed to (N, 1, 1) for consistency
        with matrix outputs.

        Returns:
            dict: Dictionary containing the total series loop impedance 'Zs_loop'
                  (complex 3D matrix, shape=(N, 1, 1)) and the high-frequency
                  resistance 'Rhf_loop' (real 3D matrix, shape=(N, 1, 1)).
        """
        if self.num_conductors != 2:
            raise ValueError("Bifilar impedance computation is only valid for 2 conductors.")

        # 1. Get the conductor data (assuming tags 0 and 1)
        conductor = self.mtl[1]
        cond_ref = self.mtl[0]

        # 2. Compute the distance (scalar)
        center_active = np.array(conductor['center_point'])
        center_ref = np.array(cond_ref['center_point'])
        D10 = np.linalg.norm(center_active - center_ref)

        # 3. Extract parameters
        w = 2 * np.pi * self.f
        ap = conductor['radius'][1]
        sigma = conductor['conductivity']

        # 4. Initialize the output vectors (CHANGE HERE)
        # Shape: (N_freq, 1, 1)
        Zs = np.zeros((self.num_freq, 1, 1), dtype=complex)
        Rhf = np.zeros((self.num_freq, 1, 1), dtype=float)

        # --- AC computation (f > 0) ---
        ac_idx = (self.f > 0)
        if np.any(ac_idx):
            with np.errstate(divide='ignore', invalid='ignore'):
                delta = np.sqrt(2 / (w[ac_idx] * self.mu[0] * sigma))
                Rs = 1 / (sigma * delta)
                Xi = np.sqrt(2) * ap / delta
                constant_term = 1 / (np.sqrt(2) * np.pi * ap * sigma * delta)

            # HF resistance (Rhf) and external inductance (Lext)
            s_2rw = D10 / (2 * ap)

            # CHANGE HERE: assign to the slice [ac_indices, 0, 0]
            Rhf[ac_idx, 0, 0] = Rs / (np.pi * ap) * s_2rw / np.sqrt(s_2rw**2 - 1)
            L_ext = self.mu[0] / np.pi * np.arccosh(s_2rw)

            # Internal impedance (Zi) with Bessel functions
            ber_bei = self.ber(Xi) + 1j * self.bei(Xi)
            beip_berp = self.bei_prime(Xi) - 1j * self.ber_prime(Xi)

            # Avoid division by zero if the denominator is null
            Zi_ac = np.full(w[ac_idx].shape, np.nan, dtype=complex) # 1D array
            valid_den = (np.abs(beip_berp) > 1e-12)
            Zi_ac[valid_den] = constant_term[valid_den] * ber_bei[valid_den] / beip_berp[valid_den]

            # Total series impedance (Zs)
            Zs[ac_idx, 0, 0] = 2 * Zi_ac + self.jw[ac_idx] * L_ext

        # --- DC computation (f = 0) ---
        dc_idx = (self.f == 0)
        if np.any(dc_idx):
            # The DC internal impedance is the DC resistance
            R_dc = 1 / (sigma * np.pi * ap**2) # scalar

            # The total loop impedance at DC is 2 * R_dc
            Zs[dc_idx, 0, 0] = 2 * R_dc
            Rhf[dc_idx, 0, 0] = 0.0
            
        return {
            'series_impedance_matrix': Zs,
            'high_frequency_limit': Rhf + self.jw[:, np.newaxis, np.newaxis] * L_ext
        }

    def bifilar_static_params(self):
        """
        Computes the static (PUL) capacitance and inductance for the
        bifilar line, using the exact and approximate formulas.
        (This method is frequency-independent).

        Returns:
            dict: A dictionary containing capacitance and inductance
                  ('exact', 'approximate').
        """
        if self.num_conductors != 2:
            raise ValueError("Bifilar static-parameter computation is only valid for 2 conductors.")

        # 2. Direct access to the conductor data (tags 0 and 1)
        c0, c1 = self.mtl[0], self.mtl[1]
        rw0, rw1 = c0['radius'][1], c1['radius'][1]

        # 3. Dynamic computation of the distance 's'
        s = np.linalg.norm(np.array(c0['center_point']) - np.array(c1['center_point']))

        # 4. Validate the external medium and define epsilon
        eps_out_0, eps_out_1 = self.epsilon_out
        pi2e = 2 * np.pi * eps_out_0
        if not np.isclose(eps_out_0, eps_out_1):
            print("Warning: the external medium must be homogeneous. Using the permittivity of conductor 0.")

        # 5. Capacitance computation
        den_approx = np.log((s**2) / (rw0 * rw1))
        capacitance_approx = pi2e / den_approx
        arg_arccosh = (s**2 - rw0**2 - rw1**2) / (2 * rw0 * rw1)
        den_exact = np.arccosh(arg_arccosh)
        capacitance_exact = pi2e / den_exact

        # 6. Inductance computation (valid for a non-magnetic medium)
        mu_0 = self.mu[0]
        inductance_approx = mu_0 * eps_out_0 / capacitance_approx
        inductance_exact = mu_0 * eps_out_0 / capacitance_exact

        # 7. Return the results
        return {
            'capacitance': {
                'exact': capacitance_exact,
                'approximate': capacitance_approx,
            },
            'inductance': {
                'exact': inductance_exact,
                'approximate': inductance_approx,
            }
        }

    def n_wires_external_inductance(self):
        """
        Computes the external inductance matrix for N conductors (N-1 active
        + 1 reference (tag=0)). (Frequency-independent).

        Returns:
            dict: Dictionary containing the inductance matrix 'L_ext' (N-1 x N-1).
        """
        # The matrix size is (N_total - 1)
        n_plus_1 = self.num_conductors
        if n_plus_1 <= 1:
            return {'L_ext': np.array([])}

        n_active = n_plus_1 - 1
        Lext = np.zeros((n_active, n_active), dtype=float)

        # Data of the reference conductor (tag=0)
        ref_center = np.array(self.mtl[0]['center_point'])
        rw0 = self.mtl[0]['radius'][1]
        mu_2pi = self.mu[0] / (2 * np.pi)

        # Active tags go from 1 to n_active
        for i in range(1, n_plus_1):
            center_i = np.array(self.mtl[i]['center_point'])
            rw_i = self.mtl[i]['radius'][1]
            di0 = np.linalg.norm(center_i - ref_center)

            for j in range(1, n_plus_1):
                center_j = np.array(self.mtl[j]['center_point'])

                # Map the tag index (1..N) to the matrix index (0..N-1)
                idx_i = i - 1
                idx_j = j - 1

                if i == j:  # Self-inductance
                    Lext[idx_i, idx_j] = mu_2pi * np.log(di0 ** 2 / (rw0 * rw_i))

                else:  # Mutual inductance
                    dj0 = np.linalg.norm(center_j - ref_center)
                    dij = np.linalg.norm(center_i - center_j)
                    Lext[idx_i, idx_j] = mu_2pi * np.log(di0 * dj0 / (rw0 * dij))

        return {'L_ext': Lext}

    def n_wires_capacitance(self, L):
        """
        Computes the per-unit-length capacitance matrix (C) from the
        inductance matrix (L) for a homogeneous medium, C = mu * eps * L^-1.
        (Frequency-independent).

        Args:
            L (np.ndarray): The inductance matrix L (n x n).

        Returns:
            dict: Dictionary containing the capacitance matrix 'C_pul' (n x n).
        """
        if L.shape[0] != L.shape[1]:
            raise ValueError("The inductance matrix must be square.")
        if L.size == 0:
            return {'C_pul': np.array([])}

        I = np.eye(L.shape[0])

        # Assume a homogeneous medium
        mu = self.mu[0]
        epsilon = self.epsilon_out[0]

        # Inversion via LU-solve (numerically stable)
        try:
            C = mu * epsilon * lu_solve(lu_factor(L), I)
            return {'C_pul': C}
        except np.linalg.LinAlgError:
            print("Warning: the inductance matrix is singular. Could not compute the capacitance.")
            return {'C_pul': np.full_like(L, np.nan)}

    def get_all_terms(self):
        """
        Computes every relevant parameter for the provided MTL configuration
        and returns them in a single dictionary.

        This is the main method to be called from outside.
        """
        results = {}

        # --- Static Parameters (N-wire, L/C matrices) ---
        # Note: assumes conductor 0 as the reference
        if self.num_conductors > 1:
            l_ext_dict = self.n_wires_external_inductance()
            results.update(l_ext_dict)

            c_pul_dict = self.n_wires_capacitance(L=l_ext_dict['L_ext'])
            results.update(c_pul_dict)

        # --- Specific Parameters (Bifilar / 2-wire) ---
        if self.num_conductors == 2:
            # Static parameters (exact and approximate)
            static_params = self.bifilar_static_params()
            results['bifilar_static'] = static_params

            # Dynamic parameters (loop impedance vs frequency)
            impedance_params = self.bifilar_series_impedance()
            results.update(impedance_params)

        return results