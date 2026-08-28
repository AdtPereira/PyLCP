"""Phase-domain propagation parameters for multiconductor transmission lines.

Given the per-unit-length series impedance ``Z'(f)`` and shunt admittance
``Y'(f)`` matrices, this module computes the phase-domain propagation-constant
matrices and the characteristic impedance / admittance matrices, following the
classic MTL theory (Paul, 2008; Ametani et al., 2015)::

    gamma_v = sqrtm(Z' Y')          gamma_i = sqrtm(Y' Z')
    Zc      = Y'^{-1} gamma_i       Yc      = Z'^{-1} gamma_v

It is the phase-domain counterpart of
``analytical_forms/modal_analysis.py`` (which performs the modal
eigen-decomposition used to reproduce Chapter 5 of Andreata's thesis).

The loop below was previously inlined in
``analytical_forms/overhead_lines.py::PerUnitParameters.pul_matrices``; it is
factored out here so both the overhead-line and single-core-cable pipelines
share a single implementation.
"""

import numpy as np
from scipy import linalg

__all__ = ["phase_domain_propagation"]


def phase_domain_propagation(Zs, Ysh):
    """Compute the phase-domain propagation matrices.

    Parameters
    ----------
    Zs, Ysh : array_like, shape ``(Nf, N, N)``, complex
        Per-unit-length series impedance (Ohm/m) and shunt admittance (S/m)
        matrices, one ``N x N`` block per frequency.

    Returns
    -------
    dict
        Keys (each ``ndarray`` of shape ``(Nf, N, N)``, complex):

        ``propagation_voltage_matrix``
            ``gamma_v = sqrtm(Zs @ Ysh)`` -- voltage propagation constant.
        ``propagation_current_matrix``
            ``gamma_i = sqrtm(Ysh @ Zs)`` -- current propagation constant.
        ``characteristic_impedance_matrix``
            ``Zc = Ysh^{-1} @ gamma_i``.
        ``characteristic_admittance_matrix``
            ``Yc = Zs^{-1} @ gamma_v``.
    """
    Zs = np.asarray(Zs, dtype=complex)
    Ysh = np.asarray(Ysh, dtype=complex)

    if Zs.ndim != 3 or Zs.shape[1] != Zs.shape[2] or Zs.shape != Ysh.shape:
        raise ValueError(
            "Zs and Ysh must both have shape (Nf, N, N); "
            f"got {Zs.shape} and {Ysh.shape}."
        )

    num_freq, n, _ = Zs.shape
    gamma_v = np.zeros((num_freq, n, n), dtype=complex)
    gamma_i = np.zeros((num_freq, n, n), dtype=complex)
    zc = np.zeros((num_freq, n, n), dtype=complex)
    yc = np.zeros((num_freq, n, n), dtype=complex)
    identity = np.identity(n)

    for i in range(num_freq):
        zs_i = Zs[i]
        ysh_i = Ysh[i]

        # Inverses via LU factorisation (matches the original overhead-line code).
        ysh_inv_i = linalg.lu_solve(linalg.lu_factor(ysh_i), identity)
        zs_inv_i = linalg.lu_solve(linalg.lu_factor(zs_i), identity)

        gamma_i[i] = linalg.sqrtm(ysh_i @ zs_i)
        gamma_v[i] = linalg.sqrtm(zs_i @ ysh_i)

        zc[i] = ysh_inv_i @ gamma_i[i]
        yc[i] = zs_inv_i @ gamma_v[i]

    return {
        "propagation_voltage_matrix": gamma_v,
        "propagation_current_matrix": gamma_i,
        "characteristic_impedance_matrix": zc,
        "characteristic_admittance_matrix": yc,
    }
