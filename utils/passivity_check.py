"""Passivity *assessment* (not enforcement) for per-unit-length line matrices.

A passive multiconductor transmission line has, at every frequency,

    eig( Re{Z'(jw)} ) >= 0      series-resistance matrix  R' = Re{Z'}  is PSD
    eig( Re{Y'(jw)} ) >= 0      shunt-conductance matrix  G' = Re{Y'}  is PSD
    eig( Re{Yc(jw)}  ) >= 0     characteristic admittance (Gustavsen 2008, eq. 3)

Negative eigenvalues beyond numerical noise point to a *formulation* problem
(an earth-return admittance model producing negative conductance, a wrong
semiconducting-layer permittivity, a bad sign in an approximation, ...).

This module only reports such violations. Repairing them would be the
FRP / FMP perturbation + Hamiltonian passivity-assessment pipeline of
Gustavsen (2008), which is deliberately out of scope for a
parameter-calculation library: `Z'` / `Y'` here are physically-grounded
closed-form samples, not a fitted rational macromodel.

REFERENCES
[1] B. Gustavsen, "Fast Passivity Enforcement for Pole-Residue Models by
    Perturbation of Residue Matrix Eigenvalues," IEEE Trans. Power Delivery,
    vol. 23, no. 4, pp. 2278-2285, Oct. 2008 -- passivity criterion, eq. (3).
"""

import numpy as np

__all__ = [
    "min_real_part_eigenvalue",
    "check_pul_passivity",
    "print_passivity_report",
]

# label -> (dict key, physical meaning of Re{M})
_MATRIX_KEYS = {
    "series_resistance": "series_impedance_matrix",         # Re{Z'} = R'
    "shunt_conductance": "shunt_admittance_matrix",          # Re{Y'} = G'
    "char_admittance": "characteristic_admittance_matrix",   # Re{Yc}
}


def min_real_part_eigenvalue(matrix):
    """Smallest eigenvalue of the (symmetrised) real part of ``M`` per frequency.

    Parameters
    ----------
    matrix : array_like, shape ``(Nf, N, N)``, complex

    Returns
    -------
    ndarray, shape ``(Nf,)``, float
    """
    matrix = np.asarray(matrix)
    real = matrix.real
    real_sym = 0.5 * (real + np.transpose(real, (0, 2, 1)))
    return np.linalg.eigvalsh(real_sym)[:, 0]


def _summarise(frequencies, matrix, rel_tol, full):
    frequencies = np.asarray(frequencies, dtype=float)
    matrix = np.asarray(matrix)

    min_eig = min_real_part_eigenvalue(matrix)
    # per-frequency scale: largest |Re{M}| entry (fallback 1.0)
    scale = np.abs(matrix.real).reshape(matrix.shape[0], -1).max(axis=1)
    scale = np.where(scale > 0, scale, 1.0)

    normalised = min_eig / scale
    violating = min_eig < (-rel_tol * scale)
    worst = int(np.argmin(normalised))

    out = {
        "passive": not bool(violating.any()),
        "n_violations": int(violating.sum()),
        "worst_min_eig": float(min_eig[worst]),
        "worst_min_eig_normalised": float(normalised[worst]),
        "worst_frequency_hz": float(frequencies[worst]),
    }
    if full:
        out["min_eig"] = min_eig
        out["violating_frequencies_hz"] = frequencies[violating]
    return out


def check_pul_passivity(frequencies, matrices, rel_tol=1e-8, full=False):
    """Assess passivity of the PUL / propagation matrices.

    Parameters
    ----------
    frequencies : array_like, shape ``(Nf,)``
    matrices : dict
        Any of ``series_impedance_matrix``, ``shunt_admittance_matrix``,
        ``characteristic_admittance_matrix`` (each ``(Nf, N, N)``). Missing
        keys are skipped -- so this accepts the dict returned by
        ``PerUnitParameters.quasi_tem_approx_matrices`` /
        ``propagation_matrices`` as-is.
    rel_tol : float
        A frequency is deemed passive for a given matrix if its minimum
        real-part eigenvalue exceeds ``-rel_tol * scale`` (scale = largest
        ``|Re{M}|`` entry at that frequency).
    full : bool
        If True, each check also carries the full ``min_eig`` array and the
        list of violating frequencies.

    Returns
    -------
    dict
        ``{'passive': bool | None, 'checks': {label: summary, ...}}``.
        ``passive`` is ``None`` when no recognised matrix was supplied.
    """
    checks = {}
    for label, key in _MATRIX_KEYS.items():
        M = matrices.get(key)
        if M is not None:
            checks[label] = _summarise(frequencies, M, rel_tol, full)
    passive = all(c["passive"] for c in checks.values()) if checks else None
    return {"passive": passive, "checks": checks}


def print_passivity_report(report, title="Passivity assessment"):
    """Pretty-print the dict returned by :func:`check_pul_passivity`."""
    line = "=" * 74
    print("\n" + line)
    print(f"  {title}  --  Gustavsen (2008), eq. 3  (assessment only)")
    print(line)
    if not report["checks"]:
        print("  no recognised matrices to check.")
        print(line + "\n")
        return
    for label, c in report["checks"].items():
        flag = "OK  " if c["passive"] else "!!  "
        print(f"  {flag}{label:18s}  min eig(Re) = {c['worst_min_eig']:+.3e}"
              f"   (norm. {c['worst_min_eig_normalised']:+.2e})"
              f"   @ {c['worst_frequency_hz']:.4g} Hz"
              f"   [{c['n_violations']} viol.]")
    print("-" * 74)
    verdict = report["passive"]
    print(f"  overall: {'PASSIVE' if verdict else 'NON-PASSIVE  -- check formulation'}")
    print(line + "\n")
