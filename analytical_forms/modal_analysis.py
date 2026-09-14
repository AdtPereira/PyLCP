"""Modal (eigen) decomposition of multiconductor transmission lines.

Reproduces the modal-domain propagation analysis of Chapter 5 of

    ANDREATA, Luis Eduardo Batista. *Analise das Caracteristicas de Propagacao
    e de Transitorios Eletromagneticos em Cabos Subterraneos Instalados em
    Tubos Nao Metalicos no Contexto de Parques Eolicos.* Dissertacao de
    Mestrado, PPGEE/UFMG, 2025.

Given the per-unit-length ``Z'(f)`` and ``Y'(f)`` matrices of an N-conductor
line, this module

1. diagonalises ``Y' Z'`` -> current transformation matrix ``T_I`` and modal
   eigenvalues ``lambda`` (eq. 5.19);
2. builds ``T_V = (T_I^T)^{-1}`` so that ``T_V^{-1} = T_I^T`` exactly
   (eq. 5.20) and ``Z' Y' = T_V lambda T_V^{-1}`` (eq. 5.18);
3. tracks the eigenvector columns across frequency so each modal curve is
   continuous -- this is the *switching-back procedure* of Gustavsen (2008),
   sec. IV-A, itself the scalar-product correlation of Wedepohl, Nguyen &
   Irwin (1996), sec. 6. It is exactly the ``intercheig`` sub-routine that
   Andreata (sec. 5.2) reused from Gustavsen's Matrix Fitting Toolbox. Here
   the sweep is anchored at the middle of the band and radiates both ways;
   the optimal assignment is solved with the Hungarian algorithm rather than
   Gustavsen's greedy row-maximum;
4. returns the modal propagation constants and characteristic impedances
   (eqs. 5.21-5.32);
5. labels the 6 modes of Configuration 1 (ground; inter-sheath 1 & 2;
   coaxial 1/2/3) by projecting the low-frequency ``T_I`` columns onto the
   reference patterns of eqs. 5.34-5.36;
6. reports (does not enforce) passivity of ``Z'`` / ``Y'`` via
   :mod:`utils.passivity_check` (Gustavsen 2008, eq. 3).

Only Configuration 1 (three directly-buried single-core cables, 6 conductors)
is in scope for now; the reference-pattern dictionary and role helpers are
written so the 7-conductor cases can be added later.

REFERENCES
[1] L. M. Wedepohl, H. V. Nguyen, G. D. Irwin, "Frequency-Dependent
    Transformation Matrices for Untransposed Transmission Lines using
    Newton-Raphson Method," IEEE Trans. Power Delivery, 11(3), 1996 -- sec. 6
    (eigenvector correlation / switch-over tracking).
[2] B. Gustavsen, "Fast Passivity Enforcement for Pole-Residue Models by
    Perturbation of Residue Matrix Eigenvalues," IEEE Trans. Power Delivery,
    23(4), 2008 -- sec. IV-A (switching-back procedure), eq. 3 (passivity).
"""

import warnings
import numpy as np
from scipy.optimize import linear_sum_assignment

__all__ = [
    "ModalDecomposition",
    "canonical_permutation",
    "default_scc_roles",
    "REFERENCE_PATTERNS_6C",
    "MODE_LABELS_6C",
]

# --- Role ordering used to bring Z'/Y' to Andreata's canonical layout ---------
# Andreata assembles Z'/Y' with all cores first, then all sheaths, then the ECC
# (sec. 5.3). pyLCP's SCC pipeline instead interleaves them per cable
# (np.kron(I_N, Z_ij) -> [c1, s1, c2, s2, ...]). ``canonical_permutation``
# converts from the pyLCP order to Andreata's.
_ROLE_RANK = {"core": 0, "sheath": 1, "armor": 2, "ecc": 3}

MODE_LABELS_6C = (
    "ground",
    "inter_sheath_1",
    "inter_sheath_2",
    "coaxial_1",
    "coaxial_2",
    "coaxial_3",
)

# Reference T_I columns in the canonical order [c1, c2, c3, s1, s2, s3].
# Values transcribed from Andreata eqs. 5.34 (1 Hz) and 5.35 (10 kHz): the
# inter-sheath / ground modes are frequency-invariant; the coaxial modes are
# taken at low frequency, where they are distinguishable (all-cores / core-core
# / core-vs-two-cores). Classification uses a *normalised projection*, so exact
# magnitudes and global sign do not matter -- only the pattern.
REFERENCE_PATTERNS_6C = {
    "ground":         np.array([0, 0, 0,  1,  1,  1], dtype=float),
    "inter_sheath_1": np.array([0, 0, 0,  1,  0, -1], dtype=float),
    "inter_sheath_2": np.array([0, 0, 0, -1,  2, -1], dtype=float),
    "coaxial_1":      np.array([1,  1,  1, 0, 0, 0], dtype=float),
    "coaxial_2":      np.array([1,  0, -1, 0, 0, 0], dtype=float),
    "coaxial_3":      np.array([-1, 2, -1, 0, 0, 0], dtype=float),
}


def default_scc_roles(num_cables, conductors_per_cable=2, num_ecc=0):
    """Conductor roles for a homogeneous SCC system in pyLCP's native order.

    Returns a list of ``(role, cable_index)`` tuples aligned with the rows of
    ``Z'`` / ``Y'`` as assembled by
    ``InternalPerUnitParameters.matrices`` (``np.kron(I_N, Z_ij)`` with the
    per-cable block ordered core, sheath[, armor]), followed by ``num_ecc``
    ``("ecc", i)`` entries for the earth continuity conductor(s) appended
    after the SCC block (Configs 3/4 -- see ``flat_scc_with_ecc_cable_model``
    in ``models/single_core_cable.py``, which places the ECC as the last
    conductor). ``canonical_permutation`` sorts ``"ecc"`` after
    ``"core"``/``"sheath"``/``"armor"`` regardless (``_ROLE_RANK``), matching
    the "by-type" conductor order the external MATLAB developer also uses for
    the ECC (see ``andreata_case3/andreata_case3.py``'s
    ``conductor_order=[0, 3, 1, 4, 2, 5, 6]`` -- the ECC stays at index 6 in
    both orderings).

    NB: for the ECC cases, do **not** derive ``num_cables``/
    ``conductors_per_cable`` from ``MulticonductorTransmissionLine.
    num_sc_cables``/``num_conductors_per_scc`` -- the strategy branch that
    builds the heterogeneous SCC+ECC model counts those attributes
    differently (one "cable" per physical conductor). Pass the known
    physical layout explicitly instead (``num_cables=3,
    conductors_per_cable=2, num_ecc=1``).
    """
    names = ["core", "sheath", "armor"][:conductors_per_cable]
    roles = [(names[k % conductors_per_cable], k // conductors_per_cable)
             for k in range(num_cables * conductors_per_cable)]
    roles += [("ecc", i) for i in range(num_ecc)]
    return roles


def canonical_permutation(conductor_roles):
    """Index array that reorders conductors to ``[cores | sheaths | ecc]``.

    ``conductor_roles`` is a list of ``(role, cable_index)`` tuples (see
    :func:`default_scc_roles`). Ties are broken by cable index, then by the
    original position, so the permutation is deterministic.
    """
    order = sorted(
        range(len(conductor_roles)),
        key=lambda i: (_ROLE_RANK.get(conductor_roles[i][0], 99),
                       conductor_roles[i][1], i),
    )
    return np.asarray(order, dtype=int)


def _refine_degenerate_clusters(lam, TI, Zs_k, rtol):
    """Re-pick eigenvectors inside near-degenerate eigenvalue clusters.

    ``np.linalg.eig`` returns an arbitrary basis of a repeated-eigenvalue
    subspace. Within such a cluster we diagonalise the restricted modal series
    impedance ``B = U^T Zs U`` (complex-symmetric, so its eigenvectors satisfy
    ``v_i^T v_j = 0``); ``U @ V`` then also diagonalises ``Z_m`` on that
    subspace, giving well-defined per-mode ``Z_cm``.
    """
    n = len(lam)
    lam = lam.copy()
    TI = TI.copy()
    used = np.zeros(n, dtype=bool)
    for i in np.argsort(-np.abs(lam)):
        if used[i]:
            continue
        cluster = [
            j for j in range(n)
            if not used[j] and abs(lam[j] - lam[i]) <= rtol * (abs(lam[i]) + 1e-30)
        ]
        for j in cluster:
            used[j] = True
        if len(cluster) < 2:
            continue
        U = TI[:, cluster]
        B = U.T @ Zs_k @ U
        _, VB = np.linalg.eig(B)
        for c in range(VB.shape[1]):
            s = VB[:, c] @ VB[:, c]
            if abs(s) > 1e-30:
                VB[:, c] = VB[:, c] / np.sqrt(s)
        new = U @ VB
        new = new / np.linalg.norm(new, axis=0, keepdims=True)
        TI[:, cluster] = new
        # keep the individual eigenvalues -- the refined basis is a better set
        # of directions within the (near-)degenerate subspace, but the small
        # eigenvalue spread carries real modal-velocity/attenuation content.
    return lam, TI


def _normalise_columns(vectors):
    """Canonicalise eigenvector columns for *comparison* (not for tracking).

    Rotates each column so its largest-magnitude entry is real-positive, then
    scales to unit infinity-norm. Returns a real-part view suitable for
    pattern projection.
    """
    out = np.array(vectors, dtype=complex)
    for j in range(out.shape[1]):
        col = out[:, j]
        k = np.argmax(np.abs(col))
        if col[k] != 0:
            col = col * np.conj(col[k]) / np.abs(col[k])
        m = np.max(np.abs(col))
        if m > 0:
            col = col / m
        out[:, j] = col
    return out


class ModalDecomposition:
    """Frequency-swept modal decomposition of an N-conductor line.

    Parameters
    ----------
    frequencies : array_like, shape (Nf,)
    Zs, Ysh : array_like, shape (Nf, N, N), complex
        Per-unit-length series impedance and shunt admittance.
    conductor_roles : list of (str, int), optional
        ``(role, cable_index)`` per conductor, in the row order of
        ``Zs``/``Ysh``. Used to reorder to the canonical
        ``[cores | sheaths | ecc]`` layout and to label modes. If omitted the
        matrices are used as-is and modes are left unlabelled.
    permutation : array_like of int, optional
        Explicit reordering (overrides ``conductor_roles``-derived one).
    """

    def __init__(self, frequencies, Zs, Ysh, conductor_roles=None, permutation=None):
        self.f = np.asarray(frequencies, dtype=float)
        self.w = 2.0 * np.pi * self.f
        Zs = np.asarray(Zs, dtype=complex)
        Ysh = np.asarray(Ysh, dtype=complex)
        if Zs.shape != Ysh.shape or Zs.ndim != 3 or Zs.shape[1] != Zs.shape[2]:
            raise ValueError(
                f"Zs and Ysh must be (Nf, N, N); got {Zs.shape} and {Ysh.shape}."
            )
        if Zs.shape[0] != self.f.size:
            raise ValueError("frequencies length must match Zs.shape[0].")

        self.n = Zs.shape[1]
        if permutation is None and conductor_roles is not None:
            permutation = canonical_permutation(conductor_roles)
        if permutation is None:
            permutation = np.arange(self.n)
        self.perm = np.asarray(permutation, dtype=int)
        self.roles = (
            [conductor_roles[i] for i in self.perm]
            if conductor_roles is not None else None
        )

        # reorder rows and columns to the canonical layout
        self.Zs = Zs[:, self.perm][:, :, self.perm]
        self.Ysh = Ysh[:, self.perm][:, :, self.perm]

        self._raw = None
        self._tracked = None
        self.diagnostics = {}

    # ------------------------------------------------------------------ #
    # 1. per-frequency eigen-decomposition                              #
    # ------------------------------------------------------------------ #
    def decompose(self, cluster_rtol=3e-2):
        """Diagonalise ``Y' Z'`` at every frequency (no cross-frequency tracking).

        ``cluster_rtol`` controls the near-degenerate-eigenvalue refinement:
        within a cluster of eigenvalues that agree to this relative tolerance,
        ``np.linalg.eig`` returns an arbitrary basis of the eigenspace, so the
        columns are re-chosen to also diagonalise the modal series impedance
        ``Z_m`` (exploiting ``v_i^T v_j = 0`` for eigenvectors of a
        complex-symmetric matrix). This removes the high-frequency glitches on
        the three near-identical coaxial modes of Configuration 1.
        """
        nf, n = self.f.size, self.n
        lam = np.zeros((nf, n), dtype=complex)
        TI = np.zeros((nf, n, n), dtype=complex)

        for k in range(nf):
            w_k, v_k = np.linalg.eig(self.Ysh[k] @ self.Zs[k])
            v_k = v_k / np.linalg.norm(v_k, axis=0, keepdims=True)
            if cluster_rtol:
                w_k, v_k = _refine_degenerate_clusters(
                    w_k, v_k, self.Zs[k], cluster_rtol)
            lam[k] = w_k
            TI[k] = v_k

        self._raw = {"lambda": lam, "TI": TI}
        return self._raw

    # ------------------------------------------------------------------ #
    # 2. cross-frequency eigenvector tracking                            #
    # ------------------------------------------------------------------ #
    def track(self, seed="mid"):
        """Reorder/rephase eigenvector columns for continuity across frequency.

        Implements the *switching-back procedure* (Gustavsen 2008, sec. IV-A;
        Wedepohl-Nguyen-Irwin 1996, sec. 6 -- the ``intercheig`` correlation
        that Andreata sec. 5.2 reused): at each step the current frame's
        columns are matched to the previous frame's by maximising
        ``sum_j |<t_j(f_k), t_j(f_{k-1})>|`` (solved optimally with the
        Hungarian algorithm), then phase-aligned.

        ``seed`` selects the anchor frame from which the sweep radiates in both
        directions:

        * ``'mid'`` -- geometric middle of the band (default; the coaxial modes
          of Config. 1 are best separated in the mid decades);
        * ``'high'`` / ``'low'`` -- highest / lowest frequency;
        * ``int`` -- explicit frame index.
        """
        if self._raw is None:
            self.decompose()

        lam = self._raw["lambda"].copy()
        TI = self._raw["TI"].copy()
        nf, n = self.f.size, self.n

        if seed == "high":
            anchor = nf - 1
        elif seed == "low":
            anchor = 0
        elif seed == "mid":
            anchor = nf // 2
        else:
            anchor = int(seed)

        def _align(prev, cur):
            sim = np.abs(TI[cur].conj().T @ TI[prev])   # (n_cur, n_prev)
            row_cur, col_prev = linear_sum_assignment(-sim)
            order = np.empty(n, dtype=int)
            order[col_prev] = row_cur                   # prev-mode j -> cur column order[j]
            TI[cur] = TI[cur][:, order]
            lam[cur] = lam[cur][order]
            for j in range(n):
                dot = np.vdot(TI[prev][:, j], TI[cur][:, j])
                if abs(dot) > 0:
                    TI[cur][:, j] *= np.conj(dot) / abs(dot)

        for cur in range(anchor + 1, nf):
            _align(cur - 1, cur)
        for cur in range(anchor - 1, -1, -1):
            _align(cur + 1, cur)

        self._tracked = {"lambda": lam, "TI": TI}
        return self._tracked

    # ------------------------------------------------------------------ #
    # 3. modal parameters                                                #
    # ------------------------------------------------------------------ #
    def modal_parameters(self, track=True, seed="mid"):
        """Return the modal propagation quantities.

        Returns
        -------
        dict with (arrays ``(Nf, N)`` unless noted):
            ``frequencies``     (Nf,)
            ``omega``           (Nf,)
            ``lambda``          Y'Z' eigenvalues
            ``gamma``           modal propagation constant Gamma_m (Re >= 0)
            ``alpha``           Re(Gamma_m)          [Np/m]
            ``beta``            Im(Gamma_m)          [rad/m]
            ``vphase``          omega / Im(Gamma_m)  [m/s]
            ``Zm``, ``Ym``      modal series impedance / shunt admittance (diag)
            ``Zcm``, ``Ycm``    modal characteristic impedance / admittance
            ``TI``, ``TV``      transformation matrices  (Nf, N, N)
            ``mode_labels``     tuple of N strings (or None if roles unknown)
            ``decoupling_freq`` f_c per sheath (eq. 5.33) or None
        """
        data = self.track(seed=seed) if track else (self._raw or self.decompose())
        lam = data["lambda"]
        TI = data["TI"]
        nf, n = self.f.size, self.n

        TV = np.zeros_like(TI)
        Zm = np.zeros((nf, n), dtype=complex)
        Ym = np.zeros((nf, n), dtype=complex)
        offdiag_z = np.zeros(nf)
        consistency = np.zeros(nf)

        for k in range(nf):
            ti = TI[k]
            ti_inv = np.linalg.inv(ti)
            TV[k] = ti_inv.T                    # T_V = (T_I^T)^-1  ->  T_V^-1 = T_I^T

            # Modal series impedance:  Z_m = T_V^-1 Z' T_I = T_I^T Z' T_I.
            # This is diagonal in exact arithmetic because  Z' t_j  is parallel
            # to the j-th column of T_V (Wedepohl). It is markedly better
            # conditioned than Y_m = T_I^-1 Y' T_V (which leans on the left
            # eigenvectors), so the modal characteristic impedance is derived
            # from Z_m and lambda rather than from Y_m directly.
            zm_full = ti.T @ self.Zs[k] @ ti
            ym_full = ti_inv @ self.Ysh[k] @ ti_inv.T
            Zm[k] = np.diag(zm_full)
            Ym[k] = np.diag(ym_full)

            dz = np.abs(zm_full - np.diag(np.diag(zm_full)))
            sz = np.max(np.abs(np.diag(zm_full))) or 1.0
            offdiag_z[k] = dz.max() / sz
            # consistency of the scalar modal relation  Z_m,j * Y_m,j = lambda_j
            cj = np.abs(Zm[k] * Ym[k] - lam[k]) / (np.abs(lam[k]) + 1e-30)
            consistency[k] = np.max(cj)

        # Physical modal propagation constant: Gamma_m = sqrt(lambda), forced to
        # the propagating branch  Re >= 0 (attenuation) and Im >= 0 (phase lag).
        g = np.sqrt(lam)
        alpha = np.abs(g.real)
        beta = np.abs(g.imag)
        gamma = alpha + 1j * beta

        with np.errstate(divide="ignore", invalid="ignore"):
            vphase = self.w[:, None] / beta
        Zcm = Zm / gamma           # modal characteristic impedance
        Ycm = 1.0 / Zcm            # modal domain is decoupled -> Y_cm = 1 / Z_cm

        from utils.passivity_check import check_pul_passivity

        passivity = check_pul_passivity(
            self.f,
            {"series_impedance_matrix": self.Zs, "shunt_admittance_matrix": self.Ysh},
        )

        self.diagnostics = {
            "max_offdiag_ratio_Zm": float(offdiag_z.max()),
            "max_scalar_relation_error": float(np.nanmax(consistency)),
            "reciprocity": "exact by construction (T_V = (T_I^T)^-1)",
            "passivity": passivity,
        }

        labels = self._classify(TI) if self.roles is not None else None

        return {
            "frequencies": self.f,
            "omega": self.w,
            "lambda": lam,
            "gamma": gamma,
            "alpha": gamma.real,
            "beta": gamma.imag,
            "vphase": vphase,
            "Zm": Zm,
            "Ym": Ym,
            "Zcm": Zcm,
            "Ycm": Ycm,
            "TI": TI,
            "TV": TV,
            "mode_labels": labels,
            "decoupling_freq": self._decoupling_frequency(),
            "diagnostics": self.diagnostics,
        }

    # ------------------------------------------------------------------ #
    # 4. mode identification                                             #
    # ------------------------------------------------------------------ #
    def _classify(self, TI):
        """Label each tracked mode by projecting its low-frequency T_I column
        onto the reference patterns (eqs. 5.34-5.36)."""
        if self.n != 6:
            warnings.warn(
                f"mode classification implemented for 6 conductors only; got {self.n}."
            )
            return None

        names_sheath = ["ground", "inter_sheath_1", "inter_sheath_2"]
        names_core = ["coaxial_1", "coaxial_2", "coaxial_3"]

        # --- stage 1: split sheath modes from coaxial modes -----------------
        # The ground / inter-sheath modes carry ~no core current at any
        # frequency (Andreata sec. 5.3, frequency-invariant); the coaxial
        # modes always involve the cores. Use the core-row energy fraction of
        # each T_I column, averaged over a mid band where every configuration
        # is decoupled (for duct configs the modes only separate above ~100 Hz).
        lo, hi = float(self.f.min()), float(self.f.max())
        band = (self.f >= max(lo, 1e2)) & (self.f <= min(hi, 1e5))
        if band.sum() < 3:
            band = np.ones_like(self.f, dtype=bool)
        kk = np.where(band)[0]

        core_frac = np.zeros(self.n)
        for k in kk:
            e = np.abs(TI[k]) ** 2
            core_frac += e[:3].sum(axis=0) / (e.sum(axis=0) + 1e-300)
        core_frac /= len(kk)

        order = np.argsort(core_frac)
        sheath_modes = sorted(order[:3].tolist())
        core_modes = sorted(order[3:].tolist())

        # --- stage 2: sub-label within each group --------------------------
        k_ref = int(kk[len(kk) // 2])
        cols = np.real(_normalise_columns(TI[k_ref]))     # (6, 6)

        labels = [None] * self.n
        quality = {}

        def _sublabel(modes, subnames, rows):
            ref = np.stack([REFERENCE_PATTERNS_6C[nm][rows] for nm in subnames], axis=1)
            ref = ref / np.linalg.norm(ref, axis=0, keepdims=True)
            block = cols[rows][:, modes]
            block = block / (np.linalg.norm(block, axis=0, keepdims=True) + 1e-30)
            sim = np.abs(ref.T @ block)                   # (3, 3)
            r_idx, c_idx = linear_sum_assignment(-sim)
            for r, c in zip(r_idx, c_idx):
                labels[modes[c]] = subnames[r]
                quality[subnames[r]] = float(sim[r, c])

        _sublabel(sheath_modes, names_sheath, slice(3, 6))
        _sublabel(core_modes, names_core, slice(0, 3))

        self.diagnostics["classification_similarity"] = quality
        self.diagnostics["classification_ref_freq_hz"] = float(self.f[k_ref])
        self.diagnostics["classification_core_fraction"] = {
            int(i): float(core_frac[i]) for i in range(self.n)
        }
        return tuple(labels)

    # ------------------------------------------------------------------ #
    # 5. decoupling frequency f_c (eq. 5.33)                             #
    # ------------------------------------------------------------------ #
    def _decoupling_frequency(self):
        """f_c = rho_sh / (pi * mu_sh * (r_out - r_in)^2) for each sheath, if
        the geometry is available on the reordered model. Returns None when the
        sheath data is not accessible from here (the caller can compute it from
        the model instead)."""
        return None
