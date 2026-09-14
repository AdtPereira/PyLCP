import sys
import os
import numpy as np
from pathlib import Path
from typing import Sequence

import scipy.io as sio

# Mode order used by the external MATLAB developer for the modal-domain
# propagation parameters (attenuation, phase, velocity, characteristic
# impedance/admittance -- see MatlabDataReader.get_modal_scenario_data).
# Confirmed (not merely assumed) by matching the per-mode |Z_cm|/alpha_m
# magnitudes against the validation table in
# testData/andreata_common/MODAL_CH5_DEVELOPMENT.md: it is Andreata's own
# thesis order (eqs. 5.34-5.36 / Figs. 5.5-5.7 legend), identical to
# analytical_forms.modal_analysis.MODE_LABELS_6C.
MODAL_MODE_ORDER_6C = (
    "ground", "inter_sheath_1", "inter_sheath_2",
    "coaxial_1", "coaxial_2", "coaxial_3",
)


class MatlabDataReader:
    """
    A dedicated class to read and parse all MATLAB .mat files
    from the 'Results' directory of a specific case, analogously to
    ComsolDataReader (utils/comsol_data.py).

    It automatically discovers and loads every .mat file, returning them
    in a structured dictionary.
    """

    def __init__(self, script_file_path: str, autoShow: bool = True):
        """
        Initializes the reader by identifying the target 'Results' directory.

        Args:
            script_file_path (str): The __file__ of the calling script
                                     (used to locate the case directory).
            autoShow (bool): If True, displays a summary of the loaded data.
        """
        project_root = Path(script_file_path).resolve().parents[2]
        if str(project_root) not in sys.path:
            sys.path.insert(0, str(project_root))

        case_name = os.path.splitext(os.path.basename(script_file_path))[0]

        self.project_root = project_root
        self.case_name = case_name
        self.results_path = project_root / 'testData' / case_name / 'Results'

        print(f"Project root configured at: {project_root}")
        print(f"\nInstantiating MatlabDataReader for case '{case_name}' ---")

        if not self.results_path.is_dir():
            print(f"  Warning: Results directory not found for case '{case_name}'. MATLAB data ignored.")

        self.data = self.load_all_results()

        if autoShow and self.data:
            self.show_summary()

    def load_all_results(self) -> dict:
        """
        Scans the 'Results' directory, loads every .mat file and returns them
        as a dictionary.

        Returns:
            A dictionary where the keys are the file names (without the .mat
            extension). The value is the variable's NumPy array when the file
            contains a single data variable (the common case), or a dictionary
            {name: array} when it contains more than one.
        """
        if not self.results_path.is_dir():
            return {}

        mat_files = list(self.results_path.glob('*.mat'))
        data = {}

        if not mat_files:
            print(f"Warning: No .mat file found in {self.results_path}")
            return {}

        print(f"Found {len(mat_files)} .mat file(s) in the Results directory of case '{self.case_name}'.")

        for file_path in mat_files:
            file_stem = file_path.stem
            try:
                data[file_stem] = self._parse_single_file(file_path)
            except Exception as e:
                print(f"Error while parsing file {file_path.name}: {e}")

        return data

    def _parse_single_file(self, file_path: Path):
        """
        Loads a single .mat file (MATLAB level 5 format, via scipy.io.loadmat)
        and returns its variables, ignoring the MATLAB metadata keys
        ('__header__', '__version__', '__globals__').

        Square 3D matrices exported in the MATLAB (N, N, num_freq) format are
        reordered to (num_freq, N, N), the convention used by the other pyLCP
        PUL matrices (see, e.g., InternalPerUnitParameters).
        """
        print(f"  -> Loading and parsing: {file_path.name}...")

        raw = sio.loadmat(file_path)
        variables = {}
        for name, value in raw.items():
            if name.startswith('__'):
                continue
            if isinstance(value, np.ndarray) and value.ndim == 3 and value.shape[0] == value.shape[1]:
                value = np.moveaxis(value, -1, 0)
            variables[name] = value

        if len(variables) == 1:
            return next(iter(variables.values()))
        return variables

    def show_summary(self):
        """
        Displays a summary of every loaded dataset (variable names, shapes
        and dtypes).
        """
        if not self.data:
            print("No data loaded to display the summary. Run 'load_all_results()' first.")
            return

        print(f"\n--- Summary of the Loaded Data for Case '{self.case_name}' ---")
        for name, value in self.data.items():
            print(f"\n==================================================")
            print(f"  File: '{name}.mat'")
            print(f"==================================================")

            variables = value if isinstance(value, dict) else {name: value}
            for var_name, arr in variables.items():
                shape = getattr(arr, 'shape', None)
                dtype = getattr(arr, 'dtype', type(arr))
                print(f"  {var_name}: shape={shape}, dtype={dtype}")

    @staticmethod
    def _reorder_conductor_matrix(matrix, conductor_order: Sequence[int]):
        """
        MATLAB exports the conductors grouped by type: [core_A, core_B, core_C,
        sheath_A, sheath_B, sheath_C]. pyLCP assembles its matrices (internal and
        quasi-TEM) grouped by cable: [core_A, sheath_A, core_B, sheath_B,
        core_C, sheath_C] (see np.kron(np.identity(N), Zij) in
        InternalPerUnitParameters.matrices()). Without this reordering, M[p,q]
        (pyLCP) and Z[p,q] (MATLAB) point to physically different conductor
        pairs for the same indices (p, q).
        """
        if matrix is None:
            return None
        return matrix[:, conductor_order, :][:, :, conductor_order]

    def get_scc_scenario_data(self, prefix: str, conductor_order: Sequence[int]) -> dict:
        """
        Assembles the MATLAB reference-data dictionary (in the format used by
        pul_data['matlab']) for an SCC case (core + sheath), from the files
        exported with the prefix `prefix` (e.g. 'andreata' ->
        'andreata_frequency_range', 'andreata_series_impedance_matrix', ...).

        `conductor_order` reorders the conductors from the MATLAB convention
        (grouped by type) to the pyLCP convention (grouped by cable); see
        `_reorder_conductor_matrix`. The two ground-return matrices (Zg, Pg)
        follow the same type-grouped layout as the others -- core and sheath of
        the same cable have redundant (identical) entries, since the ground
        return only depends on the cable position, not on which conductor inside
        it -- so they use the same `_reorder_conductor_matrix`, with no special
        slicing.

        Returns 'frequencies' as None if the corresponding file is not found --
        the fallback (e.g. to pul_data['frequencies']) is up to the caller.
        """
        frequencies = self.data.get(f'{prefix}_frequency_range')

        internal_z = self._reorder_conductor_matrix(self.data.get(f'{prefix}_internal_impedance_matrix'), conductor_order)
        series_z = self._reorder_conductor_matrix(self.data.get(f'{prefix}_series_impedance_matrix'), conductor_order)
        shunt_y = self._reorder_conductor_matrix(self.data.get(f'{prefix}_shunt_admittance_matrix'), conductor_order)
        internal_y = self._reorder_conductor_matrix(self.data.get(f'{prefix}_internal_admittance_matrix'), conductor_order)
        earth_return_pg = self._reorder_conductor_matrix(self.data.get(f'{prefix}_earth_return_potential_coefficient_matrix'), conductor_order)
        earth_return_zg = self._reorder_conductor_matrix(self.data.get(f'{prefix}_earth_return_impedance_matrix'), conductor_order)

        return {
            'frequencies': frequencies.flatten() if frequencies is not None else None,
            'scenarios': {
                'measured': {
                    'internal_impedance_matrix': internal_z,
                    'series_impedance_matrix': series_z,
                    'shunt_admittance_matrix': shunt_y,
                    'internal_admittance_matrix': internal_y,
                    'earth_return_potential_coefficient_matrix': earth_return_pg,
                    'earth_return_impedance_matrix': earth_return_zg,
                },
            },
        }

    def get_modal_scenario_data(self, prefix: str, config_index: int = 1) -> dict:
        """
        Assembles the MATLAB reference for the modal-domain propagation
        parameters of Andreata Ch. 5 (attenuation, phase constant, phase
        velocity, characteristic impedance/admittance -- one value per mode),
        from the files exported by the external developer:

            <prefix>_alpham.mat, _betam.mat, _velocm.mat, _Zcm.mat, _Ycm.mat

        Each file holds one variable per mode, named
        ``<var><mode>_<config_index>`` with ``mode`` = 1..N (e.g.
        ``alpham1_1`` .. ``alpham6_1`` for Configuration 1's 6 conductors,
        ``..._3`` for Configuration 3's 7 -- 6 + the ECC). ``N`` is inferred
        per file from the variables actually present (not hardcoded), so this
        also serves the 7-conductor ECC cases (3, 4) once their .mat files
        arrive, without silently dropping the 7th mode.

        Mode *labels* are only assigned when ``N == 6``: that is the one case
        confirmed (by matching magnitudes against
        ``MODAL_CH5_DEVELOPMENT.md``) to follow Andreata's thesis order
        ``MODAL_MODE_ORDER_6C``. For any other ``N`` (the ECC configs) mode
        classification is not implemented yet on the pyLCP side either (see
        ``analytical_forms.modal_analysis.ModalDecomposition._classify``, which
        bails out for ``n != 6``) -- there is no reference-pattern scheme to
        assign real labels to, and no thesis equations for a 7th "ECC mode" in
        the codebase yet. Labelling as ``None`` (rather than guessing e.g.
        raw index order) is deliberate: it makes downstream consumers
        (``ModalPropagationPlotter``, ``print_modal_comparison_report``) skip
        the mode-matching overlay/report instead of silently pairing up two
        arbitrarily-ordered eigenmode sets that have no verified correspondence.

        Frequencies are not stored inside these files -- the developer used
        the "standard" 91-point log grid (10 pts/decade, 1E-2..1E7 Hz), the
        same one used for the COMSOL FEM export
        (``<prefix>_frequency_range_fem.mat``), reused here if present, else
        regenerated with ``np.logspace(-2, 7, 91)``.

        Returns ``None`` if none of the 5 files are present -- the caller
        (e.g. ``ModalPropagationPlotter``) should treat that as "no MATLAB
        overlay available" rather than an error.
        """
        import re

        var_prefixes = {
            'alpha': 'alpham', 'beta': 'betam', 'vphase': 'velocm',
            'Zcm': 'Zcm', 'Ycm': 'Ycm',
        }

        out = {}
        n_modes = None
        for out_key, var_prefix in var_prefixes.items():
            raw = self.data.get(f'{prefix}_{var_prefix}')
            if not isinstance(raw, dict):
                continue

            suffix = f'_{config_index}'
            pattern = re.compile(rf'^{re.escape(var_prefix)}(\d+){re.escape(suffix)}$')
            modes_found = sorted(
                int(m.group(1)) for k in raw if (m := pattern.match(k)) is not None)
            if not modes_found or modes_found != list(range(1, len(modes_found) + 1)):
                print(f"  Warning: '{prefix}_{var_prefix}.mat' has no contiguous "
                      f"1..N '{var_prefix}<mode>{suffix}' variables -- skipping '{out_key}'.")
                continue
            n = len(modes_found)
            if n_modes is None:
                n_modes = n
            elif n != n_modes:
                print(f"  Warning: '{prefix}_{var_prefix}.mat' has {n} modes, "
                      f"but a previous file in this scenario had {n_modes} -- "
                      f"skipping '{out_key}'.")
                continue

            cols = [np.ravel(raw[f'{var_prefix}{m}{suffix}']) for m in modes_found]
            out[out_key] = np.stack(cols, axis=1)  # (Nf, n)

        if not out:
            return None

        if n_modes == len(MODAL_MODE_ORDER_6C):
            mode_labels = MODAL_MODE_ORDER_6C
        else:
            print(f"  Note: '{prefix}' modal reference has {n_modes} modes (config_index="
                  f"{config_index}); mode classification is only implemented for the "
                  f"6-conductor configs, so these modes are left unlabelled -- overlay "
                  f"plots/reports will skip them until a labelling scheme is added.")
            mode_labels = None

        freq = self.data.get(f'{prefix}_frequency_range_fem')
        if freq is not None:
            frequencies = np.ravel(freq)
        else:
            print(f"  Warning: '{prefix}_frequency_range_fem.mat' not found; assuming "
                  f"the standard grid np.logspace(-2, 7, 91) for the modal reference "
                  f"of '{prefix}' (config_index={config_index}).")
            frequencies = np.logspace(-2, 7, num=91)

        out['frequencies'] = frequencies
        out['mode_labels'] = mode_labels
        return out


# --------------------------------------------------------------------------- #
# COMSOL (FEM) -> MATLAB reference export
# --------------------------------------------------------------------------- #
#
# The MATLAB developer validates against reference .mat files
# (<prefix>_internal_impedance_matrix.mat / _internal_admittance_matrix.mat)
# whose exact structure is:
#
#   * a single variable, 'Z' (impedance) or 'Y' (admittance);
#   * shape (6, 6, Nf) -- conductor x conductor x frequency (Nf = 90 here);
#   * dtype complex128;
#   * conductors grouped BY TYPE: [core_A, core_B, core_C, sheath_A, sheath_B,
#     sheath_C] (indices 0-2 = cores, 3-5 = sheaths) -- the MATLAB convention,
#     the inverse of the pyLCP "grouped by cable" layout (see
#     MatlabDataReader._reorder_conductor_matrix and andreata_case1/BUGS_AND_FIXES.md,
#     bug 4).
#
# `export_comsol_internal_matrices_to_mat` reads the COMSOL results of the case
# (the 2x2 core+sheath matrices already parsed by ComsolPostProcessor) and
# writes '<prefix>_internal_impedance_matrix_fem.mat' /
# '<prefix>_internal_admittance_matrix_fem.mat' with that *identical* structure,
# plus '<prefix>_frequency_range_fem.mat' (the COMSOL frequency vector, same
# structure as '<prefix>_frequency_range.mat'), so the developer can swap the
# FEM data in for the reference with zero impact on notation or element ordering.
#
# NO numerical interpolation / resampling is performed: the COMSOL samples are
# written verbatim. The COMSOL frequency sweep must therefore already match the
# reference grid (same number of points, same values); otherwise the export
# aborts with a ValueError -- unless strict_reference=False, which downgrades
# the reference-grid check to a warning and exports on the COMSOL grid (the two
# COMSOL inputs must still agree with each other). Emergency use for an
# off-standard reference (e.g. 90 pts instead of the standard 91 = 10
# pts/decade over 1E-2..1E7 Hz) pending its correction by the external developer.


def _cable_block_to_type_ordered_system(block_stack: np.ndarray,
                                        conductor_order: Sequence[int]) -> np.ndarray:
    """
    Expands a per-cable 2x2 internal matrix (index 0 = core, 1 = sheath),
    identical for the 3 cables of a symmetric SCC system, into the full 6x6
    system matrix in the MATLAB *by type* ordering.

    block_stack : (Nf, 2, 2) complex
    conductor_order : the same permutation passed to
        MatlabDataReader.get_scc_scenario_data (MATLAB by-type -> pyLCP by-cable,
        e.g. [0, 3, 1, 4, 2, 5]); its inverse takes us back from by-cable to
        by-type.

    Returns (Nf, 6, 6) complex.
    """
    block_stack = np.asarray(block_stack, dtype=complex)
    nf = block_stack.shape[0]

    # Step 1: assemble in pyLCP "by cable" order [cA, sA, cB, sB, cC, sC]
    # -- block-diagonal replication of the 2x2 block (np.kron(I3, block)),
    # matching np.kron(np.identity(N), Zij) in InternalPerUnitParameters.
    by_cable = np.zeros((nf, 6, 6), dtype=complex)
    for c in range(3):
        by_cable[:, 2 * c:2 * c + 2, 2 * c:2 * c + 2] = block_stack

    # Step 2: reorder by-cable -> by-type with the inverse permutation.
    inv_perm = np.argsort(np.asarray(conductor_order))
    return by_cable[:, inv_perm][:, :, inv_perm]


def _assert_frequency_compatible(source_freq: np.ndarray, target_freq: np.ndarray,
                                 what: str, rtol: float = 1e-6) -> None:
    """
    Guards the "no interpolation" contract: `source_freq` and `target_freq`
    must match both in length and in values. `what` names the two grids being
    compared, for the error message. Raises ``ValueError`` otherwise.
    """
    source_freq = np.asarray(source_freq, dtype=float)
    target_freq = np.asarray(target_freq, dtype=float)

    if source_freq.shape != target_freq.shape:
        raise ValueError(
            f"Dimensional incompatibility ({what}): {source_freq.size} vs "
            f"{target_freq.size} frequency point(s). Interpolation is disabled, so "
            f"the two grids must have the same number of points. NB: the standard "
            f"log sweep 1E-2..1E7 Hz at 10 pts/decade has 91 points "
            f"(np.logspace(-2, 7, 91)); a 90-point reference is off-standard. "
            f"Pass strict_reference=False to export on the COMSOL grid regardless."
        )

    if not np.allclose(source_freq, target_freq, rtol=rtol, atol=0.0):
        max_rel = float(np.max(np.abs(source_freq - target_freq)
                               / np.where(target_freq != 0, np.abs(target_freq), 1.0)))
        raise ValueError(
            f"Frequency-grid value mismatch ({what}): same length "
            f"({source_freq.size}) but values differ (max relative deviation "
            f"{max_rel:.2e} > rtol {rtol:.0e}). Interpolation is disabled, so the "
            f"grids must hold the same frequencies. Pass strict_reference=False to "
            f"export on the COMSOL grid regardless."
        )


def export_comsol_internal_matrices_to_mat(
    script_file_path: str,
    prefix: str,
    conductor_order: Sequence[int] = (0, 3, 1, 4, 2, 5),
    target_frequencies: np.ndarray = None,
    output_dir: str = None,
    reference_check: bool = True,
    strict_reference: bool = True,
) -> dict:
    """
    Reads the COMSOL (FEM) internal impedance/admittance results of an SCC case
    and writes them as MATLAB .mat files structurally identical to the
    developer's reference files, so the FEM data can replace the reference with
    no impact on notation or on the ordering of the internal elements.

    Reads (from ``testData/<case>/Results/``):
      * ``cmsl_internal_impedance_matrix.txt``            (Js method, core+sheath)
      * ``cmsl_internal_admittance_charge_method.txt``    (direct charge method)
    both via :class:`utils.comsol_data.ComsolPostProcessor`, which already
    resolves every column-naming / excitation-block subtlety.

    Writes (same directory, unless ``output_dir`` is given):
      * ``<prefix>_internal_impedance_matrix_fem.mat``  -- variable ``Z``
      * ``<prefix>_internal_admittance_matrix_fem.mat`` -- variable ``Y``
        each with the exact structure of ``<prefix>_internal_impedance_matrix.mat`` /
        ``<prefix>_internal_admittance_matrix.mat``: shape ``(6, 6, Nf)``, dtype
        ``complex128``, conductors grouped **by type**
        ``[core_A, core_B, core_C, sheath_A, sheath_B, sheath_C]``.
      * ``<prefix>_frequency_range_fem.mat`` -- the COMSOL frequency vector
        (Hz), same structure as ``<prefix>_frequency_range.mat``: a single
        ``(1, Nf)`` real row vector, under the same variable name (``freq1``).

    Parameters
    ----------
    script_file_path : str
        ``__file__`` of the case script (e.g. ``testData/andreata_case1/andreata_case1.py``),
        used to locate the ``Results`` directory -- same convention as the reader classes.
    prefix : str
        File-name prefix of the reference exports (e.g. ``'andreata_case1'``).
    conductor_order : sequence of int
        MATLAB-by-type -> pyLCP-by-cable permutation, as passed to
        :meth:`MatlabDataReader.get_scc_scenario_data`. Default ``(0, 3, 1, 4, 2, 5)``.
    target_frequencies : array-like, optional
        Reference frequency grid (Hz) the COMSOL data must already match.
        Default: the vector stored in ``<prefix>_frequency_range.mat``; if that
        file is absent, the COMSOL impedance grid itself is used.
        **No interpolation is performed** -- the COMSOL samples are written
        verbatim.
    output_dir : str, optional
        Destination directory. Default: the case ``Results`` directory.
    reference_check : bool
        If True (default) and the reference .mat files exist, prints the
        FEM-vs-reference relative difference of the DC diagonal (a sanity check,
        not a pass/fail -- FEM and the analytical reference are expected to
        differ by a fraction of a percent).
    strict_reference : bool
        The two COMSOL inputs (impedance + admittance) must **always** share one
        grid -- that is enforced unconditionally. ``strict_reference`` only
        governs the check of the COMSOL grid against ``target_frequencies`` /
        ``<prefix>_frequency_range.mat``:

        * ``True`` (default) -- a mismatch raises ``ValueError``.
        * ``False`` -- a mismatch is downgraded to a warning and the export
          proceeds **on the COMSOL grid**. Emergency use for when the external
          reference .mat files are still on an off-standard grid (e.g. 90 pts
          instead of the standard 91 = 10 pts/decade over 1E-2..1E7 Hz) while a
          correct 91-point COMSOL export is already available. The companion
          ``<prefix>_frequency_range_fem.mat`` then documents the real grid; the
          reference .mat files must be regenerated on it before an
          element-by-element comparison is meaningful.

    Raises
    ------
    FileNotFoundError
        If the COMSOL ``.txt`` inputs are missing.
    ValueError
        If the two COMSOL inputs are on different grids, or (unless
        ``strict_reference=False``) if the COMSOL grid does not match the
        reference grid. Interpolation is disabled.

    Returns
    -------
    dict
        ``{'impedance': {'path', 'Z', 'frequencies'},
           'admittance': {'path', 'Y', 'frequencies'},
           'frequency_range': {'path', 'frequencies'}}``.
    """
    from utils.comsol_data import ComsolPostProcessor

    results_path = Path(script_file_path).resolve().parent / 'Results'
    out_path = Path(output_dir) if output_dir is not None else results_path
    out_path.mkdir(parents=True, exist_ok=True)

    # --- reference frequency grid (the reference .mat grid) ----------------
    freq_var_name = 'freq1'  # variable name inside <prefix>_frequency_range[_fem].mat
    if target_frequencies is not None:
        target_frequencies = np.asarray(target_frequencies, dtype=float)
    else:
        freq_file = results_path / f'{prefix}_frequency_range.mat'
        if freq_file.is_file():
            freq_var = sio.loadmat(freq_file)
            freq_var_name, ref_freq = next(
                (k, v) for k, v in freq_var.items() if not k.startswith('__'))
            target_frequencies = np.ravel(ref_freq).astype(float)
        else:
            print(f"  Warning: '{freq_file.name}' not found; the COMSOL impedance "
                  f"grid will be used as the reference for the compatibility check.")

    # --- COMSOL data (2x2 core+sheath, pyLCP by-cable convention) -----------
    cmsl = ComsolPostProcessor(script_file_path, autoShow=False)

    z_block = cmsl.get_scc_internal_impedance_matrix_combined()
    c_block = cmsl.get_scc_internal_capacitance_matrix_combined()
    if z_block is None or c_block is None:
        raise FileNotFoundError(
            "COMSOL internal impedance/admittance .txt files not found for case "
            f"'{Path(script_file_path).stem}'.")

    zi_2x2 = z_block['scenarios']['measured']['impedance_matrix']          # (Nf, 2, 2)
    ci_2x2 = c_block['scenarios']['measured']['capacitance_matrix'].astype(complex)  # (Nf, 2, 2)
    freq_z = np.asarray(z_block['frequencies'], dtype=float)
    freq_c = np.asarray(c_block['frequencies'], dtype=float)

    # --- frequency-grid compatibility (interpolation disabled) ------------
    # (1) The two COMSOL inputs must share one grid -- always enforced.
    _assert_frequency_compatible(
        freq_c, freq_z,
        "'cmsl_internal_admittance_charge_method.txt' vs 'cmsl_internal_impedance_matrix.txt'")

    # (2) The COMSOL grid vs the reference grid -- gated by strict_reference.
    if target_frequencies is None:
        target_frequencies = freq_z
    try:
        _assert_frequency_compatible(
            freq_z, target_frequencies, "COMSOL sweep vs reference grid")
    except ValueError as exc:
        if strict_reference:
            raise
        print(f"  Warning: {exc}")
        print(f"  strict_reference=False -> exporting on the COMSOL grid "
              f"({freq_z.size} pts). Regenerate the reference .mat files on this "
              f"grid (see '{prefix}_frequency_range_fem.mat') before comparing "
              f"element by element.")
        target_frequencies = freq_z

    # --- expand 2x2 -> 6x6 (by type); COMSOL samples kept verbatim ---------
    Z = _cable_block_to_type_ordered_system(zi_2x2, conductor_order)       # (Nf, 6, 6)
    C = _cable_block_to_type_ordered_system(ci_2x2, conductor_order).real  # (Nf, 6, 6)
    Y = np.zeros_like(C, dtype=complex)
    Y.imag = (2 * np.pi * freq_c)[:, None, None] * C                      # Y = jw C, exact +0 real part

    # --- lay out as (6, 6, Nf) complex128, exactly like the reference ------
    Z_out = np.ascontiguousarray(np.moveaxis(Z, 0, -1)).astype(np.complex128)
    Y_out = np.ascontiguousarray(np.moveaxis(Y, 0, -1)).astype(np.complex128)

    # COMSOL frequency vector, same structure as <prefix>_frequency_range.mat
    # (a single (1, Nf) real row vector under the same variable name).
    f_row = np.asarray(freq_z, dtype=float).reshape(1, -1)

    z_out_file = out_path / f'{prefix}_internal_impedance_matrix_fem.mat'
    y_out_file = out_path / f'{prefix}_internal_admittance_matrix_fem.mat'
    f_out_file = out_path / f'{prefix}_frequency_range_fem.mat'
    sio.savemat(z_out_file, {'Z': Z_out})
    sio.savemat(y_out_file, {'Y': Y_out})
    sio.savemat(f_out_file, {freq_var_name: f_row})

    print(f"COMSOL internal impedance : {freq_z.size} pts (verbatim, no interpolation)  ->  "
          f"{z_out_file.name}  {Z_out.shape} {Z_out.dtype}")
    print(f"COMSOL internal admittance: {freq_c.size} pts (verbatim, no interpolation)  ->  "
          f"{y_out_file.name}  {Y_out.shape} {Y_out.dtype}")
    print(f"COMSOL frequency vector   : {freq_z.size} pts  ->  "
          f"{f_out_file.name}  ('{freq_var_name}', {f_row.shape} {f_row.dtype})")

    if reference_check:
        _print_reference_check(results_path, prefix, Z_out, Y_out)

    return {
        'impedance': {'path': z_out_file, 'Z': Z_out, 'frequencies': freq_z},
        'admittance': {'path': y_out_file, 'Y': Y_out, 'frequencies': freq_c},
        'frequency_range': {'path': f_out_file, 'frequencies': freq_z},
    }


def _print_reference_check(results_path: Path, prefix: str,
                           Z_out: np.ndarray, Y_out: np.ndarray) -> None:
    """Prints the FEM-vs-reference DC-diagonal relative difference, if the
    reference .mat files are available (informational, never fatal)."""
    for var, out, tag in (('Z', Z_out, 'impedance'), ('Y', Y_out, 'admittance')):
        ref_file = results_path / f'{prefix}_internal_{tag}_matrix.mat'
        if not ref_file.is_file():
            continue
        ref = sio.loadmat(ref_file)[var]
        if ref.shape != out.shape:
            print(f"  Note: reference '{ref_file.name}' has shape {ref.shape}, "
                  f"FEM output has {out.shape} -- skipping numeric check.")
            continue
        d_ref = np.diag(ref[:, :, 0])
        d_out = np.diag(out[:, :, 0])
        denom = np.where(np.abs(d_ref) > 0, np.abs(d_ref), 1.0)
        rel = np.abs(d_out - d_ref) / denom
        print(f"  [{var}] DC-diagonal FEM vs reference: "
              f"max rel. diff = {rel.max():.3%}  (per conductor: "
              f"{', '.join(f'{r:.2%}' for r in rel)})")
