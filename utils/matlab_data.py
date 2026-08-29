import sys
import os
import numpy as np
from pathlib import Path
from typing import Sequence

import scipy.io as sio


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
