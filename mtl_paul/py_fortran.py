import re
import subprocess
import numpy as np
from pathlib import Path

class FortranRunner:
    """
    A class to control the execution of a legacy Fortran program,
    managing its input and output files.
    """
    def __init__(self, exe_path: str, input_filename: str = 'RIBBON.IN', output_filename: str = 'PUL.DAT', silent: bool = True):
        """
        Initializes the Fortran runner.

        Args:
            exe_path (str): The path to the Fortran executable file.
            input_filename (str): The name of the input file the executable expects.
            output_filename (str): The name of the output file the executable generates.
        """
        self.exe_path = Path(exe_path).resolve()
        self.input_filename = self.exe_path.parent / input_filename
        self.output_filename = self.exe_path.parent / output_filename
        self.silent = silent

        if not self.exe_path.is_file():
            raise FileNotFoundError(f"The executable file was not found at: {self.exe_path}")

        # Initialize every matrix attribute to reflect the final naming
        self.A_matrix = None
        self.B_matrix = None
        self.C_matrix = None
        self.D_matrix = None
        self.CGEN0_matrix = None
        self.IND_matrix = None
        self.CAP_matrix = None
        self.CAP0_matrix = None
        self.CGEN_matrix = None
        self.NF = None

    def _write_input_file(self, params: dict):
        """
        Writes the simulation parameters to the Fortran input file.
        """
        if not self.silent:
            print(f"Writing input file to: {self.input_filename}")

        param_order = ['N', 'NF', 'IREF', 'RW', 'TD', 'ER', 'S']

        with open(self.input_filename, 'w') as f:
            for key in param_order:
                if key not in params:
                    raise ValueError(f"Required parameter '{key}' not found in the parameter dictionary.")
                f.write(f"{params[key]}\n")

    def _run_executable(self):
        """
        Runs the compiled Fortran program.
        """
        if not self.silent:
            print(f"Running: {self.exe_path}...")
        try:
            result = subprocess.run(
                [str(self.exe_path)],
                cwd=self.exe_path.parent,
                check=True,
                capture_output=True,
                text=True,
                timeout=30
            )
            if not self.silent:
                print("Fortran execution completed successfully.")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error during execution of the Fortran program.")
            print(f"Exit code: {e.returncode}")
            print(f"Error output (stderr):\n{e.stderr}")
            return False
        except Exception as e:
            print(f"An unexpected error occurred while trying to run Fortran: {e}")
            return False

    def _parse_output_file(self):
        """
        Reads and parses the output file to extract the result matrices.
        """
        if not self.silent:
            print(f"Parsing output file: {self.output_filename}")

        if not self.output_filename.is_file():
            print(f"Error: Output file '{self.output_filename}' was not generated.")
            return

        raw_values = {}
        max_indices = {}

        pattern = re.compile(r"(\d+)\s+(\d+)\s+([0-9.E+-]+)\s+=\s*([A-Z0-9]+)\(")

        with open(self.output_filename, 'r') as f:
            for line in f:
                match = pattern.match(line.strip())
                if match:
                    i_str, j_str, val_str, key = match.groups()
                    i, j, value = int(i_str), int(j_str), float(val_str)

                    if key not in raw_values:
                        raw_values[key] = {}
                        max_indices[key] = 0

                    raw_values[key][(i, j)] = value
                    max_indices[key] = max(max_indices[key], i, j)

        # Define which matrices are symmetric based on the names in the PUL.DAT file
        symmetric_keys = ['IND', 'CAP', 'CAP0', 'CGEN', 'CGEN0']

        for key, values in raw_values.items():
            size = max_indices.get(key, 0)
            if size > 0:
                matrix = np.zeros((size, size))
                is_symmetric = key in symmetric_keys

                for (i, j), val in values.items():
                    matrix[i-1, j-1] = val
                    if is_symmetric:
                        matrix[j-1, i-1] = val

                setattr(self, f"{key}_matrix", matrix)

        if not self.silent:
            print("Output file parsing completed.")

    def run_fortran(self, params: dict):
        """
        Orchestrates the full process: writes the input, runs and parses the output.
        """
        try:
            self._write_input_file(params)
            if self._run_executable():
                self._parse_output_file()
            else:
                if not self.silent:
                    print("Simulation aborted due to an execution error.")
        except Exception as e:
            print(f"An unexpected error occurred during the simulation: {e}")
