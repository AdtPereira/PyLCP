# py_tulip.py
import json
import subprocess
import numpy as np
from pathlib import Path

class PyTulip:
    """
    A class to run pre-existing pulmtln.exe simulation cases,
    assuming a standardized directory structure.
    """
    def __init__(self, executable_path: str, case_name: str, silent: bool = True):
        """
        Initializes the PyTulip runner for a specific case.

        Args:
            executable_path (str): The path to the pulmtln.exe executable.
            case_name (str): The name of the case folder inside C:\\git\\tulip\\examples\\.
            silent (bool): If True, suppresses most print statements.
        """
        base_path = Path(r"C:\git\tulip\examples")
        
        self.executable_path = Path(executable_path).resolve()
        self.case_name = case_name
        self.working_directory = (base_path / self.case_name).resolve()
        self.input_filepath = (self.working_directory / f"{self.case_name}.pulmtln.in.json").resolve()
        self.output_filepath = (self.working_directory / "pulmtln.out.json").resolve()
        self.silent = silent

        # --- Initial Sanity Checks ---
        if not self.executable_path.is_file():
            raise FileNotFoundError(f"Executable not found at: {self.executable_path}")
        if not self.working_directory.is_dir():
            raise NotADirectoryError(f"Case directory not found at: {self.working_directory}")
        if not self.input_filepath.is_file():
            raise FileNotFoundError(f"Input JSON file not found at: {self.input_filepath}")

        # Initialize result attributes to None
        self.C_matrix = None
        self.L_matrix = None

    def _run_executable(self):
        """
        Executes pulmtln.exe and streams its output to the console in real-time.
        """
        if not self.silent:
            print(f"Executing: {self.executable_path.name}...")
            print(f"Working Directory: {self.working_directory}")
            print(f"Input File: {self.input_filepath.name}\n")
        
        command = [str(self.executable_path), "-i", self.input_filepath.name]
        
        try:
            # Use subprocess.Popen to run the command and get access to its output streams.
            # stdout=subprocess.PIPE tells Popen to capture the output.
            # stderr=subprocess.STDOUT merges the error stream into the standard output.
            # text=True and encoding='utf-8' ensures the output is decoded as text.
            with subprocess.Popen(
                command,
                cwd=self.working_directory,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding='utf-8'
            ) as process:
                # Read the output line by line from the process's stdout stream in real-time.
                if not self.silent:
                    for line in process.stdout:
                        print(line, end='', flush=True)

                # Wait for the process to complete and get the return code.
                process.wait()
                
                if process.returncode == 0:
                    if not self.silent:
                        print("\npulmtln.exe execution finished successfully.")
                    return True
                else:
                    print(f"\nError during pulmtln.exe execution for case '{self.case_name}'.")
                    print(f"Return Code: {process.returncode}")
                    return False

        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return False

    def _parse_output_json(self):
        """
        Reads and parses the JSON output file to extract the result matrices.
        """
        if not self.silent:
            print(f"\nParsing output file: {self.output_filepath.name}")

        if not self.output_filepath.is_file():
            print(f"Error: Output file '{self.output_filepath.name}' was not generated.")
            return

        try:
            with open(self.output_filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)

            if "C" in data:
                self.C_matrix = np.array(data["C"])
            if "L" in data:
                self.L_matrix = np.array(data["L"])
            
            if not self.silent:
                print("Output file parsed successfully.")
        except Exception as e:
            print(f"Error while parsing output file: {e}")
            
    def run(self):
        """
        Orchestrates the entire process: runs executable and parses output.
        """
        if self._run_executable():
            self._parse_output_json()
        else:
            if not self.silent:
                print("\nSimulation aborted due to an error during execution.")


if __name__ == '__main__':
    # --- Example Usage ---
    EXE_PATH = r"C:\git\tulip\pulmtln-build\rls\bin\Release\pulmtln.exe"
    CASE_NAME = "two_wires_open_gmsh" 

    print(f"--- Preparing to run simulation for case: '{CASE_NAME}' ---")
    
    try:
        runner = PyTulip(
            executable_path=EXE_PATH,
            case_name=CASE_NAME,
            silent=False
        )
        runner.run()

        print("\n--- Accessing Results ---")
        if runner.C_matrix is not None:
            print("\nCapacitance Matrix (C):")
            print(runner.C_matrix)
        
        if runner.L_matrix is not None:
            print("\nInductance Matrix (L):")
            print(runner.L_matrix)

    except (FileNotFoundError, NotADirectoryError) as e:
        print(f"\nConfiguration Error: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")