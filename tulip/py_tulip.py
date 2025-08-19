# py_tulip.py
import sys
import json
import gmsh
import subprocess
import numpy as np
from pathlib import Path


class GmshMeshGenerator:
    """
    A dedicated class to generate .msh files using Gmsh.
    """
    def __init__(self, case_name: str, output_dir: Path, silent: bool = True, overwrite: bool = True):
        """
        Initializes the mesh generator.

        Args:
            case_name (str): The name of the case, used for the model and output file.
            output_dir (Path): The directory where the .msh file will be saved.
            silent (bool): If True, suppresses the Gmsh GUI.
            overwrite (bool): If True, will generate a new mesh even if one exists.
        """
        self.case_name = case_name
        self.output_dir = output_dir
        self.output_filepath = (output_dir / f"{case_name}.msh").resolve()
        self.silent = silent
        self.overwrite = overwrite

    def create_two_wires_mesh(self, rw: float = 0.010, s: float = 0.021, msize_conductors: float = 0.0005, msize_boundary: float = 0.05):
        """
        Creates the specific mesh for the 'two_wires_open' case.
        This method is an adaptation of the provided script.
        Args:
            rw (float): Radius of the wires (mm).
            s (float): Distance between the wires (mm).
            msize_conductors (float): Mesh size for the conductors.
            msize_boundary (float): Mesh size for the boundary.
        """
         # If the mesh file already exists and overwrite is False, do nothing.
        if self.output_filepath.is_file() and not self.overwrite:
            if not self.silent:
                print(f"Mesh file '{self.output_filepath.name}' already exists. Skipping generation.")
            return

        if not self.silent:
            print(f"Generating mesh file: {self.output_filepath.name}")

        gmsh.initialize()
        try:
            gmsh.model.add(self.case_name)
            conductor_1 = gmsh.model.occ.addCircle(+s/2, 0, 0, rw)
            conductor_2 = gmsh.model.occ.addCircle(-s/2, 0, 0, rw)
            boundary = gmsh.model.occ.addCircle(0, 0, 0, 8*s)

            cl1 = gmsh.model.occ.addCurveLoop([conductor_1])
            cl2 = gmsh.model.occ.addCurveLoop([conductor_2])
            cl_boundary = gmsh.model.occ.addCurveLoop([boundary])

            s_vacuum = gmsh.model.occ.addPlaneSurface([cl_boundary, cl1, cl2])
            gmsh.model.occ.synchronize()

            # -- Define Physical Groups --
            gmsh.model.addPhysicalGroup(1, [conductor_1], name="Conductor_1")
            gmsh.model.addPhysicalGroup(1, [conductor_2], name="Conductor_2")
            gmsh.model.addPhysicalGroup(1, [boundary], name="OpenBoundary_0")
            gmsh.model.addPhysicalGroup(2, [s_vacuum], name="Vacuum")

            # -- Define Mesh Size Fields --
            field_dist_1 = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(field_dist_1, "CurvesList", [conductor_1])
            field_dist_2 = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(field_dist_2, "CurvesList", [conductor_2])
            field_min = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(field_min, "FieldsList", [field_dist_1, field_dist_2])
            field_thresh = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(field_thresh, "IField", field_min)
            gmsh.model.mesh.field.setNumber(field_thresh, "SizeMin", msize_conductors)
            gmsh.model.mesh.field.setNumber(field_thresh, "SizeMax", msize_boundary)
            gmsh.model.mesh.field.setNumber(field_thresh, "DistMin", 0.01)
            gmsh.model.mesh.field.setNumber(field_thresh, "DistMax", 0.2)            
            gmsh.model.mesh.field.setAsBackgroundMesh(field_thresh)

            # -- Generate Mesh and Write to File --
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.setOrder(3)

            gmsh.write(str(self.output_filepath))
            
            if not self.silent:
                print("Mesh file generated successfully.")
                # The GUI will only run if silent is explicitly set to False
                # and 'close' is not in the system arguments.
                if 'close' not in sys.argv:
                    gmsh.fltk.run()

        finally:
            # Ensure gmsh is always finalized to clean up resources.
            gmsh.finalize()


class PyTulip:
    """
    A class to run pre-existing pulmtln.exe simulation cases,
    assuming a standardized directory structure.
    This class is intended to be run from the PyLCP project root (C:\\git\\PyLCP).
    """
    def __init__(self, executable_path: str, case_name: str, silent: bool = True, overwrite_mesh: bool = True):
        """
        Initializes the PyTulip runner for a specific case.

        Args:
            executable_path (str): The path to the pulmtln.exe executable.
            case_name (str): The name of the case folder inside C:\\git\\tulip\\testData\\.
            silent (bool): If True, suppresses most print statements.
            overwrite_mesh (bool): If True, forces regeneration of the mesh file.
        """
        # --- Path Definitions ---
        base_path = Path(r"C:\git\tulip\testData")
        
        self.executable_path = Path(executable_path).resolve()
        self.case_name = case_name
        self.working_directory = (base_path / self.case_name).resolve()
        self.input_filepath = (self.working_directory / f"{self.case_name}.pulmtln.in.json").resolve()
        self.mesh_filepath = (self.working_directory / f"{self.case_name}.msh").resolve()
        self.output_filepath = (self.working_directory / "pulmtln.out.json").resolve()
        self.silent = silent
        self.overwrite_mesh = overwrite_mesh

        # --- Initial Sanity Checks and Directory/File Creation ---
        if not self.working_directory.exists():
            if not self.silent:
                print(f"Case directory not found. Creating: {self.working_directory}")
            try:
                self.working_directory.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                raise OSError(f"Failed to create directory at {self.working_directory}: {e}")

        if not self.executable_path.is_file():
            raise FileNotFoundError(f"Executable not found at: {self.executable_path}")
        
        # Check if input JSON exists. If not, generate a default one.
        if not self.input_filepath.is_file():
            if not self.silent:
                print(f"Input file not found at {self.input_filepath}.")
            self._create_input_json()

        # REVISED BEHAVIOR: Create mesh if it doesn't exist OR if overwrite is True.
        if self.overwrite_mesh or not self.mesh_filepath.is_file():
            self._create_mesh_file()

        # Initialize result attributes to None
        self.C_matrix = None
        self.L_matrix = None

    def _create_mesh_file(self):
        """
        Creates the required .msh file by calling the appropriate generator.
        """
        mesh_generator = GmshMeshGenerator(
            case_name=self.case_name,
            output_dir=self.working_directory,
            silent=self.silent,
            overwrite=self.overwrite_mesh
        )
        # Add logic here for different cases in the future.
        # For now, we only have the 'two_wires_open' case.
        if self.case_name == "two_wires_open":
            mesh_generator.create_two_wires_mesh()
        else:
            # If we don't know how to generate the mesh, raise an error.
            raise FileNotFoundError(
                f"Mesh file '{self.mesh_filepath.name}' not found, "
                f"and no generation method exists for case '{self.case_name}'."
            )

    def _create_input_json(self):
        """
        Generates a default .pulmtln.in.json file for the given case name.
        """
        if not self.silent:
            print(f"Generating default input file: {self.input_filepath.name}")

        json_data = {
          "analysis": {
            "order": 3,
            "exportParaViewSolution": False,
            "exportFolder": f"Results/{self.case_name}/"
          },
          "model": {
            "materials": {
              "Conductor_0":    {"type": "PEC", "tag": 1 },
              "Conductor_1":    {"type": "PEC", "tag": 2 },
              "OpenBoundary_0": {"type": "OpenBoundary", "tag": 3},
              "Vacuum_0":       {"type": "Vacuum", "tag": 4}
            },  
            "gmshFile": f"{self.case_name}.msh"
          }
        }
        try:
            with open(self.input_filepath, 'w', encoding='utf-8') as f:
                json.dump(json_data, f, indent=2)
            if not self.silent:
                print("Default input file created successfully.")
        except IOError as e:
            print(f"Error: Could not write input file to {self.input_filepath}: {e}")
            raise

    def _run_executable(self):
        """
        Executes pulmtln.exe and streams its output to the console in real-time.
        """
        if not self.silent:
            print(f"\nExecuting: {self.executable_path.name}...")
            print(f"Working Directory: {self.working_directory}")
            print(f"Input File: {self.input_filepath.name}\n")
        
        command = [str(self.executable_path), "-i", self.input_filepath.name]
        
        try:
            with subprocess.Popen(
                command, cwd=self.working_directory,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding='utf-8'
            ) as process:
                if not self.silent:
                    for line in process.stdout:
                        print(line, end='', flush=True)
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


# if __name__ == '__main__':
#     # --- Example Usage ---
#     EXE_PATH = r"C:\git\tulip\pulmtln-build\rls\bin\Release\pulmtln.exe"
#     CASE_NAME = "two_wires_open" 

#     print(f"--- Preparing to run simulation for case: '{CASE_NAME}' ---")
    
#     try:
#         runner = PyTulip(
#             executable_path=EXE_PATH,
#             case_name=CASE_NAME,
#             silent=False,                # Set to True to suppress output
#             overwrite_mesh=True          # Set to True to force mesh regeneration
#         )
#         runner.run()

#         print("\n--- Accessing Results ---")
#         if runner.C_matrix is not None:
#             print("\nCapacitance Matrix (C):")
#             print(runner.C_matrix)
        
#         if runner.L_matrix is not None:
#             print("\nInductance Matrix (L):")
#             print(runner.L_matrix)

#     except (FileNotFoundError, NotADirectoryError, OSError) as e:
#         print(f"\nConfiguration Error: {e}")
#     except Exception as e:
#         print(f"\nAn unexpected error occurred: {e}")