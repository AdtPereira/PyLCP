import subprocess
import os
import json

def print_matrix(matrix_name, matrix_data):
    """
    Helper function to print a matrix in a readable format.
    """
    print(f"\n--- {matrix_name} Matrix ---")
    if not matrix_data or not isinstance(matrix_data, list):
        print("Invalid or empty matrix.")
        return
    
    for row in matrix_data:
        # Format each number in scientific notation with 6 decimal places
        formatted_row = [f"{value:.6e}" for value in row]
        print("  ".join(formatted_row))
    print("-" * (len(matrix_name) + 9))


def read_and_print_results(output_filepath):
    """
    Reads the JSON output file, extracts the L and C matrices, and prints them.
    """
    print("\nReading the results file...")
    try:
        with open(output_filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Extract matrices from the Python dictionary created from the JSON
        c_matrix = data.get("C")
        l_matrix = data.get("L")
        
        if c_matrix is not None:
            print_matrix("Capacitance (C)", c_matrix)
        else:
            print("Capacitance (C) matrix not found in the output file.")

        if l_matrix is not None:
            print_matrix("Inductance (L)", l_matrix)
        else:
            print("Inductance (L) matrix not found in the output file.")

    except FileNotFoundError:
        print(f"ERROR: The output file was not found at '{output_filepath}'")
    except json.JSONDecodeError:
        print(f"ERROR: The file '{output_filepath}' is not a valid JSON.")
    except Exception as e:
        print(f"An unexpected error occurred while reading the results: {e}")


def run_tulip_simulation():
    """
    Executes the pulmtln simulation and, upon success, reads and displays the results.
    """
    # --- Path Configuration ---
    executable_path = r"C:\git\tulip\pulmtln-build\rls\bin\Release\pulmtln.exe"
    working_directory = r"C:\git\tulip\examples\two_wires_open"
    input_filename = "two_wires_open.pulmtln.in.json"
    output_filename = "pulmtln.out.json" # Name of the output file
    
    full_output_path = os.path.join(working_directory, output_filename)

    # --- Command Construction and Execution ---
    command = [executable_path, "-i", input_filename]
    
    print("Starting pulmtln execution...")
    # ... (existence checks can be kept here) ...
    
    try:
        # Execute the simulation
        simulation_result = subprocess.run(
            command,
            cwd=working_directory,
            check=True,
            capture_output=True,
            text=True,
            encoding='utf-8'
        )
        
        print("--- Program Output ---")
        print(simulation_result.stdout)
        print("----------------------")
        print("\nExecution completed successfully!")
        
        # --- Read Results ---
        # If the simulation was successful, call the function to read the results.
        read_and_print_results(full_output_path)

    except subprocess.CalledProcessError as e:
        print("\nERROR: The program returned an error code.")
        print(f"Return code: {e.returncode}")
        print("\n--- Standard Error (stderr) ---")
        print(e.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    run_tulip_simulation()