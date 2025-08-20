# File: run_all.py
# Placed at the root of the PyLCP project.
# Purpose: To sequentially run all example scripts to verify they work correctly.

import subprocess
import sys
import os

# --- Configuration ---
SCRIPTS_TO_RUN = [
    os.path.join("examples", "buried_cables", "one_scc.py"),
]

def main():
    """
    Executes a list of Python scripts sequentially.
    Stops execution if any script fails.
    """
    print("--- Starting all example scripts ---")
    
    for script_path in SCRIPTS_TO_RUN:
        print(f"\n>>> Running: {script_path}")
        
        if not os.path.exists(script_path):
            print(f"[ERROR] Script not found: {script_path}", file=sys.stderr)
            print("--- Sequence aborted ---")
            sys.exit(1)

        try:
            # **CHANGE 1: Adjusted encoding for Windows default in Portuguese (cp1252).**
            # The 'errors' parameter will replace any problematic character, preventing crashes.
            result = subprocess.run(
                [sys.executable, script_path], 
                check=True, 
                capture_output=True,
                text=True,
                encoding='cp1252', # Changed from 'utf-8'
                errors='replace'   # Added for more robustness
            )
            print(f"--- Successfully finished: {script_path} ---")
            if result.stdout:
                print("Output:\n", result.stdout)

        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Script failed with a non-zero exit code: {script_path}", file=sys.stderr)
            print(f"Return Code: {e.returncode}", file=sys.stderr)
            print(f"\n--- Output (stdout) ---\n{e.stdout}", file=sys.stderr)
            print(f"\n--- Error Output (stderr) ---\n{e.stderr}", file=sys.stderr)
            print("--- Sequence aborted due to error ---")
            sys.exit(1)
        
        # **CHANGE 2: Added a block to catch other exceptions, like encoding errors.**
        except Exception as e:
            print(f"[UNEXPECTED ERROR] An error occurred while running {script_path}: {e}", file=sys.stderr)
            print("This could be an encoding issue or another problem.", file=sys.stderr)
            print("--- Sequence aborted due to error ---")
            sys.exit(1)

    print("\n\n--- All example scripts executed successfully! ---")

if __name__ == "__main__":
    main()