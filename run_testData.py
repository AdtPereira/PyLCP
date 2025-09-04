# File: run_all.py
# Placed at the root of the PyLCP project.
# Purpose: To sequentially run all example scripts to verify they work correctly.

import subprocess
import sys
import os

# --- CHANGE 1: Use module notation instead of file paths ---
SCRIPTS_TO_RUN = [
    "testData.pipe_trefoil_patel.pipe_trefoil_patel",
    "testData.scc_flat_xue.scc_flat_xue",
]

def main():
    """
    Executes a list of Python modules sequentially.
    Stops execution if any module fails.
    """
    print("--- Starting all example scripts ---")
    
    for script_module in SCRIPTS_TO_RUN:
        print(f"\n>>> Running module: {script_module}")
        
        # We can't easily check if a module exists without trying to import it,
        # so we'll let subprocess handle the error if it's not found.

        try:
            # --- CHANGE 2: Add the '-m' flag to run as a module ---
            result = subprocess.run(
                [sys.executable, "-m", script_module], 
                check=True, 
                capture_output=True,
                text=True,
                encoding='utf-8', # It's safer to standardize on utf-8
                errors='replace'
            )
            print(f"--- Successfully finished: {script_module} ---")
            if result.stdout:
                print("Output:\n", result.stdout)

        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Module failed with a non-zero exit code: {script_module}", file=sys.stderr)
            print(f"Return Code: {e.returncode}", file=sys.stderr)
            print(f"\n--- Output (stdout) ---\n{e.stdout}", file=sys.stderr)
            print(f"\n--- Error Output (stderr) ---\n{e.stderr}", file=sys.stderr)
            print("--- Sequence aborted due to error ---")
            sys.exit(1)
        
        except Exception as e:
            print(f"[UNEXPECTED ERROR] An error occurred while running {script_module}: {e}", file=sys.stderr)
            print("--- Sequence aborted due to error ---")
            sys.exit(1)

    print("\n\n--- All example scripts executed successfully! ---")

if __name__ == "__main__":
    main()