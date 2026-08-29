# debug_runner.py
import runpy
import sys
from pathlib import Path

if __name__ == "__main__":
    # Get the file path passed as an argument by launch.json
    file_to_debug = Path(sys.argv[1])

    # Assume the project root is the current directory
    project_root = Path.cwd()

    # Convert the file path into a module name
    # e.g. "TestData/bare_bifilar_s100/bare_bifilar_s100.py" -> "TestData.bare_bifilar_s100.bare_bifilar_s100"
    module_path = file_to_debug.relative_to(project_root).with_suffix('').as_posix().replace('/', '.')

    print(f"--- Debug Runner: Running module '{module_path}' ---")

    # Use runpy to run the module in the correct context,
    # which makes relative imports work.
    runpy.run_module(module_path, run_name='__main__')