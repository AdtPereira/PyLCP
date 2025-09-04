# C:\git\PyLCP\utils\case_utils.py
import importlib.util
from pathlib import Path

def load_model(script_file_path):
    """
    Dynamically loads the 'MODEL' variable from a '.in.py' file.
    It assumes the '.in.py' file has the same base name as the calling script
    and is located in the same directory (is a sibling file).

    Args:
        script_file_path (str): The __file__ attribute from the calling script.

    Returns:
        The 'MODEL' data structure from the corresponding '.in.py' file.
    """
    script_path = Path(script_file_path)
    script_dir = script_path.resolve().parent
    case_name = script_path.stem  # Gets the filename without the .py extension
    model_file_name = f"{case_name}.in.py"
    model_path = script_dir / model_file_name
    
    if not model_path.exists():
        raise FileNotFoundError(f"The model file could not be found at: {model_path}")

    spec = importlib.util.spec_from_file_location(case_name, model_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create module spec from file: {model_path}")
        
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    if not hasattr(module, 'MODEL'):
        raise AttributeError(f"The model file '{model_path}' must contain a 'MODEL' variable.")
    
    print(f"Successfully loaded model from: {model_path}")
    return module.MODEL