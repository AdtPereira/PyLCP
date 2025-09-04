# root_folder/model/case_utils.py

import json
from pathlib import Path

def load_json_parameters(script_file_path, show_content=False):
    """
    Dynamically loads parameters from a '.json' file.
    It assumes the '.in.json' file has the same base name as the
    calling script and is located in the same directory.

    Args:
        script_file_path (str): The __file__ attribute from the calling script.
        show_content (bool): If True, prints the content of the loaded
                             dictionary to the console. Defaults to False.

    Returns:
        dict: A dictionary with the parameters loaded from the JSON file.
    """
    script_path = Path(script_file_path)
    script_dir = script_path.resolve().parent
    case_name = script_path.stem
    model_file_name = f"{case_name}.json"
    model_path = script_dir / model_file_name
    
    if not model_path.exists():
        raise FileNotFoundError(f"The parameter file could not be found at: {model_path}")

    with open(model_path, 'r') as f:
        parameters = json.load(f)
    
    print(f"Successfully loaded parameters from: {model_path}")
    
    # --- NOVO BLOCO DE CÓDIGO ---
    # Se o parâmetro show_content for True, imprime o dicionário
    if show_content:
        print(f"\n--- Content of {model_file_name} ---")
        # Usa json.dumps para uma impressão "pretty print"
        print(json.dumps(parameters, indent=2))
        print("--------------------------------------\n")
        
    return parameters