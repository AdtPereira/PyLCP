# debug_runner.py
import runpy
import sys
from pathlib import Path

if __name__ == "__main__":
    # Pega o caminho do arquivo passado como argumento pelo launch.json
    file_to_debug = Path(sys.argv[1])

    # Assume que a raiz do projeto é o diretório atual
    project_root = Path.cwd()

    # Converte o caminho do arquivo em um nome de módulo
    # Ex: "TestData/bare_bifilar_s100/bare_bifilar_s100.py" -> "TestData.bare_bifilar_s100.bare_bifilar_s100"
    module_path = file_to_debug.relative_to(project_root).with_suffix('').as_posix().replace('/', '.')

    print(f"--- Debug Runner: Executando o módulo '{module_path}' ---")

    # Usa runpy para executar o módulo no contexto correto,
    # o que faz com que as importações relativas funcionem.
    runpy.run_module(module_path, run_name='__main__')