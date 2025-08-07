import os
import sys
import time
from pathlib import Path
import matplotlib.pyplot as plt

# RAIZ DO PROJETO E DIRETÓRIOS
os.system('cls' if os.name == 'nt' else 'clear')
try:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parents[1]
    sys.path.append(str(project_root))        
    
    print("Caminhos do projeto configurados com sucesso.")
except IndexError:
    raise FileNotFoundError(
        "Não foi possível encontrar a raiz do projeto. "
        "Certifique-se de que o script está em 'examples/coated_wires'."
    )

# IMPORTAÇÕES DOS MÓDULOS E MODELO DE DADOS
try:
    from analyzer.compare_methods import PULparameters
    from mtl_data.models import BIFILAR_BARE_WIRE_S21
    from mtl_data.graphics import MTLRepresentation
    print("Módulos e modelo de dados importados com sucesso.")
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

if __name__ == "__main__":
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    start_time = time.time()
    analyzer = PULparameters(project_root, BIFILAR_BARE_WIRE_S21)
    analyzer.show_header()
    MTLRepresentation(BIFILAR_BARE_WIRE_S21).bare_and_coated_wires()
    
    analyzer.run_fortran()
    analyzer.run_mom(autoPlots=False)
    analyzer.run_analytical()
    analyzer.run_mom_so()

    print(f"\nRotinas de cálculo finalizadas em {(time.time()-start_time):.2f} segundos.")
    analyzer.plot_resistance_results()
    analyzer.plot_inductance_results()
    analyzer.plot_capacitance_results()
    plt.show()