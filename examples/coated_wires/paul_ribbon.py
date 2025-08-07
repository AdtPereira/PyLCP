import os
import sys
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
    from analyzer.convergence import ConvergenceAnalyzer
    from mtl_data.models import PAUL_RIBBON_CABLE
    from mom.coated_wire_systems import TwoCoatedWireSystem
    
    print("Módulos e modelo de dados importados com sucesso.")
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

if __name__ == '__main__':
    analyzer = ConvergenceAnalyzer(project_root, PAUL_RIBBON_CABLE, NF_MAX=10)
    analyzer.run_single_simulation()
    analyzer.run_simulations()
    analyzer.plot_results()
    plt.show()
