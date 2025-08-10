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
    from analyzer.convergence import BifilarConvergenceAnalyzer
    from mtl_data.graphics import MTLRepresentation
    from mtl_data.models import PAUL_RIBBON_CABLE
    print("Módulos e modelo de dados importados com sucesso.")
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

if __name__ == '__main__':
    st = time.time()
    analyzer = BifilarConvergenceAnalyzer(project_root, PAUL_RIBBON_CABLE, SUM_MAX=10)
    analyzer.run_single_fortran_simulation()
    analyzer.run_convergence()
    
    print(f"\nRotinas de cálculo finalizadas em {(time.time()-st):.2f} segundos.")
    analyzer.plot_paul_fig514a()
    analyzer.plot_paul_fig514b()
    analyzer.plot_paul_fig514c()
    analyzer.plot_paul_fig514d()
    MTLRepresentation(PAUL_RIBBON_CABLE).bare_and_coated_wires()
    plt.show()
