import os
import sys
import time
from pathlib import Path
import matplotlib.pyplot as plt

# RAIZ DO PROJETO E DIRETÓRIOS
os.system('cls' if os.name == 'nt' else 'clear')
try:
    script_dir = Path(__file__).resolve().parent
    print(f"Script directory: {script_dir}")
    project_root = script_dir.parents[1]
    print(f"Project root: {project_root}")
    sys.path.append(str(project_root))
    print("Caminhos do projeto configurados com sucesso.")
except IndexError:
    raise FileNotFoundError(
        "Não foi possível encontrar a raiz do projeto. "
        "Certifique-se de que o script está em 'examples/coated_wires'."
    )

# IMPORTAÇÕES DOS MÓDULOS E MODELO DE DADOS
try:
    from mtl_data.wire_models import BIFILAR_BARE_WIRE_S21_K10 as MTL
    from mtl_data.graphics import MTLRepresentation
    from analyzer.bifilar_bare_wires import BifilarBareWirePULParameters as PUL
    print("Módulos e modelo de dados importados com sucesso.") 
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)

if __name__ == "__main__":
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    st = time.time()
    pul = PUL(project_root, MTL, SUM_MAX=18)    
    pul.run_single_fortran()
    pul.run_analytical()
    pul.run_fortran()
    pul.run_py_mom()
    pul.run_mom_so()
    pul.run_srw_rates()
    pul.run_convergence()

    print(f"\nRotinas de cálculo finalizadas em {(time.time()-st):.2f} segundos.")
    pul.show_header()    
    pul.plot_resistance_results()
    pul.plot_inductance_results()
    pul.plot_capacitance_results()
    pul.plot_srw_rates()
    pul.plot_generalized_capacitance_convergence()    
    pul.plot_free_space_capacitance_convergence()
    MTLRepresentation(MTL, units='millimeter').bared_and_coated_wires()
    plt.show()