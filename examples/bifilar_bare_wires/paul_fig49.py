"""
Reprodução da Figura 4.9 do livro "Introduction to Electromagnetic Compatibility".

Este script calcula e plota a capacitância por unidade de comprimento de uma
linha de transmissão bifilar usando as fórmulas exata e aproximada. O objetivo
é recriar a Figura 4.9, que compara essas duas formulações em função da razão
entre a separação dos condutores e seu raio (s/r_w).

Estrutura:
- Define-se uma faixa de valores para a razão s/r_w.
- Para cada valor, uma geometria de linha bifilar é criada.
- A classe Bifilar é usada para calcular as capacitâncias.
- Os resultados são plotados usando Matplotlib.
"""

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
    from mtl_data.models import BIFILAR_BARE_WIRE_S21
    from analyzer.bifilar_wires import BifilarPULParameters
    print("Módulos e modelo de dados importados com sucesso.")
except ImportError as e:
    print(f"Erro ao importar módulos: {e}")
    sys.exit(1)


if __name__ == "__main__":
    """ Função principal para orquestrar a análise, cálculo e visualização dos resultados. """
    start_time = time.time()
    analyzer = BifilarPULParameters(project_root, BIFILAR_BARE_WIRE_S21)
    analyzer.srw_rates_analytical()
    analyzer.srw_rates_py_mom()
    
    print(f"\nRotinas de cálculo finalizadas! Tempo de simulação: {(time.time() - start_time):.2f} segundos.")
    analyzer.plot_paul_fig49()
    plt.show()
