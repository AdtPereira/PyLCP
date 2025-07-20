import numpy as np
from scipy.constants import mu_0 as MUO
from scc_user_functions import MatrixDef as md
from scc_user_functions import Trigonometric as tg

from scc_data import DC_M02_3SCC_1C as Model

R = Model()

class ClassTest:

    def __init__(self):
        pass  

    def RetonarRaio():
        print(f"Pela classe ClassTest, chamo o método RetornarRaio. O raio ra = {R.a} m")

print(f"O raio ra = {R.a} m\n")