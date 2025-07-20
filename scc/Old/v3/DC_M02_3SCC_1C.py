import os
# Clear the terminal window
os.system('cls' if os.name == 'nt' else 'clear')

from scc_data import scc_models_list
from scc_systm import SystemType
from scc_result import Result
from scc_impedance import GroundReturnImpedance as gri
from scc_impedance import SeriesImpedanceMatrix as sim

# Define system configuration and soil parameters
# Model Name: DC_M02_3SCC_1C
scc = scc_models_list[1]
syst = SystemType(Model=scc, rhog=100, erg=1, Syst_id='#1')

# Angular Frequency
jw = 9.010913347279287e+03 + 6.283185307179586e+03j 

# Ground Return Impedance
Gri = gri(jw, syst)
Zg = Gri.DeConti(syst)

# SCC Internal Impedance
Zc = sim(jw).CoreImpedanceAprox(syst.ra, syst.rho_c)
Zeins = sim(jw).InsulationImpedance(syst.rf, syst.re)
Zi = {"Zc": Zc, "Zeins": Zeins}

# Loop Impedance Matrix (NODA,2008)
zL = sim(jw).LoopImpedance(syst, Zg[0][0])
A, ZL = sim(jw).MatrixTransformation(syst, zL, Zg)

# Series Impedance Matrix
Zs = sim(jw).SeriesImpedance(A, ZL)

# Log Result
Result.show_results(syst, Gri, Zg, Zi, A, ZL, Zs)