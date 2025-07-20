import os
# Clear the terminal window
os.system('cls' if os.name == 'nt' else 'clear')

from scc_data import PRY_M01_1SCC_1C as model, SystemType
from scc_result import Result
from scc_impedance_ground_return import GroundReturnImpedance
from scc_impedance import InternalImpedanceAprox, SeriesImpedanceMatrix

# Define single core cable and soil parameters
scc = model()

# Define system configuration
syst = SystemType(SystId=1, SCC=scc, rhog=100, erg=1)

# Angular Frequency
jw = 9.010913347279287e+03 + 6.283185307179586e+03j 

# Ground Return Impedance
Gri = GroundReturnImpedance(jw, scc, syst)
Zg = Gri.DeConti()

# SCC Internal Impedance
zI = InternalImpedanceAprox(jw, scc)

# Loop Impedance Matrix (NODA,2008)
A, ZL = SeriesImpedanceMatrix().MatrixTransformation(zI, Zg, syst)

# Series Impedance Matrix
Zs = SeriesImpedanceMatrix().SeriesImpedance(zI, Zg, syst)

# Log Result
Result.show_results(scc, syst, Gri, Zg, zI, A, ZL, Zs, 'PRY_M01_1SCC_1C')