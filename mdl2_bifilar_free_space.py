"""
This script analyzes the behavior of a two-wire transmission line using the method of moments (MoM).
The code is structured into classes and functions, facilitating a modular approach to the problem. 

REFERENCES:
[1] 

"""

import os
import time
import numpy as np
from data import systems, multiconductor as mtl
from free_space.pul_parameters import Bifilar
from data.multiconductor import GraphicRepresentation as graph
from mom_so import green, patel


# Multiconductor Transmission Line choices
# 1 - Bifilar Transmission Line (TWT) with 100 mm spacing and Np=Nq=0
MTL = systems.MULTICONDUCTOR_TRANSMISSION_LINE[1]


def main():
    """ Main function to perform the calculations and display the results."""

    # 1. Clears the console screen and starts the timer
    os.system('cls' if os.name == 'nt' else 'clear')
    print("Calculations were started ... ...")
    start_time = time.time()

    # 2. Define the frequency range
    frequency = [np.logspace(0, 7, num=200), np.logspace(0, 7, num=30)]

    # 3. Analytical routine
    Zs = []  # pylint: disable=invalid-name
    Rhf = []  # pylint: disable=invalid-name
    Lext = []  # pylint: disable=invalid-name
    for f in frequency[0]:
        Zs.append(Bifilar(MTL).series_impedance(f)[0])
        Rhf.append(Bifilar(MTL).series_impedance(f)[1])
        Lext.append(Bifilar(MTL).series_impedance(f)[2])

    # 4. MoM-SO routine
    Zs_mom = []  # pylint: disable=invalid-name
    green_list = ['Analytically', 'Numerically']
    green_matrix = green.QuasiStatic(MTL).g_tanaka(mode=green_list[0])
    post_processing = patel.LosslessPostProcessing(MTL)
    for f in frequency[1]:
        mom_so = patel.HomogeneousLosslessMedium(MTL, f)
        z_partial = mom_so.z_partial(green_matrix)
        Zs_mom.append(post_processing.z_total(z_partial))

    # 5. End the timer
    elapsed_time = time.time() - start_time
    print(f"End of the routine! Time spent on simulation: {
          elapsed_time:.2f} seconds.\n")

    # 6. Display the geometry of the transmission line
    graph(MTL).wires_and_cables(line_type='bifilar')

    # 7. Plot the series resistance as a function of frequency
    data = [0, MTL['data'][0]['fourier_order'],
            mtl.MulticonductorTransmissionLine(MTL).D[0]]
    Bifilar(MTL).plot_series_resistance(frequency, [Zs, Zs_mom], Rhf, data)
    Bifilar(MTL).plot_series_inductance(frequency, [Zs, Zs_mom], Lext, data)


if __name__ == "__main__":
    main()
