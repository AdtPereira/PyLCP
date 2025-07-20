# clean screen
import os
os.system('cls' if os.name == 'nt' else 'clear')

# Python library 
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

# User library
from scc_data import DC_M02_3SCC_1C as model, SystemType
from scc_impedance_ground_return import GroundReturnImpedance

# Define single core cable 
scc = model()

# Define system configuration and soil parameters
syst_A = SystemType(1, scc, 1E2, 10)
syst_B = SystemType(2, scc, 1E3, 10)
syst_C = SystemType(3, scc, 1E4, 10)

# Generate the logarithmic range
decade_start = 4        # Initial decade (1E4)
decade_end = 7          # Final decade (1E7)
num_points = 300        # Number of points
f = np.logspace(decade_start, decade_end, num=num_points)
s = 1j * (2 * np.pi * f)

# Element-31 Ground-Return Impedance Matrices
Z31A = []
Z31B = []
Z31C = []

for k in range (len(s)):
    Zg_A = GroundReturnImpedance(s[k], scc, syst_A).DeConti()
    Zg_B = GroundReturnImpedance(s[k], scc, syst_B).DeConti()
    Zg_C = GroundReturnImpedance(s[k], scc, syst_C).DeConti()

    Z31A.append(Zg_A[2][0])
    Z31B.append(Zg_B[2][0])
    Z31C.append(Zg_C[2][0])

# Import the .mat files
MatData = sio.loadmat(r'C:\Users\adilt\OneDrive\01 ACADEMIA\05 PPGEE\EEE003 MODELOS\1.DECONTI\MatLabData.mat')
Mat_f = MatData['f']['A1'][0][0][0]
Mat_Z31A = MatData['Z31']['A5'][0][0][0]
Mat_Z31B = MatData['Z31']['B5'][0][0][0]
Mat_Z31C = MatData['Z31']['C5'][0][0][0]

# Create the semilogx graphs
plt.figure(figsize=(8, 6))
plt.semilogx(f, np.abs(Z31A), color='red', label='Python Data')
plt.semilogx(f, np.abs(Z31B), color='red')
plt.semilogx(f, np.abs(Z31C), color='red')

plt.semilogx(Mat_f, np.abs(Mat_Z31A), color='black', linestyle='--', label=r'MATLAB Data')
plt.semilogx(Mat_f, np.abs(Mat_Z31B), color='black', linestyle='--')
plt.semilogx(Mat_f, np.abs(Mat_Z31C), color='black', linestyle='--')

# Add labels and title
plt.xlabel(r'$Frequency\;(Hz)$')
plt.ylabel(r'$|Z_g|\;(\Omega/m)$')
plt.legend(loc='upper left')
plt.xlim(1E4, 1E7)
plt.ylim(0, 25)
plt.grid(False)
plt.show()