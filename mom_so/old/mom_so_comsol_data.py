""" This script reads the data from a COMSOL simulation and plots it. """

import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the file
file_path = 'C:\\Users\\adilt\\OneDrive\\01 ACADEMIA\\06 MODELOS\\7.MoM-SO\\data\\comsol_resistance_coax.txt'
data = pd.read_csv(file_path, sep='\s+', comment='%')

# Rename the columns
data.columns = ['freq (Hz)', 'Analytic (DC)', 'Analytic (HF)', 'COMSOL (mf/ec)']

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(data['freq (Hz)'], data['Analytic (DC)'], label='Analytic (DC)')
plt.plot(data['freq (Hz)'], data['Analytic (HF)'], label='Analytic (HF)')
plt.plot(data['freq (Hz)'], data['COMSOL (mf/ec)'], label='COMSOL (mf/ec)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Resistance')
plt.title('Resistance vs Frequency')
plt.legend()
plt.grid(True)
plt.show()
