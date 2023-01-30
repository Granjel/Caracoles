"""
author: J. W. Spaak
plot NFD_distribution
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

try:
    NFD_data
except NameError:
    NFD_data = pd.read_csv("NFD_uncertanty.csv")


ND_cols = [True if col[:2] == "ND" else False for col in NFD_data.columns]
ND_cols = NFD_data.columns[ND_cols]
FD_cols = [True if col[:2] == "FD" else False for col in NFD_data.columns]
FD_cols = NFD_data.columns[FD_cols]
equi_cols = [True if col[:4] == "equi" else False for col in NFD_data.columns]
equi_cols = NFD_data.columns[equi_cols]    

ND = NFD_data[ND_cols].values.flatten()
FD = NFD_data[FD_cols].values.flatten()
equi = NFD_data[equi_cols].values.flatten()
equi[np.isnan(equi)] = -1

color = 1*((equi>0) == (ND > FD))
color = np.array(["red", "green"])[color]

fig = plt.figure()
perc = 5
perc = [perc, 100-perc]
plt.scatter(ND, FD, s = 2, alpha = 0.1, c = color)
plt.xlim(np.nanpercentile(ND, perc))
plt.ylim([np.nanpercentile(FD, perc)[0], 1])
x = np.linspace(*np.nanpercentile(ND, perc), 101)
plt.plot(x,x, 'k')

fig.savefig("NFD_distribution")

fig, ax = plt.subplots(5,9, figsize = (18,15), sharex = True, sharey = True)
ax[0,0].set_xlim(np.nanpercentile(ND, perc))
ax[0,0].set_ylim([np.nanpercentile(FD, perc)[0], 1])
ax = ax.flatten()

for i, case in enumerate(sorted(set(NFD_data.case))):
    data = NFD_data[NFD_data.case == case]
    ND = data[ND_cols].values.flatten()
    FD = data[FD_cols].values.flatten()
    equi = data[equi_cols].values.flatten()
    equi[np.isnan(equi)] = -1
    
    color = 1*((equi>0) == (ND > FD))
    color = np.array(["red", "green"])[color]
    ax[i].scatter(ND, FD, s = 2, alpha = 0.2, c = color)
    ax[i].set_title(case)
    ax[i].plot(x,x, 'k')
    
fig.savefig("NFD_per_plot.pdf")