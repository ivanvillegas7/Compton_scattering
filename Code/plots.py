#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 12:43:31 2023

@author: ivanvi
"""

import numpy as np

import matplotlib.pyplot as plt

def plot_spectra(file: str):
    
    data = np.loadtxt(f'../Data/{file}.asc', skiprows=2)
    
    E: np.array(float) = data[:, 0]
    
    I: np.array(float) = data[:, 1]
    
    
    plt.plot(E, I, label=file)
#%% 
from typing import List
   
files: List[str] = ['0 deg', '30 deg', '60 deg', '90 deg', '120 deg',\
                    '90 deg block', '120 deg block']

for file in files:
    
    plt.figure()
    plot_spectra(file)
    plt.xlabel(r'$E$ [keV]')
    plt.ylabel(r"Intensity")
    plt.title(r"Compton scattering")
    plt.legend()
    plt.grid()
    plt.savefig(f'../Plots/spectra {file}.pdf')
#%%
plt.figure()
plot_spectra('sodium cal')
plt.xlabel(r'$E$ [keV]')
plt.ylabel(r"Intensity")
plt.title(r"Compton scattering")
plt.legend()
plt.grid()
plt.savefig('../Plots/spectra sodium cal.pdf')

plt.figure()
plot_spectra('cesium cal')
plt.xlabel(r'$E$ [keV]')
plt.ylabel(r"Intensity")
plt.title(r"Compton scattering")
plt.legend()
plt.grid()
plt.savefig('../Plots/spectra cesium cal.pdf')