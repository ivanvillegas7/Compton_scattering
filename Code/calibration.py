#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 10:34:57 2023

@author: Iván
"""

import numpy as np

import matplotlib.pyplot as plt

from typing import List

def calibrate():
    
    data = np.loadtxt('../Data/sodium cal.asc', skiprows=2)
    
    E: np.array(float) = data[:, 0]
    
    I: np.array(float) = data[:, 1]
    
    peak_1I: List[float] = []
    
    peak_1E: List[float] = []
    
    peak_2I: List[float] = []
    
    peak_2E: List[float] = []
    
    for i in range(len(E)):
    
        if E[i]>200 and E[i]<400:
            
            peak_1I.append(I[i])
            
            peak_1E.append(E[i])
            
        elif E[i]>600 and E[i]<800:
            
            peak_2I.append(I[i])
            
            peak_2E.append(E[i])
            
    peak1E: np.array(float) = np.array(peak_1E)
    
    peak2E: np.array(float) = np.array(peak_2E)
    
    peak1I: np.array(float) = np.array(peak_1I)
    
    peak2I: np.array(float) = np.array(peak_2I)
    
    energies: np.array(float) = np.array([511, 1274.5])
    
    x: np.array(float) = np.array([peak1E[np.argmax(peak1I)],\
                                   peak2E[np.argmax(peak2I)]])
    
    params: np.array(float)
    
    cov: List[np.array(float), np.array(float)]
    
    params, cov = np.polyfit(x, energies, 1, w=np.ones(2)/0.5, cov='unscaled')
    
    m: float = params[0]
    
    n:float = params[1]
    
    plt.figure()
    plt.plot(x, m*x+n, label='Fit')
    plt.plot(x, energies, label='Data points', ls='none', marker='.')
    plt.ylabel(r'$E$ [keV]')
    plt.xlabel(r"Channel")
    plt.title(r"Calibration points")
    plt.legend()
    plt.grid()
    
    return (m, n, np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1]))
    
    