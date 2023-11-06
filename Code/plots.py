#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 12:43:31 2023

@author: Iván
"""

import numpy as np

import matplotlib.pyplot as plt

from typing import List

import calibration

def plot_spectra(file: str, m: float, n: float, option: str):
    
    if option.upper()=='TEST':
        
        data = np.loadtxt(f'../Test/{file}.asc')
        
    else:
    
        data = np.loadtxt(f'../Data/{file}.asc', skiprows=2)
        
    I: np.array(float) = data[:, 1]
    
    E: np.array(float) = data[:, 0]*m+n*np.ones(len(I))
    
    plt.plot(E, I, label=file, marker='.', ls='none')
    
def plot(option: str):
    
    if option.upper()=='TEST':
        
        files: List[str] = ['guessed_zero', '35_degrees', '65_degrees',
                            '95_degrees', '120_degrees']
            
        m: float = 1.862
        
        n: float = -17.446
        
        for file in files:
            
            plt.figure()
            plot_spectra(file, m, n, option)
            plt.xlabel(r'$E$ [keV]')
            plt.ylabel(r"Counts")
            plt.title(r"Compton scattering")
            plt.xlim(0, 800)
            plt.legend()
            plt.grid()
            plt.savefig(f'../Test/spectra {file}.pdf')

    else:    

        files: List[str] = ['0 deg', '30 deg', '90 deg', '120 deg',\
                            '90 deg block', '120 deg block', 'sodium cal',\
                            'cesium cal']
        
        params: np.array(float) = np.array(calibration.calibrate())
        
        m: float = params[0]
        
        n: float = params[1]
            
        for file in files:
            
            plt.figure()
            plot_spectra(file, m, n, option)
            plt.xlabel(r'$E$ [keV]')
            plt.ylabel(r"Counts")
            plt.title(r"Compton scattering")
            plt.xlim(0, 800)
            plt.legend()
            plt.grid()
            plt.savefig(f'../Plots/spectra {file}.pdf')
