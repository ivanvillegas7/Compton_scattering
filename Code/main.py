# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 15:04:05 2021

@author: Iván
"""

import numpy as np

import matplotlib.pyplot as plt

from typing import List

import scipy as sc

import plots

def peaks(option: bool):
    
    if option.upper()=='N':
    
        files: List[str] = ['0 deg', '30 deg', '60 deg', '90 deg', '120 deg']
    
    else:
        
        files: List[str] = ['0 deg', '30 deg', '60 deg', '90 deg block',\
                            '120 deg block']
    
    peak: List[float] = []
    
    for file in files:
        
        data = np.loadtxt(f'../Data/{file}.asc', skiprows=2)
        
        E: np.array(float) = data[:, 0]
        
        I: np.array(float) = data[:, 1]
    
        i: int = np.argmax(I)
        
        print(f'\nPosition for the highest peak for {file}: {i}.')
        
        print(f'\nScattered energy peak for {file}: {E[i]} keV.')
        
        peak.append(E[i])
        
    return peak

def curve(x:np.array(float), m: float, n: float):
    
    return m*x+n

def main():

    m_e: float = 511#keV/c²
    
    c: float = 1#Natural units
    
    E: float = 661.7#keV
    
    theta_th: np.array(float) = np.linspace(0, np.pi)
    
    E_prime_th: np.array(float) = 1/((1/E) + (1-np.cos(theta_th))/(m_e*c**2))
    
    theta: np.array(float) = np.array([0, 30, 60, 90, 120])*np.pi/180
    
    print('')
    
    option: bool = input('Do you want to use the corrected data? [Y/N]: ')
    
    while option.upper()!='Y' and option.upper()!='N':
        
        print(f'\n{option} is not a valid character. Type Y or N.')
        
        option = input('Do you want to use the corrected data? [Y/N]: ')
    
    E_prime: np.array(float) = np.array(peaks(option))
        
    err_E_prime: np.array(float) = 1/E_prime**2
    
    err_theta: np.array(float) = 0.1*(np.sin(theta))
    
    params: np.array(float)
    
    cov: List[np.array(float), np.array(float)]
              
    params, cov = sc.optimize.curve_fit(curve, 1-np.cos(theta), 1/E_prime,\
                                        sigma=err_E_prime)
    
    m: float = params[0]
    
    n: float = params[1]
    
    x: np.array(float) = 1-np.cos(theta)
    
    x_th: np.array(float) = 1-np.cos(theta_th)
    
    plt.figure()
    plt.errorbar(x, 1/E_prime, xerr=err_theta, yerr=err_E_prime,\
                 label='Experimental data', marker= '.', linestyle='none')
    plt.plot(x_th, m*(1-np.cos(theta_th))+n,\
             label=f'Fit: {m:.3f}x+({n:.3f})')
    plt.plot(x_th, 1/E_prime_th, label='Theoretical curve')
    plt.xlabel(r'$1-\cos{(\theta)}$')
    plt.ylabel(r"$E'^{-1}$ [keV$^{-1}$]")
    plt.title(r"$E'^{-1}$ against $1-\cos{(\theta)}$ in the Compton scattering")
    plt.xlim(right=2)
    plt.legend()
    plt.grid()
    plt.savefig(f'../Plots/plot_{option}_fit.pdf')
    
    print(f'\nThe determined mass of the electron is m_e = ({1/m}±{cov[0, 0]/m**2}) keV/c².')
    
    print(f'\nThe determined incident energy is E = ({1/n}±{cov[1, 1]/n**2}) keV.')
    
main()
