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

import calibration

def peaks(option: bool):
    
    peak: List[float] = []
    
    if option.upper()=='N':
    
        files: List[str] = ['0 deg', '30 deg', '90 deg', '120 deg']
        
    elif option.upper()=='TEST':
        
        files: List[str] = ['guessed_zero', '35_degrees', '65_degrees',\
                            '95_degrees', '120_degrees']
            
        for file in files:
            
            data = np.loadtxt(f'../Test/{file}.asc')
            
            Ch: np.array(float) = data[:, 0]
            
            I: np.array(float) = data[:, 1]
            
            peak.append(Ch[np.argmax(I)])
            
        return np.array(peak)
    
    else:
        
        files: List[str] = ['0 deg', '30 deg', '90 deg block', '120 deg block']
       
    for file in files:
        
        data = np.loadtxt(f'../Data/{file}.asc', skiprows=2)
        
        Ch: np.array(float) = data[:, 0]
        
        I: np.array(float) = data[:, 1]
        
        peak.append(Ch[np.argmax(I)])
        
    return np.array(peak)

def curve(x:np.array(float), m: float, n: float):
    
    return m*x+n

def main():

    m_e: float = 511#keV/c²
    
    c: float = 1#Natural units
    
    E: float = 661.7#keV
    
    print('')
    
    option: bool = input('Do you want to use the corrected data? [Y/N]: ')
    
    while option.upper()!='Y' and option.upper()!='N' and option.upper()!='TEST':
        
        print(f'\n{option} is not a valid character. Type Y or N.')
        
        option = input('Do you want to use the corrected data? [Y/N]: ')
        
    print('')
    
    m: float
    
    n: float
    
    sigma_m: float
    
    sigma_n: float
    
    if option.upper()=='TEST':
        
        m = 1.862
        
        n = -17.446
        
        sigma_m = 0
        
        sigma_n = 0
        
    else:
    
        m, n, sigma_m, sigma_n = calibration.calibrate()
        
    plots.plot(option)
    
    theta_th: np.array(float) = np.linspace(0, np.pi)
    
    E_prime_th: np.array(float) = 1/((1/E) + (1-np.cos(theta_th))/(m_e*c**2))
    
    theta: np.array(float) = np.array([0, 30, 60, 90, 120])*np.pi/180
        
    E_prime: np.array(float) = peaks(option)*m+n*np.ones(len(theta))
           
    sigma_E: np.array(float)
    
    sigma_E = np.sqrt(m**2*0.5**2+peaks(option)**2*sigma_m**2+sigma_n**2)
    
    err_E_prime: np.array(float) = sigma_E/E_prime**2
    
    err_theta: np.array(float) = 0.1*(np.sin(theta))
    
    err_theta[0] = np.max(err_theta)
    
    params: np.array(float)
    
    cov: List[np.array(float), np.array(float)]
              
    params, cov = sc.optimize.curve_fit(curve, 1-np.cos(theta), 1/E_prime,\
                                        sigma=err_E_prime)
    
    a: float = params[0]
    
    b: float = params[1]
    
    x: np.array(float) = 1-np.cos(theta)
    
    x_th: np.array(float) = 1-np.cos(theta_th)
    
    plt.figure()
    plt.errorbar(x, 1/E_prime, xerr=err_theta, yerr=err_E_prime,\
                 label='Experimental data', marker= '.', linestyle='none')
    plt.plot(x_th, a*(1-np.cos(theta_th))+b,\
             label=f'Fit: {a:.3f}x+({b:.3f})')
    plt.plot(x_th, 1/E_prime_th, label='Theoretical curve')
    plt.xlabel(r'$1-\cos{(\theta)}$')
    plt.ylabel(r"$E'^{-1}$ [keV$^{-1}$]")
    plt.title(r"$E'^{-1}$ against $1-\cos{(\theta)}$ in the Compton scattering")
    plt.xlim(right=2)
    plt.legend()
    plt.grid()
    
    if option.upper()=='TEST':
        plt.savefig('../Test/plot.pdf')
    else:
        plt.savefig(f'../Plots/plot {option}.pdf')
    
    print(f'\nThe determined mass of the electron is m_e = ({1/a}±{cov[0, 0]/a**2}) keV/c².')
    
    print(f'\nThe determined incident energy is E = ({1/b}±{cov[1, 1]/b**2}) keV.')
    
main()
