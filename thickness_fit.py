#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 14:14:08 2023

@author: young
"""

import os
import happi
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from cycler import cycler
plt.rcParams["mathtext.fontset"]='cm'
plt.rcParams["font.size"]=12

def get_thickness_and_peak_pressure(sim_string):
    S = happi.Open(sim_string)
    mass_ratio = S.namelist.mass_ratio
    omega_ratio = S.namelist.omega_ratio
    T = S.namelist.T
    dx = S.namelist.dx
    beta = S.namelist.beta
    xmax = S.namelist.xmax
    x = np.arange(0,xmax)*dx
    x -= x[int(len(x)/2)]
    l = S.namelist.l
    B0 = mass_ratio**0.5*omega_ratio
    every = S.namelist.every
    By_Field = S.Field(0,"By")
    n_i_Bin = S.ParticleBinning(0)
    U_i_Bin = [S.ParticleBinning(i) for i in range(1,4)]
    P_i_Bin = [S.ParticleBinning(i) for i in range(4,7)]
    n_e_Bin = S.ParticleBinning(10)
    U_e_Bin = [S.ParticleBinning(i) for i in range(11,14)]
    P_e_Bin = [S.ParticleBinning(i) for i in range(14,17)]
    timesteps = By_Field.getTimesteps()
    def fit_tanh(x,B,l):
        return B * np.tanh(x/l)
    def fit_sech2(x,P,l,P0):
        return P - P0 * np.tanh(x/l)**2
    
    def savgol(x):
        return x
        # return savgol_filter(x, 20, 4)
    #%% Find Thickness
    thickness = [curve_fit(fit_tanh,x,savgol(By_Field.getData(t)[0][1:]),
                           p0=[By_Field.getData(0)[0][10],l/2.5])[0][1]
                 for t in timesteps] ### time-dependent thickness
    
    # plt.plot(thickness)
    
        
    #%% Pressure
    
    def get_pressure(timestep):
        n_i = savgol(n_i_Bin.getData(timestep)[0]) + .1  ### add small number to prevent NaNs during division
        n_e = savgol(n_e_Bin.getData(timestep)[0]) + .1
        return savgol(np.mean([P_i_Bin[i].getData(timestep)[0] 
                        - U_i_Bin[i].getData(timestep)[0]**2/n_i/mass_ratio
                        + P_e_Bin[i].getData(timestep)[0]
                        - U_e_Bin[i].getData(timestep)[0]**2/n_e
                        for i in range(1)],axis=0))
    peak_pressure = [curve_fit(fit_sech2,x,get_pressure(t),
                               p0=[T,l/4,0],
                               maxfev = 8000,
                               )[0][0]
                     # /(2 * T / beta * mass_ratio)
                     /(B0**2/2)
                     ### factors to convert to beta ####
                     for t in timesteps]
    return np.array([thickness,peak_pressure])

sim_list = ["1Dpinch%02i"%(i) for i in range(1,19)]
metadata = np.array([get_thickness_and_peak_pressure(sim_name)
                     for sim_name in sim_list])

#%% Plot
n=np.size(metadata,0) ### number of simulation points
timesteps=np.size(metadata,-1) ### number of timesteps
max_time=int(timesteps*1/2) ###maximum timestep to plot
plot_range = (3,n) #### range of simulations to plot
plt.figure(figsize=(6.4,4))
plasma_cmap = [plt.get_cmap('viridis')(1. * i/n) for i in range(*plot_range)]
plt.rc('axes', prop_cycle=cycler('color',plasma_cmap))
[plt.loglog(metadata[i,0,:max_time],metadata[i,1,:max_time],alpha=0.25) for i in range(*plot_range)]
start = [plt.scatter(metadata[i,0,0],metadata[i,1,0],s=50) for i in range(*plot_range)]
[plt.scatter(metadata[i,0,max_time],metadata[i,1,max_time],marker='D') for i in range(*plot_range)]
plt.xlabel(r'$\lambda/d_i$')
plt.ylabel(r'$\hat{\beta}_{xx}$')


"""
Draw arrows to represent simulations that are presented later in 2D
"""
origins = [[0.54,0.06],
           [0.53,0.65],
           [0.75,0.5],
           [0.75,0.65]]
diff    = [[-0.39,0.68],
           [-0.10,0.18],
           [-0.20,0.31],
           [-0.12,0.19]]
[plt.arrow(origins[0],origins[1],diff[0],diff[1],
            width = 0.007,
            ec = 'w',
            fc = 'k',
            alpha = 0.5,
            transform = plt.gca().transAxes) for origins,diff in zip(origins,diff)]
[plt.text(origins[0]+0.01,origins[1],i,
          transform = plt.gca().transAxes) for origins,i in zip(origins,['(d)','(c)','(b)','(a)'])]

"""
Plot theoretical predictions of where the current sheets will end up

The theory is based on double adiabatic theory, i.e., P/nB = const.

P_f / P_i = 1 + 1/beta
n_f / n_i = l_i / l_f
B_f / B_i = l_i / l_f

So, l_f = l_i * (1 + 1/beta)**-1/gamma

"""
gamma = 1.1 #adiabatic index
[plt.scatter(metadata[i,0,0] * (1 + 1/metadata[i,1,0])**(-1/gamma),
             (1 + metadata[i,1,0]),
             marker="*",
             edgecolors='k',
             linewidths=0.5) for i in range(*plot_range)]
lambda_ = np.linspace(1.01,50,1000)
lambda2 = np.linspace(0.1,1.01,1000)
lambda3 = np.linspace(0.0001,0.101,1000)
alpha = 0.1
plt.fill_between(lambda_, 1/(lambda_**gamma-1),30,alpha=alpha,color='r',linewidth=0)
plt.fill_between(lambda_, 1/((lambda_/0.1)**gamma-1),1/(lambda_**gamma-1),alpha=alpha,color='y',linewidth=0)
plt.fill_between(lambda2, 1/((lambda2/0.1)**gamma-1),30,alpha=alpha,color='y',linewidth=0)
plt.fill_between(lambda_,1e-9, 1/((lambda_/0.1)**gamma-1),alpha=alpha,color='g',linewidth=0)
plt.fill_between(lambda2,1e-9, 1/((lambda2/0.1)**gamma-1),alpha=alpha,color='g',linewidth=0)
plt.fill_between(lambda3,1e-9,30,alpha=alpha,color='g',linewidth=0)
plt.xlim([0.03,20])
plt.ylim([0.025,2.5])
# select_sim_num = np.array([13,17,10,11])
# [plt.arrow(metadata[i,0,0],
#            metadata[i,1,0],
#            metadata[i,0,-1]-metadata[i,0,0],
#            metadata[i,1,-1]-metadata[i,1,0],
#            width = 0.01) for i in select_sim_num-1]
# gamma_1 = 5/3.
# gamma_2 = 1.
# Lambda_1 = np.linspace(0.1,0.5,100)
# Lambda_2 = np.linspace(0.05,0.1,100)
# beta_1   = 0.02*Lambda_1**(-gamma_1)
# beta_2   = 0.1*Lambda_2**(-gamma_2)
# plt.plot(Lambda_1,beta_1,c='k',linestyle='--',alpha=0.35)
# plt.plot(Lambda_2,beta_2,c='k',linestyle='--',alpha=0.35)
# plt.text(0.05,0.5,'$\gamma=5/3$')
# plt.text(0.04,1.2,'$\gamma=1$')

# plt.savefig('../../manuscript/pinch_pathways_xxpressure.pdf',dpi=300,bbox_inches='tight')
