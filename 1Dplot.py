#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:27:33 2024

@author: young
"""

import happi
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["mathtext.fontset"]='cm'
plt.rcParams["font.size"]=12
S           = happi.Open('1Dpinch16')
Tmax        = S.namelist.Main.simulation_time
dt          = S.namelist.Main.timestep
mime        = S.namelist.mass_ratio
wcewpe      = S.namelist.omega_ratio
Nx          = S.namelist.xmax
dx          = S.namelist.dx
Lx          = Nx * dx


B0          = mime**.5 * wcewpe
N0          = mime
J0          = N0
PB          = B0**2/2 #Asymptotic magnetic pressure
#%% Plot
fig,ax = plt.subplots(1,3,sharex=1,sharey=1,figsize=[7,1.7])
By_plot = ax[0].imshow(np.array(S.Field(0,'By',subset={'x':[Lx/4,3*Lx/4]}).getData())/B0,
                       origin='lower',
                       aspect='auto',
                       extent=[-Lx/4,Lx/4,0,Tmax],
                       cmap='PiYG')
plt.colorbar(By_plot,ax=ax[0])
Jz_plot = ax[1].imshow(np.array(S.Field(0,'Jz',subset={'x':[Lx/4,3*Lx/4]}).getData())/J0,
                       origin='lower',
                       aspect='auto',
                       extent=[-Lx/4,Lx/4,0,Tmax],
                       cmap='plasma')
plt.colorbar(Jz_plot,ax=ax[1])
P_plot = ax[2].imshow(np.array(S.ParticleBinning('(#4-#1*#1/#0/1600)+(#14-#11*#11/#10)',
                                                 subset={'x':[Lx/4,3*Lx/4]}).getData())/PB,
                       origin='lower',
                       aspect='auto',
                       extent=[-Lx/4,Lx/4,0,Tmax],
                       cmap='hot')
plt.colorbar(P_plot,ax=ax[2])

[ax.set_xlabel('$x/d_i$')  for ax in ax]
ax[0].set_ylabel('$t\omega_{pi}$')
for i,str_ in enumerate([r'(a) $B_y$',r'(b) $J_z$',r'(c) $\beta_{xx}$']):
    ax[i].text(0.1,0.8,str_,transform=ax[i].transAxes,
               bbox=dict(facecolor='w',alpha=0.75,edgecolor='None',boxstyle='round'))

# plt.savefig('../../manuscript/1Dpinch.pdf',dpi=300,bbox_inches='tight')
