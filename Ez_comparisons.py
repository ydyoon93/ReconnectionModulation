#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 11:01:40 2024

@author: young
"""

import happi
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from itertools import cycle
from string import ascii_lowercase as alc
plt.rcParams["mathtext.fontset"]='cm'

Sims = [happi.Open("2Dpinch{:02d}".format(i+1)) for i in range(3,18)]
n=np.size(Sims,0) ### number of simulation points
plot_range = (0,n) #### range of simulations to plot
plasma_cmap = [plt.get_cmap('viridis')(1. * i/n) for i in range(*plot_range)]
prop_cycle = cycler('color',plasma_cmap)
colors = cycle(prop_cycle.by_key()['color'])
plt.rc('axes', prop_cycle=prop_cycle)
#%% Calculate Max Ez
Ez_max = []
for S in Sims:
    Lx = S.namelist.Main.grid_length[0]
    mime = S.namelist.mass_ratio
    wcewpe = S.namelist.omega_ratio
    B0 = wcewpe*mime**.5
    vA = wcewpe/mime**.5
    Ez_max.append(np.max(np.array(S.Field(0,'Ez',subset={'x':Lx/2}).getData()),axis=1)/vA/B0)
    
Ez_max = np.array(Ez_max)
#%% Plot
fig,axes = plt.subplots(6,3,sharex='col',sharey=1,figsize=(6.4,4))
t = np.linspace(0,1000,Ez_max.shape[1])
for ax,Ez in zip(axes[-4:-7:-1,2],Ez_max[0:3]):
    ax.fill_between(t,Ez,color=next(colors))

for ax,Ez in zip(axes[-1:-7:-1,1],Ez_max[3:9]):
    ax.fill_between(t,Ez,color=next(colors))

for ax,Ez in zip(axes[-1:-7:-1,0],Ez_max[9:15]):
    ax.fill_between(t,Ez,color=next(colors))

for ax in axes.flatten():
    ax.set_ylim([0,0.9])
    ax.set_xlim([t[0],t[-1]])
    ax.hlines(0.2,xmin=t[0],xmax=t[-1],color='k',lw=0.3,linestyle='--')
    ax.hlines(0.5,xmin=t[0],xmax=t[-1],color='m',lw=0.3,linestyle='--')
[ax.set_xlabel(r'$t c /\lambda_i$') for ax in axes[-1,0:2]]
[ax.set_ylabel(r'$E_{z,\rm{max}}$') for ax in axes[:,0]]
axes[2,2].set_xlabel(r'$t c /\lambda_i$')
axes[2,2].tick_params(axis='x',reset=1,top=0)
### delete unused subplots ####
[fig.delaxes(axes[i,2]) for i in [3,4,5]]


### subplot labels #####
for ax,text in zip(axes.transpose().flatten()[0:15],
                   alc[0:15],
                   ):
    ax.text(0.07,0.66,'('+text+')',transform=ax.transAxes)
# plt.savefig('/Users/young/Dropbox/APCTP/Research/2023/2023_ReconnectionOnset/manuscript/Ez_max.pdf',dpi=300,bbox_inches='tight')
