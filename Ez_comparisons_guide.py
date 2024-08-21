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
plt.rcParams["font.size"]=12

Sims = [happi.Open("2DpinchGuide{:02d}".format(i+1)) for i in range(6,12)]
n=np.size(Sims,0) ### number of simulation points
plot_range = (0,n) #### range of simulations to plot
plasma_cmap = [plt.get_cmap('viridis')(1. * i/15) for i in range(3,9)]
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
fig,axes = plt.subplots(6,3,sharex='col',sharey=0,squeeze=0,figsize=[8,4])
t = np.linspace(0,4000,Ez_max.shape[1])
for ax,Ez in zip(axes[-1:-7:-1,0],Ez_max):
    ax.fill_between(t,Ez,color=next(colors))

for ax in axes[:,0]:
    ax.set_ylim([0,0.9])
    ax.set_xlim([t[0],t[-1]])
    ax.hlines(0.2,xmin=t[0],xmax=t[-1],color='k',lw=0.3,linestyle='--')
    ax.hlines(0.5,xmin=t[0],xmax=t[-1],color='m',lw=0.3,linestyle='--')
    ax.set_ylabel(r'$E_{z,\rm{max}}$')
axes[-1,0].set_xlabel(r'$t \omega_{pi}$')

### subplot labels #####
for ax,text in zip(axes.transpose().flatten()[0:6],
                   alc[0:6],
                   ):
    ax.text(0.07,0.7,'('+text+')',transform=ax.transAxes)
    
### plot three timesteps for the simulation 07 and 09 ###
gs1 = [axes[i,1].get_gridspec() for i in np.arange(3)*2]
gs2 = [axes[i,2].get_gridspec() for i in np.arange(3)*2]
### delete unused subplots ####
[fig.delaxes(axes.transpose().flatten()[i]) for i in range(6,18)]
ax1 = [fig.add_subplot(gs1[0][i:i+2,1]) for i in np.arange(3)*2] 
ax2 = [fig.add_subplot(gs1[0][i:i+2,2]) for i in np.arange(3)*2] 
[ax1[i].set_xticklabels([]) for i in [0,1]]
[ax2[i].set_xticklabels([]) for i in [0,1]]
[ax2[i].set_yticklabels([]) for i in range(3)]

### plot time-dependent current denisty for the guide field case ####
for simName, axCol,time_idx,label in zip(['2DpinchGuide08','2DpinchGuide09'],
                                         [ax1,ax2],
                                         [[0,80,150],[0,80,199]],
                                         ['(g)','(h)']):
    S = happi.Open(simName)
    Jz = S.Field(0,'Jz_ions+Jz_electrons')
    tsteps = Jz.getTimesteps()
    times  = Jz.getTimes()
    xaxis = Jz.getAxis('x')
    yaxis = Jz.getAxis('y')
    n0    = S.namelist.mass_ratio
    for axes,t in zip(axCol,time_idx):
        c = axes.imshow(Jz.getData(tsteps[t])[0]/n0,aspect='equal',
                        extent = [yaxis[0],yaxis[-1],xaxis[0],xaxis[-1]],
                        origin = 'lower',
                        cmap = 'Reds',
                        vmin = 0,
                        )
        ###### Calculate flux functions and plot their contours
        By = S.Field(0,'By')
        Bx = S.Field(0,'Bx')
        dx = S.namelist.dx
        dy = S.namelist.dy
        Nx = S.namelist.xmax
        Ny = S.namelist.ymax
        Az = -np.cumsum(By.getData(tsteps[t])[0],axis=0)*dx
        f  = np.cumsum(Bx.getData(tsteps[t])[0][0],axis=0)*dy
        f  = np.tile(f,(int(Nx+1),1))
        Az += f 
        axes.contour(yaxis,xaxis,Az,levels=10,
                     colors = 'k',
                     linewidths=0.5,
                     linestyles='solid')
        #######
        # axes.set_ylim([xaxis[-1]/4,xaxis[-1]*3/4])
        plt.colorbar(c,ax=axes)
        axes.text(0.01,1.1,'$t\omega_{pi}=%.2f$'%times[t],transform=axes.transAxes,c='k')
    axCol[0].text(0.01,1.35,label,transform=axCol[0].transAxes)
# [ax.set_ylabel(r'$x/d_i$') for ax in ax1]
[ax.set_xlabel(r'$y/d_i$') for ax in [ax1[-1],ax2[-1]]]

# plt.tight_layout()
# plt.savefig('/Users/young/Dropbox/APCTP/Research/2023/2023_ReconnectionOnset/manuscript/Ez_max_guide.pdf',dpi=300,bbox_inches='tight')
