#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:17:25 2020

@author: young

smilei namelist for runaway scattering simulation
"""
import numpy as np
#%% Smilei Parameters
mass_ratio = 100. #m_i/m_e
omega_ratio = 1./2.  #omega_ce/omega_pe
beta = 1./32.
l = 4. #sheath length scale in terms of d_i
dx = 0.0125*l
dy = 0.0125*l
xwindow = 6.4*l
ywindow = 12.8*l
xmax = xwindow/dx
ymax = ywindow/dy
B0 = np.sqrt(mass_ratio)*omega_ratio
b  = 0. * B0 #perturbation strength... No perturbation here (pinching)
n0 = mass_ratio
nb = 0
vA = omega_ratio/mass_ratio**.5 #Alfven speed / c
ppc = 200 #particles per cell
CFL = dx*dy/(dx**2+dy**2)**.5
dt = 0.95 * CFL
every = int(l*1000/200*(1/dt))
save_every = every
#%% Smilei namelist
Main(
    geometry = "2Dcartesian",
    number_of_cells = [xmax,ymax],
    cell_length = [dx,dy],
    simulation_time = l*1000,
    timestep = dt,
    number_of_patches = [xmax/8,ymax/8],
    EM_boundary_conditions = [["silver-muller"],["periodic"]],
    print_every=every
)
LoadBalancing(
    initial_balance = True,
)
Vectorization(
    mode = "on",
)
#%% Background Field
def By(x,y):
    return B0*np.tanh((x-1/2*xwindow)/l)
def By_(x,y): 
    xx = (x - 1/2 * xwindow)/l
    yy = (y - 1/2 * ywindow)/l
    return By(x,y) - b * xx  * np.exp(-(xx**2+yy**2)/2)
def Bx_(x,y):
    xx = (x - 1/2 * xwindow)/l
    yy = (y - 1/2 * ywindow)/l
    return b * yy * np.exp(-(xx**2+yy**2)/2)
ExternalField(
    field = "By",
    profile = By_,
)
ExternalField(
    field = "Bx",
    profile = Bx_,
)
#%% Species
def n(x,y):
    return n0
def j(x,y):
    return B0/l*np.cosh((x-1/2*xwindow)/l)**(-2)
def V_i(x,y):
    return j(x,y)/n(x,y)/2
def V_e(x,y):
    return -j(x,y)/n(x,y)/2
T = omega_ratio**2/4*beta
Species(
    name = "ions",
    position_initialization = "random",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = ppc,
    mass = mass_ratio,
    number_density = n,
    charge = 1,
    mean_velocity = [0,0,V_i],
    temperature = [T],
    boundary_conditions = [["remove"],["periodic"]],
)
Species(
    name = "electrons",
    position_initialization = "random",
    momentum_initialization = "maxwell-juettner",
    particles_per_cell = ppc,
    mass = 1,
    number_density = n,
    charge = -1,
    mean_velocity = [0,0,V_e],
    temperature = [T],
    boundary_conditions = [["remove"],["periodic"]],
)
ParticleInjector(
        name = "ions1",
        species = "ions",
        box_side = "xmin",
)
ParticleInjector(
        name = "ions2",
        species = "ions",
        box_side = "xmax",
)
ParticleInjector(
        name = "electrons1",
        species = "electrons",
        box_side = "xmin",
)
ParticleInjector(
        name = "electrons2",
        species = "electrons",
        box_side = "xmax",
)
#%% Diagnostics
DiagScalar(
    every = save_every,
    )
DiagFields(
    every = save_every,
    time_average = 100,
    fields = ["Ex","Ey","Ez",
            "Bx","By","Bz",
            "Jx_ions","Jy_ions","Jz_ions",
            "Jx_electrons","Jy_electrons","Jz_electrons",
            "Rho_ions","Rho_electrons"
            ]
    )
DiagParticleBinning(
        name = "ion_n",
        deposited_quantity = "weight",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Ux",
        deposited_quantity = "weight_px",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Uy",
        deposited_quantity = "weight_py",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Uz",
        deposited_quantity = "weight_pz",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Pxx",
        deposited_quantity = "weight_vx_px",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Pyy",
        deposited_quantity = "weight_vy_py",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Pzz",
        deposited_quantity = "weight_vz_pz",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Pxy",
        deposited_quantity = "weight_vx_py",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Pyz",
        deposited_quantity = "weight_vy_pz",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "ion_Pzx",
        deposited_quantity = "weight_vx_pz",
        every = save_every,
        time_average = 100,
        species = ["ions"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )



DiagParticleBinning(
        name = "electron_n",
        deposited_quantity = "weight",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Ux",
        deposited_quantity = "weight_px",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Uy",
        deposited_quantity = "weight_py",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Uz",
        deposited_quantity = "weight_pz",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Pxx",
        deposited_quantity = "weight_vx_px",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Pyy",
        deposited_quantity = "weight_vy_py",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Pzz",
        deposited_quantity = "weight_vz_pz",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Pxy",
        deposited_quantity = "weight_vx_py",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Pyz",
        deposited_quantity = "weight_vy_pz",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )
DiagParticleBinning(
        name = "electron_Pzx",
        deposited_quantity = "weight_vx_pz",
        every = save_every,
        time_average = 100,
        species = ["electrons"],
        axes = [
            ["x",0,xwindow,xmax+1],
            ["y",0,ywindow,ymax+1],
            ]
        )

