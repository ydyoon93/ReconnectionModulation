#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:17:25 2020

@author: young

smilei namelist for runaway scattering simulation
"""
import numpy as np
#%% Smilei Parameters
mass_ratio = 1600. #m_i/m_e
omega_ratio = 5.  #omega_ce/omega_pe
l = 16. #sheath length scale in terms of d_i
window = 10. #window length in terms of l
beta = 1./1. #peak beta
xmax = 2**12
dx = l*window/xmax
every = 500
save_every = 500
#%% Smilei namelist
Main(
    geometry = "1Dcartesian",
    number_of_cells = [xmax],
    cell_length = [dx],
    simulation_time = 1600,
    timestep_over_CFL = 0.95,
    number_of_patches = [xmax/8],
    EM_boundary_conditions = [["silver-muller"]],
    print_every=every
)
LoadBalancing(
    initial_balance = True,
)
#%% Background Field
def By(x):
    return np.sqrt(mass_ratio)*omega_ratio*np.tanh(x/l-1/2*window)
ExternalField(
    field = "By",
    profile = By
)
#%% Species
def n(x):
    return mass_ratio
def j(x):
    return np.sqrt(mass_ratio)*omega_ratio/l*np.cosh(x/l-1/2*window)**(-2)
def V_i(x):
    return j(x)/n(x)/2
def V_e(x):
    return -j(x)/n(x)/2
def ppc(x):
    #return 1000*np.cosh(x/l-1/2*window)**(-2)+100
    return 400
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
    boundary_conditions = [["remove"]],
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
    boundary_conditions = [["remove"]],
)
ParticleInjector(
    name = "i_xmin",
    species = "ions",
    box_side = "xmin", 
)
ParticleInjector(
    name = "i_xmax",
    species = "ions",
    box_side = "xmax", 
)
ParticleInjector(
    name = "e_xmin",
    species = "electrons",
    box_side = "xmin", 
)
ParticleInjector(
    name = "e_xmax",
    species = "electrons",
    box_side = "xmax", 
)
#%% Diagnostics
DiagScalar(
    every = save_every,
    )
DiagFields(
    every = save_every,
    )
DiagParticleBinning(
    deposited_quantity = "weight",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_px",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_py",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_pz",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vx_px",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vy_py",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vz_pz",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vx_py",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vx_pz",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vy_pz",
    every = save_every,
    species = ["ions"],
    axes = [["x",0,xmax*dx,xmax]]
)


DiagParticleBinning(
    deposited_quantity = "weight",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_px",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_py",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_pz",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vx_px",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vy_py",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vz_pz",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vx_py",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vx_pz",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)
DiagParticleBinning(
    deposited_quantity = "weight_vy_pz",
    every = save_every,
    species = ["electrons"],
    axes = [["x",0,xmax*dx,xmax]]
)



