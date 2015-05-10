#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import LBDyn.plot as LBPlot
import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Parameters of the simulation
nx = 101
ny = 21
deltav = 1e-5
tau = 1.0
omega = 1.0/tau
t_max = 10000

#####################################
model = LatticeBoltzmann(ny, nx)

# Initial conditions
v = np.zeros((model.d, ny, nx))
rho = np.ones((ny,nx))
dens = model.calc_eq_dens(rho, v)

# Define tringular obstacle
h = ny/2.0
l = nx/8.0
c = nx/2.0
triangle_func = np.fromfunction(lambda y, x: y - (h - abs(h/l*(x-c))), (ny, nx))
obstacle = (triangle_func < 0)

# Define boundaries (boolean array)
wall = np.zeros((ny, nx)).astype(bool)
wall[0,:] = np.ones(nx).astype(bool)
wall[ny-1,:] = np.ones(nx).astype(bool)
wall = np.logical_or(wall, obstacle) 
fluid = np.logical_not(wall)

# Dynamics
for t in np.arange(t_max):
    
    dens = model.stream(dens)
    v, rho = model.calc_v_rho(dens)
    
    dens = model.bounce_back(dens, wall)

    # Mean velocity
    v_mean = np.mean(v[0,:,:], axis=1)
    
    # Add velocity due to pressure gradient
    v[0,fluid] = v[0,fluid] + deltav

    # Calculate equilibrium and collide
    dens_eq = model.calc_eq_dens(rho, v)
    dens[:,fluid] = (1-omega)*dens[:,fluid] + omega*dens_eq[:,fluid]

    print t, v_mean[ny/2]

# Plot a stream with color and a vector field
LBPlot.streamline_plot_2D(nx, ny, v, wall)

plt.show()
