#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import LBDyn.plot as LBPlot
import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Parameters of the simulation
nx = 51
ny = 27
deltav = 1e-2
tau = 5.0
omega = 1.0/tau
t_max = 2000

#####################################
model = LatticeBoltzmann(ny, nx)

# Initial conditions
v = np.zeros((model.d, ny, nx))
rho = np.ones((ny,nx))
dens = model.calc_eq_dens(rho, v)

# Define obstacle
# h = ny/2.0
# l = nx/8.0
# c = nx/2.0
# triangle_func = np.fromfunction(lambda y, x: y - (h - abs(h/l*(x-c))), (ny, nx))
# obstacle = (triangle_func < 0)

cx = nx/4.0
cy = ny/2.0
r = ny/4.0
obstacle = np.fromfunction(lambda y, x: (y-cy)**2 + (x-cx)**2 < r**2, (ny, nx))

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

    #v, rho = model.calc_v_rho(dens)

    # Mean velocity
    v_mean = np.mean(v[0,:,:], axis=1)
    
    # Add velocity due to pressure gradient
    v[0,fluid] = v[0,fluid] + deltav

    # Calculate equilibrium and collide
    dens_eq = model.calc_eq_dens(rho, v)
    dens[:,fluid] = (1-omega)*dens[:,fluid] + omega*dens_eq[:,fluid]

    print t, v_mean[ny/2]-v_mean[0]

# Plot a stream with color and a vector field
LBPlot.streamline_plot_2D(nx, ny, v, wall)
#LBPlot.poiseuille_profile_with_exact(ny, v_mean, tau, deltav)

plt.show()
