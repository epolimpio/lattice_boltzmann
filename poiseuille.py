#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import LBDyn.plot as LBPlot
import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Parameters of the simulation
nx = 50
ny = 27
deltav = 1e-2
tau = 7.0
omega = 1.0/tau
t_max = 5000

#####################################
model = LatticeBoltzmann(ny, nx)

# Initial conditions
v = np.zeros((model.d, ny, nx))
rho = np.ones((ny,nx))
dens = model.calc_eq_dens(rho, v)

# Define boundaries (boolean array)
wall = np.zeros((ny, nx)).astype(bool)
wall[0,:] = np.ones(nx).astype(bool)
wall[ny-1,:] = np.ones(nx).astype(bool)
fluid = np.logical_not(wall)

# Dynamics
for t in np.arange(t_max):
    
    dens = model.stream(dens)
    dens = model.bounce_back(dens, wall)

    v, rho = model.calc_v_rho(dens)

    # Mean velocity
    v_mean = np.mean(v[0,:,:], axis=1)
    
    # Add velocity due to pressure gradient
    v[0,fluid] = v[0,fluid] + deltav

    # Calculate equilibrium and collide
    dens_eq = model.calc_eq_dens(rho, v)
    dens = (1-omega)*dens + omega*dens_eq

    print t, v_mean[ny/2]-v_mean[0]

# Printing in a file
data_file = open("data_tau_{}.txt".format(tau), "w+")
for i in range(ny):
     data_file.write('{!r}\t{!r}\n'.format(i, v_mean[i]))
data_file.close()

# Plot a stream with color and a vector field
LBPlot.streamline_plot_2D(nx, ny, v, wall)
LBPlot.poiseuille_profile_with_exact(ny, v_mean, tau, deltav)

plt.show()
