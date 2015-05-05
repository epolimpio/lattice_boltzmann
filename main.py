#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Parameters of the simulation
nx = 100
ny = 50
PGrad = 1e-5
tau = 0.9
omega = 1.0/tau
t_max = 10000

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

# Initiate graphics
plt.ion()
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

# Dynamics
for t in np.arange(t_max):
    
    dens = model.stream(dens)
    dens = model.bounce_back(dens, wall)

    v, rho = model.calc_v_rho(dens)
    
    # Plot the velocity profile
    ax1.clear()
    ax1.plot(np.arange(ny), np.mean(v[0,:,:], axis=1))
    #ax1.pcolor(v[0,:,:])
    plt.draw()
    
    # Add velocity due to pressure gradient
    v[0,fluid] = v[0,fluid] + PGrad

    # Calculate equilibrium and collide
    dens_eq = model.calc_eq_dens(rho, v)
    dens[:,fluid] = (1-omega)*dens[:,fluid] + omega*dens_eq[:,fluid]

    print t