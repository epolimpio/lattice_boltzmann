#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Parameters of the simulation
nx = 30
ny = 11
deltav = 1e-5
tau = 1
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

# Exact Poiseuille flow
visc = (tau-0.5)/3
peak = deltav/tau/(2*visc)
L = (ny-2)/2
dist = np.linspace(-(L-0.5), (L-0.5), num = ny-2)
exact = -peak*(dist**2-L**2)

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
    v_mean = np.mean(v[0,:,:], axis=1)
    ax1.plot(dist, v_mean[1:ny-1])
    ax1.plot(dist, exact)
    #ax1.pcolor(v[0,:,:])
    plt.draw()
    
    # Add velocity due to pressure gradient
    v[0,:,:] = v[0,:,:] + deltav

    # Calculate equilibrium and collide
    dens_eq = model.calc_eq_dens(rho, v)
    dens = (1-omega)*dens + omega*dens_eq

    print t