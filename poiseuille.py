#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Parameters of the simulation
nx = 100
ny = 32
deltav = 1e-2
tau = 1.0
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
visc = (tau-0.5)/3.0
peak = deltav/tau/visc/2.0
L = ny/2.0 - 1.0
dist = np.linspace(-(L-0.5), L-0.5, num = ny-2)
exact = -peak*(dist**2-(L+0.5)**2)

# Initiate graphics
fig = plt.figure()
fig2 = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
ax3 = fig2.add_subplot(2,1,1)
ax4 = fig2.add_subplot(2,1,2)

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

    print t, exact[ny/2-1], v_mean[ny/2]-v_mean[0]

# Plot a stream with color and a vector field
x = np.linspace(0,nx,nx)
y = np.linspace(-(L+1.0),(L+1.0),ny)
o = v[0,:,:]
p = v[1,:,:]
speed = np.sqrt(o*o + p*p)
lw = 5*speed/speed.max()

ax3.streamplot(x, y, o, p, density=0.6, color = 'k', linewidth=lw)
ax3.imshow(~wall, extent=(0, nx, -(L+1.0), (L+1.0)), alpha=0.5, interpolation='nearest', cmap=plt.cm.gray)
ax4.quiver(x, y, o, p, angles = 'xy')
ax4.imshow(~wall, extent=(0, nx, -(L+1.0), (L+1.0)), alpha=0.5, interpolation='nearest', cmap=plt.cm.gray)
plt.draw()

# Plot the velocity profile
ax1.plot(dist, v_mean[1:ny-1])
ax1.plot(dist, exact)
v_mean = v_mean - np.min(v_mean[1:ny-1])
exact = exact - np.min(exact)
ax2.plot(dist, v_mean[1:ny-1])
ax2.plot(dist, exact)

plt.show()
