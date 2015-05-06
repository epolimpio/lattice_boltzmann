#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Parameters of the simulation
nx = 30
ny = 17
<<<<<<< HEAD:main.py
deltav = 1e-2
tau = 1.0
=======
deltav = 1e-5
tau = 5
>>>>>>> 69be368235e15b70524839b2cd8bdfdf358287c8:poiseuille.py
omega = 1.0/tau
t_max = 100

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
<<<<<<< HEAD:main.py
visc = (tau-0.5)/3
peak = deltav/tau/visc
L = ny/2
dist = np.linspace(-L, L, num = ny)
exact = -peak*(dist**2-(L-0.5)**2)

####################################
Y, X = np.mgrid[-L:L, 0:nx]
o = -1 - X**2 + Y
p = 1 + X - Y**2
speed = np.sqrt(o*o + p*p)
lw = 5*speed/speed.max()

plt.streamplot(X, Y, o, p, density=0.6, color='k', linewidth=2)
plt.show()

##############################################
=======
visc = (tau-0.5)/3.0
peak = deltav/tau/visc/2.0
L = ny/2.0 - 1.0
dist = np.linspace(-(L-0.5), L-0.5, num = ny-2)
exact = -peak*(dist**2-(L+0.5)**2)
>>>>>>> 69be368235e15b70524839b2cd8bdfdf358287c8:poiseuille.py

# Initiate graphics
fig = plt.figure()
fig2 = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

# Dynamics
for t in np.arange(t_max):
    
    dens = model.stream(dens)
    dens = model.bounce_back(dens, wall)

    v, rho = model.calc_v_rho(dens)

<<<<<<< HEAD:main.py
=======
    # Mean velocity
    v_mean = np.mean(v[0,:,:], axis=1)
    
>>>>>>> 69be368235e15b70524839b2cd8bdfdf358287c8:poiseuille.py
    # Add velocity due to pressure gradient
    v[0,fluid] = v[0,fluid] + deltav

    # Calculate equilibrium and collide
    dens_eq = model.calc_eq_dens(rho, v)
    dens = (1-omega)*dens + omega*dens_eq

<<<<<<< HEAD:main.py
    #print t, exact[ny/2], v_mean[ny/2]


###################################
Y,X = np.mgrid[0:nx,0:(ny-1)]
prov = v.transpose(0,2,1)
o = prov[1,:,:]
p = prov[0,:,:]
print X, p
speed = np.sqrt(o*o + p*p)
lw = 5*speed/speed.max()
=======
    print t, exact[ny/2-1], v_mean[ny/2]-v_mean[0]

# Plot the velocity profile
ax1.plot(dist, v_mean[1:ny-1])
ax1.plot(dist, exact)
v_mean = v_mean - np.min(v_mean[1:ny-1])
exact = exact - np.min(exact)
ax2.plot(dist, v_mean[1:ny-1])
ax2.plot(dist, exact)

plt.show()
>>>>>>> 69be368235e15b70524839b2cd8bdfdf358287c8:poiseuille.py

plt.streamplot(X, Y, o, p, density=0.6, color='k', linewidth=1)
plt.draw()
