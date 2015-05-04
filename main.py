#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Parameters of the simulation
nx = 100
ny = 50
PGrad = 0.00002
tau = 0.5
omega = 1.0/tau
t_max = 1000

#####################################
model = LatticeBoltzmann(ny, nx)
# initial conditions
v = np.zeros((model.d, ny, nx))
rho = np.ones((ny,nx))
dens = model.calc_eq_dens(rho, v)

# define boundaries (boolean array)

wall = np.zeros((ny, nx)).astype(bool)
wall[0,:] = np.ones(nx).astype(bool)
wall[ny-1,:] = np.ones(nx).astype(bool)
fluid = np.logical_not(wall)
#model.wall = wall

# Initiate graphics
plt.ion()
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

# evolute
for t in np.arange(t_max):
	dens = model.move(dens)
	dens = model.bounce_back(dens, wall)
	v, rho = model.calc_v_rho(dens)
	ax1.clear()
	ax1.plot(np.arange(ny), np.mean(v[0,:,:], axis=1))
	#ax1.pcolor(v[0,:,:])
	plt.draw()
	v[0,:,:] = v[0,:,:] + PGrad
	dens_eq = model.calc_eq_dens(rho, v)
	#dens[:,fluid] = (1-omega)*dens[:,fluid] + omega*dens_eq[:,fluid]
	dens = (1-omega)*dens + omega*dens_eq