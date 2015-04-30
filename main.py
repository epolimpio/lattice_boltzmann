#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import numpy as np

# Parameters of the simulation

nx = 20
ny = 7
PGrad = 0.0001
tau = 0.5
omega = 1.0/tau
t_max = 10

#####################################
model = LatticeBoltzmann(ny, nx)
# initial conditions
v = np.zeros((model.d, ny, nx))
rho = np.ones((ny,nx))
dens_new = np.ones ((model.q, ny, nx))
dens_old = np.ones ((model.q, ny, nx))
dens_eq = np.ones ((model.q, ny, nx))

# define boundaries (boolean array)

wall = np.zeros((ny, nx))
wall[0,:] = np.ones(nx)
wall[ny-1,:] = np.ones(nx)
model.wall = wall

# evolute
for t in np.linspace(0, t_max, 1):
	#dens_old = model.move(dens_new)
	v, rho = model.calc_v_rho(dens_old)
	v[0,:] = v[0,:] + PGrad
	dens_eq = model.calc_eq_dens(rho, v)
	dens_new = (1-omega)*dens_old + omega*dens_eq

