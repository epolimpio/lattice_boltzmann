#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann
import numpy as np

# Parameters of the simulation

nx = 20
ny = 7
PGrad = 0.1/nx
tau = 0.5

#####################################
model = LatticeBoltzmann(nx, ny)
# initial conditions
v = np.zero(...)
rho = np.ones(...)
dens = ...

# define boundaries (boolean array)

model.wall = wall

# evolute
for t in np.linspace(0, t_max, 1):
	dens_old = model.move(dens_new)
	v, rho = model.calc_v_rho(v, rho)
	v(0,:) = v(0,:) + 0.0001
	dens_eq = model.calc_eq_dens(rho, v)
	dens_new = (1-omega)*dens_old + omega*dens_eq

