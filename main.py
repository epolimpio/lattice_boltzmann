#!/usr/bin/python

from LBDyn.dynamics import LatticeBoltzmann

# Parameters of the simulation

nx = 20
ny = 7
PGrad = 0.1/nx
tau = 0.5

test = LatticeBoltzmann(nx, ny, tau, PGrad)
