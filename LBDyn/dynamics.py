from numpy import *
import numpy as np

class LatticeBoltzmann:

    def __init__(self, ny, nx):

        self.q = 9 # number of velocities
        self.d = 2 # dimension
        self.init_grid_var()

        # Define class constants for use in its functions
        self.nx = nx
        self.ny = ny

    def init_grid_var(self):

        # Define the d2q9 vectors for each nodal point
        self.ei = array([(0,0),
                         (1,0), (1,1), (0,1), (-1,1),
                         (-1,0), (-1,-1), (0,-1), (1,-1)])

        # Define the weights for the vectors
        self.w = 1.0/9.0 * ones(self.q)
        self.w[0] = 4.0/9.0
        self.w[[2, 4, 6, 8]] = 10./36.0

        return True

    def calc_eq_dens(self, rho, v):

        # Density is in format (q, ny, nx) and ei is (q, dimension)
        # v is the velocity vector, with dimension (dimension, ny, nx)
        # rho is the density at each node, rho(ny, nx)

        dens_eq = np.ones ((self.q, self.ny, self.nx))

        # Make the dot product of ei and v.
        # Dimensions -> edotv(q, ny, nx)
        edotv = 3.0*dot(self.ei, v.transpose(1,0,2))

        # Calculate the square velocity at each node, usqr(ny, nx)
        vsqr = sum(v**2, axis = 0)

        # Calculate the equilibrium density
        for i in range(self.q):
            dens_eq[i,:,:] = self.w[i]*rho*(1+edotv[i] + 0.5*edotv[i]**2.0 - 1.5*vsqr)

        return dens_eq

    def calc_v_rho(self, dens):

        # Densities are dens(q, ny, nx) and rho(ny, nx)
        rho = sum(dens, axis = 0)

        # ei(q, d) and dens(q, ny, nx). We need v(d, ny, nx).
        v = dot(self.ei.transpose(), dens.transpose(1,0,2))/rho

        return v, rho

    def move(self, dens):

        dens_new[0,:,:] = dens[0,:,:]

        dens_new[1,:,1:self.nx-1] = dens[1,:,0:self.nx-2]  
        dens_new[1,:,0] = dens[1,:,self.nx-1]

        dens_new[3,1:self.ny-1,:] = dens[3,0:self.ny-2,:]  
        dens_new[3,0,:] = dens[3,self.ny-1,:]

        dens_new[5,:,0:self.nx-2] = dens[5,:,1:self.nx-1]  
        dens_new[5,:,self.nx-1] = dens[5,:,0]

        dens_new[7,0:self.ny-2,:] = dens[7,1:self.ny-1,:]  
        dens_new[7,self.ny-1,:] = dens[7,0,:]

        dens_new[2,1:self.ny-1,1:self.nx-1] = dens[2,0:self.ny-2,0:self.nx-2]
        dens_new[2,:,0] = dens[2,:,self.nx-1]  
        dens_new[2,0,:] = dens[2,self.ny-1,:]

        dens_new[4,1:self.ny-1,0:self.nx-2] = dens[4,0:self.ny-2,1:self.nx-1]
        dens_new[4,:,self.nx-1] = dens[4,:,0]  
        dens_new[4,0,:] = dens[4,self.ny-1,:]

        dens_new[6,0:self.ny-2,0:self.nx-2] = dens[6,1:self.ny-1,1:self.nx-1]
        dens_new[6,:,self.nx-1] = dens[6,:,0]  
        dens_new[6,self.ny-1,:] = dens[6,0,:]

        dens_new[8,0:self.ny-2,1:self.nx-1] = dens[8,1:self.ny-1,0:self.nx-2]
        dens_new[8,:,0] = dens[8,:,self.nx-1]  
        dens_new[8,self.ny-1,:] = dens[8,0,:]

        return dens_new


