from numpy import *

class LatticeBoltzmann:

    def __init__(self, nx, ny):

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

        # Density is in format (q, nx, ny) and ei is (q, dimension)
        # v is the velocity vector, with dimension (dimension, nx, ny)
        # rho is the density at each node, rho(nx, ny)

        # Make the dot product of ei and v.
        # Dimensions -> edotv(q, nx, ny)
        edotv = 3.0*dot(self.ei, v.transpose(1,0,2))

        # Calculate the square velocity at each node, usqr(nx, ny)
        vsqr = sum(v**2, axis = 0)

        # Calculate the equilibrium density
        for i in range(self.q):
            dens_eq[i,:,:] = self.w[i]*rho*(1+edotv[i] + 
                            0.5*edotv[i]**2.0 - 1.5*vsqr)

        return dens_eq

    def calc_v_rho(self, dens):

        # Densities are dens(q, nx, ny) and rho(nx, ny)
        rho = sum(dens, axis = 0)

        # ei(q, d) and dens(q, nx, ny). We need v(d, nx, ny).
        v = dot(self.ei.transpose(), dens.transpose(1,0,2))/rho

        return v, rho


