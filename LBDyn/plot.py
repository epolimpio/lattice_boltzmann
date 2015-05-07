import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

def streamline_plot_2D(nx, ny, v, wall):
    
    # Initialize figure
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    # Calculate protting parameters
    x = np.linspace(0,nx,nx)
    L = ny/2.0
    y = np.linspace(-L,L,ny)
    speed = np.sqrt(np.sum(v**2, axis = 0))
    lw = 5*speed/speed.max()
    p = np.ma.array(v[0,:,:], mask=wall)
    q = np.ma.array(v[1,:,:], mask=wall)
    lw = np.ma.array(lw, mask=wall)

    # Plot streamline and wall
    ax1.streamplot(x, y, p, q, density=0.6, color = '#004C99', linewidth=lw)
    ax1.imshow(~wall, extent=(0, nx, -L, L), alpha=1, interpolation='nearest', cmap=plt.cm.gray)
    c = ax2.contourf(x,y,lw, cmap = 'summer')
    b = plt.colorbar(c, orientation='vertical')
    ax2.quiver(x, y, p, q, angles = 'xy')
    ax2.imshow(~wall, extent=(0, nx, -L, L), alpha=1, interpolation='nearest', cmap=plt.cm.gray)

    plt.draw()

def poiseuille_profile_with_exact(ny, v_mean, tau, deltav):

    # Initiate graphics
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    # Exact Poiseuille flow
    visc = (tau-0.5)/3.0
    peak = deltav/tau/visc/2.0
    L = ny/2.0 - 1.0
    dist = np.linspace(-(L-0.5), L-0.5, num = ny-2)
    exact = -peak*(dist**2-(L+0.5)**2)

    # Plot the velocity profile
    ax1.plot(v_mean[1:ny-1], dist, 'k', linewidth = 2.0, label = "Calculated")
    ax1.plot(exact, dist, 'r', linewidth = 2.0, label = "Exact")
    ax1.set_title("Abs values")
    ax1.legend()

    # Plot the unbiased variables
    v_mean = v_mean - np.min(v_mean[1:ny-1])
    exact = exact - np.min(exact)
    ax2.set_title("Unbiased")
    ax2.scatter(v_mean[1:ny-1], dist, color = 'k', marker = 'x', label = "Calculated")
    ax2.plot(exact, dist, color = 'r', linewidth = 2.0, label = "Exact")
    ax2.legend()
