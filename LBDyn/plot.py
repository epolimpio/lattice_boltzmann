import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

def streamline_plot_2D(nx, ny, v, wall):
    
    # Initialize figure
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.tick_params(axis='both', which='major', labelsize=16)

    # Calculate protting parameters
    x = np.linspace(0,nx-1,nx)
    L = ny/2.0
    y = np.linspace(0,ny-1,ny)

    # Calculate norm
    speed = np.sqrt(np.sum(v**2, axis = 0))

    # Masking the data to plot only the fluid
    p = np.ma.array(v[0,:,:], mask=wall)
    q = np.ma.array(v[1,:,:], mask=wall)
    lw = np.ma.array(speed, mask=wall)
    
    # Plot colourmap
    c = ax1.contourf(x,y, lw, cmap = 'summer')
    b = plt.colorbar(c, orientation='vertical')
    b.set_label('Speed', size = 18, rotation = 270, labelpad = 20)
    b.ax.tick_params(axis='y', which='major', labelsize=16)
    ax1.set_xlabel('X', size = 18)
    ax1.set_ylabel('Y', size = 18)

    ax1.quiver(x, y, p, q, angles = 'xy')
    ax1.patch.set_facecolor('black')

    plt.draw()

def poiseuille_profile_with_exact(ny, v_mean, tau, deltav):

    # Initiate graphics
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.tick_params(axis='both', which='major', labelsize=16)

    # Exact Poiseuille flow
    visc = (tau-0.5)/3.0
    peak = deltav/tau/visc/2.0
    L = ny/2.0 - 1.0
    dist = np.linspace(-(L+0.5), L+0.5, num = ny)
    exact = -peak*(dist**2-(L+0.5)**2)

    ax1.scatter(v_mean, dist, color = 'k', marker = 'x', label = "Calculated")
    ax1.plot(exact, dist, color = 'r', linewidth = 2.0, label = "Exact")
    ax1.set_xlabel('Velocity', size = 16)
    ax1.set_ylabel('Dist. from the center of the pipe', size = 16)
    ax1.legend()

def poiseuille_profile_with_exact_unbias(ny, v_mean, tau, deltav):

    # Initiate graphics
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.tick_params(axis='both', which='major', labelsize=16)

    # Exact Poiseuille flow
    visc = (tau-0.5)/3.0
    peak = deltav/tau/visc/2.0
    L = ny/2.0 - 1.0
    dist = np.linspace(-(L-0.5), L-0.5, num = ny-2)
    exact = -peak*(dist**2-(L+0.5)**2)

    exact = exact - np.min(exact)
    v_mean = v_mean - np.min(v_mean[1:ny-1])

    ax1.scatter(v_mean[1:ny-1], dist, color = 'k', marker = 'x', label = "Calculated")
    ax1.plot(exact, dist, color = 'r', linewidth = 2.0, label = "Exact")
    ax1.set_xlabel('Velocity - Bias', size = 16)
    ax1.set_ylabel('Dist. from the center of the pipe', size = 16)
    ax1.legend()