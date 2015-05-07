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

	# Plot streamline and wall
	ax1.streamplot(x, y, v[0,:,:], v[1,:,:], density=0.6, color = 'k', linewidth=lw)
	ax1.imshow(~wall, extent=(0, nx, -L, L), alpha=0.5, interpolation='nearest', cmap=plt.cm.gray)
	ax2.quiver(x, y, v[0,:,:], v[1,:,:], angles = 'xy')
	ax2.imshow(~wall, extent=(0, nx, -L, L), alpha=0.5, interpolation='nearest', cmap=plt.cm.gray)
	plt.draw()