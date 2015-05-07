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


        c = ax2.contourf(x,y,lw, cmap = 'summer')#, np.linspace(-1,1,11))
        b = plt.colorbar(c, orientation='vertical')


	# Plot streamline and wall
	ax1.streamplot(x, y, p, q, density=0.6, color = '#004C99', linewidth=lw)
	ax1.imshow(~wall, extent=(0, nx, -L, L), alpha=1, interpolation='nearest', cmap=plt.cm.gray)
	ax2.quiver(x, y, p, q, angles = 'xy')
	ax2.imshow(~wall, extent=(0, nx, -L, L), alpha=1, interpolation='nearest', cmap=plt.cm.gray)
	plt.draw()
