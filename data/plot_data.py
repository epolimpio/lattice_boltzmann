import numpy as np
import matplotlib as matl
import matplotlib.pyplot as plt

# Initiate graphics
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=16)

deltav = 1e-2

all_tau = [1.0, 1.5, 2.0]
colors = 'krbgm'
i = 0

for tau in all_tau:

	# Read from file
	data = np.loadtxt('data_tau_' + str(tau) + '.txt')

	ny = len(data[:,0])
	v_mean = data[:,1]

	L = ny/2.0 - 1.0

	# Exact Poiseuille flow
	visc = (tau-0.5)/3.0
	peak = deltav/tau/visc/2.0
	dist = np.linspace(-(L+0.5), L+0.5, num = ny)
	exact = -peak*(dist**2-(L+0.5)**2)

	ax.scatter(v_mean, dist, color = colors[i], marker = 'x', label = r'$\tau$ = ' + str(tau))
	ax.plot(exact, dist, color = colors[i], linewidth = 2.0)

	i += 1

ax.legend()
ax.set_xlabel('Velocity', size = 16)
ax.set_xlim([0,6])
ax.set_ylabel('Dist. from the center of the pipe', size = 16)

plt.show()