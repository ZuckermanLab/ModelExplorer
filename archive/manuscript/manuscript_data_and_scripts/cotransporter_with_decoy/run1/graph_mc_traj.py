import numpy as np 
from matplotlib import pyplot as plt 
import matplotlib.ticker as mtick

# import data
D = np.genfromtxt("run2_cluster_data.dat", delimiter=',', invalid_raise=False)
x = D[:,0]
y = D[:,1]

# print(x)
# print(y)

fig = plt.figure(figsize=(15, 15))

ax = fig.add_subplot('111')  # (122) for top right 
ax.set_title("ModelExplorer trajectory in model space (Run 2)\n", size=36)

ax.set_ylabel('Monte Carlo \'energy\' [arb. flux units]', size=36)
ax.set_xlabel('Monte Carlo iteration number [n]', size=36)
ax.set_xlim(0,1e6)
ax.set_ylim(-1e-3, 1e-3)
plt.grid(b=True, which='both', axis='both')
ax.tick_params(axis='both', which='major', labelsize=24)
ax.tick_params(axis='both', which='minor', labelsize=20)
#ax.yaxis.set_minor_locator(mtick.MultipleLocator(5))
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
ax.plot(x, y, linewidth=1, color='Blue')

plt.savefig("test.png")
plt.show()
