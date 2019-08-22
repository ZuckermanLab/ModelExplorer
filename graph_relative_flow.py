import numpy as np 
from matplotlib import pyplot as plt 
import matplotlib.ticker as mtick

D = np.genfromtxt("antiporter_flows", skip_header=1)

fig = plt.figure(figsize=(10, 10))
#plt.suptitle("Model 1800 flux:", size = 30)
# ax1 = fig.add_subplot('111')  # top left 
# ax1.plot(D[:, 0], np.abs(D[:, 2]/D[:, 1, ]), label="|Substrate/Ion flux| [into cell]", linewidth=3, linestyle='--')
# #ax1.plot(D[:, 0], D[:, 2]/D[:, 3], label="Substrate/Decoy flux [into cell]", linewidth=3)  # proofreading only
# #ax1.set_ylim([1e-2,1e2])
# #ax1.set_yscale('log')  # proofreading only
# ax1.set_title("|Relative flux|", size=28)
# ax1.set_ylabel('|Relative flux|', size = 26)
# fig.gca().set_xlabel(
#     r'd$\mu$ ion (change in chemical potential of driving ion) [kT]', size=26)
# ax1.tick_params(axis='both', which='major', labelsize=20)
# ax1.tick_params(axis='both', which='minor', labelsize=18)
# plt.grid(b=True, which='both', axis='both')
# plt.legend(fontsize='large')

ax2 = fig.add_subplot('111')  # (122) for top right 
ax2.set_title("Model 845000 flux:", size=38)
ax2.plot(D[:, 0], D[:, 1]/D[0, 1], label=r'$J_{ion}$',
         linewidth=40, linestyle='--', color='Blue', alpha=0.7)
ax2.plot(D[:, 0], D[:, 2]/D[0, 1], label=r'$J_{substrate}$',
         linewidth=20, linestyle='-', color='Green')
#ax2.plot(D[:, 0], D[:, 3]/D[0, 1], label=r'$J_{decoy}$', linewidth=5, color='Red')  # proofreading only

ax2.set_ylabel(
    r'$J$ (flux) [arb. units]', size=38)
#ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
fig.gca().set_xlabel(
    r'$\Delta\mu_{ion}$ (chemical potential difference of driving ion) [$k_{B}T$]', size=38)
#ax2.set_yscale('log')
#ax1.set_ylim([-2,2])
plt.grid(b=True, which='both', axis='both')
ax2.tick_params(axis='both', which='major', labelsize=28)
ax2.tick_params(axis='both', which='minor', labelsize=24)
plt.legend(prop={'size': 32})
plt.tight_layout()
plt.show()
