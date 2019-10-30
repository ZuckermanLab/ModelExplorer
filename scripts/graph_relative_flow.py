import numpy as np 
from matplotlib import pyplot as plt 
import matplotlib.ticker as mtick


### to plot flux (stoichiometry)
# D = np.genfromtxt("flux_829000", skip_header=1)

# fig = plt.figure(figsize=(15, 15))

# ax2 = fig.add_subplot('111')  # (122) for top right 
# ax2.set_title("Model D cost\n", size=36)

### for selectivity (flux ratio) plot
#D2 = np.genfromtxt("flux_29000_no_leak", skip_header=1)  # for no-leak comparison
#ax2.plot(D[:, 0], D[:, 2]/D[:, 3],  
#          linewidth=10, linestyle='-', color='Blue', alpha=0.75)
# ax2.plot(D2[:, 0], D2[:, 2]/D2[:, 3], label="Model w/o ion leak",  # for no-leak comparison
#          linewidth=5, linestyle='-', color='Green', alpha=0.75)
#ax2.axhline(y=2.71, xmin=0, xmax=1, linestyle='--',linewidth=20,
 #           label="exp(1) ~ 2.7", color='Red', alpha=0.75)
#ax2.set_ylabel(r'$J_{substrate}$/$J_{decoy}$', size=48)
#ax2.set_yscale('log')
#ax2.set_ylim(1e-2, 1e4)

### for cost (flux ratio) plot
# ax2.plot(D[:, 0], D[:, 1]/D[:, 2],
#          linewidth=10, linestyle='-', color='Blue', alpha=0.75)
# ax2.set_ylabel(r'$J_{ion}$/$J_{substrate}$', size=48)
# ax2.set_ylim(-10, 10)


### for flux plot
# ax2.plot(D[:, 0], D[:, 1]/D[0, 1], label=r'$J_{ion}$',  # ion flux
#          linewidth=40, linestyle='--', color='Blue', alpha=0.7)
# ax2.plot(D[:, 0], D[:, 2]/D[0, 1], label=r'$J_{substrate}$',  # substrate flux
#          linewidth=20, linestyle='-', color='Green')
# ax2.plot(D[:, 0], D[:, 3]/D[0, 1], label=r'$J_{decoy}$', linewidth=5, color='Red')  # decoy (proofreading only)
# ax2.set_ylabel(r'$J$ (flux) [arb. units]', size=36)


#fig.gca().set_xlabel(
#    r'$\Delta\mu_{ion}$ (chemical potential difference of driving ion) [$k_{B}T$]', size=36)
#plt.grid(b=True, which='both', axis='both')
#ax2.tick_params(axis='both', which='major', labelsize=28)
#ax2.tick_params(axis='both', which='minor', labelsize=24)
#ax2.yaxis.set_minor_locator(mtick.MultipleLocator(5))
#plt.legend(prop={'size': 28})  # 42

### for inset image...
# ax = plt.axes([0.15, 0.50, 0.45, 0.45], frameon=True)  # inset image location (L,B,dX,dY)
# ax.set_title("Kinetic pathway", size=32)
# image = plt.imread(
#     'C:\\Users\\georgeau\\Desktop\\preprint\\preprint_figures\\manuscript\\fig5b_inset.png')
# ax.tick_params(axis='both', left='off', top='off', right='off', bottom='off',
#                 labelleft='off', labeltop='off', labelright='off', labelbottom='off')
# ax.imshow(image)


# vertical line on main axis
# ax2.axvline(-4, linewidth=5, color='gray', linestyle='--')

### for multiplot
fig = plt.figure(figsize=(15, 15))

D1 = np.genfromtxt("flux_3500", skip_header=1)
D2 = np.genfromtxt("flux_29000", skip_header=1)
D3 = np.genfromtxt("flux_3000", skip_header=1)
D4 = np.genfromtxt("flux_829000", skip_header=1)
D_raw = [D1,D2,D3,D4]
model_list = ["A","B","C","D"]

for i in range(0, 4):
    plot_n = 221+i
    
    D = D_raw[i]
    ax2 = plt.subplot(plot_n)
    ax2.plot(D[:, 0], D[:, 1]/D[0, 1], label=r'$J_{ion}$',  # ion flux
             linewidth=20, linestyle='--', color='Blue', alpha=0.7)
    ax2.plot(D[:, 0], D[:, 2]/D[0, 1], label=r'$J_{substrate}$',  # substrate flux
             linewidth=10, linestyle='-', color='Green')
    ax2.plot(D[:, 0], D[:, 3]/D[0, 1], label=r'$J_{decoy}$', linewidth=5, color='Red')  # decoy (proofreading only)

    if i == 0:
        ax2.set_ylabel(r'$J$ (flux) [arb. units]', size=36)
    if i == 2:
        ax2.set_ylabel(r'$J$ (flux) [arb. units]', size=36)
        fig.gca().set_xlabel(
            r'$\Delta\mu_{ion}$ [$k_{B}T$]', size=36)  
    if i ==3:
        fig.gca().set_xlabel(
            r'$\Delta\mu_{ion}$ [$k_{B}T$]', size=36)  
    ax2.set_ylim(-0.2,1.0)
    #plt.suptitle("Ion, substrate, and decoy fluxes\n", size = 38)
    ax2.set_title("Model %s flux" % (model_list[i]), size=36)  # flux

    plt.grid(b=True, which='both', axis='both')
    ax2.tick_params(axis='both', which='major', labelsize=28)
    ax2.tick_params(axis='both', which='minor', labelsize=24)
    plt.legend(prop={'size': 28})  # 42


plt.savefig("test.png")
plt.show()
