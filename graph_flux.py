import numpy as np 
from matplotlib import pyplot as plt 
import matplotlib.ticker as mtick
import argparse


def analysis_plots(datafile, y1 =[-0.5,1.25], y2=[-30,50], y3=[-15,15]):
    D = np.genfromtxt(datafile, skip_header=1)  # col 0,1,2,3 = ddmu_now, N,S,W flux, row 0, 1-n: labels, ddmu_now 

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1,figsize=(15, 10), tight_layout=True)
    #fig.suptitle('Model flux analysis', size=24)

    # absolute flux
    ax1.set_ylim(y1[0],y1[1])
    ax1.set_title("Flux", size=18)
    ax1.axvline(-4, linewidth=2, color='gray', linestyle='--')
    ax1.set_ylabel(r'$J$ (scaled by max)', size=18)
    ax1.plot(D[:, 0], D[:, 1]/D[0, 1], label=r'$J_{ion}$', linewidth=2, linestyle='--', color='Blue', alpha=0.7)  # ion flux
    ax1.plot(D[:, 0], D[:, 2]/D[0, 1], label=r'$J_{substrate}$', linewidth=2, linestyle='-', color='Green')  # substrate flux
    ax1.plot(D[:, 0], D[:, 3]/D[0, 1], label=r'$J_{decoy}$', linewidth=2, color='Red')  # decoy (proofreading only)
    ax1.grid(b=True, which='both', axis='both')
    ax1.legend()

    # selectivity (s/w flux)
    ax2.set_ylim(y2[0],y2[1])
    ax2.set_title("Selectivity",size=18)
    ax2.axvline(-4, linewidth=2, color='gray', linestyle='--')
    ax2.set_ylabel(r'$J_{substrate}$/$J_{decoy}$',size=22)
    ax2.plot(D[:, 0], D[:, 2]/D[:, 3], linewidth=2, linestyle='-', color='Blue', alpha=0.75, label="model")
    ax2.axhline(y=10*2.71, xmin=0, xmax=1, linestyle='--',linewidth=2,
           label="(ref) filter thresh. ~ 27", color='Red', alpha=0.75)
    ax2.grid(b=True, which='both', axis='both')
    ax2.legend()

    # cost (s/n flux)
    ax3.set_ylim(y3[0],y3[1])
    ax3.set_title("Cost",size=18)
    ax3.axvline(-4, linewidth=2, color='gray', linestyle='--')
    ax3.set_xlabel(r'$\Delta\mu_{ion}$ (chemical potential difference of driving ion) [$k_{B}T$]', size=18)
    ax3.set_ylabel(r'$J_{ion}$/$J_{substrate}$', size=22)
    ax3.plot(D[:, 0], D[:, 1]/D[:, 2], linewidth=2, linestyle='-', color='Blue', alpha=0.75, label = "model")
    ax3.axhline(y=10, xmin=0, xmax=1, linestyle='--',linewidth=2,
           label="(ref) filter thresh. = 10", color='Red', alpha=0.75)
    ax3.grid(b=True, which='both', axis='both')
    ax3.legend()


### main

# parse input from command line
parser = argparse.ArgumentParser(description='Graph the flux, selectivity, and cost')
parser.add_argument("--f_name", required=True, type=str, help="flux datafile name. Example: 'name.dat' ")  # filename required
parser.add_argument("--s_name", required=True, type=str, help="save datafile name (w/ .png extension) Example: 'graph.png'.")  # save filename required 
parser.add_argument("--y_lim1", default=[-1,1], type=list, help="relative flux graph y limits [y_min,y_max]. Default: [-1,1].")
parser.add_argument("--y_lim2", default=[-100,100], type=list, help="selectivity graph y limits [y_min,y_max]. Default: [-100,100]")
parser.add_argument("--y_lim3", default=[-10,10], type=list, help="cost graph y limits [y_min,y_max]. Default: [-10,10]")
args = parser.parse_args()
f_name = args.f_name
s_name = args.s_name
y_lim1 = args.y_lim1
y_lim2 = args.y_lim2
y_lim3 = args.y_lim3

analysis_plots(f_name, y_lim1, y_lim2, y_lim3)
plt.savefig(s_name, bbox_inches='tight')
plt.show()


# try:
#     analysis_plots(f_name, y_lim1, y_lim2, y_lim3)
#     plt.show()
#     plt.save(s_name)
# except:
#     print("Error: could not graph the flux")

