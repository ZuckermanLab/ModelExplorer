import numpy as np 
from matplotlib import pyplot as plt 
import matplotlib.ticker as mtick
import argparse


### main
datafile = 'flows_3090_00_51'
n_ATP = 1
n_H = 3
dmu_ATP = -10
dmu_H = -0.1

D = np.genfromtxt(datafile, skip_header=1)

### thermodynamics
# ideal equilibrium: ideal H+ stoich. coeff. * dmu_H+ =  ideal ATP stoich. coeff. * dmu_ATP

# equilibrum when net flow = 0 <=> flow_H+ = flow_ATP
# ideal case when difference in concentrations are equal

# find when H flow ~ ATP flow
# compare each pair of flows for a given dmu, return pair with min distance in values. 
dist_list = [np.inf]
eq_idx = 0
for i in range(np.shape(D)[0]):
       x1 = D[i][1]
       x2 = D[i][2]
       dist = np.abs(x2-x1)
       if dist < dist_list[-1]:
              dist_list.append(dist)
              eq_idx = i

dmu_H_eq_idl = (n_ATP*dmu_ATP/n_H)
dmu_H_eq_act = D[int(eq_idx)][0]


print("ideal equilibrium conditions (vary dmu_H): dmu_H = %.1f n_H = %s dmu_ATP = %s n_ATP = %s" %(dmu_H_eq_idl,n_H, dmu_ATP, n_ATP))
print("actual equilibrium dmu_H = %.1f" % dmu_H_eq_act)

# unit test (thermodynamics at equilibrium)
# actual equilibrium point cannot be less than the ideal equilibrium point (2nd law)  
if np.abs(dmu_H_eq_act) < np.abs(dmu_H_eq_idl):
       print("Failed thermodynamic reality check")
       exit()

flux_labels: ["ATP->ADP+Pi", "ADP+Pi->ATP (against gradient)", "H+ -> out (against gradient)", "H+ -> in"]

fig,  ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2,figsize=(10, 10), tight_layout=True)

# absolute flux (top left)
ax1.set_title("Flux", size=18)
ax1.axvline(dmu_H, linewidth=2, color='gray', linestyle='--', label='sim. conditions')
ax1.axvline(dmu_H_eq_idl, linewidth=2, color='black', linestyle='--', label='ideal eq. point')
ax1.set_ylabel(r'Flux (scaled by max)', size=14)
ax1.plot(D[:, 0], D[:, 1]/np.max([np.max(np.abs(D[:,1])),np.max(np.abs(D[:,2]))]), label=r'ATP', linewidth=2, linestyle='--', color='Blue', alpha=0.7)  # ion flux
ax1.plot(D[:, 0], D[:, 2]/np.max([np.max(np.abs(D[:,1])),np.max(np.abs(D[:,2]))]), label=r'H+', linewidth=2, linestyle='-', color='Green')  # substrate flux
ax1.grid(b=True, which='both', axis='both')
ax1.axvspan(dmu_H_eq_act, 0, color='red', alpha=0.3, label='region of interest')
ax1.set_ylim(-1,1)
ax1.set_xlim(-10,0)
ax1.set_xlabel(r'$\Delta\mu_{H^+}$ (chemical potential difference of $H^+$) [$k_{B}T$]', size=14)
ax1.legend()

# coupling ratio (bottom left)
ax2.set_title("Coupling ratio (H+/ATP)",size=18)
#ax2.set_ylim(-1*(np.abs(n_H)+1),1*(np.abs(n_H)+1))
ax2.set_ylim(-3.5,0)
ax2.set_xlim(-10,0)
ax2.axvline(dmu_H, linewidth=2, color='gray', linestyle='--', label='sim. conditions')
ax2.axvline(dmu_H_eq_idl, linewidth=2, color='black', linestyle='--', label='ideal eq. point')
ax2.set_xlabel(r'$\Delta\mu_{H^+}$ (chemical potential difference of $H^+$) [$k_{B}T$]', size=14)
ax2.set_ylabel('Coupling ratio (H+/ATP) flux', size=14)
ax2.plot(D[:, 0], D[:, 2]/D[:, 1], linewidth=2, linestyle='-', color='Blue', alpha=0.75, label='H+:ATP flux')
ax2.grid(b=True, which='both', axis='both')
ax2.axvspan(dmu_H_eq_act, 0, color='red', alpha=0.3, label='region of interest')
ax2.axhline(y=-1*n_H/n_ATP, color = 'red', linestyle= '--', label='ideal coupling ratio')
ax2.legend()

# zoomed absolute flux (top right)
ax3.set_title("Flux [region of interest]", size=18)
ax3.axvline(dmu_H, linewidth=2, color='gray', linestyle='--', label='sim. conditions')
ax3.axvline(dmu_H_eq_idl, linewidth=2, color='black', linestyle='--', label='ideal eq. point')
ax3.plot(D[:, 0], D[:, 1]/np.max([np.max(np.abs(D[:,1])),np.max(np.abs(D[:,2]))]), label=r'ATP', linewidth=2, linestyle='--', color='Blue', alpha=0.7)  # ion flux
ax3.plot(D[:, 0], D[:, 2]/np.max([np.max(np.abs(D[:,1])),np.max(np.abs(D[:,2]))]), label=r'H+', linewidth=2, linestyle='-', color='Green')  # substrate flux
ax3.grid(b=True, which='both', axis='both')
ax3.set_ylim(-1,1)
ax3.set_xlim(dmu_H_eq_act,0)
ax3.set_xlabel(r'$\Delta\mu_{H^+}$ (chemical potential difference of $H^+$) [$k_{B}T$]', size=14)
ax1.set_ylabel(r'Flux (scaled by max)', size=14)
ax3.legend()

# zoomed coupling ratio (bottom right)
ax4.set_title("Coupling ratio (H+/ATP) [region of interest]",size=18)
#ax4.set_ylim(-1*(np.abs(n_H)+1),0)
ax4.set_ylim(-3.5,0)
ax4.set_xlim(dmu_H_eq_act,0)
ax4.axvline(dmu_H, linewidth=2, color='gray', linestyle='--', label='sim. conditions')
ax4.axvline(dmu_H_eq_idl, linewidth=2, color='black', linestyle='--', label='ideal eq. point')
ax4.set_xlabel(r'$\Delta\mu_{H^+}$ (chemical potential difference of $H^+$) [$k_{B}T$]', size=14)
ax4.set_ylabel('Coupling ratio (H+/ATP) flux', size=14)
ax4.plot(D[:, 0], D[:, 2]/D[:, 1], linewidth=2, linestyle='-', color='Blue', alpha=0.75, label='H+:ATP flux')
ax4.grid(b=True, which='both', axis='both')
ax4.axhline(y=-1*n_H/n_ATP, color = 'red', linestyle= '--', label='ideal coupling ratio')
ax4.legend()



#plt.figtext(.02, .02, "dmu < 0: [ATP]/[ADP] > [ATP_eq]/[ADP_eq], [H_in] > [H_out]")
plt.show()




## filters

#print(D[70][2]/D[70][1])  # coupling ratio at dmu_H_now = -1 (simulation settings)

