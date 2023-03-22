import numpy as np
import matplotlib.pyplot as plt
from XFOIL2 import xfoil
from XFOIL_ID import xfoil_id

# %% INPUTS - PART I #
NACA = "2614"
re = 1.2e06
visc = "on"
n_crit = 9

# %% EXTRACTING DATA FROM X-FOIL %% #
alpha = np.arange(-2, 9, 1)
cl_polar = []
cd_polar = []
bot_xtr_pol = []
top_xtr_pol = []
for i in alpha:
    cl, cd, bot_xtr, top_xtr, cf_list, x_list = xfoil(NACA, i, re, visc)
    cl_polar.append(cl)
    cd_polar.append(cd)
    bot_xtr_pol.append(bot_xtr)
    top_xtr_pol.append(top_xtr)
    if i == 0:
        cf_a0 = cf_list
        x_a0 = x_list
    if i == 4:
        cf_a4 = cf_list
        x_a4 = x_list

# Part II
cl_ii = 0.4
cp_org, x_org = xfoil_id(2614, re, n_crit, cl_ii)
cp_opt, x_opt = xfoil_id("Optimized_Airfoil.txt", re, n_crit, cl_ii)

# Part III
cp_list = []
for R in [2E05, 4E05, 6E05]:
    cp_re, x_nouse = xfoil_id(2614, R, 12, 0.7)
    cp_list.append(cp_re)

# %% PLOTTING %% #

# Plotting lift polar
fig = plt.figure()
plt.plot(alpha, cl_polar, color='tab:blue', marker='o', markersize=2)
plt.xlabel(chr(945))
plt.ylabel('Cl')
plt.title("Lift Polar")
plt.grid()
plt.show()

# Plotting drag polar
fig2 = plt.figure()
plt.plot(cd_polar, cl_polar, color='tab:blue', marker='o', markersize=2)
plt.xlabel('Cd')
plt.ylabel('Cl')
plt.grid()
plt.show()

# Plotting location of transition point
fig3 = plt.figure()
plt.plot(bot_xtr_pol, alpha, label='Bottom', color='tab:blue', marker='o', markersize=2)
plt.plot(top_xtr_pol, alpha, label='Top', color='tab:orange', marker='o', markersize=2)
plt.xlabel('Transition point (x/c)')
plt.ylabel(chr(945))
plt.legend()
plt.grid()
plt.show()

# Plotting skin friction
fig3 = plt.figure()
plt.plot(x_a0, cf_a0, label='a = 0', color='tab:blue')
plt.plot(x_a4, cf_a4, label='a = 4', color='tab:orange')
plt.xlabel('Chordwise position (x/c)')
plt.ylabel('Skin friction (Cf)')
plt.legend()
plt.grid()
plt.show()

# Plotting cp distributions
fig4 = plt.figure()
plt.plot(x_org, cp_org, label='Original', color='tab:blue')
plt.plot(x_opt, cp_opt, label='Optimized', color='tab:orange')
plt.gca().invert_yaxis()
plt.xlabel('Chordwise position (x/c)')
plt.ylabel('Pressure coefficient (Cp)')
plt.legend()
plt.grid()
plt.show()

# Plotting cp distributions for different Re numbers
fig5 = plt.figure()
plt.plot(x_org, cp_list[0], label=r'$Re = 2,0*10^5$', color='tab:blue')
plt.plot(x_org, cp_list[1], label=r'$Re = 4,0*10^5$', color='tab:orange')
plt.plot(x_org, cp_list[2], label=r'$Re = 6,0*10^5$', color='tab:green')
plt.gca().invert_yaxis()
plt.xlabel('Chordwise position (x/c)')
plt.ylabel('Pressure coefficient (Cp)')
plt.legend()
plt.grid()
plt.show()



