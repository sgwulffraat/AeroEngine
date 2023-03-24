import numpy as np
import matplotlib.pyplot as plt
from XFOIL2 import xfoil
from XFOIL_ID import xfoil_id, xfoil_trs

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

cd_list = []
trans = np.arange(0.4, 0.75, 0.05)
for xtr in trans:
    cd, cp, x = xfoil_trs(2614, 2E05, 12, 0.7, xtr)
    cd_list.append(cd)

cd_init, cp_init, x_init = xfoil_trs(2614, 2E05, 12, 0.7, 1)
cd_new, cp_new, x_new = xfoil_trs(2614, 2E05, 12, 0.7, 0.65)



#Original xtr = 0.6969 and  CD =  0.01120
print()
print("Forced transition at: ", trans)
print("Respective drag: ", cd_list)
# %% PLOTTING %% #
#Plotting flag
flag = [
    1,      # Set to 1 to show Lift polar plot
    1,      # Set to 1 to show Drag polar plot
    1,      # Set to 1 to show Transition point plot
    1,      # Set to 1 to show Skin friction plot
    1,      # Set to 1 to show Optimized cp plot
    1,      # Set to 1 to show Cp for different Re plot
    1,      # Set to 1 to show Cp for fixed transition plot
]

# Plotting lift polar
if flag[0] == 1:
    fig1 = plt.figure()
    plt.plot(alpha, cl_polar, color='tab:blue', marker='o', markersize=2)
    plt.xlabel(chr(945))
    plt.ylabel('Cl')
    plt.grid()
    plt.show()

# Plotting drag polar
if flag[1] == 1:
    fig2 = plt.figure()
    plt.plot(cd_polar, cl_polar, color='tab:blue', marker='o', markersize=2)
    plt.xlabel('Cd')
    plt.ylabel('Cl')
    plt.grid()
    plt.show()

# Plotting location of transition point
if flag[2] == 1:
    fig3 = plt.figure()
    plt.plot(bot_xtr_pol, alpha, label='Bottom', color='tab:blue', marker='o', markersize=2)
    plt.plot(top_xtr_pol, alpha, label='Top', color='tab:orange', marker='o', markersize=2)
    plt.xlabel('Transition point (x/c)')
    plt.ylabel(chr(945))
    plt.legend()
    plt.grid()
    plt.show()

# Plotting skin friction
if flag[3] == 1:
    fig4 = plt.figure()
    plt.plot(x_a0, cf_a0, label=r'α = 0', color='tab:blue')
    plt.plot(x_a4, cf_a4, label=r'α = 4', color='tab:orange')
    plt.xlabel('Chordwise position (x/c)')
    plt.ylabel('Skin friction (Cf)')
    plt.legend()
    plt.grid()
    plt.show()

# Plotting cp distributions
if flag[4] == 1:
    fig5 = plt.figure()
    plt.plot(x_org, cp_org, label='Original', color='tab:blue')
    plt.plot(x_opt, cp_opt, label='Optimized', color='tab:orange')
    plt.gca().invert_yaxis()
    plt.xlabel('Chordwise position (x/c)')
    plt.ylabel('Pressure coefficient (Cp)')
    plt.legend()
    plt.grid()
    plt.show()

# Plotting cp distributions for different Re numbers
if flag[5] == 1:
    fig6 = plt.figure()
    plt.plot(x_org, cp_list[0], label=r'$Re = 2,0*10^5$', color='tab:blue')
    plt.plot(x_org, cp_list[1], label=r'$Re = 4,0*10^5$', color='tab:orange')
    plt.plot(x_org, cp_list[2], label=r'$Re = 6,0*10^5$', color='tab:green')
    plt.gca().invert_yaxis()
    plt.xlabel('Chordwise position (x/c)')
    plt.ylabel('Pressure coefficient (Cp)')
    plt.legend()
    plt.grid()
    plt.show()


# Plotting cp distributions for different Re numbers
if flag[6] == 1:
    fig7 = plt.figure()
    plt.plot(x_init, cp_init, label='Free transition', color='tab:blue')
    plt.plot(x_new, cp_new, label='Fixed transition at x/c = 0.65', color='tab:orange')
    plt.gca().invert_yaxis()
    plt.xlabel('Chordwise position (x/c)')
    plt.ylabel('Pressure coefficient (Cp)')
    plt.legend()
    plt.grid()
    plt.show()


