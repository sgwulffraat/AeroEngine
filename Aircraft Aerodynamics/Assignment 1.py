from Airfoilselector import NACAcalculator
import numpy as np
import matplotlib.pyplot as plt
from Airfoilimporter import import_airfoil
from Airfoilgeometrycalculator import airfoil_geometery
from Cp_calculator import Cp_calculatorVPSP
from XFOIL import xfoil
# %% INPUTS ##

# Flight conditions parameters
AoA = 15     #deg
V_fs = 1    #m/s


# Initiate airfoil selector
Naca = NACAcalculator()


# %% COMPUTATIONS ##

# Transfer Angle of attack to radians
AoAr = AoA * np.pi/180

# Importing selected airfoil geometry from DAT file
AF1 = import_airfoil(Naca[2])
AF2 = import_airfoil(Naca[3])
X1 = AF1[0]     # x coordinates of the first airfoil
Y1 = -AF1[1]    # y coordinates of the first airfoil
X2 = AF2[0]     # x coordinates of the second airfoil
Y2 = -AF2[1]    # y coordinates of the first airfoil

# Split airfoil into (U)pper and (L)ower
X1_U = X1[Y1 >= 0]
X1_L = X1[Y1 < 0]
Y1_U = Y1[Y1 >= 0]
Y1_L = Y1[Y1 < 0]

X2_U = X2[Y2 >= 0]
X2_L = X2[Y2 < 0]
Y2_U = Y2[Y2 >= 0]
Y2_L = Y2[Y2 < 0]

# Using DAT file to compute geometery parameters
XC1, YC1, S1, phi1, beta1, N_pan1 = airfoil_geometery(X1, Y1, AoA)
XC2, YC2, S2, phi2, beta2, N_pan2 = airfoil_geometery(X2, Y2, AoA)
XC1U = XC1[YC1 >= 0]
XC1L = XC1[YC1 < 0]

# Computing results using Vortex Panel and Source Panel methods drawn from [1]
Cp1, Cl1, Cm1 = Cp_calculatorVPSP(X1, Y1, XC1, YC1, S1, phi1, beta1, V_fs, AoAr)
Cp2, Cl2, Cm2 = Cp_calculatorVPSP(X2, Y2, XC2, YC2, S2, phi2, beta2, V_fs, AoAr)

# %% XFOIL INTERGRATION ###

# Getting naca code in the form of e.g. '0012'
NACAcode1 = Naca[2][4:]
NACAcode2 = Naca[3][4:]
Numnodes1 = len(X1)
Numnodes2 = len(X2)

# Computing CP array from xfoil
xfoil_cp_u_1, xfoil_cp_l_1, xfoil_x_u_1, xfoil_x_l_1, xfoil_cl_1 = xfoil(NACAcode1, AoA, Numnodes1)
xfoil_cp_u_2, xfoil_cp_l_2, xfoil_x_u_2, xfoil_x_l_2, xfoil_cl_2 = xfoil(NACAcode2, AoA, Numnodes2)

# %% LIFT POLAR %% ###
alpha = [-10, 0, 10]
cl_polar1 = []
cl_polar2 = []
cl_polar1_xfoil = []
cl_polar2_xfoil = []
for i in alpha:
    i_r = i * np.pi/180
    xc1, yc1, s1, phi_1, beta_1, N_pan_1 = airfoil_geometery(X1, Y1, i)
    xc2, yc2, s2, phi_2, beta_2, N_pan_2 = airfoil_geometery(X2, Y2, i)
    cp_u_1, cp_l_1, x_u_1, x_l_1, cl_1 = xfoil(NACAcode1, i, Numnodes1)
    cp_u_2, cp_l_2, x_u_2, x_l_2, cl_2 = xfoil(NACAcode2, i, Numnodes2)
    Cp_1, Cl_1, Cm_1 = Cp_calculatorVPSP(X1, Y1, XC1, YC1, S1, phi_1, beta_1, V_fs, i_r)
    Cp_2, Cl_2, Cm_2 = Cp_calculatorVPSP(X2, Y2, XC2, YC2, S2, phi_2, beta_2, V_fs, i_r)
    cl_polar1.append(Cl_1)
    cl_polar2.append(Cl_2)
    cl_polar1_xfoil.append(cl_1)
    cl_polar2_xfoil.append(cl_2)

# Reference Data
abbot = np.loadtxt("Lift_polar_abbot.dat", skiprows=6)
abbot_alpha = abbot[:, 0]
abbot_cl = abbot[:, 1]
gregory = np.loadtxt("Lift_polar_gregory.dat", skiprows=3)
gregory_alpha = gregory[:, 0]
gregory_cl = gregory[:, 1]
ladson_cp_6free_0 = np.loadtxt("CP_Ladson_6free_alpha_0.dat", skiprows=5)
ladson_cp_6free_10 = np.loadtxt("CP_Ladson_6free_alpha_10.dat", skiprows=1)
ladson_cp_6free_15 = np.loadtxt("CP_Ladson_6free_alpha_15.dat", skiprows=1)
ladson_cp_3fixed_0 = np.loadtxt("CP_Ladson_3fixed_alpha_0.dat", skiprows=5)
ladson_cp_3fixed_10 = np.loadtxt("CP_Ladson_3fixed_alpha_10.dat", skiprows=1)
ladson_cp_3fixed_15 = np.loadtxt("CP_Ladson_3fixed_alpha_15.dat", skiprows=1)
gregory_cp_0 =np.loadtxt("CP_Gregory_alpha_0.dat", skiprows=6)
gregory_cp_10 =np.loadtxt("CP_Gregory_alpha_10.dat", skiprows=1)
gregory_cp_15 =np.loadtxt("CP_Gregory_alpha_15.dat", skiprows=1)



# %% PLOTTING ###

# Plotting the paneled geometry
# fig = plt.figure(1)
# plt.plot(X1, Y1, color='orange', marker='o', markersize=3, label=Naca[2])
# plt.plot(X2, Y2, color='blue', marker='o', markersize=3, label=Naca[3])
# plt.plot(0, 0, color='w', label = "AoA = "+str(AoA)+"\N{DEGREE SIGN}")
# plt.xlabel('X-Axis')
# plt.ylabel('Y-Axis')
# plt.title('Panel Geometry')
# plt.axis('equal')
# plt.legend()
# plt.show()

# Plotting the pressure distribution Airfoil 1
fig1 = plt.figure() # Airfoil middle point of VPM data
midpoint = int(np.floor(len(Cp1)/2)) # Separating top and bottom side of airfoil
plt.plot(XC1[midpoint + 1:len(XC1)], Cp1[midpoint + 1:len(XC1)],  color = 'green', marker='^', label='Model Upper', markersize=2)
plt.plot(XC1[0:midpoint], Cp1[0:midpoint], color='orange', marker='^', label='Model Lower', markersize=2)
plt.plot(xfoil_x_u_1, xfoil_cp_u_1, linestyle=':', color = 'red', marker='o', label= 'XFOIL Upper', markersize=2)
plt.plot(xfoil_x_l_1, xfoil_cp_l_1, linestyle=':', color = 'blue', marker='o', label='XFOIL Lower', markersize=2)
plt.plot([], [], ' ', label=chr(945) + " = " + str(AoA))
plt.gca().invert_yaxis()
plt.xlim([0, 1])
plt.xlabel('x/c')
plt.ylabel('C_p')
plt.title("Cp distribution "+str(Naca[2]))
plt.legend()
plt.grid()
plt.show()

# Plotting the pressure distribution reference Airfoil
if AoA == 0 and NACAcode1 == '0012':
    ladsonref1 = ladson_cp_6free_0
    ladsonref2 = ladson_cp_3fixed_0
    gregoryref = gregory_cp_0
elif AoA == 10 and NACAcode1 == '0012':
    ladsonref1 = ladson_cp_6free_10
    ladsonref2 = ladson_cp_3fixed_10
    gregoryref = gregory_cp_10
elif AoA == 15 and NACAcode1 == '0012':
    ladsonref1 = ladson_cp_6free_15
    ladsonref2 = ladson_cp_3fixed_15
    gregoryref = gregory_cp_15
else:
    ladsonref1 = []
    ladsonref2 = []
    gregoryref = []

if ladsonref1 != []:
    fig2 = plt.figure() # Airfoil middle point of VPM data
    midpoint = int(np.floor(len(Cp1)/2)) # Separating top and bottom side of airfoil
    plt.plot(xfoil_x_u_1, xfoil_cp_u_1, color='blue', marker='o', label= 'XFOIL', markersize=2)
    plt.plot(xfoil_x_l_1, xfoil_cp_l_1, color='blue', marker='o', markersize=2)
    plt.plot(XC1, Cp1,  color='orange', marker='.', label='Model')
    plt.scatter(ladsonref1[:, 0], ladsonref1[:, 1], color='red', marker='^', s=4, label='Ladson Re=6E6, free transition')
    plt.scatter(ladsonref2[:, 0], ladsonref2[:, 1], color='green', marker='D', s=4, label='Ladson Re=3E6, fixed transition')
    plt.scatter(gregoryref[:, 0], gregoryref[:, 1], color='purple', marker='s', s=4, label='Gregory Re=3E6, free transition')
    plt.plot([], [], ' ', label=chr(945) + " = " + str(AoA))
    plt.gca().invert_yaxis()
    plt.xlim([0, 1])
    plt.xlabel('x/c')
    plt.ylabel('Cp')
    plt.title("Reference Cp distribution "+str(Naca[2]))
    plt.legend(loc='lower right')
    plt.grid()
    plt.show()


# Plotting the pressure distribution Airfoil 2
fig3 = plt.figure() # Airfoil middle point of VPM data
midpoint = int(np.floor(len(Cp2)/2)) # Separating top and bottom side of airfoil
plt.plot(XC2[midpoint + 1:len(XC2)], Cp2[midpoint + 1:len(XC2)], color = 'green', marker='^', label='Model Upper', markersize=2)
plt.plot(XC2[0:midpoint], Cp2[0:midpoint], color='orange', marker='^', label='Model Lower', markersize=2)
plt.plot(xfoil_x_u_2, xfoil_cp_u_2, linestyle=':', color='red', marker='o', label='XFOIL Upper', markersize=2)
plt.plot(xfoil_x_l_2, xfoil_cp_l_2, linestyle=':', color='blue', marker='o', label='XFOIL Lower', markersize=2)
plt.plot([], [], ' ', label=chr(945) + " = " + str(AoA))
plt.gca().invert_yaxis()
plt.xlim([0, 1])
plt.xlabel('x/c')
plt.ylabel('Cp')
plt.title("Cp distribution " + str(Naca[3]))
plt.legend()
plt.grid()
plt.show()

# Plotting lift polar airfoil 1
fig4 = plt.figure()
if NACAcode1 == '0012':
    plt.plot(abbot_alpha, abbot_cl, label='Abbot', linestyle="", marker='o', markersize=3, color='red')
    plt.plot(gregory_alpha, gregory_cl, label='Gregory', linestyle="", marker='o', markersize=3, color='green')
plt.plot(alpha, cl_polar1_xfoil, label='Xfoil', color='blue')
plt.plot(alpha, cl_polar1, label='Model', color='orange')
plt.xlabel(chr(945))
plt.ylabel('Cl')
plt.title("Lift polar " + str(Naca[2]))
plt.legend()
plt.grid()
plt.show()



# Plotting lift polar airfoil 2
fig5 = plt.figure()
plt.plot(alpha, cl_polar2_xfoil, label='Xfoil', color='blue')
plt.plot(alpha, cl_polar2, label='Model', color='orange')
plt.xlabel(chr(945))
plt.ylabel('Cl')
plt.title("Lift polar " + str(Naca[3]))
plt.legend()
plt.grid()
plt.show()

