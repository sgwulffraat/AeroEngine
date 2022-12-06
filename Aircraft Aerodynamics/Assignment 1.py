from Airfoilselector import NACAcalculator
import numpy as np
import matplotlib.pyplot as plt
from Airfoilimporter import import_airfoil
from Airfoilgeometrycalculator import airfoil_geometery
from Cp_calculator import Cp_calculatorVPSP
from XFOIL import xfoil
# %% INPUTS ##

# Flight conditions parameters
AoA = 5     #deg
V_fs = 1    #m/s


# Initiate airfoil selector
Naca = NACAcalculator()
flag1 = Naca[4]
Naca2 = NACAcalculator()
flag2 = Naca2[4]


# %% COMPUTATIONS ##

# Transfer Angle of attack to radians
AoAr = AoA * np.pi/180

# Importing selected airfoil geometry from DAT file
AF1 = import_airfoil(Naca[2])
AF2 = import_airfoil(Naca[3])
AF3 = import_airfoil(Naca2[2])
AF4 = import_airfoil(Naca2[3])
X1 = AF1[0]     # x coordinates of the first airfoil
Y1 = AF1[1]    # y coordinates of the first airfoil
X2 = AF2[0]     # x coordinates of the second airfoil
Y2 = AF2[1]    # y coordinates of the second airfoil
X3 = AF3[0]     # x coordinates of the third airfoil
Y3 = AF3[1]    # y coordinates of the third airfoil
X4 = AF4[0]     # x coordinates of the fourth airfoil
Y4 = AF4[1]    # y coordinates of the fourth airfoil

# Using DAT file to compute geometery parameters
XC1, YC1, S1, phi1, beta1, N_pan1, Xnew_1, Ynew_1 = airfoil_geometery(X1, Y1, AoA)
XC2, YC2, S2, phi2, beta2, N_pan2, Xnew_2, Ynew_2 = airfoil_geometery(X2, Y2, AoA)
XC3, YC3, S3, phi3, beta3, N_pan3, Xnew_3, Ynew_3 = airfoil_geometery(X3, Y3, AoA)
XC4, YC4, S4, phi4, beta4, N_pan4, Xnew_4, Ynew_4 = airfoil_geometery(X4, Y4, AoA)


# Computing results using Vortex Panel and Source Panel methods drawn from [1]
Cp1, Cl1, Cm1 = Cp_calculatorVPSP(Xnew_1, Ynew_1, XC1, YC1, S1, phi1, beta1, V_fs, AoAr)
Cp2, Cl2, Cm2 = Cp_calculatorVPSP(Xnew_2, Ynew_2, XC2, YC2, S2, phi2, beta2, V_fs, AoAr)
Cp3, Cl3, Cm3 = Cp_calculatorVPSP(Xnew_3, Ynew_3, XC3, YC3, S3, phi3, beta3, V_fs, AoAr)
Cp4, Cl4, Cm4 = Cp_calculatorVPSP(Xnew_4, Ynew_4, XC4, YC4, S4, phi4, beta4, V_fs, AoAr)

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
cl_polar3 = []
cl_polar4 = []
cl_polar1_xfoil = []
cl_polar2_xfoil = []
for i in alpha:
    i_r = i * np.pi/180
    xc1, yc1, s1, phi_1, beta_1, N_pan_1, xnew1, ynew1 = airfoil_geometery(X1, Y1, i)
    xc2, yc2, s2, phi_2, beta_2, N_pan_2, xnew2, ynew2 = airfoil_geometery(X2, Y2, i)
    xc3, yc3, s3, phi_3, beta_3, N_pan_3, xnew3, ynew3 = airfoil_geometery(X3, Y3, i)
    xc4, yc4, s4, phi_4, beta_4, N_pan_4, xnew4, ynew4 = airfoil_geometery(X4, Y4, i)
    cp_u_1, cp_l_1, x_u_1, x_l_1, cl_1 = xfoil(NACAcode1, i, Numnodes1)
    cp_u_2, cp_l_2, x_u_2, x_l_2, cl_2 = xfoil(NACAcode2, i, Numnodes2)
    Cp_1, Cl_1, Cm_1 = Cp_calculatorVPSP(xnew1, ynew1, xc1, yc1, s1, phi_1, beta_1, V_fs, i_r)
    Cp_2, Cl_2, Cm_2 = Cp_calculatorVPSP(xnew2, ynew2, xc2, yc2, s2, phi_2, beta_2, V_fs, i_r)
    Cp_3, Cl_3, Cm_3 = Cp_calculatorVPSP(xnew3, ynew3, xc3, yc3, s3, phi_3, beta_3, V_fs, i_r)
    Cp_4, Cl_4, Cm_4 = Cp_calculatorVPSP(xnew4, ynew4, xc4, yc4, s4, phi_4, beta_4, V_fs, i_r)
    cl_polar1.append(Cl_1)
    cl_polar2.append(Cl_2)
    cl_polar3.append(Cl_3)
    cl_polar4.append(Cl_4)
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



### %% PLOTTING ###

### GEOMETERY ###
### Plotting the paneled geometry
if flag1 != "Reference":
    fig = plt.figure(1)
    plt.plot(X1, Y1, color='orange', marker='o', markersize=3, label=Naca[2])
    plt.plot(X2, Y2, color='blue', marker='o', markersize=3, label=Naca[3])
    plt.plot(0, 0, color='w', label = "AoA = "+str(AoA)+"\N{DEGREE SIGN}")
    plt.xlabel('X-Axis')
    plt.ylabel('Y-Axis')
    plt.title('Panel Geometry')
    plt.axis('equal')
    plt.legend()
    plt.show()

### REFERENCE ###

# Plotting the pressure distribution reference Airfoil
if flag1 == "Reference":
    for i in [0, 10, 15]:
        i_r = i * np.pi / 180
        XC1r, YC1r, S1r, phi1r, beta1r, N_pan1r, Xnew_1r, Ynew_1r = airfoil_geometery(X1, Y1, i)
        Cp1r, Cl1r, Cm1r = Cp_calculatorVPSP(Xnew_1r, Ynew_1r, XC1r, YC1r, S1r, phi1r, beta1r, V_fs, i_r)
        cp_u_1r, cp_l_1r, x_u_1r, x_l_1r, cl_1r = xfoil(NACAcode1, i, Numnodes1)

        if i == 0:
            ladsonref1 = ladson_cp_6free_0
            ladsonref2 = ladson_cp_3fixed_0
            gregoryref = gregory_cp_0
        elif i == 10:
            ladsonref1 = ladson_cp_6free_10
            ladsonref2 = ladson_cp_3fixed_10
            gregoryref = gregory_cp_10
        elif i == 15:
            ladsonref1 = ladson_cp_6free_15
            ladsonref2 = ladson_cp_3fixed_15
            gregoryref = gregory_cp_15
        else:
            ladsonref1 = []
            ladsonref2 = []
            gregoryref = []

        fig2 = plt.figure() # Airfoil middle point of VPM data
        plt.plot(x_u_1r, cp_u_1r, color='blue', marker='o', label= 'XFOIL', markersize=2)
        plt.plot(x_l_1r, cp_l_1r, color='blue', marker='o', markersize=2)
        plt.plot(XC1r, Cp1r,  color='orange', marker='.', label='Model')
        plt.scatter(ladsonref1[:, 0], ladsonref1[:, 1], color='red', marker='^', s=4, label='Ladson Re=6E6, free transition')
        plt.scatter(ladsonref2[:, 0], ladsonref2[:, 1], color='green', marker='D', s=4, label='Ladson Re=3E6, fixed transition')
        plt.scatter(gregoryref[:, 0], gregoryref[:, 1], color='purple', marker='s', s=4, label='Gregory Re=3E6, free transition')
        plt.plot([], [], ' ', label=chr(945) + " = " + str(i))
        plt.gca().invert_yaxis()
        plt.xlim([0, 1])
        plt.xlabel('x/c')
        plt.ylabel('Cp')
        plt.title("Reference Cp distribution "+str(Naca[2]))
        plt.legend(loc='lower right')
        plt.grid()
        plt.show()


### THICK VS THIN ###

if flag1 == "None":
    fig3 = plt.figure()
    plt.plot(XC1, Cp1, color='orange', marker='^', label=str(Naca[2]), markersize=2)
    plt.plot(XC3, Cp3, color='blue', marker='^', label=str(Naca2[2]), markersize=2)
    plt.plot([], [], ' ', label=chr(945) + " = " + str(AoA))
    plt.gca().invert_yaxis()
    plt.xlim([0, 1])
    plt.xlabel('x/c')
    plt.ylabel('C_p')
    plt.title("Cp distribution " + str(Naca[2]) + " & " + str(Naca2[2]))
    plt.legend()
    plt.grid()
    plt.show()

    fig4 = plt.figure()
    plt.plot(XC2, Cp2, color='orange', marker='^', label=str(Naca[3]), markersize=2)
    plt.plot(XC4, Cp4, color='blue', marker='^', label=str(Naca2[3]), markersize=2)
    plt.plot([], [], ' ', label=chr(945) + " = " + str(AoA))
    plt.gca().invert_yaxis()
    plt.xlim([0, 1])
    plt.xlabel('x/c')
    plt.ylabel('C_p')
    plt.title("Cp distribution " + str(Naca[3]) + " & " + str(Naca2[3]))
    plt.legend()
    plt.grid()
    plt.show()


### LIFT POLAR ###

# Plotting lift polar airfoil 1
fig5 = plt.figure()
if flag1 == "Reference":
    plt.plot(abbot_alpha, abbot_cl, label='Abbot', linestyle="", marker='o', markersize=3, color='red')
    plt.plot(gregory_alpha, gregory_cl, label='Gregory', linestyle="", marker='o', markersize=3, color='green')
    plt.plot(alpha, cl_polar1_xfoil, label='Xfoil', color='blue')
    plt.plot(alpha, cl_polar1, label='Model', color='orange')
if flag1 == "None":
    plt.plot(alpha, cl_polar1, label=str(Naca[2]), color='orange')
    plt.plot(alpha, cl_polar3, label=str(Naca2[2]), color='blue')
    plt.plot(alpha, cl_polar2, label=str(Naca[3]), color='green')
    plt.plot(alpha, cl_polar4, label=str(Naca2[3]), color='red')
if flag1 == "Camber":
    plt.plot(alpha, cl_polar1, label=str(Naca[2]), color='orange')
    plt.plot(alpha, cl_polar2, label=str(Naca[3]), color='blue')
plt.xlabel(chr(945))
plt.ylabel('Cl')
plt.title("Lift polars")
plt.legend()
plt.grid()
plt.show()

### CAMBERED ###
# Plot to compare CP cambered vs uncambered
if flag1 == "Camber":
    fig6 = plt.figure()
    plt.plot(XC1, Cp1, color='orange', marker='^', label=str(Naca[2]), markersize=2)
    plt.plot(XC2, Cp2, color='blue', marker='^', label=str(Naca[3]), markersize=2)
    plt.plot([], [], ' ', label=chr(945) + " = " + str(AoA))
    plt.gca().invert_yaxis()
    plt.xlim([0, 1])
    plt.xlabel('x/c')
    plt.ylabel('C_p')
    plt.title("Cp distribution " + str(Naca[2]) + " & " + str(Naca[3]))
    plt.legend()
    plt.grid()
    plt.show()

### Panels ###
# if flag1 == "Panels 1":

