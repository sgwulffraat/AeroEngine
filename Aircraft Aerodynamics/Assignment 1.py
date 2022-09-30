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

# Computing results using Vortex Panel and Source Panel methods drawn from [1]
Results1 = Cp_calculatorVPSP(X1, Y1, XC1, YC1, S1, phi1, beta1, V_fs, AoAr)
Results2 = Cp_calculatorVPSP(X2, Y2, XC2, YC2, S2, phi2, beta2, V_fs, AoAr)

# %% XFOIL INTERGRATION ###

# Getting naca code in the form of e.g. '0012'
NACAcode1 = Naca[2][4:]
NACAcode2 = Naca[3][4:]
Numnodes1 = len(X1)
Numnodes2 = len(X2)

# Computing CP array from xfoil
xfoil_cp_u_1, xfoil_cp_l_1 = xfoil(NACAcode1, AoA, Numnodes1, Y1)
xfoil_cp_u_2, xfoil_cp_l_2 = xfoil(NACAcode2, AoA, Numnodes2, Y2)


# %% PLOTTING ###

# Plotting the paneled geometry
fig = plt.figure(1)
plt.plot(X1, Y1, color='orange', marker='o', markersize = 3, label = Naca[2])
plt.plot(X2, Y2, color='blue', marker='o', markersize = 3, label = Naca[3])
plt.plot(0, 0, color='w', label = "AoA = "+str(AoA)+"\N{DEGREE SIGN}")
plt.xlabel('X-Axis')
plt.ylabel('Y-Axis')
plt.title('Panel Geometry')
plt.axis('equal')
plt.legend()
#plt.show()

# Plotting the pressure distribution
fig1 = plt.figure() # Airfoil middle point of VPM data
midpoint = int(np.floor(len(Results1)/2)) # Separating top and bottom side of airfoil
plt.plot(XC1[midpoint + 1:len(XC1)], Results1[midpoint + 1:len(XC1)], markerfacecolor='b', label='Upper')
plt.plot(XC1[0:midpoint], Results1[0:midpoint], markerfacecolor='r', label='Lower')
plt.plot(X1_U, xfoil_cp_u_1, label= 'XFOIL Upper')
plt.plot(X1_L, xfoil_cp_l_1, label='XFOIL Lower')
plt.gca().invert_yaxis()
plt.xlim([0,1])
plt.xlabel('X-Axis')
plt.ylabel('C_p')
plt.title('C_p distribution')
plt.legend()
plt.grid()
plt.show()