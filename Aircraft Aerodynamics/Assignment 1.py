from Airfoilselector import NACAcalculator
import math as m
import numpy as np
import matplotlib.pyplot as plt
from Airfoilimporter import import_airfoil
from Airfoilgeometrycalculator import airfoil_geometery
from Cp_calculator import Cp_calculatorSP, Cp_calculatorVP, Cp_calculatorVPSP


# Flight conditions parameters
AoA = 5     #deg
V_fs = 1    #m/s
M = 0.3

# Initiate airfoil selector
Naca = NACAcalculator()

# Importing Airfoil Geometry
AF1 = import_airfoil(Naca[2])
AF2 = import_airfoil(Naca[3])
X1 = AF1[0]
Y1 = -AF1[1]
X2 = AF2[0]
Y2 = -AF2[1]


#Geometery parameters
Geom1 = airfoil_geometery(X1, Y1, AoA)
Geom2 = airfoil_geometery(X2, Y2, AoA)


# %% PLOTTING

# Plot the paneled geometry
fig = plt.figure(1)
plt.plot(X1, Y1, color='orange', marker='o', markersize = 3, label = Naca[2])
plt.plot(X2, Y2, color='blue', marker='o', markersize = 3, label = Naca[3])
plt.plot(0, 0, color='w', label = "AoA = "+str(AoA)+"\N{DEGREE SIGN}")
plt.plot(0, 0, color='w', label = "M = "+str(M))
plt.xlabel('X-Axis')
plt.ylabel('Y-Axis')
plt.title('Panel Geometry')
plt.axis('equal')
plt.legend()
#plt.show()

AoAr = AoA * np.pi/180
Results1 = Cp_calculatorVPSP(X1, Y1, Geom1[0], Geom1[1], Geom1[2], Geom1[3], Geom1[4], V_fs, AoAr)
Results2 = Cp_calculatorVPSP(X2, Y2, Geom2[0], Geom2[1], Geom2[2], Geom2[3], Geom2[4], V_fs, AoAr)

fig1 = plt.figure()

midIndS = int(np.floor(len(Results1)/2))                    # Airfoil middle index for VPM data
plt.plot(Geom1[0][midIndS + 1:len(Geom1[0])], Results1[midIndS + 1:len(Geom1[0])], markerfacecolor='b', label='VPM Upper')
plt.plot(Geom1[0][0:midIndS], Results1[0:midIndS], markerfacecolor='r', label='VPM Lower')
plt.gca().invert_yaxis()
plt.xlim([0,1])
plt.xlabel('X-Axis')
plt.ylabel('C_p')
plt.title('C_p distribution')
plt.legend()
plt.grid()
plt.show()

#plt.legend()
# print(Results1[1])
# plot2 = plt.figure(3)
# plt.cla()
# plt.contourf(Results1[2],Results1[3],Results1[1],500,cmap='jet')            # Plot contour
# plt.fill(X1,Y1,'k')                                                         # Plot airfoil as black polygon
# plt.xlabel('X Units')                                                       # Set X-label
# plt.ylabel('Y Units')                                                       # Set Y-label
# plt.gca().set_aspect('equal')                                               # Set axes equal
# plt.show()                                                                  # Display plot