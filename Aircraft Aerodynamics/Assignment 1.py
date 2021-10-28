from Airfoilselector import NACAcalculator
import math as m
import numpy as np
import matplotlib.pyplot as plt
from Airfoilimporter import import_airfoil
from Airfoilgeometrycalculator import airfoil_geometery
from Cp_calculator import Cp_calculatorSP, Cp_calculatorVP

# Flight conditions parameters
AoA = 5
V_fs = 50 #m/s
M = 0.2

# Initiate airfoil selector
Naca = NACAcalculator()

# Importing Airfoil Geometry
AF1 = import_airfoil(Naca[2])
AF2 = import_airfoil(Naca[3])
X1 = AF1[0]
Y1 = AF1[1]
X2 = AF2[0]
Y2 = AF2[1]

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
plt.show()

AoAR = AoA * np.pi/180
Results1 = Cp_calculatorVP(X1, Y1, Geom1[0], Geom1[1], Geom1[2], Geom1[3], Geom1[4], Geom1[5], V_fs, AoAR)
Results2 = Cp_calculatorVP(X2, Y2, Geom2[0], Geom2[1], Geom2[2], Geom2[3], Geom2[4], Geom2[5], V_fs, AoAR)

print(Results2)

fig = plt.figure(2)
plt.plot(Geom1[0][:34], Results1[:34], color = 'red')
plt.plot(Geom1[0][34:], Results1[34:], color = 'blue')
plt.gca().invert_yaxis()
plt.xlabel('X-Axis')
plt.ylabel('C_p')
plt.title('C_p distribution')
plt.axis('equal')
#plt.legend()
plt.show()