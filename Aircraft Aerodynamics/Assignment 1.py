from Airfoilselector import NACAcalculator
import math as m
import numpy as np
import matplotlib.pyplot as plt
from Airfoilimporter import import_airfoil
from Airfoilgeometrycalculator import airfoil_geometery

# Flight conditions parameters
AoA = 0
M = 0.8

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