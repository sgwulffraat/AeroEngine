from Airfoilselector import NACAcalculator
from math import *
import numpy as np
import matplotlib.pyplot as plt
from Airfoilimporter import import_airfoil

# Airfoil thicknesses
Naca = NACAcalculator()
N1 = Naca[0]
N2 = Naca[1]


#Importing Airfoil Geometry
AF1 = import_airfoil(Naca[2])
AF2 = import_airfoil(Naca[3])
X1 = AF1[0]
Y1 = AF1[1]
X2 = AF2[0]
Y2 = AF2[1]

x1,y1 = AF1
x2,y2 = AF2
# plt.figure()
# plt.plot(X1,Y1)
# plt.show()
plt.figure()
plt.plot(x1,y1)
plt.show()