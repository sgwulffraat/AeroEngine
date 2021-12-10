
import numpy as np
import math as math
import matplotlib.pyplot as plt
from Airfoilselector import NACAcalculator
from Airfoilimporter import import_airfoil
from COMPUTE_IJ_SPM import COMPUTE_IJ_SPM
from COMPUTE_KL_VPM import COMPUTE_KL_VPM


# User-defined knowns
Vinf = 1  # Freestream velocity [] (just leave this at 1)
AoA = 0  # Angle of attack [deg]
NACA = '0012'  # NACA airfoil to load [####]

# Convert angle of attack to radians
AoAR = AoA * (np.pi / 180)  # Angle of attack [rad]

# Plotting flags
flagPlot = [0,  # Airfoil with panel normal vectors
            0,  # Geometry boundary pts, control pts, first panel, second panel
            0,  # Cp vectors at airfoil surface panels
            1]  # Pressure coefficient comparison (XFOIL vs. VPM)

# Initiate airfoil selector
Naca = NACAcalculator()

# Importing Airfoil Geometry
AF1 = import_airfoil(Naca[2])
XB = AF1[0]
YB = AF1[1]

# Number of boundary points and panels
numPts = len(XB)  # Number of boundary points
numPan = numPts - 1  # Number of panels (control points)


# %% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

# Check for direction of points
edge = np.zeros(numPan)  # Initialize edge value array
for i in range(numPan):  # Loop over all panels
    edge[i] = (XB[i + 1] - XB[i]) * (YB[i + 1] + YB[i])  # Compute edge values

sumEdge = np.sum(edge)  # Sum all edge values

# If panels are CCW, flip them (don't if CW)
if (sumEdge < 0):  # If panels are CCW
    XB = np.flipud(XB)  # Flip the X-data array
    YB = np.flipud(YB)  # Flip the Y-data array


# %% PANEL METHOD GEOMETRY
# Initialize variables
XC = np.zeros(numPan)  # Initialize control point X-coordinate array
YC = np.zeros(numPan)  # Initialize control point Y-coordinate array
S = np.zeros(numPan)  # Initialize panel length array
phi = np.zeros(numPan)  # Initialize panel orientation angle array [deg]

# Find geometric quantities of the airfoil
for i in range(numPan):  # Loop over all panels
    XC[i] = 0.5 * (XB[i] + XB[i + 1])  # X-value of control point
    YC[i] = 0.5 * (YB[i] + YB[i + 1])  # Y-value of control point
    dx = XB[i + 1] - XB[i]  # Change in X between boundary points
    dy = YB[i + 1] - YB[i]  # Change in Y between boundary points
    S[i] = (dx ** 2 + dy ** 2) ** 0.5  # Length of the panel
    phi[i] = math.atan2(dy, dx)  # Angle of panel (positive X-axis to inside face)
    if (phi[i] < 0):  # Make all panel angles positive [rad]
        phi[i] = phi[i] + 2 * np.pi

# Compute angle of panel normal w.r.t. horizontal and include AoA
delta = phi + (np.pi / 2)  # Angle from positive X-axis to outward normal vector [rad]
beta = delta - AoAR  # Angle between freestream vector and outward normal vector [rad]
beta[beta > 2 * np.pi] = beta[beta > 2 * np.pi] - 2 * np.pi  # Make all panel angles between 0 and 2pi [rad]


# %% COMPUTE SOURCE AND VORTEX PANEL STRENGTHS - REF [10]

# Geometric integrals for SPM and VPM (normal [I,K] and tangential [J,L])
I, J = COMPUTE_IJ_SPM(XC, YC, XB, YB, phi, S)  # Call COMPUTE_IJ_SPM function (Refs [2] and [3])
K, L = COMPUTE_KL_VPM(XC, YC, XB, YB, phi, S)  # Call COMPUTE_KL_VPM function (Refs [6] and [7])

# Populate A matrix
# - Simpler option: A = I + np.pi*np.eye(numPan,numPan)
A = np.zeros([numPan, numPan])  # Initialize the A matrix
for i in range(numPan):  # Loop over all i panels
    for j in range(numPan):  # Loop over all j panels
        if (i == j):  # If the panels are the same
            A[i, j] = np.pi  # Set A equal to pi
        else:  # If panels are not the same
            A[i, j] = I[i, j]  # Set A equal to I

# Right column of A matrix
newAV = np.zeros((numPan, 1))  # Used to enlarge the A matrix to account for gamma column
A = np.hstack((A, newAV))  # Horizontally stack the A matrix with newAV to get enlarged matrix
for i in range(numPan):  # Loop over all i panels (rows)
    A[i, numPan] = -sum(K[i, :])  # Add gamma term to right-most column of A matrix

# Bottom row of A matrix
newAH = np.zeros((1, numPan + 1))  # Used to enlarge the A matrix to account for Kutta condition equation
A = np.vstack((A, newAH))  # Vertically stack the A matrix with newAH to get enlarged matrix
for j in range(numPan):  # Loop over all j panels (columns)
    A[numPan, j] = J[0, j] + J[numPan - 1, j]  # Source contribution of Kutta condition equation
A[numPan, numPan] = -(sum(L[0, :] + L[numPan - 1, :])) + 2 * np.pi  # Vortex contribution of Kutta condition equation

# Populate b array
# - Simpler option: b = -Vinf*2*np.pi*np.cos(beta)
b = np.zeros(numPan)  # Initialize the b array
for i in range(numPan):  # Loop over all i panels (rows)
    b[i] = -Vinf * 2 * np.pi * np.cos(beta[i])  # Compute RHS array

# Last element of b array (Kutta condition)
b = np.append(b, -Vinf * 2 * np.pi * (
            np.sin(beta[0]) + np.sin(beta[numPan - 1])))  # Add Kutta condition equation RHS to b array

# Compute result array
resArr = np.linalg.solve(A, b)  # Solve system of equation for all source strengths and single vortex strength

# Separate lam and gamma values from result
lam = resArr[0:len(resArr) - 1]  # All panel source strengths
gamma = resArr[len(resArr) - 1]  # Constant vortex strength


# %% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

# Compute velocities
Vt = np.zeros(numPan)  # Initialize tangential velocity
Cp = np.zeros(numPan)  # Initialize pressure coefficient
for i in range(numPan):  # Loop over all panels
    term1 = Vinf * np.sin(beta[i])  # Uniform flow term
    term2 = (1 / (2 * np.pi)) * sum(lam * J[i, :])  # Source panel terms when j is not equal to i
    term3 = gamma / 2  # Vortex panel term when j is equal to i
    term4 = -(gamma / (2 * np.pi)) * sum(L[i, :])  # Vortex panel terms when j is not equal to i

    Vt[i] = term1 + term2 + term3 + term4  # Compute tangential velocity on panel i
    Cp[i] = 1 - (Vt[i] / Vinf) ** 2  # Compute pressure coefficient on panel i

# %% COMPUTE LIFT AND MOMENT COEFFICIENTS

# Compute normal and axial force coefficients
CN = -Cp * S * np.sin(beta)  # Normal force coefficient []
CA = -Cp * S * np.cos(beta)  # Axial force coefficient []

# Compute lift and moment coefficients
CL = sum(CN * np.cos(AoAR)) - sum(CA * np.sin(AoAR))  # Decompose axial and normal to lift coefficient []
CM = sum(Cp * (XC - 0.25) * S * np.cos(phi))  # Moment coefficient []

# Print the results to the Console
print("======= RESULTS =======")
print("Lift Coefficient (CL)")
print("  SPVP : %2.8f" % CL)  # From this SPVP code
print("  K-J  : %2.8f" % (2 * sum(gamma * S)))  # From Kutta-Joukowski lift equation
print("Moment Coefficient (CM)")
print("  SPVP : %2.8f" % CM)  # From this SPVP code


# %% PLOTTING

# FIGURE: Airfoil with panel normal vectors
if (flagPlot[0] == 1):
    fig = plt.figure(1)  # Create the figure
    plt.cla()  # Clear the axes
    plt.fill(XB, YB, 'k')  # Plot the airfoil
    X = np.zeros(2)  # Initialize 'X'
    Y = np.zeros(2)  # Initialize 'Y'
    for i in range(numPan):  # Loop over all panels
        X[0] = XC[i]  # Set X start of panel orientation vector
        X[1] = XC[i] + S[i] * np.cos(delta[i])  # Set X end of panel orientation vector
        Y[0] = YC[i]  # Set Y start of panel orientation vector
        Y[1] = YC[i] + S[i] * np.sin(delta[i])  # Set Y end of panel orientation vector
        if (i == 0):  # If it's the first panel index
            plt.plot(X, Y, 'b-', label='First Panel')  # Plot normal vector for first panel
        elif (i == 1):  # If it's the second panel index
            plt.plot(X, Y, 'g-', label='Second Panel')  # Plot normal vector for second panel
        else:  # If it's neither the first nor second panel index
            plt.plot(X, Y, 'r-')  # Plot normal vector for all other panels
    plt.xlabel('X Units')  # Set X-label
    plt.ylabel('Y Units')  # Set Y-label
    plt.title('Panel Geometry')  # Set title
    plt.axis('equal')  # Set axes equal
    plt.legend()  # Display legend
    plt.show()  # Display plot

# FIGURE: Geometry with the following indicated:
# - Boundary points, control points, first panel, second panel
if (flagPlot[1] == 1):
    fig = plt.figure(2)  # Create figure
    plt.cla()  # Get ready for plotting
    plt.plot(XB, YB, 'k-')  # Plot airfoil panels
    plt.plot([XB[0], XB[1]], [YB[0], YB[1]], 'b-', label='First Panel')  # Plot first panel
    plt.plot([XB[1], XB[2]], [YB[1], YB[2]], 'g-', label='Second Panel')  # Plot second panel
    plt.plot(XB, YB, 'ko', markerfacecolor='k', label='Boundary Pts')  # Plot boundary points (black circles)
    plt.plot(XC, YC, 'ko', markerfacecolor='r', label='Control Pts')  # Plot control points (red circles)
    plt.xlabel('X Units')  # Set X-label
    plt.ylabel('Y Units')  # Set Y-label
    plt.axis('equal')  # Set axes equal
    plt.legend()  # Display legend
    plt.show()  # Display plot

# FIGURE: Cp vectors at airfoil control points
if (flagPlot[2] == 1):
    fig = plt.figure(3)  # Create figure
    plt.cla()  # Get ready for plotting
    Cps = np.absolute(Cp * 0.15)  # Scale and make positive all Cp values
    X = np.zeros(2)  # Initialize X values
    Y = np.zeros(2)  # Initialize Y values
    for i in range(len(Cps)):  # Loop over all panels
        X[0] = XC[i]  # Control point X-coordinate
        X[1] = XC[i] + Cps[i] * np.cos(delta[i])  # Ending X-value based on Cp magnitude
        Y[0] = YC[i]  # Control point Y-coordinate
        Y[1] = YC[i] + Cps[i] * np.sin(delta[i])  # Ending Y-value based on Cp magnitude

        if (Cp[i] < 0):  # If pressure coefficient is negative
            plt.plot(X, Y, 'r-')  # Plot as a red line
        elif (Cp[i] >= 0):  # If pressure coefficient is zero or positive
            plt.plot(X, Y, 'b-')  # Plot as a blue line
    plt.fill(XB, YB, 'k')  # Plot the airfoil as black polygon
    plt.xlabel('X Units')  # Set X-label
    plt.ylabel('Y Units')  # Set Y-label
    plt.gca().set_aspect('equal')  # Set aspect ratio equal
    plt.show()  # Show the plot

# FIGURE: Pressure coefficient comparison (XFOIL vs. VPM)
if (flagPlot[3] == 1):
    fig = plt.figure(4)  # Create figure
    plt.cla()  # Get ready for plotting
    midIndS = int(np.floor(len(Cp) / 2))  # Airfoil middle index for VPM data
    plt.plot(XC[midIndS + 1:len(XC)], Cp[midIndS + 1:len(XC)],  # Plot Cp for upper surface of airfoil from panel method
             'ks', markerfacecolor='b', label='VPM Upper')
    plt.plot(XC[0:midIndS], Cp[0:midIndS],  # Plot Cp for lower surface of airfoil from panel method
             'ks', markerfacecolor='r', label='VPM Lower')
    plt.xlim(0, 1)  # Set X-limits
    plt.xlabel('X Coordinate')  # Set X-label
    plt.ylabel('Cp')  # Set Y-label
    plt.title('Pressure Coefficient')  # Set title
    plt.legend()  # Display legend
    plt.gca().invert_yaxis()  # Invert Cp (Y) axis
    plt.grid()  #Turn on grid
    plt.show()  # Display plot
