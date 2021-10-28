import numpy as np
import math as m

def airfoil_geometery(Xdata, Ydata, AoA):
    # Number of panels
    N_pan= len(Xdata) - 1

    # Check for direction of points
    edge = np.zeros(N_pan)
    for i in range(N_pan):
        edge[i] = (Xdata[i+1]-Xdata[i])*(Ydata[i+1]+Ydata[i])
    sumEdge = np.sum(edge)

    # If panels are CCW, flip them (don't if CW)
    if (sumEdge < 0):
        Xdata = np.flipud(Xdata)
        Ydata = np.flipud(Ydata)

    #COMPUTE GEOMETRIC VARIABLES

    # Initialize variables
    XC = np.zeros(N_pan)
    YC = np.zeros(N_pan)
    S = np.zeros(N_pan)
    phi = np.zeros(N_pan)

    # Find geometric quantities of the airfoil
    for i in range(N_pan):
        XC[i] = 0.5 * (Xdata[i] + Xdata[i + 1])
        YC[i] = 0.5 * (Ydata[i] + Ydata[i + 1])
        dx = Xdata[i + 1] - Xdata[i]
        dy = Ydata[i + 1] - Ydata[i]
        S[i] = (dx ** 2 + dy ** 2) ** 0.5
        phi[i] = m.atan2(dy, dx)
        if (phi[i] < 0):
            phi[i] = phi[i] + 2 * np.pi

    # Compute angle of panel normal w.r.t. horizontal and include AoA
    AoAr = AoA*np.pi/180
    delta = phi + (np.pi / 2)  # Panel normal angle [rad]
    beta = delta - (AoAr * (np.pi / 180))  # Angle between freestream and panel normal [rad]
    beta[beta > 2 * np.pi] = beta[beta > 2 * np.pi] - 2 * np.pi

    return XC, YC, S, phi, N_pan, beta