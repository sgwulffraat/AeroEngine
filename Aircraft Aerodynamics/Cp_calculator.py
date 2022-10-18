import numpy as np
import math as math

def Cp_calculatorVPSP(Xdata, Ydata, XC, YC, S, phi, beta, V_fs, AoAr):
    # Number of panels
    N_pan = len(XC)  # Number of panels/control points

    # %% Source Panel Method
    # Initialize arrays
    I = np.zeros([N_pan, N_pan])  # Initialize I integral matrix
    J = np.zeros([N_pan, N_pan])  # Initialize J integral matrix

    # Compute integral
    for i in range(N_pan):      # Loop over i panels
        for j in range(N_pan):  # Loop over j panels
            if (j != i):        # If the i and j panels are not the same
                # Compute intermediate values
                A = -(XC[i] - Xdata[j]) * np.cos(phi[j]) - (YC[i] - Ydata[j]) * np.sin(phi[j])  # A term
                B = (XC[i] - Xdata[j]) ** 2 + (YC[i] - Ydata[j]) ** 2  # B term
                Cn = np.sin(phi[i] - phi[j])  # C term (normal)
                Dn = -(XC[i] - Xdata[j]) * np.sin(phi[i]) + (YC[i] - Ydata[j]) * np.cos(phi[i])  # D term (normal)
                Ct = -np.cos(phi[i] - phi[j])  # C term (tangential)
                Dt = (XC[i] - Xdata[j]) * np.cos(phi[i]) + (YC[i] - Ydata[j]) * np.sin(phi[i])  # D term (tangential)
                E = np.sqrt(B - A ** 2)  # E term
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(
                        E)):  # If E term is 0 or complex or a NAN or an INF
                    I[i, j] = 0  # Set I value equal to zero
                    J[i, j] = 0  # Set J value equal to zero
                else:
                    # Compute I (needed for normal velocity), Ref [1]
                    term1 = 0.5 * Cn * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)  # First term in I equation
                    term2 = ((Dn - A * Cn) / E) * (
                                math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in I equation
                    I[i, j] = term1 + term2  # Compute I integral

                    # Compute J (needed for tangential velocity), Ref [2]
                    term1 = 0.5 * Ct * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)  # First term in I equation
                    term2 = ((Dt - A * Ct) / E) * (
                                math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in I equation
                    J[i, j] = term1 + term2  # Compute J integral

            # Zero out any problem values
            if (np.iscomplex(I[i, j]) or np.isnan(I[i, j]) or np.isinf(
                    I[i, j])):  # If I term is complex or a NAN or an INF
                I[i, j] = 0  # Set I value equal to zero
            if (np.iscomplex(J[i, j]) or np.isnan(J[i, j]) or np.isinf(
                    J[i, j])):  # If J term is complex or a NAN or an INF
                J[i, j] = 0  # Set J value equal to zero


    # %% Vortex Panel Method
    # Initialize arrays
    K = np.zeros([N_pan, N_pan])  # Initialize K integral matrix
    L = np.zeros([N_pan, N_pan])  # Initialize L integral matrix

    # Compute integral
    for i in range(N_pan):  # Loop over i panels
        for j in range(N_pan):  # Loop over j panels
            if (j != i):  # If panel j is not the same as panel i
                # Compute intermediate values
                A = -(XC[i] - Xdata[j]) * np.cos(phi[j]) - (YC[i] - Ydata[j]) * np.sin(phi[j])  # A term
                B = (XC[i] - Xdata[j]) ** 2 + (YC[i] - Ydata[j]) ** 2  # B term
                Cn = -np.cos(phi[i] - phi[j])  # C term (normal)
                Dn = (XC[i] - Xdata[j]) * np.cos(phi[i]) + (YC[i] - Ydata[j]) * np.sin(phi[i])  # D term (normal)
                Ct = np.sin(phi[j] - phi[i])  # C term (tangential)
                Dt = (XC[i] - Xdata[j]) * np.sin(phi[i]) - (YC[i] - Ydata[j]) * np.cos(phi[i])  # D term (tangential)
                E = np.sqrt(B - A ** 2)  # E term
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(
                        E)):  # If E term is 0 or complex or a NAN or an INF
                    K[i, j] = 0  # Set K value equal to zero
                    L[i, j] = 0  # Set L value equal to zero
                else:
                    # Compute K
                    term1 = 0.5 * Cn * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)  # First term in K equation
                    term2 = ((Dn - A * Cn) / E) * (
                            math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in K equation
                    K[i, j] = term1 + term2  # Compute K integral

                    # Compute L
                    term1 = 0.5 * Ct * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)  # First term in L equation
                    term2 = ((Dt - A * Ct) / E) * (
                            math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in L equation
                    L[i, j] = term1 + term2  # Compute L integral

            # Zero out any problem values
            if (np.iscomplex(K[i, j]) or np.isnan(K[i, j]) or np.isinf(
                    K[i, j])):  # If K term is complex or a NAN or an INF
                K[i, j] = 0  # Set K value equal to zero
            if (np.iscomplex(L[i, j]) or np.isnan(L[i, j]) or np.isinf(
                    L[i, j])):  # If L term is complex or a NAN or an INF
                L[i, j] = 0  # Set L value equal to zero

    # Populate A matrix
    A = np.zeros([N_pan, N_pan])  # Initialize the A matrix
    for i in range(N_pan):  # Loop over all i panels
        for j in range(N_pan):  # Loop over all j panels
            if (i == j):  # If the panels are the same
                A[i, j] = np.pi  # Set A equal to pi
            else:  # If panels are not the same
                A[i, j] = I[i, j]  # Set A equal to I

    # Right column of A matrix
    newAV = np.zeros((N_pan, 1))  # Used to enlarge the A matrix to account for gamma column
    A = np.hstack((A, newAV))  # Horizontally stack the A matrix with newAV to get enlarged matrix
    for i in range(N_pan):  # Loop over all i panels (rows)
        A[i, N_pan] = -sum(K[i, :])  # Add gamma term to right-most column of A matrix

    # Bottom row of A matrix
    newAH = np.zeros((1, N_pan + 1))  # Used to enlarge the A matrix to account for Kutta condition equation
    A = np.vstack((A, newAH))  # Vertically stack the A matrix with newAH to get enlarged matrix
    for j in range(N_pan):  # Loop over all j panels (columns)
        A[N_pan, j] = J[0, j] + J[N_pan - 1, j]  # Source contribution of Kutta condition equation
    A[N_pan, N_pan] = -(sum(L[0, :] + L[N_pan - 1, :])) + 2 * np.pi  # Vortex contribution of Kutta condition equation

    # Populate b array
    b = np.zeros(N_pan)  # Initialize the b array
    for i in range(N_pan):  # Loop over all i panels (rows)
        b[i] = -V_fs * 2 * np.pi * np.cos(beta[i])  # Compute RHS array

    # Last element of b array (Kutta condition)
    b = np.append(b, -V_fs * 2 * np.pi * (
                np.sin(beta[0]) + np.sin(beta[N_pan - 1])))  # Add Kutta condition equation RHS to b array

    # Compute result array
    resArr = np.linalg.solve(A, b)  # Solve system of equation for all source strengths and single vortex strength

    # Separate lam and gamma values from result
    lam = resArr[0:len(resArr) - 1]  # All panel source strengths
    gamma = resArr[len(resArr) - 1]  # Constant vortex strength


    # %% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS ###

    # Compute velocities
    Vt = np.zeros(N_pan)  # Initialize tangential velocity
    Cp = np.zeros(N_pan)  # Initialize pressure coefficient
    for i in range(N_pan):  # Loop over all panels
        term1 = V_fs * np.sin(beta[i])  # Uniform flow term
        term2 = (1 / (2 * np.pi)) * sum(lam * J[i, :])  # Source panel terms when j is not equal to i
        term3 = gamma / 2  # Vortex panel term when j is equal to i
        term4 = -(gamma / (2 * np.pi)) * sum(L[i, :])  # Vortex panel terms when j is not equal to i

        Vt[i] = term1 + term2 + term3 + term4  # Compute tangential velocity on panel i
        Cp[i] = 1 - (Vt[i] / V_fs) ** 2  # Compute pressure coefficient on panel i

    # %% COMPUTE LIFT AND MOMENT COEFFICIENTS

    # Compute normal and axial force coefficients
    Cn = -Cp * S * np.sin(beta)  # Normal force coefficient []
    CA = -Cp * S * np.cos(beta)  # Axial force coefficient []

    # Compute lift and moment coefficients
    Cl = sum(Cn * np.cos(AoAr)) - sum(
        CA * np.sin(AoAr))  # Decompose axial and normal to lift coefficient []
    Cm = sum(Cp * (XC - 0.25) * S * np.cos(phi))  # Moment coefficient []

    return Cp, Cl, Cm
