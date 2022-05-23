import numpy as np
import math as math


def Cp_calculatorSP(Xdata, Ydata, XC, YC, S, phi, beta, V_fs, AoAr):
    # Number of panels
    N_pan = len(XC)  # Number of panels/control points

    # Initialize arrays
    I = np.zeros([N_pan, N_pan])  # Initialize I integral matrix
    J = np.zeros([N_pan, N_pan])  # Initialize J integral matrix

    # Compute integral
    for i in range(N_pan):  # Loop over i panels
        for j in range(N_pan):  # Loop over j panels
            if (j != i):  # If the i and j panels are not the same
                # Compute intermediate values
                A = -(XC[i] - Xdata[j]) * np.cos(phi[j]) - (YC[i] - Ydata[j]) * np.sin(phi[j])  # A term
                B = (XC[i] - Xdata[j])**2 + (YC[i] - Ydata[j])**2  # B term
                Cn = np.sin(phi[i] - phi[j])  # C term (normal)
                Dn = -(XC[i] - Xdata[j]) * np.sin(phi[i]) + (YC[i] - Ydata[j]) * np.cos(phi[i])  # D term (normal)
                Ct = -np.cos(phi[i] - phi[j])  # C term (tangential)
                Dt = (XC[i] - Xdata[j]) * np.cos(phi[i]) + (YC[i] - Ydata[j]) * np.sin(phi[i])  # D term (tangential)
                E = np.sqrt(B - A**2)  # E term
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
                    I[i, j] = 0  # Set I value equal to zero
                    J[i, j] = 0  # Set J value equal to zero
                else:
                    # Compute I (needed for normal velocity), Ref [1]
                    term1 = 0.5 * Cn * np.log((S[j]**2 + 2 * A * S[j] + B) / B)  # First term in I equation
                    term2 = ((Dn - A * Cn) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in I equation
                    I[i, j] = term1 + term2  # Compute I integral

                    # Compute J (needed for tangential velocity), Ref [2]
                    term1 = 0.5 * Ct * np.log((S[j]**2 + 2 * A * S[j] + B) / B)  # First term in I equation
                    term2 = ((Dt - A * Ct) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in I equation
                    J[i, j] = term1 + term2  # Compute J integral

            # Zero out any problem values
            if (np.iscomplex(I[i, j]) or np.isnan(I[i, j]) or np.isinf(I[i, j])):  # If I term is complex or a NAN or an INF
                I[i, j] = 0  # Set I value equal to zero
            if (np.iscomplex(J[i, j]) or np.isnan(J[i, j]) or np.isinf(J[i, j])):  # If J term is complex or a NAN or an INF
                J[i, j] = 0  # Set J value equal to zero

    # Populate A matrix
    # - Simpler option: A = I + np.pi*np.eye(N_pan,N_pan)
    A = np.zeros([N_pan, N_pan])  # Initialize the A matrix
    for i in range(N_pan):  # Loop over all i panels
        for j in range(N_pan):  # Loop over all j panels
            if (i == j):  # If the panels are the same
                A[i, j] = np.pi  # Set A equal to pi
            else:  # If panels are not the same
                A[i, j] = I[i, j]  # Set A equal to I

    # Populate b array
    # - Simpler option: b = -V_fs*2*np.pi*np.cos(beta)
    b = np.zeros(N_pan)  # Initialize the b array
    for i in range(N_pan):  # Loop over all i panels (rows)
        b[i] = -V_fs * 2 * np.pi * np.cos(beta[i])  # Compute RHS array

    #solve linear system of equations
    lam = np.linalg.solve(A, b)

    # Compute velocities
    # - Simpler method: Vt = V_fs*np.sin(beta) + np.dot(J,lam)/(2*np.pi)
    #                   C_p = 1 - (Vt/V_fs)**2
    Vt = np.zeros(N_pan)  # Initialize tangential velocity array
    C_p = np.zeros(N_pan)  # Initialize pressure coefficient array
    for i in range(N_pan):  # Loop over all i panels
        addVal = 0  # Reset the summation value to zero
        for j in range(N_pan):  # Loop over all j panels
            addVal = addVal + (lam[j] / (2 * np.pi)) * J[i, j]  # Sum all tangential source panel terms

        Vt[i] = V_fs * np.sin(beta[i]) + addVal  # Compute tangential velocity by adding uniform flow term
        C_p[i] = 1 - (Vt[i] / V_fs) ** 2  # Compute pressure coefficient
    # %% COMPUTE LIFT AND DRAG

    # Compute normal and axial force coefficients
    # CN = -C_p * S * np.sin(beta)  # Normal force coefficient []
    # CA = -C_p * S * np.cos(beta)  # Axial force coefficient []

    # Compute lift, drag, and moment coefficients
    # C_l = sum(CN * np.cos(AoA)) - sum(CA * np.sin(AoA))  # Decompose axial and normal to lift coefficient []
    # C_d = sum(CN * np.sin(AoA)) + sum(CA * np.cos(AoA))  # Decompose axial and normal to drag coefficient []
    # C_m = sum(C_p * (XC - 0.25) * S * np.cos(phi))

    # # Compute grid point velocity magnitude and pressure coefficient
    # Vxy  = np.sqrt(Vx**2 + Vy**2)                                               # Compute magnitude of velocity vector
    # C_pXY = 1 - (Vxy/V_fs)**2                                                    # Pressure coefficient []
    #
    # plt.figure()  # Create figure
    # plt.cla()  # Get ready for plotting
    # plt.contourf(XX, YY, C_pXY, 500, cmap='jet')  # Plot contour
    # plt.fill(Xdata, Ydata, 'k')  # Plot airfoil as black polygon
    # plt.xlabel('X-Axis')  # Set X-label
    # plt.ylabel('Y-Axis')  # Set Y-label
    # plt.gca().set_aspect('equal')  # Set axes equal
    # plt.xlim(xVals)  # Set X-limits
    # plt.ylim(yVals)  # Set Y-limits
    # plt.show()  # Display plot

    return C_p#, C_l, C_d, C_m


def Cp_calculatorVP(Xdata, Ydata, XC, YC, S, phi, beta, V_fs, AoAr):

    # Number of panels
    N_pan = len(XC)

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
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):  # If E term is 0 or complex or a NAN or an INF
                    K[i, j] = 0  # Set K value equal to zero
                    L[i, j] = 0  # Set L value equal to zero
                else:
                    # Compute K
                    term1 = 0.5 * Cn * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)  # First term in K equation
                    term2 = ((Dn - A * Cn) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in K equation
                    K[i, j] = term1 + term2  # Compute K integral

                    # Compute L
                    term1 = 0.5 * Ct * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)  # First term in L equation
                    term2 = ((Dt - A * Ct) / E) * (math.atan2((S[j] + A), E) - math.atan2(A, E))  # Second term in L equation
                    L[i, j] = term1 + term2  # Compute L integral

            # Zero out any problem values
            if (np.iscomplex(K[i, j]) or np.isnan(K[i, j]) or np.isinf(K[i, j])):  # If K term is complex or a NAN or an INF
                K[i, j] = 0  # Set K value equal to zero
            if (np.iscomplex(L[i, j]) or np.isnan(L[i, j]) or np.isinf(L[i, j])):  # If L term is complex or a NAN or an INF
                L[i, j] = 0  # Set L value equal to zero

    # Populate A matrix
    A = np.zeros([N_pan, N_pan])  # Initialize the A matrix
    for i in range(N_pan):  # Loop over all i panels
        for j in range(N_pan):  # Loop over all j panels
            if (i == j):  # If the panels are the same
                A[i, j] = 0  # Set A equal to zero
            else:  # If panels are not the same
                A[i, j] = -K[i, j]  # Set A equal to negative geometric integral

    # Populate b array
    b = np.zeros(N_pan)  # Initialize the b array
    for i in range(N_pan):  # Loop over all panels
        b[i] = -V_fs * 2 * np.pi * np.cos(beta[i])  # Compute RHS array

    # Satisfy the Kutta condition
    pct = 100  # Panel replacement percentage
    panRep = int((pct / 100) * N_pan) - 1  # Replace this panel with Kutta condition equation
    if (panRep >= N_pan):  # If we specify the last panel
        panRep = N_pan - 1  # Set appropriate replacement panel index
    A[panRep, :] = 0  # Set all colums of the replaced panel equation to zero
    A[panRep, 0] = 1  # Set first column of replaced panel equal to 1
    A[panRep, N_pan - 1] = 1  # Set last column of replaced panel equal to 1
    b[panRep] = 0  # Set replaced panel value in b array equal to zero

    # Compute gamma values
    gamma = np.linalg.solve(A, b)  # Compute all vortex strength values

    # Compute velocities
    Vt = np.zeros(N_pan)  # Initialize tangential velocity array
    C_p = np.zeros(N_pan)  # Initialize pressure coefficient array
    for i in range(N_pan):  # Loop over all i panels
        addVal = 0  # Reset summation value to zero
        for j in range(N_pan):  # Loop over all j panels
            addVal = addVal - (gamma[j] / (2 * np.pi)) * L[i, j]  # Sum all tangential vortex panel terms

        Vt[i] = V_fs * np.sin(beta[i]) + addVal + gamma[i] / 2  # Compute tangential velocity by adding uniform flow and i=j terms
        C_p[i] = 1 - (Vt[i] / V_fs) ** 2  # Compute pressure coefficient

    # # %% COMPUTE LIFT AND MOMENT COEFFICIENTS
    #
    # # Compute normal and axial force coefficients
    # CN = -C_p * S * np.sin(beta)  # Normal force coefficient []
    # CA = -C_p * S * np.cos(beta)  # Axial force coefficient []
    #
    # # Compute lift, drag, and moment coefficients
    # CL = sum(CN * np.cos(AoAr)) - sum(CA * np.sin(AoAr))  # Decompose axial and normal to lift coefficient []
    # CM = sum(C_p * (XC - 0.25) * S * np.cos(phi))  # Moment coefficient []
    #
    # # Print the results to the Console
    # print("======= RESULTS =======")
    # print("Lift Coefficient (CL)")
    # print("  K-J  : %2.8f" % (2 * sum(gamma * S)))  # From Kutta-Joukowski lift equation
    # print("  VPM  : %2.8f" % CL)  # From this VPM code
    # print("Moment Coefficient (CM)")
    # print("  VPM  : %2.8f" % CM)  # From this VPM code
    #
    # # %% COMPUTE STREAMLINES - REF [4]
    #
    #
    # nGridX = 150  # X-grid for streamlines and contours
    # nGridY = 150  # Y-grid for streamlines and contours
    # xVals = [-0.5, 1.5]  # X-grid extents [min, max]
    # yVals = [-0.3, 0.3]  # Y-grid extents [min, max]
    #
    # # Streamline parameters
    # slPct = 25  # Percentage of streamlines of the grid
    # Ysl = np.linspace(yVals[0], yVals[1], int((slPct / 100) * nGridY))  # Create array of Y streamline starting points
    # Xsl = xVals[0] * np.ones(len(Ysl))  # Create array of X streamline starting points
    # XYsl = np.vstack((Xsl.T, Ysl.T)).T  # Concatenate X and Y streamline starting points
    #
    # # Generate the grid points
    # Xgrid = np.linspace(xVals[0], xVals[1], nGridX)  # X-values in evenly spaced grid
    # Ygrid = np.linspace(yVals[0], yVals[1], nGridY)  # Y-values in evenly spaced grid
    # XX, YY = np.meshgrid(Xgrid, Ygrid)  # Create meshgrid from X and Y grid arrays
    #
    # # Initialize velocities
    # Vx = np.zeros([nGridX, nGridY])  # Initialize X velocity matrix
    # Vy = np.zeros([nGridX, nGridY])  # Initialize Y velocity matrix
    #
    # # Path to figure out if grid point is inside polygon or not
    # AF = np.vstack((Xdata.T, Ydata.T)).T  # Concatenate Xdata and Ydata geometry points
    # afPath = path.Path(AF)  # Create a path for the geometry
    #
    # # Solve for grid point X and Y velocities
    # for m in range(nGridX):  # Loop over X-grid points
    #     for n in range(nGridY):  # Loop over Y-grid points
    #         XP = XX[m, n]  # Current iteration's X grid point
    #         YP = YY[m, n]  # Current iteration's Y grid point
    #         Nx, Ny = STREAMLINE_VPM(XP, YP, Xdata, Ydata, phi, S)  # Compute Nx and Ny geometric integrals
    #
    #         # Check if grid points are in object
    #         # - If they are, assign a velocity of zero
    #         if afPath.contains_points([(XP, YP)]):  # If (XP,YP) is in the body
    #             Vx[m, n] = 0  # Set X-velocity equal to zero
    #             Vy[m, n] = 0  # Set Y-velocity equal to zero
    #         else:
    #             Vx[m, n] = V_fs * np.cos(AoAr) + sum(-gamma * Nx / (2 * np.pi))  # Compute X-velocity
    #             Vy[m, n] = V_fs * np.sin(AoAr) + sum(-gamma * Ny / (2 * np.pi))  # Compute Y-velocity
    #
    # # Compute grid point velocity magnitude and pressure coefficient
    # Vxy = np.sqrt(Vx ** 2 + Vy ** 2)  # Compute magnitude of velocity vector []
    # CpXY = 1 - (Vxy / V_fs) ** 2  # Pressure coefficient []

    return C_p#, CpXY, XX, YY

def Cp_calculatorVPSP(Xdata, Ydata, XC, YC, S, phi, beta, V_fs, AoAr):
    # Number of panels
    N_pan = len(XC)  # Number of panels/control points
    #SP Method
    # Initialize arrays
    I = np.zeros([N_pan, N_pan])  # Initialize I integral matrix
    J = np.zeros([N_pan, N_pan])  # Initialize J integral matrix

    # Compute integral
    # Compute integral
    for i in range(N_pan):  # Loop over i panels
        for j in range(N_pan):  # Loop over j panels
            if (j != i):  # If the i and j panels are not the same
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


    #VP Method
    # Initialize arrays
    K = np.zeros([N_pan, N_pan])  # Initialize K integral matrix
    L = np.zeros([N_pan, N_pan])  # Initialize L integral matrix

    # Compute integral
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
    # - Simpler option: A = I + np.pi*np.eye(N_pan,N_pan)
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
    # - Simpler option: b = -V_fs*2*np.pi*np.cos(beta)
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


    # %% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

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
    CN = -Cp * S * np.sin(beta)  # Normal force coefficient []
    CA = -Cp * S * np.cos(beta)  # Axial force coefficient []

    # Compute lift and moment coefficients
    CL = sum(CN * np.cos(AoAr)) - sum(CA * np.sin(AoAr))  # Decompose axial and normal to lift coefficient []
    CM = sum(Cp * (XC - 0.25) * S * np.cos(phi))  # Moment coefficient []

    # Print the results to the Console
    print("======= RESULTS =======")
    print("Lift Coefficient (CL)")
    print("  SPVP : %2.8f" % CL)  # From this SPVP code
    print("  K-J  : %2.8f" % (2 * sum(gamma * S)))  # From Kutta-Joukowski lift equation
    print("Moment Coefficient (CM)")
    print("  SPVP : %2.8f" % CM)  # From this SPVP code

    return Cp

