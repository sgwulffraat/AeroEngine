import os
import numpy as np
import matplotlib.pyplot as plt


def xfoil(NACA, AoA, numNodes, YB):

    # %% CREATE LOADING FILE
    # Knowns
    NACA = NACA
    AoA = str(AoA)
    numNodes = str(numNodes)
    saveFlnmAF = 'Save_Airfoil.txt'
    saveFlnmCp = 'Save_Cp.txt'
    xfoilFlnm = 'xfoil_input.txt'

    # Delete files if they exist
    if os.path.exists(saveFlnmAF):
        os.remove(saveFlnmAF)

    if os.path.exists(saveFlnmCp):
        os.remove(saveFlnmCp)

    # Create the airfoil
    fid = open(xfoilFlnm, "w")
    fid.write("NACA " + NACA + "\n")
    fid.write("PPAR\n")
    fid.write("N " + numNodes + "\n")
    fid.write("\n\n")
    fid.write("PSAV " + saveFlnmAF + "\n")
    fid.write("OPER\n")
    fid.write("ALFA " + AoA + "\n")
    fid.write("CPWR " + saveFlnmCp + "\n")
    fid.close()

    # Run the XFoil calling command
    os.system("xfoil.exe < xfoil_input.txt")

    # Delete file after running
    if os.path.exists(xfoilFlnm):
        os.remove(xfoilFlnm)



    # %% READ DATA FILE: PRESSURE COEFFICIENT

    # Load the data from the text file
    dataBuffer = np.loadtxt(saveFlnmCp, skiprows=3)

    # Extract data from the loaded dataBuffer array
    X_0 = dataBuffer[:, 0]
    Y_0 = dataBuffer[:, 1]
    Cp_0 = dataBuffer[:, 2]

    # Delete file after loading
    if os.path.exists(saveFlnmCp):
        os.remove(saveFlnmCp)

    # %% EXTRACT UPPER AND LOWER AIRFOIL DATA

    # Split XFoil results into (U)pper and (L)ower
    print(Cp_0)
    Cp_U = Cp_0[YB >= 0]
    Cp_L = Cp_0[YB < 0]
    X_U = X_0[YB >= 0]
    X_L = X_0[YB < 0]

    return Cp_U, Cp_L, X_U, X_L

