import os
import numpy as np
import matplotlib.pyplot as plt


def xfoil(NACA, AoA, numNodes):

    # %% CREATE LOADING FILE
    # Knowns
    NACA = NACA
    AoA = str(AoA)
    numNodes = str(numNodes-1)
    saveFlnmAF = 'Save_Airfoil.txt'
    saveFlnmCp = 'Save_Cp.txt'
    saveFlnmPol = 'Save_Pol.txt'
    xfoilFlnm = 'xfoil_input.txt'

    # Delete files if they exist
    if os.path.exists(saveFlnmAF):
        os.remove(saveFlnmAF)



    # Create the airfoil
    fid = open(xfoilFlnm, "w")
    fid.write("NACA " + NACA + "\n")
    fid.write("PPAR\n")
    fid.write("N " + numNodes + "\n")
    fid.write("\n\n")
    fid.write("PSAV " + saveFlnmAF + "\n")
    if os.path.exists(saveFlnmAF):
        fid.write("y \n")
    fid.write("OPER\n")
    fid.write("Pacc 1 \n")
    fid.write("\n\n")
    fid.write("ALFA " + AoA + "\n")
    fid.write("CPWR " + saveFlnmCp + "\n")
    if os.path.exists(saveFlnmCp):
        fid.write("y \n")
    fid.write("PWRT\n")
    fid.write(saveFlnmPol + "\n")
    if os.path.exists(saveFlnmPol):
        fid.write("y \n")

    fid.close()

    # Run the XFoil calling command
    os.system("xfoil.exe < xfoil_input.txt")

    # Delete file after running
    if os.path.exists(xfoilFlnm):
        os.remove(xfoilFlnm)



    # %% READ DATA FILE: PRESSURE COEFFICIENT

    # Load the data from the text file
    fin = open(saveFlnmCp, 'r+')
    fout = open("Save_Cpcheck.txt", 'w+')
    lines = fin.readlines()
    for line in lines:
        fout.write(line.replace('-', ' -'))
    fin.close()
    fout.close()
    dataBuffer = np.loadtxt('Save_Cpcheck.txt', skiprows=3)


    # Extract data from the loaded dataBuffer array
    X_0 = dataBuffer[:, 0]
    Y_0 = dataBuffer[:, 1]
    Cp_0 = dataBuffer[:, 2]

    # Polar
    dataBuffer2 = np.loadtxt(saveFlnmPol, skiprows=12)
    Cl = dataBuffer2[1]


    # %% EXTRACT UPPER AND LOWER AIRFOIL DATA

    # Split XFoil results into (U)pper and (L)ower
    Cp_U = Cp_0[Y_0 >= 0]
    Cp_L = Cp_0[Y_0 < 0]
    X_U = X_0[Y_0 >= 0]
    X_L = X_0[Y_0 < 0]

    return Cp_U, Cp_L, X_U, X_L, Cl

