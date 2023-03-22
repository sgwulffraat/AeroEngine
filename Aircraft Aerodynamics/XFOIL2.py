import os
import numpy as np


def xfoil(NACA, AoA, re, visc):
    # %% CREATE LOADING FILE
    # Knowns
    AoA = str(AoA)
    re = str(re)
    N = 200
    saveFlnmAF = 'Save_Airfoil.txt'
    saveFlnmCp = 'Save_Cp.txt'
    saveFlnmPol = 'Save_Pol.txt'
    saveFlnmDump = 'Save_Dump.txt'
    xfoilFlnm = 'xfoil_input.txt'
    iter = "20"

    # Delete files if they exist
    if os.path.exists(saveFlnmAF):
        os.remove(saveFlnmAF)
    if os.path.exists(saveFlnmDump):
        os.remove(saveFlnmDump)

    # Create the airfoil
    fid = open(xfoilFlnm, "w")
    fid.write("NACA " + NACA + "\n")


    fid.write("PSAV " + saveFlnmAF + "\n")
    if os.path.exists(saveFlnmAF):
        fid.write("y \n")
    fid.write("PPAR\n")
    fid.write("N " + str(N) + "\n")
    fid.write("\n\n")
    fid.write("OPER\n")
    fid.write("iter\n")
    fid.write(iter + "\n")
    if visc == "on":
        fid.write("visc\n")
        fid.write(re + "\n")
    fid.write("Pacc 1 \n")
    fid.write("\n\n")
    fid.write("ALFA " + AoA + "\n")
    #if AoA == "3":
    #    fid.write("")
    fid.write("CPWR " + saveFlnmCp + "\n")
    if os.path.exists(saveFlnmCp):
        fid.write("y \n")
    fid.write("PWRT\n")
    fid.write(saveFlnmPol + "\n")
    if os.path.exists(saveFlnmPol):
        fid.write("y \n")
    fid.write("DUMP\n")
    fid.write(saveFlnmDump + "\n")
    if os.path.exists(saveFlnmDump):
        fid.write("y \n")

    fid.close()

    # Run the XFoil calling command
    os.system("xfoil.exe < xfoil_input.txt")

    # Delete file after running
    #if os.path.exists(xfoilFlnm):
     #   os.remove(xfoilFlnm)

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

    # Polar
    dataBuffer2 = np.loadtxt(saveFlnmPol, skiprows=12)
    Cl = dataBuffer2[1]
    Cd = dataBuffer2[2]
    bot_xtr = dataBuffer2[6]
    top_xtr = dataBuffer2[5]

    dataBuffer3 = np.loadtxt(saveFlnmDump, skiprows=1)
    Cf_list = []
    x_list = []
    for i in np.arange(1, N+1, 1):
        if dataBuffer3[i, 1] <= 1.0 and dataBuffer3[i, 2] > 0:
            Cf_list.append(dataBuffer3[i, 6])
            x_list.append(dataBuffer3[i, 1])



    return Cl, Cd, bot_xtr, top_xtr, Cf_list, x_list

