import os
import numpy as np


def xfoil_id(filename, re, n_crit, cl):
    # %% CREATE LOADING FILE
    # Knowns
    CL = cl
    loadFlnmAF = filename
    saveFlnmAF = 'Save_Airfoil.txt'
    saveFlnmCp = 'Save_Cp.txt'
    xfoilFlnm = 'xfoil_input.txt'
    iter = "30"
    ncrit = str(n_crit)

    # Delete files if they exist
    #if os.path.exists(saveFlnmAF):
        #os.remove(saveFlnmAF)

    # Create the airfoil
    fid = open(xfoilFlnm, "w")
    if isinstance(filename, str) == True:
        fid.write("load \n")
        fid.write(loadFlnmAF + "\n")
    else:
        fid.write("NACA " + str(filename) + "\n")
    fid.write("PSAV " + saveFlnmAF + "\n")
    if os.path.exists(saveFlnmAF):
        fid.write("y \n")
    fid.write("OPER\n")
    fid.write("iter\n")
    fid.write(iter + "\n")
    fid.write("visc\n")
    fid.write(str(re) + "\n")
    if n_crit != 9:
        fid.write("VPAR\n")
        fid.write("n " + str(n_crit) + "\n\n")
    fid.write("cl \n")
    fid.write(str(CL)+ "\n")
    fid.write("CPWR " + saveFlnmCp + "\n")
    if os.path.exists(saveFlnmCp):
        fid.write("y \n")
    fid.close()

    # Run the XFoil calling command
    os.system("xfoil.exe < xfoil_input.txt")

    # Delete file after running
    #if os.path.exists(xfoilFlnm):
       # os.remove(xfoilFlnm)

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
    X = dataBuffer[:, 0]
    Cp = dataBuffer[:, 2]

    return Cp, X


def xfoil_trs(naca, re, n_crit, cl, trs_point):
    # %% CREATE LOADING FILE
    # Knowns
    CL = cl
    saveFlnmAF = 'Save_Airfoil.txt'
    saveFlnmCp = 'Save_Cp.txt'
    xfoilFlnm = 'xfoil_input.txt'
    saveFlnmPol = 'Save_Pol.txt'
    iter = "30"
    ncrit = str(n_crit)

    # Delete files if they exist
    #if os.path.exists(saveFlnmAF):
        #os.remove(saveFlnmAF)

    # Create the airfoil
    fid = open(xfoilFlnm, "w")
    fid.write("NACA " + str(naca) + "\n")
    fid.write("PSAV " + saveFlnmAF + "\n")
    if os.path.exists(saveFlnmAF):
        fid.write("y \n")
    fid.write("OPER\n")
    fid.write("iter\n")
    fid.write(iter + "\n")
    fid.write("visc\n")
    fid.write(str(re) + "\n")
    fid.write("VPAR\n")
    if n_crit != 9:
        fid.write("n " + ncrit + "\n")
    fid.write("xtr" + "\n")
    fid.write(str(trs_point) + "\n")
    fid.write("1 \n \n")
    fid.write("Pacc 1 \n")
    fid.write("\n\n")
    fid.write("cl \n")
    fid.write(str(CL)+ "\n")
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
    #if os.path.exists(xfoilFlnm):
       # os.remove(xfoilFlnm)

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
    Cd = dataBuffer2[2]

    # Extract data from the loaded dataBuffer array
    X = dataBuffer[:, 0]
    Cp = dataBuffer[:, 2]

    return Cd, Cp, X


