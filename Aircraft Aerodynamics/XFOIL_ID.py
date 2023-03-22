import os
import numpy as np


def xfoil_id(filename):
    # %% CREATE LOADING FILE
    # Knowns
    re = "1.2E06"
    CL = 0.4
    loadFlnmAF = filename
    saveFlnmAF = 'Save_Airfoil.txt'
    saveFlnmCp = 'Save_Cp.txt'
    xfoilFlnm = 'xfoil_input.txt'
    iter = "20"

    # Delete files if they exist
    if os.path.exists(saveFlnmAF):
        os.remove(saveFlnmAF)

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
    fid.write(re + "\n")
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


