# This is a script to determine the thickness of two 4-digit NACA symmetrical airfoils.

# Use statement like:
# Naca_12 = NACAcalculator()

# Thickness of the first airfoil:
# N1 = Naca_12[0]

# Thickness of the second airfoil:
# N2 = Naca_12[1]

### Airfoil Parameter calculation ###
def check_user_input(input):
    if input == "Test camber":
        print("Entry confirmed.")
        return True
    elif len(input) != 3:
        print("Invalid input detected...")
        print()
        return False

    else:
        try:
            # Convert it into integer
            val = int(input)
            print("Entry confirmed.")
            return True
        except ValueError:
            print("Invalid input detected...")
            print()
            return False


def Loop():
    print("Would you like to try again?")
    r = input()
    if r == "yes" or r == "y":
        return True
    if r == "n" or r == "no":
        print("Exiting Code..")
        exit()


def NACAcalculator(): #Ouput: (N1,N2)
    #Input
    print("Enter first 3 digits of your student number: ...")
    sn = input()
    if check_user_input(sn) == True:
        if sn == "Test camber":
            NC = (12,12,"naca0012","naca4412")

        else:
            #Sum of digits
            N = int(str(sn)[0])+int(str(sn)[1])+int(str(sn)[2])

            #Checking for oddness
            if (N % 2) ==0:
                N = N
            else:
                N = N + 1

            #Calculating Airfoil thicknesses
            if N<=12:
                N1 = N
                N2 = int(N*2)
            else:
                N1 = int(N/2)
                N2 = N


            #Determining NACA code
            if N1 < 10:
                N_1 = "naca000" + str(N1)
            else:
                N_1 = "naca00" + str(N1)

            if N2 < 10:
                N_2 = "naca000" + str(N2)
            else:
                N_2 = "naca00" + str(N2)
            NC = (N1, N2, N_1, N_2)
        print()
        print("Selected Airfoils:", NC[2], "and", NC[3])
        print()
        return NC
    else:
        if Loop() == True:
            print()
            NC2 = NACAcalculator()
            return NC2
        else:
            exit()
