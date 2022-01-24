from math import *
import numpy as np
import pandas as pd

def adiabatic_flame_temperature(phi,Tad_guess,T3):
    "CxHy"
    x = 12
    y = 23
    Xo2 = 0.209476
    Xco2 = 0.000319
    Xn2 = 0.780840
    Xar = 0.009365
    Mo2 = 32
    Mco2 = 44.01
    Mn2 = 28.01
    Mar = 39.948
    Mc12h23 = 167.31102
    hf0_o2 = 0e3
    hf0_co2 = -393.5224e3
    hf0_h2o = -241.8264e3
    hf0_n2 = 0e3
    hf0_Ar = 0e3
    epsilon = (x + y / 4) / Xo2
    nco2 = x + (epsilon/phi) * Xco2
    nh2o = y / 2
    nn2 = (epsilon/phi) * Xn2
    nar = (epsilon/phi) * Xar
    no2 = (epsilon/phi) * Xo2 - y / 4 - x
    Hco2 = H_calculator_CO2(T3)
    Hn2 = H_calculator_N2(T3)
    Har = H_calculator_Ar(T3)
    Ho2 = H_calculator_O2(T3)
    Hkerosene = -211.46
    Hreac = (Hkerosene + (epsilon/phi) * (Xco2*Hco2 + Xn2*Hn2 + Xar*Har + Xo2*Ho2)) * 1000
    Tad = 0
    while abs((Tad/Tad_guess - 1)) > 1e-6:
        T_in = (Tad_guess + T3)/2
        Tad = Tad_guess
        Tad_guess = 298.15 + ((Hreac - nco2*hf0_co2 - nh2o*hf0_h2o)/(no2*CP_calculator_O2(T_in)+nar*CP_calculator_Ar(T_in)+nco2 *
                        CP_calculator_CO2(T_in)+nn2*CP_calculator_N2(T_in)+nh2o*CP_calculator_H20(T_in)))
    return Tad

"functions to calculate standard enthalpy and heat capacity for Argon for a given temperature"
def CP_calculator_Ar(T):
    R = T/1000
    if T < 298:
        print("No value for Cp Ar available")
        Cp = None
    elif 298 <= T < 6000:
        Cp = 20.786 + 2.825911e-7 * R - 1.464191e-7 * R**2 + 1.092131e-8 * R**3 - 3.661371e-8 * R**-2
    else:
        print("No value for Cp Ar available")
        Cp = None
    return Cp

def H_calculator_Ar(T):
    R = T/1000
    if T<298:
        print("No value for H Ar available")
        H = None
    elif 298 <= T < 6000:
        H = 20.786 * R + 2.825911e-7 * R*R/2 - 1.464191e-7 *R*R*R/3 + 1.092131e-8 * R*R*R*R/4 + 3.661371e-8 * R**-1 - 6.197350
    else:
        print("No value for H Ar available")
        H = None
    return H

"functions to calculate standard enthalpy and heat capacity for carbondioxide for a given temperature"
def CP_calculator_CO2(T):
    R = T/1000
    if T < 298:
        print("No value for Cp CO2 available")
        Cp = None
    elif 298 <= T < 1200:
        Cp = 24.99735 + 55.18696 * R - 33.69137 * R**2 + 7.948387 * R**3 - 0.136638 * R**-2
    elif 1200 <= T < 6000:
        Cp = 58.16639 + 2.720074 * R - 0.492289 * R**2 + 0.038844 * R**3 - 6.447293 * R**-2
    else:
        print("No value for Cp CO2 available")
        Cp = None
    return Cp

def H_calculator_CO2(T):
    R = T/1000
    if T < 298:
        print("No value for H CO2 available")
        H = None
    elif 298 <= T < 1200:
        H = 24.99735 * R + 55.18696 * R**2/2 - 33.69137 * R**3/3 + 7.948387 * R**4/4 + 0.136638 * R**-1 - 403.6075
    elif 1200 <= T < 6000:
        H = 58.16639 * R + 2.720074 * R**2/2 - 0.492289 * R**3/3 + 0.038844 * R**4/4 + 6.447293 * R**-1 - 425.9186
    else:
        print("No value for Cp CO2 available")
        H = None
    return H

"functions to calculate standard enthalpy and heat capacity for water vapour for a given temperature"
def CP_calculator_H20(T):
    R = T/1000
    if T < 500:
        print("No value for Cp H2O available")
        Cp = None
    elif 500 <= T < 1700:
        Cp = 30.09200 + 6.832514 * R + 6.793435 * R**2 - 2.534480 * R**3 + 0.082139 * R**-2
    elif 1700 <= T < 6000:
        Cp = 41.96426 + 8.622053 * R - 1.499780 * R**2 + 0.098119 * R**3 - 11.15764 * R**-2
    else:
        print("No value for Cp H2O available")
        Cp = None
    return Cp

def H_calculator_H20(T):
    R = T/1000
    if T < 500:
        print("No value for H H2O available")
        H = None
    elif 500 <= T < 1700:
        H = 30.09200*R + 6.832514 * R**2/2 + 6.793435 * R**3/3 - 2.534480 * R**4/4 - 0.082139 * R**-1 - 250.8810
    elif 1700 <= T < 6000:
        H = 41.96426*R + 8.622053 * R**2/2 - 1.499780 * R**3/3 + 0.098119 * R**4/4 + 11.15764 * R**-1 - 272.1797
    else:
        print("No value for H H2O available")
        H = None
    return H

"functions to calculate standard enthalpy and heat capacity for nitrogen for a given temperature"
def CP_calculator_N2(T):
    R = T/1000
    if T < 100:
        print("No value for Cp N2 available")
        Cp = None
    elif 100 <= T < 500:
        Cp = 28.98641 + 1.853978 * R - 9.647459 * R**2 + 16.63537 * R**3 + 0.000117 * R**-2
    elif 500 <= T < 2000:
        Cp = 19.50583 + 19.88705 * R - 8.598535 * R**2 + 1.369784 * R**3 + 0.527601 * R**-2
    elif 2000 <= T < 6000:
        Cp = 35.51872 + 1.128728 * R - 0.196103 * R ** 2 + 0.014662 * R ** 3 - 4.553760 * R ** -2
    else:
        print("No value for Cp N2 available")
        Cp = None
    return Cp

def H_calculator_N2(T):
    R = T/1000
    if T < 100:
        print("No value for Cp N2 available")
        H = None
    elif 100 <= T < 500:
        H = 28.98641 * R + 1.853978 * R**2/2 - 9.647459 * R**3/3 + 16.63537 * R**4/4 - 0.000117 * R**-1 -8.671914
    elif 500 <= T < 2000:
        H = 19.50583 * R + 19.88705 * R**2/2 - 8.598535 * R**3/3 + 1.369784 * R**4/4 - 0.527601 * R**-1 - 4.935202
    elif 2000 <= T < 6000:
        H = 35.51872 * R + 1.128728 * R**2/2 - 0.196103 * R ** 3/3 + 0.014662 * R ** 4/4 + 4.553760 * R ** -1 - 18.97091
    else:
        print("No value for Cp N2 available")
        H = None
    return H

"functions to calculate standard enthalpy and heat capacity for oxygen for a given temperature"
def CP_calculator_O2(T):
    R = T/1000
    if T < 100:
        print("No value for Cp O2 available")
        Cp = None
    elif 100 <= T < 700:
        Cp = 31.32234 - 20.23531 * R + 57.86644 * R**2 - 36.50624 * R**3 - 0.007374 * R**-2
    elif 700 <= T < 2000:
        Cp = 30.03235 + 8.772972 * R - 3.988133 * R**2 + 0.788313 * R**3 - 0.741599 * R**-2
    elif 2000 <= T < 6000:
        Cp = 20.91111 + 10.72071 * R - 2.020498 * R**2 + 0.146449 * R ** 3 + 9.245722 * R ** -2
    else:
        print("No value for Cp O2 available")
        Cp = None
    return Cp

def H_calculator_O2(T):
    R = T/1000
    if T < 100:
        print("No value for Cp O2 available")
        H = None
    elif 100 <= T < 700:
        H = 31.32234*R - 20.23531 * R**2/2 + 57.86644 * R**3/3 - 36.50624 * R**4/4 + 0.007374 * R**-1 - 8.903471
    elif 700 <= T < 2000:
        H = 30.03235*R + 8.772972 * R**2/2 - 3.988133 * R**3/3 + 0.788313 * R**4/4 + 0.741599 * R**-1 - 11.32468
    elif 2000 <= T < 6000:
        H = 20.91111*R + 10.72071 * R**2/2 - 2.020498 * R**3/3 + 0.146449 * R ** 4/4 - 9.245722 * R ** -1 + 5.337651
    else:
        print("No value for Cp O2 available")
        H = None
    return H

n_comb = 0.99   #[-]
c_p_g = 1150    #[-]
LHV = ((48*393.5+46*241.826)-211.46*4)/(167.31102*4)    #[MJ/kg]
HHV = ((48*393.5+46*285.83)-211.46*4)/(167.31102*4)    #[MJ/kg]
T_3 = np.array([612, 723, 765, 820])    #[K]
T_4 = np.array([1234, 1423, 1523, 1631])#[K]
T_mean = (T_3 + T_4)/2                  #[K]
P_3 = np.array([1166,2188,2700,3100]) * 1000    #[kPa]
V_can = 0.02        #[m^3]
gas_constant = 287      #[J/mol/K]
m_dot_air = np.array([17.12, 29.17, 34.58, 39.58])  #[kg/s]
m_dot_f = (m_dot_air * c_p_g * (T_4 - T_3)) / (n_comb * LHV *10**6) #[kg/s]
print("the lhv is", LHV)
"CxHy"
x = 12
y = 23
Xo2 = 0.209476
Xco2 = 0.000319
Xn2 = 0.780840
Xar = 0.009365
Mo2 = 32
Mco2 = 44.01
Mn2 = 28.01
Mar = 39.948
Mc12h23 = 167.31102
epsilon = (x+y/4)/Xo2
Ma = Xo2 * Mo2 + Xn2 * Mn2 + Xco2 * Mco2 + Xar * Mar
FAR_stoich = 1/epsilon*Mc12h23/Ma
phi_overall = m_dot_f / m_dot_air / FAR_stoich

AF_stoch = (71*3.76*28.01+71*32)/(4*167.31102)  #[-]
AF = m_dot_air/m_dot_f
AF_PZ = 0.32*AF
AF_SZ = 0.52*AF
AF_DZ = 0.80*AF
ER = AF_stoch/AF    #[-]
ER_PZ = AF_stoch/AF_PZ
ER_SZ = AF_stoch/AF_SZ
ER_DZ = AF_stoch/AF_DZ

print(ER_PZ)
print(T_3)
print("Adiabatic flame temperature pz, sz, dz F = 40kN")
print(adiabatic_flame_temperature(ER_PZ[0],2000,T_3[0]))
print(adiabatic_flame_temperature(ER_SZ[0],2000,T_3[0]))
print(adiabatic_flame_temperature(ER_DZ[0],2000,T_3[0]))

print("Adiabatic flame temperature pz, sz, dz F = 80kN")
print(adiabatic_flame_temperature(ER_PZ[1],2000,T_3[1]))
print(adiabatic_flame_temperature(ER_SZ[1],2000,T_3[1]))
print(adiabatic_flame_temperature(ER_DZ[1],2000,T_3[1]))

print("Adiabatic flame temperature pz, sz, dz F = 100kN")
print(adiabatic_flame_temperature(ER_PZ[2],2000,T_3[2]))
print(adiabatic_flame_temperature(ER_SZ[2],2000,T_3[2]))
print(adiabatic_flame_temperature(ER_DZ[2],2000,T_3[2]))

print("Adiabatic flame temperature pz, sz, dz F = 120kN")
print(adiabatic_flame_temperature(ER_PZ[3],2000,T_3[3]))
print(adiabatic_flame_temperature(ER_SZ[3],2000,T_3[3]))
print(adiabatic_flame_temperature(ER_DZ[3],2000,T_3[3]))


tau = (V_can*P_3)/(m_dot_air*gas_constant*T_mean)
print(tau)

heat_density  = (m_dot_f*LHV*1000000*n_comb)/(V_can*P_3)
print(heat_density/1000000)

data = {'air mass flow [kg/s] ':[m_dot_air[0], m_dot_air[1], m_dot_air[2], m_dot_air[3]],
        'Fuel mass flow [kg/s]':[m_dot_f[0], m_dot_f[1], m_dot_f[2], m_dot_f[3]],
        'HPC exit Pressure [Pa]':[P_3[0], P_3[1], P_3[2], P_3[3]],
        'HPC exit Temperature [K]':[T_3[0], T_3[1], T_3[2], T_3[3]],
        'Burner exit Temperature [K]':[T_4[0], T_4[1], T_4[2], T_4[3]],
        'Overall Equivalence ratio [-]': [phi_overall[0], phi_overall[1], phi_overall[2], phi_overall[3]],
        'Equivalence ratio PZ [-]': [ER_PZ[0], ER_PZ[1], ER_PZ[2], ER_PZ[3]],
        'Equivalence ratio SZ [-]': [ER_SZ[0], ER_SZ[1], ER_SZ[2], ER_SZ[3]],
        'Equivalence ratio DZ [-]': [ER_DZ[0], ER_DZ[1], ER_DZ[2], ER_DZ[3]],
        'Residence Time [s]':[tau[0], tau[1], tau[2], tau[3]],
        'Heat Density [MW/m^3]':[heat_density[0], heat_density[1], heat_density[2], heat_density[3]]}
df = pd.DataFrame(data, index=['40', '80', '100', '120'])

with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(df)