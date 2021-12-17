from math import *
import numpy as np

n_comb = 0.99   #[-]
c_p_g = 1150    #[-]
LHV = ((48*393.5+46*241.826)-211.46*4)/(167.31102*4)    #[MJ/kg]
T_3 = np.array([612, 723, 765, 820])    #[K]
T_4 = np.array([1234, 1423, 1523, 1631])    #[K]
m_dot_air = np.array([17.12, 29.17, 34.58, 39.58])  #[kg/s]
m_dot_f = (m_dot_air * c_p_g * (T_4 - T_3)) / (n_comb * LHV *10**6) #[kg/s]
print("mfuel is", m_dot_f)
AF_stoch = (71*3.76*28.01+71*32)/(4*167.31102)  #[-]
AF = m_dot_air/m_dot_f#[-]
AF_1 = 0.32*AF
AF_2 = 0.72*AF

ER = AF_stoch/AF    #[-]
ER_1 = AF_stoch/AF_1
ER_2 = AF_stoch/AF_2
print("equivalence ratio in PZ", ER_1)
print("equivalence ratio in SZ", ER_2)
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