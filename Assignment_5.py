import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import *

#Constants
nozz_eff = 0.98
A7 = 0.52           #m^2
k_g = 1.44
k_a = 1.3
Cp_g = 1150         #J/kg*K
R = 287

# For convergent nozzle CD = False, For CD nozzle CD = True
CD = True

## Conditions defined as classes
class TO_Rh:
    def __init__(self):
        self.h = 0              #m
        self.M = 0
        self.Reheat = True
        self.P7 = 303000      #Pa
        self.T7 = 1900        #K
        self.m_dot = 103.5      #kg/s
        self.T_amb = 288.15     #K
        self.P_amb = 101325     #kPa

class TO:
    def __init__(self):
        self.h = 0              #m
        self.M = 0
        self.Reheat = False
        self.P7 = 330000      #Pa
        self.T7 = 800         #K
        self.m_dot = 100        #kg/s
        self.T_amb = 288.15     #K
        self.P_amb = 101325     #Pa

class Cruise:
    def __init__(self):
        self.h = 10000          #m
        self.M = 0.8
        self.Reheat = False
        self.P7 = 147000      #Pa
        self.T7 = 800         #K
        self.m_dot = 44         #kg/s
        self.T_amb = 223.150    #K
        self.P_amb = 26436.3    #Pa

class Cruise_Rh:
    def __init__(self):
        self.h = 12000          #m
        self.M = 1.4
        self.Reheat = True
        self.P7 = 173000      #Pa
        self.T7 = 1800        #K
        self.m_dot = 57.4       #kg/s
        self.T_amb = 223.150    #K
        self.P_amb = 19330.4    #Pa


#Definitions
def critPR(n_j, k):
    e_c = (1/(1-1/n_j*((k-1)/(k+1)))**(k/(k-1)))
    return e_c

def Cyclecalculator(P7, T7, T_amb, P_amb, m_dot, M, R, k_g, k_a, nozz_eff, Cp_g, CD):
    P8 = P7/PR_crit
    T8 = T7 * (2/(k_g+1))
    rho8 = P8/(R * T8)
    V8 = sqrt(k_g * R * T8)
    A8 = m_dot/(rho8 * V8)
    Vflight = M * sqrt(k_a * R * T_amb)
    if CD == True:
        dT = T7 * nozz_eff * (1 - (P_amb/P7)**((k_g - 1)/k_g))
        T9 = T7 - dT
        P9 = P_amb
        V9 = sqrt(2 * Cp_g * dT)
        M9 = V9/sqrt(k_g * R * T9)
        A9 = A8 * ((k_g + 1)/2)**(-1*(k_g + 1)/(2 * (k_g -1))) * ((1 + (k_g -1)/2 * M9**2)**((k_g + 1)/(2 * (k_g - 1))))/M9
        r8 = sqrt(A8/np.pi)
        r9 = sqrt(A9/np.pi)
        L_nozz = (r9-r8)/np.tan(10*(np.pi/180))
        F_N = m_dot * (V9 - Vflight)
        return P8, T8, A8, V8, P9, T9, A9, V9, M9, F_N, L_nozz
    else:
        F_m = m_dot * (V8 - Vflight)
        F_p = A8 * (P8 - P_amb)
        F_N = F_m + F_p
        return P8, T8, A8, V8, F_m, F_p, F_N



#Cycle calculations for the various conditions
PR_crit = critPR(nozz_eff, k_g)

if CD == True:
    TO_P8, TO_T8, TO_A8, TO_V8, TO_P9, TO_T9, TO_A9, TO_V9, TO_M9, TO_F_N, TO_L_nozz = Cyclecalculator(
        TO().P7, TO().T7,TO().T_amb, TO().P_amb,TO().m_dot, TO().M, R, k_g, k_a, nozz_eff, Cp_g, CD)
    TO_Rh_P8, TO_Rh_T8, TO_Rh_A8, TO_Rh_V8, TO_Rh_P9, TO_Rh_T9, TO_Rh_A9, TO_Rh_V9, TO_Rh_M9, TO_Rh_F_N, TO_Rh_L_nozz = Cyclecalculator(
        TO_Rh().P7, TO_Rh().T7, TO_Rh().T_amb, TO_Rh().P_amb, TO_Rh().m_dot, TO_Rh().M, R, k_g, k_a, nozz_eff, Cp_g, CD)
    Cruise_P8, Cruise_T8, Cruise_A8, Cruise_V8, Cruise_P9, Cruise_T9, Cruise_A9, Cruise_V9, Cruise_M9, Cruise_F_N, Cruise_L_nozz = Cyclecalculator(
        Cruise().P7, Cruise().T7, Cruise().T_amb, Cruise().P_amb, Cruise().m_dot, Cruise().M, R, k_g, k_a, nozz_eff, Cp_g, CD)
    Cruise_Rh_P8, Cruise_Rh_T8, Cruise_Rh_A8, Cruise_Rh_V8, Cruise_Rh_P9, Cruise_Rh_T9, Cruise_Rh_A9, Cruise_Rh_V9, Cruise_Rh_M9, Cruise_Rh_F_N, Cruise_Rh_L_nozz = Cyclecalculator(
        Cruise_Rh().P7, Cruise_Rh().T7, Cruise_Rh().T_amb, Cruise_Rh().P_amb, Cruise_Rh().m_dot, Cruise_Rh().M, R, k_g, k_a, nozz_eff, Cp_g, CD)
    data = {'P_8 [Pa]': [TO_P8, TO_Rh_P8, Cruise_P8, Cruise_Rh_P8],
            'T_8 [K]': [TO_T8, TO_Rh_T8, Cruise_T8, Cruise_Rh_T8],
            'A_8 [m^2]': [TO_A8, TO_Rh_A8, Cruise_A8, Cruise_Rh_A8],
            'V_8 [m/s]': [TO_V8, TO_Rh_V8, Cruise_V8, Cruise_Rh_V8],
            'M_8 [-]': [1, 1, 1, 1],
            'P_9 [Pa]': [TO_P9, TO_Rh_P9, Cruise_P9, Cruise_Rh_P9],
            'T_9 [K]': [TO_T9, TO_Rh_T9, Cruise_T9, Cruise_Rh_T9],
            'A_9 [m^2]': [TO_A9, TO_Rh_A9, Cruise_A9, Cruise_Rh_A9],
            'V_9 [m/s]': [TO_V9, TO_Rh_V9, Cruise_V9, Cruise_Rh_V9],
            'M_9 [-]': [TO_M9, TO_Rh_M9, Cruise_M9, Cruise_Rh_M9],
            'F_N [N]': [TO_F_N, TO_Rh_F_N, Cruise_F_N, Cruise_Rh_F_N],
            'Lenght Nozzle [m]': [TO_L_nozz, TO_Rh_L_nozz, Cruise_L_nozz, Cruise_Rh_L_nozz]}
else:
    TO_P8, TO_T8, TO_A8, TO_V8, TO_F_m, TO_F_p, TO_F_N = Cyclecalculator(TO().P7, TO().T7, TO().T_amb,
        TO().P_amb, TO().m_dot, TO().M, R, k_g, k_a, nozz_eff, Cp_g, CD)
    TO_Rh_P8, TO_Rh_T8, TO_Rh_A8, TO_Rh_V8, TO_Rh_F_m, TO_Rh_F_p ,TO_Rh_F_N = Cyclecalculator(TO_Rh().P7,
        TO_Rh().T7, TO_Rh().T_amb, TO_Rh().P_amb, TO_Rh().m_dot, TO_Rh().M, R, k_g, k_a, nozz_eff, Cp_g, CD)
    Cruise_P8, Cruise_T8, Cruise_A8, Cruise_V8, Cruise_F_m, Cruise_F_p, Cruise_F_N = Cyclecalculator(
        Cruise().P7, Cruise().T7, Cruise().T_amb, Cruise().P_amb, Cruise().m_dot, Cruise().M, R, k_g, k_a, nozz_eff, Cp_g, CD)
    Cruise_Rh_P8, Cruise_Rh_T8, Cruise_Rh_A8, Cruise_Rh_V8, Cruise_Rh_F_m, Cruise_Rh_F_p , Cruise_Rh_F_N = Cyclecalculator(
        Cruise_Rh().P7, Cruise_Rh().T7, Cruise_Rh().T_amb, Cruise_Rh().P_amb, Cruise_Rh().m_dot, Cruise_Rh().M, R, k_g, k_a, nozz_eff, Cp_g, CD)
    data = {'P_8 [Pa]': [TO_P8, TO_Rh_P8, Cruise_P8, Cruise_Rh_P8],
            'T_8 [K]': [TO_T8, TO_Rh_T8, Cruise_T8, Cruise_Rh_T8],
            'A_8 [m^2]': [TO_A8, TO_Rh_A8, Cruise_A8, Cruise_Rh_A8],
            'V_8 [m/s]': [TO_V8, TO_Rh_V8, Cruise_V8, Cruise_Rh_V8],
            'M_8 [-]': [1, 1, 1, 1],
            'F_m [N]': [TO_F_m, TO_Rh_F_m, Cruise_F_m, Cruise_Rh_F_m],
            'F_p [N]': [TO_F_p, TO_Rh_F_p, Cruise_F_p, Cruise_Rh_F_p],
            'F_N [N]': [TO_F_N, TO_Rh_F_N, Cruise_F_N, Cruise_Rh_F_N]}
pd.set_option("display.max_columns", 6)
pd.set_option("display.precision", 3)
df = pd.DataFrame(data, index=['Take-Off', 'Take-Off (Afterburner)', 'Cruise', 'Cruise (Afterburner)'])
df_t = df.T

print(df_t)

delta_F_m = np.array([52057.975, 83034.663,	12748.783,21634.573]) - np.array([73901.321, 114234.681, 27583.282, 57589.544])
print("delta F_m = ",delta_F_m)

delta_F_n = np.array([66614.235, 103183.144, 23090.084,	45994.613]) - np.array([73901.321, 114234.681, 27583.282, 57589.544])
print("delta F_n = ",delta_F_n)

delta_A = np.array([])

