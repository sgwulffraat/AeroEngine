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


## Conditions defined as classes
class TO_Rh:
    def __init__(self):
        self.CD = True
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
        self.CD = True
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
        self.CD = False
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
        self.CD = False
        self.h = 12000          #m
        self.M = 1.4
        self.Reheat = True
        self.P7 = 173000      #Pa
        self.T7 = 1800        #K
        self.m_dot = 57.4       #kg/s
        self.T_amb = 223.150    #K
        self.P_amb = 19330.4    #Pa


### Functions ###
def isentropcomp_T(T1, n_is, p1, p2, k):
    T2 = (1 + 1 / n_is * ((p2 / p1) ** ((k - 1) / k) - 1)) * T1
    return T2

def isentropexp_T(T1, n_is, p1, p2, k):
    T2 = (-T1 * n_is * (1 - (p2 / p1) ** ((k - 1) / k) )) + T1
    return T2

def isentropcomp_p(p1, n_is, M, k):
    p2 = (1 + n_is * ( (k-1) / 2 ) * M**2 ) ** ( k / (k - 1)) * p1
    return p2

def totalT(T, M, k):
    Tt = T*(1+(k-1)/2*M**2)
    return Tt

def totalp(p, Tt, T, k):
    pt = p*(Tt/T)**(k/(k-1))
    return pt

def critPR(n_j, k):
    e_c = (1/(1-1/n_j*((k-1)/(k+1)))**(k/(k-1)))
    return e_c

PR_crit = critPR(nozz_eff, k_g)

def Cyclecalculator(P7, T7, T_amb, P_amb, m_dot, M, CD, R, k_g, k_a, nozz_eff, Cp_g):
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
        rho9 = P9/(R * T9)
        V9 = sqrt(2 * Cp_g * dT)
        A9 = m_dot/(rho9 * V9)
        M9 = V9/sqrt(k_g * R * T9)
        F_N = m_dot * (V9 - Vflight)
        return P8, T8, A8, V8, P9, T9, A9, V9, M9, F_N
    else:
        F_m = m_dot * (V8 - Vflight)
        F_p = A8 * (P8 - P_amb)
        F_N = F_m + F_p
        return P8, T8, A8, V8, F_m, F_p, F_N



#Cycle calculations for the various conditions
TO_P8, TO_T8, TO_A8, TO_V8, TO_P9, TO_T9, TO_A9, TO_V9, TO_M9, TO_F_N = Cyclecalculator(TO().P7, TO().T7, TO().T_amb, TO().P_amb, TO().m_dot, TO().M, TO().CD, R, k_g, k_a, nozz_eff, Cp_g)
TO_Rh_P8, TO_Rh_T8, TO_Rh_A8, TO_Rh_V8, TO_Rh_P9, TO_Rh_T9, TO_Rh_A9, TO_Rh_V9, TO_Rh_M9, TO_Rh_F_N = Cyclecalculator(TO_Rh().P7, TO_Rh().T7, TO_Rh().T_amb, TO_Rh().P_amb, TO_Rh().m_dot, TO_Rh().M, TO_Rh().CD, R, k_g, k_a, nozz_eff, Cp_g)
data = {'P_8 [Pa]': [TO_P8, TO_Rh_P8],
        'T_8 [K]': [TO_T8, TO_Rh_T8],
        'A_8 [m^2]': [TO_A8, TO_Rh_A8],
        'V_8 [m/s]' [TO_V8, TO_Rh_V8]}
df = pd.DataFrame(data, index=['Take-Off (No afterburner)', 'Take-Off (Afterburner)'])
df_t = df.T
print(df_t)