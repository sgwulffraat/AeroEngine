from math import *
import numpy as np

### Inputs ####
M = 0.78
h = 10668           #m
m_dot = 23.81       #kg/s
m_fuel = 0.4267     #kg/s
PR_comp = 5.5
T_combexit = 1400   #K
n_is_C = 0.92
n_is_T = 0.9
n_mech = 0.99
n_comb = 1
PR_comb = 0.96
n_nozzle = 1
T_a = 218.8         #K
p_a = 23842         #Pa
R = 287             #J/kg K
LHV = 43            #43MJ/kg
c_p_a = 1000        #J/kg K
c_p_g = 1150        #J/kg K
k_a = 1.4
k_g = 1.33


### Functions ###
def isentropcomp_T(T1, n_is, p1, p2, k):
    T2 = (1 + 1 / n_is * ((p2 / p1) ** ((k - 1) / k) - 1)) * T1
    return T2

def isentropexp_T(T1, n_is, p1, p2, k):
    T2 = T1 * n_is * (1 - (p2 / p1) ** ((k - 1) / k) )
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
    e_c = (1-1/n_j*((k-1)/(k+1)))**(-k/(k-1))
    return e_c

### Inlet conditions ###
T_t = totalT(T_a, M, k_a)
p_t = totalp(p_a, T_t, T_a, k_a)

### Compressor conditions ###
p_t2 = PR_comp * p_t
T_t2 = isentropcomp_T(T_t, n_is_C, p_t, p_t2, k_a)
W_req_C = m_dot * c_p_a * (T_t2 - T_t)

### Combustion conditions ###
T_combexit = T_t2 + (m_dot_f * n_comb * LHV*10**6)/(m_dot*c_p_g)
p_t4 = PR_comb * p_t3


### Turbine conditions ###
T_t5 = T_combexit - (W_req_C / (n_mech * m_dot * c_p_g))
p_t5 = p_t4*(1 - 1/n_is_T*(1- T_t45/T_combexit)) ** (k_g/(k_g - 1))

### nozzle conditions ###
e_c = critPR(n_nozzle, k_g)
PR_noz = p_t5/p_a
v_fs = M * sqrt(k_a*R*T_a)
if PR_noz > e_c:
    T_8 = T_t5 * (2/(k_g+1))
    p_8 = p_t5/e_c
    v_8 = sqrt(k_g*R*T_8)
    rho_8 = p_8/(R*T_8)
    A_8 = m_dot/(rho_8*v_8)
    F_core = m_dot*(v_8-v_fs) + A_8 * (p_8-p_a)
    v_9eff = F_core / m_dot + v_fs
    #print("Nozzle is choked")
else:
    p_8 = p_a
    T_8 = isentropexp_T(T_t5, n_nozzle, p_t5, p_a, k_g)
    if T_t5 > T_8:
        v_8 = sqrt(2 * c_p_g * (T_t5 - T_8))
    else:
        v_8 = sqrt(2 * c_p_g * (T_8 - T_t5))
    F_core = m_dot * (v_8 - v_fs)
    v_9eff = v_8
    #print("Nozzle is not choked")


