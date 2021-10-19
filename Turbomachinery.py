from math import *
import numpy as np

### Inputs ####
M = 0.78
h = 10668           #m
PR_intake = 0.98
m_dot = 173         #kg/s
PR_fan = 1.4
PR_LPC = 1.7
PR_HPC = 12.5
T_combexit = 1400   #K
n_is_fan = 0.9
n_is_C = 0.92
n_is_T = 0.9
n_mech = 0.99
n_comb = 0.995
PR_comb = 0.96
n_nozzle = 0.98
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
p_t2 = p_t * PR_intake

### LPC conditions ###
p_t25 = PR_LPC * p_t21
T_t25 = isentropcomp_T(T_t21, n_is_C, p_t21, p_t25, k_a)
W_req_LPC = m_dot_core * c_p_a * (T_t25 - T_t21)

### HPC conditions ###
p_t3 = p_t25 * PR_HPC
T_t3 = isentropcomp_T(T_t25, n_is_C, p_t25, p_t3, k_a)
W_req_HPC = m_dot_core * c_p_a * (T_t3- T_t25)

### Combustion conditions ###
m_dot_f = (m_dot_core * c_p_g * (T_combexit - T_t3)) / (n_comb * LHV *10**6)
m_dot_4 = m_dot_core + m_dot_f
p_t4 = PR_comb * p_t3

### HPT conditions ###
T_t45 = T_combexit - (W_req_HPC / (n_mech * m_dot_4 * c_p_g))
p_t45 = p_t4*(1 - 1/n_is_T*(1- T_t45/T_combexit)) ** (k_g/(k_g - 1))

### LPT conditions ###
T_gg = T_t45 - ((W_req_LPC + 1/(BPR+1)* W_req_fan) / (n_mech * m_dot_4 * c_p_g))
p_gg = p_t45 * (1 - 1/n_is_T*(1- T_gg/T_t45)) ** (k_g/(k_g - 1))
T_t5 = T_t45 - ((W_req_LPC +  W_req_fan) / (n_mech * m_dot_4 * c_p_g))
p_t5 = p_t45 * (1 - 1/n_is_T*(1- T_t5/T_t45)) ** (k_g/(k_g - 1))



### nozzle conditions ###
e_c = critPR(n_nozzle, k_g)
PR_noz = p_t5/p_a
v_fs = M * sqrt(k_a*R*T_a)
if PR_noz > e_c:
    T_8 = T_t5 * (2/(k_g+1))
    p_8 = p_t5/e_c
    v_8 = sqrt(k_g*R*T_8)
    rho_8 = p_8/(R*T_8)
    A_8 = m_dot_4/(rho_8*v_8)
    F_core = m_dot_4*(v_8-v_fs) + A_8 * (p_8-p_a)
    v_9eff = F_core / m_dot_4 + v_fs
    #print("Core nozzle is choked")
else:
    p_8 = p_a
    T_8 = isentropexp_T(T_t5, n_nozzle, p_t5, p_a, k_g)
    if T_t5 > T_8:
        v_8 = sqrt(2 * c_p_g * (T_t5 - T_8))
    else:
        v_8 = sqrt(2 * c_p_g * (T_8 - T_t5))
    F_core = m_dot_4 * (v_8 - v_fs)
    v_9eff = v_8
    #print("Core nozzle is not choked")


