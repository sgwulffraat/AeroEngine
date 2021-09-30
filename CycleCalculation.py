from math import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyromat as pm
from scipy import interpolate
#test 2

#test
### Inputs ####
M = 0.78
h = 10668           #m
PR_intake = 0.98
m_dot = 173         #kg/s
BPR = 12
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
    e_c = (1-1/n_j*((k-1)/(k+1)))**(-k/(k-1))
    return e_c

### Inlet conditions ###
T_t = totalT(T_a, M, k_a)
p_t = totalp(p_a, T_t, T_a, k_a)
p_t2 = p_t * PR_intake
m_dot_core = m_dot/(BPR+1)
m_dot_bypass = m_dot_core * BPR

### Fan conditions ###
p_t21 = p_t2 * PR_fan
T_t21 = isentropcomp_T(T_t, n_is_fan, p_t2, p_t21, k_a)
W_req_fan = m_dot*c_p_a*(T_t21-T_t)

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
T_gg = T_t45 - ((W_req_LPC + 1/(BPR+1)*W_req_fan) / (n_mech * m_dot_4 * c_p_g))
p_gg = p_t45 * (1 - 1/n_is_T*(1- T_gg/T_t45)) ** (k_g/(k_g - 1))
T_t5 = T_t45 - ((W_req_LPC + W_req_fan) / (n_mech * m_dot_4 * c_p_g))
p_t5 = p_t45 * (1 - 1/n_is_T*(1- T_t5/T_t45)) ** (k_g/(k_g - 1))



### Core nozzle conditions ###
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

### Bypass nozzle conditions ###
e_c_bypass = critPR(n_nozzle, k_a)
PR_bypass = p_t21/p_a
if PR_bypass > e_c_bypass:
    T_18 = T_t21*(2/(k_a+1))
    p_18 = p_t21/e_c_bypass
    v_18 = sqrt(k_a * R * T_18)
    rho_18 = p_18 / (R * T_18)
    A_18 = m_dot_bypass / (rho_18 * v_18)
    F_bypass = m_dot_bypass * (v_18 - v_fs) + A_18 * (p_18 - p_a)
    v_19eff = F_bypass / m_dot_bypass + v_fs
    #print("Bypass nozzle is choked")
else:
    p_18 = p_a
    T_18 = isentropexp_T(T_t21, n_nozzle, p_t21, p_a, k_a)
    v_18 = sqrt(2*c_p_a*(T_t21-T_18))
    F_bypass = m_dot_bypass * (v_18 - v_fs)
    v_19eff = v_18
    print("Bypasss nozzle is not choked")

### Overall Performance ###
F_N = F_core + F_bypass
TSFC = m_dot_f/F_N *1000000 #g/kN.S



#Propulsive efficiency#
n_prop = (m_dot_4*(v_9eff-v_fs)+m_dot_bypass*(v_19eff-v_fs))*v_fs/(((0.5*m_dot_4)*(v_9eff**2-v_fs**2)+(0.5*m_dot_bypass)*(v_19eff**2-v_fs**2)))

#Thermodynamic efficiency #
P_gg = m_dot_4*c_p_g*T_gg*(1-(p_t/p_gg)**((k_g-1)/k_g))-0.5*m_dot_4*v_fs**2
n_thdy = P_gg / (m_dot_core*c_p_g*(T_combexit-T_t3))

#Gas generation efficiency #
n_gas = ((0.5*m_dot_4)*(v_9eff**2-v_fs**2)+(0.5*m_dot_bypass)*(v_19eff**2-v_fs**2))/P_gg

#Thermal efficiency #
n_th = ((0.5 * m_dot_4)*(v_9eff**2 - v_fs**2) + (0.5 * m_dot_bypass)*(v_19eff**2 - v_fs**2))/(m_dot_f * LHV*1000000)


#Overall efficiency#
n_tot = (v_fs*F_N)/(LHV*m_dot_f*1000000)
n_totcheck = n_th*n_prop
OPR = PR_LPC*PR_HPC

data = {'Temperature [K] ':[T_t, T_t21, T_t21,T_18 ,T_t25,T_t3, T_combexit, T_t45, T_gg, T_t5, T_t5, T_8],
        'Pressure [Pa]':[p_t2, p_t21, p_t21, p_18, p_t25, p_t3, p_t4, p_t45, p_gg, p_t5, p_t5, p_8 ],
        'Mass flow [kg/s]':[m_dot, m_dot_core, m_dot_bypass, m_dot_bypass, m_dot_core, m_dot_core, m_dot_4, m_dot_4, m_dot_4, m_dot_4, m_dot_4, m_dot_4]}
df = pd.DataFrame(data, index=['2', '21', '13', '18', '25', '3', '4', '45', 'gg', '5', '7', '8'])
print()
print("Part One-cycle calculation")
print()
print(df)
print()
print("At cruise conditions:")
print("OPR = ", OPR)
print("Thrust = ", F_N, "N")
print("TSFC =", TSFC, "g/kN.S")
print()
print("Efficiencies:")
print("Thermodynamic efficiendcy = ", n_thdy)
print("Gas generation efficiendcy = ", n_gas)
print("Thermal efficiendcy = ", n_th)
print("Propulsive efficiendcy = ", n_prop)
print("------------------------------------------")
print("Total efficiency = ", n_tot)
print("Combined efficiencies total check =", n_th*n_prop)
print("Combined efficiencies thermal check =", n_gas*n_thdy*n_comb)
print("Combined efficiencies prop check =", n_gas*n_thdy*n_comb*n_prop)









###### T-s plot ######

air = pm.get('ig.air')
pm.config['unit_pressure'] = 'Pa'
pm.config['unit_energy'] = 'J'


## ideal cycle ##
PR_tot = PR_fan*PR_HPC*PR_LPC*PR_intake
p1 = p_t
T1 = T_t
s1 = air.s(T1, p1)
p2 = p1*PR_tot
T2 = air.T_s(s=s1, p=p2)
T3 = T_combexit
p3 = p2
s3 = air.s(T3, p3)
s4 = s3
p4 = p1
T4 = air.T_s(s=s4, p=p4)

plt.figure()
plt.plot([s1, s1], [T1, T2], 'r', linewidth=1)
T = np.linspace(T2, T3, 20)
plt.plot(air.s(T=T, p=p2), T, 'r', linewidth=1)
plt.plot([s3, s3], [T3, T4], 'r', linewidth=1)
T = np.linspace(T1, T4, 20)
plt.plot(air.s(T=T, p=p1), T, 'r--', linewidth=1)

#points
plt.plot(s1, T2, 'ro', markersize=3)
plt.plot(s3, T3, 'ro', markersize=3)
plt.plot(s4, T4, 'ro', markersize=3)
plt.annotate("3'", (s1, T2), textcoords="offset points", xytext=(-5, 0), ha='right')
plt.annotate("4'", (s3, T3), textcoords="offset points", xytext=(-5, 0), ha='right')
plt.annotate("8'", (s4, T4), textcoords="offset points", xytext=(5, 5), ha='left')

## real cycle ##
st = air.s(T_t, p_t)             #c_p_a*np.log(T_t/T_a) - R*np.log(p_t/p_a)
s2 = air.s(T_t, p_t2)            #c_p_a*np.log(T_t/T_t) - R*np.log(p_t2/p_t) + st
s21 = air.s(T_t21, p_t21)        #c_p_a*np.log(T_t21/T_t) - R*np.log(p_t21/p_t2) + s2
s25 = air.s(T_t25, p_t25)        #c_p_a*np.log(T_t25/T_t21) - R*np.log(p_t25/p_t21) + s21
s3 = air.s(T_t3, p_t3)           #c_p_a*np.log(T_t3/T_t25) - R*np.log(p_t3/p_t25) + s25
s4 = air.s(T_combexit, p_t4)     #c_p_g*np.log(T_combexit/T_t3) - R*np.log(p_t4/p_t3) + s3
s45 = air.s(T_t45, p_t45)        #c_p_g*np.log(T_t45/T_combexit) - R*np.log(p_t45/p_t4) + s4
s5 = air.s(T_t5, p_t5)           #c_p_g*np.log(T_t5/T_t45) - R*np.log(p_t5/p_t45) + s45
sgg = air.s(T_gg, p_gg) + 5.795
s7 = s5
s8 = air.s(T_8, p_8)             #c_p_g*np.log(T_8/T_t5) - R*np.log(p_8/p_t5) + s7

plt.plot([st, s2], [T_t, T_t], 'b', linewidth=1)
plt.plot([s21, s2], [T_t21, T_t], 'b', linewidth=1)
plt.plot([s25, s21], [T_t25, T_t21], 'b', linewidth=1)
plt.plot([s3, s25], [T_t3, T_t25], 'b', linewidth=1)
T = np.linspace(T_t3, T_combexit, 20)
p = np.linspace(p_t3, p_t4, 20)
plt.plot(air.s(T=T, p=p), T, 'b', linewidth=1)
plt.plot([s45, s4], [T_t45, T_combexit], 'b', linewidth=1)
plt.plot([s5, s45], [T_t5, T_t45], 'b', linewidth=1)
plt.plot([s8, s7], [T_8, T_t5], 'b', linewidth=1)
T = np.linspace(T_t, T_8, 20)
p = np.linspace(p_t, p_8, 20)
plt.plot(air.s(T=T, p=p), T, 'b--', linewidth=1)

#points
plt.plot(st, T_t, 'bo', markersize=3)
plt.annotate(0, (st, T_t), textcoords="offset points", xytext=(-5, -5), ha='right')
plt.plot(s21, T_t21, 'bo', markersize=3)
plt.annotate(21, (s21, T_t21), textcoords="offset points", xytext=(-5, 0), ha='right')
plt.plot(s25, T_t25, 'bo', markersize=3)
plt.annotate(25, (s25, T_t25), textcoords="offset points", xytext=(5, 0), ha='left')
plt.plot(s3, T_t3, 'bo', markersize=3)
plt.annotate(3, (s3, T_t3), textcoords="offset points", xytext=(-5, 0), ha='right')
plt.plot(s4, T_combexit, 'bo', markersize=3)
plt.annotate(4, (s4, T_combexit), textcoords="offset points", xytext=(5, 0), ha='left')
plt.plot(s45, T_t45, 'bo', markersize=3)
plt.annotate(45, (s45, T_t45), textcoords="offset points", xytext=(5, 0), ha='left')
plt.plot(sgg, T_gg, 'bo', markersize=3)
plt.annotate("gg", (sgg, T_gg), textcoords="offset points", xytext=(5, 0), ha='left')
plt.plot(s5, T_t5, 'bo', markersize=3)
plt.annotate("5,7", (s5, T_t5), textcoords="offset points", xytext=(5, 0), ha='left')
plt.plot(s8, T_8, 'bo', markersize=3)
plt.annotate(8, (s8, T_8), textcoords="offset points", xytext=(5, 0), ha='left')


## plot ##
ax = plt.gca()
ax.set_ylim([200, 1500])
plt.xlabel('Entropy, s (J/kg/K)')
plt.ylabel('Temperature, T (K)')
plt.grid('on')
plt.title('LEAP-1A: T-s Diagram')
plt.show()
