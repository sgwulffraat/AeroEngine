from math import *
import numpy as np
import pandas as pd

### Inputs ####
M = 0.78
h = 0               #m
m_dot = 23.81       #kg/s
# m_fuel = 0.4267     #kg/s
T_combexit = 1150    #K
PR_comp = 5.5
n_is_C = 0.83
n_is_T = 0.82
n_mech = 0.96
n_comb = 1
PR_comb = 0.96
n_nozzle = 1
#T_a = 288           #K
#p_a = 100000        #Pa
T_a = 218.808       #K @ h = 35000ft
p_a = 23842.3       #pa @ h = 35000ft
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
     #T_combexit = T_t2 + (m_fuel * n_comb * LHV*10**6)/(m_dot*c_p_g)
m_fuel = (m_dot * c_p_g * (T_combexit - T_t2)) / (n_comb * LHV *10**6)
p_t4 = PR_comb * p_t2
m_dot_4 = m_fuel + m_dot

### Turbine conditions ###
T_t5 = T_combexit - (W_req_C / (n_mech * m_dot_4 * c_p_g))
p_t5 = p_t4*(1 - 1/n_is_T*(1- T_t5/T_combexit)) ** (k_g/(k_g - 1))



### Core nozzle conditions ###
e_c = critPR(n_nozzle, k_g)
PR_noz = p_t5/p_a
v_fs = M * sqrt(k_a*R*T_a)
if PR_noz > e_c:
    A_ratio = (p_t5 / p_t4) ** ((2 * k_g - n_is_T * (k_g - 1)) / (2 * k_g))
    T_8 = T_t5 * (2/(k_g+1))
    p_8 = p_t5/e_c
    v_8 = sqrt(k_g*R*T_8)
    rho_8 = p_8/(R*T_8)
    A_8 = m_dot_4/(rho_8*v_8)
    F_N = m_dot_4*(v_8-v_fs) + A_8 * (p_8-p_a)
    F_gross = m_dot_4*v_8 + A_8 * (p_8-p_a)
    v_9eff = F_N / m_dot_4 + v_fs
    nzchoked = True
    print("Core nozzle is choked")
else:
    p_8 = p_a
    T_78 = isentropexp_T(T_t5, n_nozzle, p_t5, p_a, k_g)
    if T_78 > 0:
        v_8 = sqrt(2 * c_p_g * (T_78))
        v_9eff = v_8
        F_N = m_dot_4 * (v_8 - v_fs)
        F_gross = m_dot_4 * v_8
    else:
        F_N = 0
        F_gross = 0
    nzchoked = False
    print("Core nozzle is not choked")

data = {'Temperature [K] ':[T_t, T_t2, T_combexit,T_t5 ,T_8],
        'Pressure [Pa]':[p_t, p_t2, p_t4, p_t5, p_8],
        'Mass flow [kg/s]':[m_dot, m_dot, m_dot_4, m_dot_4, m_dot_4]}
df = pd.DataFrame(data, index=['Inlet', '3', '4', '5', '8'])

print(T_combexit)
print(F_gross)
print(A_ratio)

print(df)

#Design parameters
RPM = 13800
omega = RPM*(2*pi)/(60)
psi = 0.25
phi = 0.35
stages = 3
r_c = 0.5

power = W_req_C

w =  power/(m_dot*stages)
U_s = sqrt((w)/(psi))
lambda_s = (2*w)/(U_s**2)
v_m = phi*U_s
r_out = U_s/omega

#velocity triangles

A = np.array([[-phi, 0, 0, phi], [-phi, 0, 0, 0], [-1, 0, 1, 0], [0, 1, 0 ,-1]])
b = np.array([[psi - 1], [r_c - 1 + psi/2], [-1/phi], [1/phi]])
x = np.linalg.solve(A, b)
xr = np.arctan(x)


#Calculate Dimensional Velocities
Vm1 = phi * U_s
V1 = Vm1 / np.cos(xr[0])
W1 = Vm1 / np.cos(xr[2])
Vm2 = Vm1
V2 = Vm2 / np.cos(xr[1])
W2 = Vm2 / np.cos(xr[3])
Vt1 = Vm1 * np.tan(xr[0])
Vt2 = Vm2 * np.tan(xr[1])
dVt = Vt2 - Vt1
w = U_s * dVt





#interstage thermodynamics properties

#after first stage
n_is_Cs = 0.8536
W_req_Cs = W_req_C
T_t01 = (W_req_Cs)/(stages*m_dot*c_p_a) + T_t
p_t01 = p_t*((n_is_Cs*(T_t01-T_t)/T_t)+1)**(k_a/(k_a - 1))
Ts_01 = T_t01 - (k_a-1)/2*(Vm1*Vm1/k_a/R)
M_01 = Vm1 / np.sqrt(k_a*R*Ts_01)
ps_01 = p_t01 * (1+((k_a-1)/2)*M_01*M_01)**(-k_a/(k_a-1))
rho_01 = ps_01 / R / Ts_01
pr_ratio_stage = p_t01/p_t

#after second stage
T_t02 = (W_req_Cs)/(stages*m_dot*c_p_a) + T_t01
p_t02 = p_t01*((n_is_Cs*(T_t02-T_t01)/T_t01)+1)**(k_a/(k_a - 1))
Ts_02 = T_t02 - (k_a-1)/2*(Vm1*Vm1/k_a/R)
M_02 = Vm1 / np.sqrt(k_a*R*Ts_02)
ps_02 = p_t02 * (1+((k_a-1)/2)*M_02*M_02)**(-k_a/(k_a-1))
rho_02 = ps_02 / R / Ts_02
pr_ratio_stage2 = p_t02/p_t01

#after third stage
T_t03 = (W_req_Cs)/(stages*m_dot*c_p_a) + T_t02
p_t03 = p_t02*((n_is_Cs*(T_t03-T_t02)/T_t02)+1)**(k_a/(k_a - 1))
Ts_03 = T_t03 - (k_a-1)/2*(Vm1*Vm1/k_a/R)
M_03 = Vm1 / np.sqrt(k_a*R*Ts_03)
ps_03 = p_t03 * (1+((k_a-1)/2)*M_03*M_03)**(-k_a/(k_a-1))
rho_03 = ps_03 / (R * Ts_03)
pr_ratio_stage3 = p_t03/p_t02

# print('de waarde voor t01 en p01 =', T_t01, p_t01)
# print('de waarde voor t02 en p02 =', T_t02, p_t02)
# print('de waarde voor t03 en p03 =', T_t03, p_t03)
# print('de waarde voor de pressure ratio tussen de stages= ', pr_ratio_stage*pr_ratio_stage2*pr_ratio_stage3)

#Area calculations first stage
#Ts_ref = T_t - (k_a-1)/2*(Vm1*Vm1/k_a/R)
#M_ref = Vm1 / np.sqrt(k_a*R*Ts_ref)
#print('Mref=', M_ref)
Ts_00 = T_t*(1+(k_a-1)/2 * M**2)**-1
ps_00 = p_t*(1+(k_a-1)/2 * M**2)**(-k_a/(k_a-1))
rho00 = ps_00/(R*Ts_00)
A00 = m_dot/(rho00*M*sqrt(R*k_a*Ts_00))
r00 = sqrt((np.pi*r_out**2 - A00)/np.pi)

A1 = m_dot/(rho_01*Vm1)
r1 = sqrt((np.pi*r_out**2 - A1)/np.pi)

A2 = m_dot/(rho_02*Vm1)
r2 = sqrt((np.pi*r_out**2 - A2)/np.pi)

A3 = m_dot/(rho_03*Vm1)
r3 = sqrt((np.pi*r_out**2 - A3)/np.pi)

print("r:")
print(r_out)
print(r00)
print(r1)
print(r2)
print(r3)
