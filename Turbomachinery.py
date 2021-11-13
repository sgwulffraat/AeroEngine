from math import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import pyromat as pyro
import matplotlib.pyplot as plt

### Inputs ####
M = 0.78
h = 0               #m
m_dot = 23.81       #kg/s
#m_fuel = 0.4267     #kg/s
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
        T_8 = T_t5 + T_78
        v_8 = sqrt(2 * c_p_g * (T_8))
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
#print(A_ratio)

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

V1vector = [float(V1*np.sin(xr[0])),float(-V1*np.cos(xr[0]))]
W1vector = [float(W1*np.sin(xr[2])),float(-W1*np.cos(xr[2]))]

ax = plt.axes()
plt.axvline(0,linestyle='--',color='grey')
ax.arrow(0,0,V1vector[0],V1vector[1],length_includes_head=True,width=3,head_width=8,head_length=20,fc='red',ec='red',label='V1')
ax.arrow(0,0,W1vector[0],W1vector[1],length_includes_head=True,width=3,head_width=8,head_length=20,fc='blue',ec='blue',label='W1')
ax.arrow(W1vector[0],W1vector[1],(V1vector[0]-W1vector[0]),0,length_includes_head=True,width=3,head_width=8,head_length=20,fc='green',ec='green',label='U')
ax.text(7,-30,r'$\alpha_1$')
ax.text(-27,-30,r'$\beta_1$')
plt.title('Velocity triangle in front of the rotor')
plt.legend()
plt.ylabel('Axial velocity [m/s]')
plt.xlabel('Tangential velocity [m/s]')
plt.show()


V2vector = [float(V2*np.sin(xr[1])),float(-V2*np.cos(xr[1]))]
W2vector = [float(W2*np.sin(xr[3])),float(-W2*np.cos(xr[3]))]

ax = plt.axes()
ax.arrow(0,0,V2vector[0],V2vector[1],length_includes_head=True,width=3,head_width=8,head_length=20,fc='red',ec='red',label='V2')
ax.arrow(0,0,W2vector[0],W2vector[1],length_includes_head=True,width=3,head_width=8,head_length=20,fc='blue',ec='blue',label='W2')
ax.arrow(W2vector[0],W2vector[1],(V2vector[0]-W2vector[0]),0,length_includes_head=True,width=3,head_width=8,head_length=20,fc='green',ec='green',label='U')
ax.text(7,-30,r'$\alpha_2$')
ax.text(-23,-30,r'$\beta_2$')
plt.axvline(0,linestyle='--',color='grey')
plt.title('Velocity triangle after the rotor')
plt.legend()
plt.ylabel('Axial velocity [m/s]')
plt.xlabel('Tangential velocity [m/s]')
plt.show()



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
Ts_ref = T_t - (k_a-1)/2*(Vm1*Vm1/k_a/R)
M_ref = Vm1 / np.sqrt(k_a*R*Ts_ref)
#print('Mref=', M_ref)
#Ts_00 = T_t*(1+(k_a-1)/2 * M**2)**-1
ps_00 = p_t*(1+(k_a-1)/2 * M_ref**2)**(-k_a/(k_a-1))
rho00 = ps_00/(R*Ts_ref)
A00 = m_dot/(rho00*M_ref*sqrt(R*k_a*Ts_ref))
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
""""
#hs diagram for 3 stages
air = pyro.get('ig.air')
p1 = p_t
T1 = T_t
p2 = p_t01
T2 = T_t01
p3 = p_t02
T3 = T_t02
p4 = p_t03
T4 = T_t03
pyro.config['unit_temperature'] = 'K'
pyro.config['unit_pressure'] = 'Pa'

s1 = air.s(T1,p1)
s2 = air.s(T2,p2)
s3 = air.s(T3,p3)
s4 = air.s(T4,p4)

h1 = air.h(T1,p1) + 60
h2 = air.h(T2,p2) + 60
h3 = air.h(T3,p3) + 60
h4 = air.h(T4,p4) + 60

plt.plot([s1,s2],[h1,h2],'r',linewidth=1.5)
plt.plot([s2,s3],[h2,h3],'b',linewidth=1.5)
plt.plot([s3,s4],[h3,h4],'y',linewidth=1.5)
plt.annotate("Stage 1", (s1, h1), textcoords="offset points", xytext=(55, 0), ha='right')
plt.annotate("Stage 2", (s2, h2), textcoords="offset points", xytext=(-5, 0), ha='right')
plt.annotate("Stage 3", (s3, h3), textcoords="offset points", xytext=(-5, 0), ha='right')



plt.axis()

plt.xlabel('Entropy, s (kJ/kg/K)')
plt.ylabel('Enthalpy, h (J/kg)')
plt.grid('on')

plt.title('Interstage Compressor h-s diagram ')

plt.show()
"""
#Plot aero-thermal flow properties
Stages = [0,1,2,3]
Rotor = [1,2,3]

Tt_ratio00 = T_t / T_t
Tt_ratio01 = T_t01 / T_t
Tt_ratio02 = T_t02 / T_t
Tt_ratio03 = T_t03 / T_t
Tt_ratio = [float(Tt_ratio00),float(Tt_ratio01),float(Tt_ratio02),float(Tt_ratio03)]


Pt_ratio00 = p_t / p_t
Pt_ratio01 = p_t01 / p_t
Pt_ratio02 = p_t02 / p_t
Pt_ratio03 = p_t03 / p_t
Pt_ratio = [float(Pt_ratio00),float(Pt_ratio01),float(Pt_ratio02),float(Pt_ratio03)]

P_ratio00 = ps_00 / p_t
P_ratio01 = ps_01 / p_t
P_ratio02 = ps_02 / p_t
P_ratio03 = ps_03 / p_t
P_ratio = [float(P_ratio00),float(P_ratio01),float(P_ratio02),float(P_ratio03)]

beta1 = T_t01 / T_t
beta2 = T_t02 / T_t01
beta3 = T_t03 / T_t02
beta = [beta1,beta2,beta3]

fig, axs = plt.subplots(2,2)
fig.suptitle('Plots of aero-thermal flow properties')
axs[0, 0].plot(Stages,Tt_ratio)
plt.setp(axs[0,0],xlabel="# stage",ylabel=r'$\frac{T_{t}}{T_{t,0}}$')
axs[0,0].set_title("Total temperature over reference total temperature")
axs[0, 1].plot(Stages,Pt_ratio)
axs[0,1].set_title("Total pressure over reference total pressure")
plt.setp(axs[0,1],xlabel="# stage",ylabel=r'$\frac{p_{t}}{p_{t,0}}$')
axs[1, 0].plot(Stages,P_ratio)
axs[1,0].set_title("Static pressure over reference total pressure")
plt.setp(axs[1,0],xlabel="# stage",ylabel=r'$\frac{p_{s}}{p_{t,0}}$')
axs[1, 1].plot(Rotor,beta)
axs[1,1].set_title(r"$\beta_{tt}$ over rotor stages")
plt.setp(axs[1,1],xlabel="# stage",ylabel=r'$\beta_{tt}$')
plt.show()

## Meriodial Gas Path ##
## Geometry ##
wR1 = (r_out-r1)/1.5
wR2 = (r_out-r2)/1.5
wR3 = (r_out-r3)/1.5
ds = 0.01

## Stage 1 ##
xR1 = [0, 0, wR1, wR1, 0]
yR1 = [r1-(r2-r1)/2, r_out, r_out, r1, r1-(r2-r1)/2]
xS1 = [wR1 + ds, wR1 + ds, 2 * wR1 + ds, 2 * wR1 + ds, wR1 + ds]
yS1 = [r1, r_out, r_out, r1+(r2-r1)/2, r1]

## Stage 2 ##
xR2 = [2*wR1 + 3* ds, 2*wR1 + 3*ds, 2*wR1 + wR2 + 3*ds, 2*wR1 + wR2 + 3*ds, 2*wR1 + 3*ds]
yR2 = [r2-(r2-r1)/3, r_out, r_out, r2, r2-(r2-r1)/3]
xS2 = [2*wR1 + wR2 + 4*ds, 2*wR1 + wR2 + 4*ds, 2*wR1 + 2*wR2 + 4*ds, 2*wR1 + 2*wR2 + 4*ds, 2*wR1 + wR2 + 4*ds]
yS2 = [r2, r_out, r_out, r2+(r2-r1)/3, r2]

## Stage 3 ##
xR3 = [2*wR1 + 2*wR2 + 6*ds, 2*wR1 + 2*wR2 + 6*ds, 2*wR1 + 2*wR2 + wR3 + 6*ds, 2*wR1 + 2*wR2 + wR3 + 6*ds, 2*wR1 + 2*wR2 + 6*ds]
yR3 = [r3-(r3-r2)/3, r_out, r_out, r3, r3-(r3-r2)/3]
xS3 = [2*wR1 + 2*wR2 + wR3 + 7*ds, 2*wR1 + 2*wR2 + wR3 + 7*ds, 2*wR1 + 2*wR2 + 2*wR3 + 7*ds, 2*wR1 + 2*wR2 + 2*wR3 + 7*ds, 2*wR1 + 2*wR2 + wR3 + 7*ds]
yS3 = [r3, r_out, r_out, r3+(r3-r2)/3, r3]

plt.figure(0)
# plt.axhline(y = r_out+0.005, color = 'black', linestyle = '-')
# plt.axhline(y = r_out+0.02, color = 'black', linestyle = '-')
# plt.axhline(y = 0, color = 'black', linestyle = '--')
plt.plot(xR1, yR1, marker = 'o', color = 'tab:red', label = 'Rotor')
plt.plot(xS1, yS1, marker = 'o', color = 'tab:blue', label = 'Stator')
plt.plot(xR2, yR2, marker = 'o', color = 'tab:red')
plt.plot(xS2, yS2, marker = 'o', color = 'tab:blue')
plt.plot(xR3, yR3, marker = 'o', color = 'tab:red')
plt.plot(xS3, yS3, marker = 'o', color = 'tab:blue')
plt.xlabel('Axial length [m]')
plt.ylabel('Radius [m]')
plt.ylim(-0.05, r_out + 0.1)
plt.xlim(-0.05, r_out + 0.1)
plt.grid()
plt.legend()
plt.show()



