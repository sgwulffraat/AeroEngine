import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### Read file ###
df = pd.read_excel (r'220208_Scenarios.xlsx', skiprows=5, header=None)
df.columns = "Time", "Flown_km_inc", "Fuel_use_inc", "Flown_km_cnst", "Fuel_use_cnst", "Flown_km_decr", "Fuel_use_decr", "Flown_km_bc", "Fuel_use_bc"

### Part A ###

Ebc = (df.Fuel_use_bc * 1000/167.3110 * 12 * 44.01)/1E12      #Tg
Einc = (df.Fuel_use_inc * 1000/167.3110 * 12 * 44.01)/1E12    #Tg
Ecnst = (df.Fuel_use_cnst * 1000/167.3110 * 12 * 44.01)/1E12  #Tg
Edecr = (df.Fuel_use_decr * 1000/167.3110 * 12 * 44.01)/1E12  #Tg

plt.figure(0)
plt.plot(df.Time, Einc, label = 'Increasing')
plt.plot(df.Time, Edecr, label = 'Decreasing')
plt.plot(df.Time, Ecnst, label = 'Constant')
plt.plot(df.Time, Ebc, label = 'Base Case')
plt.xlim(1940, 2120)
plt.ylim(0, max(Einc)*1.05)
plt.title("CO2 emissions per year")
plt.xlabel('Time (years)')
plt.ylabel('CO2 Emissions (Tg/year)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'small')
plt.show()

### Part B ###

#S&S Method
def Gc(t):
    Gc = (0.067 + 0.1135*np.exp(-t/313.8) + 0.152*np.exp(-t/79.8) + 0.0970*np.exp(-t/18.8) + 0.041*np.exp(-t/1.7))/1000
    return Gc

def convsolver(E):
    dCO2 = []
    for i in range(0,len(df.Time)):
        if i == 0:
            G = Gc(df.Time[i]-1940)
            dCO2.append(G * E[i])
        else:
            x = []
            for k in range(0,i+1):
                x.append(Gc(df.Time[i]-df.Time[k])*E[k])
            G = np.sum(x)
            dCO2.append(G)
    return dCO2

#B&R Method
def f(t):
    f = (0.217 + 0.259*np.exp(-t/172.9) + 0.338*np.exp(-t/18.51) + 0.186*np.exp(-t/1.186))/2.123
    return f

def BRconvsolver(E):
    dCO2 = []
    for i in range(0,len(df.Time)):
        if i == 0:
            F = f(df.Time[i]-1940)
            dCO2.append(F * E[i]/1000)
        else:
            y = []
            for k in range(0,i+1):
                y.append(f(df.Time[i]-df.Time[k])*E[k]/1000)
            F = np.sum(y)
            dCO2.append(F)
    return dCO2

SS_CO2_inc = convsolver(Einc)
SS_CO2_decr = convsolver(Edecr)
SS_CO2_cnst = convsolver(Ecnst)
SS_CO2_bc = convsolver(Ebc)

BR_CO2_inc = BRconvsolver(Einc)
BR_CO2_decr = BRconvsolver(Edecr)
BR_CO2_cnst = BRconvsolver(Ecnst)
BR_CO2_bc = BRconvsolver(Ebc)

plt.figure(1)
plt.plot(df.Time, SS_CO2_inc, label= 'S&S Increasing', color = 'tab:blue', linewidth=1.5)
plt.plot(df.Time, SS_CO2_decr, label= 'S&S Decreasing', color = 'tab:orange', linewidth=1.5)
plt.plot(df.Time, SS_CO2_cnst, label= 'S&S Constant', color = 'tab:green', linewidth=1.5)
plt.plot(df.Time, SS_CO2_bc, label= 'S&S Base Case', color = 'tab:red', linewidth=1.5)
plt.plot(df.Time, BR_CO2_inc, label= 'B&R Increasing', color = 'tab:blue', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_CO2_decr, label= 'B&R Decreasing', color = 'tab:orange', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_CO2_cnst, label= 'B&R Constant', color = 'tab:green', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_CO2_bc, label= 'B&R Base Case', color = 'tab:red', linestyle='dashed', linewidth=1.5)
plt.xlim(1940, 2120)
plt.ylim(0, max(SS_CO2_inc)*1.05)
plt.title("CO2 concentration in atmosphere")
plt.xlabel('Time (years)')
plt.ylabel('Delta CO2 (ppmv)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'medium')
plt.show()


### Part C ###
## Question A ##
C0 = np.zeros(len(SS_CO2_inc)) + 380.

# IPCC Method
def RFcalc(dC):
    C = dC + C0
    RF = 5.35 * np.log(C/C0)
    return RF

SS_RF_inc = RFcalc(SS_CO2_inc)
SS_RF_decr = RFcalc(SS_CO2_decr)
SS_RF_cnst = RFcalc(SS_CO2_cnst)
SS_RF_bc = RFcalc(SS_CO2_bc)

BR_RF_inc = RFcalc(BR_CO2_inc)
BR_RF_decr = RFcalc(BR_CO2_decr)
BR_RF_cnst = RFcalc(BR_CO2_cnst)
BR_RF_bc = RFcalc(BR_CO2_bc)

plt.figure(2)
plt.plot(df.Time, SS_RF_inc, label= 'S&S Increasing', color = 'tab:blue', linewidth=1.5)
plt.plot(df.Time, SS_RF_decr, label= 'S&S Decreasing', color = 'tab:orange', linewidth=1.5)
plt.plot(df.Time, SS_RF_cnst, label= 'S&S Constant', color = 'tab:green', linewidth=1.5)
plt.plot(df.Time, SS_RF_bc, label= 'S&S Base Case', color = 'tab:red', linewidth=1.5)
plt.plot(df.Time, BR_RF_inc, label= 'B&R Increasing', color = 'tab:blue', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_RF_decr, label= 'B&R Decreasing', color = 'tab:orange', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_RF_cnst, label= 'B&R Constant', color = 'tab:green', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_RF_bc, label= 'B&R Base Case', color = 'tab:red', linestyle='dashed', linewidth=1.5)
plt.xlim(1940, 2120)
plt.ylim(0, max(SS_RF_inc)*1.05)
plt.title("Radiative Forcing due to CO2 emissions")
plt.xlabel('Time (years)')
plt.ylabel('Radiative Forcing (W/m^2)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'medium')


## Question B ##
def RFcalccontrail(km):
    RF = (1.82E-9  * km)/1000
    return RF

contrailRF_inc = RFcalccontrail(df.Flown_km_inc)
contrailRF_decr = RFcalccontrail(df.Flown_km_decr)
contrailRF_cnst = RFcalccontrail(df.Flown_km_cnst)
contrailRF_bc = RFcalccontrail(df.Flown_km_bc)

plt.figure(3)
plt.plot(df.Time, contrailRF_inc, label= 'S&S Increasing', color = 'tab:blue', linewidth=1.5)
plt.plot(df.Time, contrailRF_decr, label= 'S&S Decreasing', color = 'tab:orange', linewidth=1.5)
plt.plot(df.Time, contrailRF_cnst, label= 'S&S Constant', color = 'tab:green', linewidth=1.5)
plt.plot(df.Time, contrailRF_bc, label= 'S&S Base Case', color = 'tab:red', linewidth=1.5)
plt.xlim(1940, 2120)
plt.ylim(0, max(contrailRF_inc)*1.05)
plt.title("Radiative Forcing due to contrails")
plt.xlabel('Time (years)')
plt.ylabel('Radiative Forcing (W/m^2)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'medium')
plt.show()

### Part D ###

# S&S Temperature calculations
def Gt(t):
    Gt = 2.246/36.8 * np.exp(-t/36.8)
    return Gt

def Tconv_SS(RF):
    dT = []
    beta = 0.45
    for i in range(0,len(df.Time)):
        if i == 0:
            G = Gt(df.Time[i]-1940)
            dT.append(beta * G * RF[i])
        else:
            x = []
            for k in range(0,i+1):
                x.append(Gt(df.Time[i]-df.Time[k])*RF[k]*beta)
            G = np.sum(x)
            dT.append(G)
    return dT

# B&R Temperature calculations
def deltaT(t):
    deltaT = (0.631/8.4 * np.exp(-t/8.4) + 0.429/409.5 *np.exp(-t/409.5))
    return deltaT

def Tconv_BR(RF):
    dT = []
    for i in range(0,len(df.Time)):
        if i == 0:
            F = deltaT(df.Time[i]-1940)
            dT.append(F * RF[i])
        else:
            y = []
            for k in range(0,i+1):
                y.append(deltaT(df.Time[i]-df.Time[k])*RF[k])
            F = np.sum(y)
            dT.append(F)
    return dT

#Delta T due to CO2
SS_dT_inc = Tconv_SS(SS_RF_inc)
SS_dT_decr = Tconv_SS(SS_RF_decr)
SS_dT_cnst = Tconv_SS(SS_RF_cnst)
SS_dT_bc = Tconv_SS(SS_RF_bc)

BR_dT_inc = Tconv_BR(BR_RF_inc)
BR_dT_decr = Tconv_BR(BR_RF_decr)
BR_dT_cnst = Tconv_BR(BR_RF_cnst)
BR_dT_bc = Tconv_BR(BR_RF_bc)

plt.figure(4)
plt.plot(df.Time, SS_dT_inc, label= 'S&S Increasing', color = 'tab:blue', linewidth=1.5)
plt.plot(df.Time, SS_dT_decr, label= 'S&S Decreasing', color = 'tab:orange', linewidth=1.5)
plt.plot(df.Time, SS_dT_cnst, label= 'S&S Constant', color = 'tab:green', linewidth=1.5)
plt.plot(df.Time, SS_dT_bc, label= 'S&S Base Case', color = 'tab:red', linewidth=1.5)
plt.plot(df.Time, BR_dT_inc, label= 'B&R Increasing', color = 'tab:blue', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_dT_decr, label= 'B&R Decreasing', color = 'tab:orange', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_dT_cnst, label= 'B&R Constant', color = 'tab:green', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BR_dT_bc, label= 'B&R Base Case', color = 'tab:red', linestyle='dashed', linewidth=1.5)
plt.xlim(1940, 2120)
plt.ylim(0, max(SS_dT_inc)*1.05)
plt.title("Temperature change over time due to RF (CO2)")
plt.xlabel('Time (years)')
plt.ylabel('Delta T (K)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'medium')

#delta T due to contrails
SSc_dT_inc = Tconv_SS(contrailRF_inc)
SSc_dT_decr = Tconv_SS(contrailRF_decr)
SSc_dT_cnst = Tconv_SS(contrailRF_cnst)
SSc_dT_bc = Tconv_SS(contrailRF_bc)

BRc_dT_inc = Tconv_BR(contrailRF_inc)
BRc_dT_decr = Tconv_BR(contrailRF_decr)
BRc_dT_cnst = Tconv_BR(contrailRF_cnst)
BRc_dT_bc = Tconv_BR(contrailRF_bc)

plt.figure(5)
plt.plot(df.Time, SSc_dT_inc, label= 'S&S Increasing', color = 'tab:blue', linewidth=1.5)
plt.plot(df.Time, SSc_dT_decr, label= 'S&S Decreasing', color = 'tab:orange', linewidth=1.5)
plt.plot(df.Time, SSc_dT_cnst, label= 'S&S Constant', color = 'tab:green', linewidth=1.5)
plt.plot(df.Time, SSc_dT_bc, label= 'S&S Base Case', color = 'tab:red', linewidth=1.5)
plt.plot(df.Time, BRc_dT_inc, label= 'B&R Increasing', color = 'tab:blue', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BRc_dT_decr, label= 'B&R Decreasing', color = 'tab:orange', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BRc_dT_cnst, label= 'B&R Constant', color = 'tab:green', linestyle='dashed', linewidth=1.5)
plt.plot(df.Time, BRc_dT_bc, label= 'B&R Base Case', color = 'tab:red', linestyle='dashed', linewidth=1.5)
plt.xlim(1940, 2120)
plt.ylim(0, max(SSc_dT_inc)*1.05)
plt.title("Temperature change over time due RF (contrails)")
plt.xlabel('Time (years)')
plt.ylabel('Delta T (K)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'medium')
plt.show()

### Part E ###
#ERF
cERF_inc = contrailRF_inc*0.42
cERF_decr = contrailRF_decr*0.42
cERF_cnst = contrailRF_cnst*0.42
cERF_bc = contrailRF_bc*0.42


#GWP function
def GWP(cRFx, cRFbc, Th, RFx, RFbc):
    GWP = (sum(cRFx[80:81+Th])-sum(cRFbc[80:81+Th]))/(sum(RFx[80:81+Th])-sum(RFbc[80:81+Th]))
    return GWP

def GTP(dTcx, dTcbc, Th, dTx, dTbc):
    GTP = (dTcx[80+Th]-dTcbc[80+Th])/(dTx[80+Th]-dTbc[80+Th])
    return GTP

#ERFdeltaT
SSc2_dT_inc = Tconv_SS(cERF_inc)
SSc2_dT_decr = Tconv_SS(cERF_decr)
SSc2_dT_cnst = Tconv_SS(cERF_cnst)
SSc2_dT_bc = Tconv_SS(cERF_bc)
BRc2_dT_inc = Tconv_BR(cERF_inc)
BRc2_dT_decr = Tconv_BR(cERF_decr)
BRc2_dT_cnst = Tconv_BR(cERF_cnst)
BRc2_dT_bc = Tconv_BR(cERF_bc)

#GWP20
GWPinc20SS = GWP(cERF_inc, cERF_bc, 20, SS_RF_inc, SS_RF_bc)
GWPinc20BR = GWP(cERF_inc, cERF_bc, 20, BR_RF_inc, BR_RF_bc)
GWPcnst20SS = GWP(cERF_cnst, cERF_bc, 20, SS_RF_cnst, SS_RF_bc)
GWPcnst20BR = GWP(cERF_cnst, cERF_bc, 20, BR_RF_cnst, BR_RF_bc)
GWPdecr20SS = GWP(cERF_decr, cERF_bc, 20, SS_RF_decr, SS_RF_bc)
GWPdecr20BR = GWP(cERF_decr, cERF_bc, 20, BR_RF_decr, BR_RF_bc)
#GWP100
GWPinc100SS = GWP(cERF_inc, cERF_bc, 100, SS_RF_inc, SS_RF_bc)
GWPinc100BR = GWP(cERF_inc, cERF_bc, 100, BR_RF_inc, BR_RF_bc)
GWPcnst100SS = GWP(cERF_cnst, cERF_bc, 100, SS_RF_cnst, SS_RF_bc)
GWPcnst100BR = GWP(cERF_cnst, cERF_bc, 100, BR_RF_cnst, BR_RF_bc)
GWPdecr100SS = GWP(cERF_decr, cERF_bc, 100, SS_RF_decr, SS_RF_bc)
GWPdecr100BR = GWP(cERF_decr, cERF_bc, 100, BR_RF_decr, BR_RF_bc)
#GTP20
GTPinc20SS = GTP(SSc2_dT_inc, SSc2_dT_bc, 20, SS_dT_inc, SS_dT_bc)
GTPinc20BR = GTP(BRc2_dT_inc, BRc2_dT_bc, 20, BR_dT_inc, BR_dT_bc)
GTPcnst20SS = GTP(SSc2_dT_cnst, SSc2_dT_bc, 20, SS_dT_cnst, SS_dT_bc)
GTPcnst20BR = GTP(BRc2_dT_cnst, BRc2_dT_bc, 20, BR_dT_cnst, BR_dT_bc)
GTPdecr20SS = GTP(SSc2_dT_decr, SSc2_dT_bc, 20, SS_dT_decr, SS_dT_bc)
GTPdecr20BR = GTP(BRc2_dT_decr, BRc2_dT_bc, 20, BR_dT_decr, BR_dT_bc)
#GTP100
GTPinc100SS = GTP(SSc2_dT_inc, SSc2_dT_bc, 100, SS_dT_inc, SS_dT_bc)
GTPinc100BR = GTP(BRc2_dT_inc, BRc2_dT_bc, 100, BR_dT_inc, BR_dT_bc)
GTPcnst100SS = GTP(SSc2_dT_cnst, SSc2_dT_bc, 100, SS_dT_cnst, SS_dT_bc)
GTPcnst100BR = GTP(BRc2_dT_cnst, BRc2_dT_bc, 100, BR_dT_cnst, BR_dT_bc)
GTPdecr100SS = GTP(SSc2_dT_decr, SSc2_dT_bc, 100, SS_dT_decr, SS_dT_bc)
GTPdecr100BR = GTP(BRc2_dT_decr, BRc2_dT_bc, 100, BR_dT_decr, BR_dT_bc)


print("######## GWP 20 ########")
print("GWP 20 SS inc =", GWPinc20SS)
print("GWP 20 BR inc =", GWPinc20BR)
print("GWP 20 SS cnst =", GWPcnst20SS)
print("GWP 20 BR cnst =", GWPcnst20BR)
print("GWP 20 SS decr =", GWPdecr20SS)
print("GWP 20 BR decr =", GWPdecr20BR)
print("######## GWP 100 ########")
print("GWP 100 SS inc =", GWPinc100SS)
print("GWP 100 BR inc =", GWPinc100BR)
print("GWP 100 SS cnst =", GWPcnst100SS)
print("GWP 100 BR cnst =", GWPcnst100BR)
print("GWP 100 SS decr =", GWPdecr100SS)
print("GWP 100 BR decr =", GWPdecr100BR)
print("######## GTP 20 ########")
print("GTP 20 SS inc =", GTPinc20SS)
print("GTP 20 BR inc =", GTPinc20BR)
print("GTP 20 SS cnst =", GTPcnst20SS)
print("GTP 20 BR cnst =", GTPcnst20BR)
print("GTP 20 SS decr =", GTPdecr20SS)
print("GTP 20 BR decr =", GTPdecr20BR)
print("######## GTP 100 ########")
print("GTP 100 SS inc =", GTPinc100SS)
print("GTP 100 BR inc =", GTPinc100BR)
print("GTP 100 SS cnst =", GTPcnst100SS)
print("GTP 100 BR cnst =", GTPcnst100BR)
print("GTP 100 SS decr =", GTPdecr100SS)
print("GTP 100 BR decr =", GTPdecr100BR)