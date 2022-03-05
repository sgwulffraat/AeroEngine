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

# plt.figure(0)
# plt.plot(df.Time, Ecnst, label = 'Constant')
# plt.plot(df.Time, Einc, label = 'Increasing')
# plt.plot(df.Time, Edecr, label = 'Decreasing')
# plt.plot(df.Time, Ebc, label = 'Base Case')
# plt.xlim(1940, 2120)
# plt.ylim(0, max(Einc)*1.05)
# plt.title("CO2 emissions per year")
# plt.xlabel('Time (years)')
# plt.ylabel('CO2 Emissions (Tg/year)')
# plt.grid(visible=True, which='major', color='#666666', linestyle='-')
# plt.minorticks_on()
# plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
# plt.legend(loc = 'best', fontsize = 'small')
# plt.show()

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
plt.plot(df.Time, SS_CO2_inc, label= 'S&S Increasing', color = 'tab:blue', linewidth=1)
plt.plot(df.Time, SS_CO2_decr, label= 'S&S Decreasing', color = 'tab:orange', linewidth=1)
plt.plot(df.Time, SS_CO2_cnst, label= 'S&S Constant', color = 'tab:green', linewidth=1)
plt.plot(df.Time, SS_CO2_bc, label= 'S&S Base Case', color = 'tab:red', linewidth=1)
plt.plot(df.Time, BR_CO2_inc, label= 'B&R Increasing', color = 'tab:blue', linestyle='dotted', linewidth=1)
plt.plot(df.Time, BR_CO2_decr, label= 'B&R Decreasing', color = 'tab:orange', linestyle='dotted', linewidth=1)
plt.plot(df.Time, BR_CO2_cnst, label= 'B&R Constant', color = 'tab:green', linestyle='dotted', linewidth=1)
plt.plot(df.Time, BR_CO2_bc, label= 'B&R Base Case', color = 'tab:red', linestyle='dotted', linewidth=1)
plt.xlim(1940, 2120)
plt.ylim(0, max(SS_CO2_inc)*1.05)
plt.title("CO2 concentration evolution")
plt.xlabel('Time (years)')
plt.ylabel('delta CO2 Emissions (ppmv)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'small')
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
plt.plot(df.Time, SS_RF_inc, label= 'S&S Increasing', color = 'tab:blue', linewidth=1)
plt.plot(df.Time, SS_RF_decr, label= 'S&S Decreasing', color = 'tab:orange', linewidth=1)
plt.plot(df.Time, SS_RF_cnst, label= 'S&S Constant', color = 'tab:green', linewidth=1)
plt.plot(df.Time, SS_RF_bc, label= 'S&S Base Case', color = 'tab:red', linewidth=1)
plt.plot(df.Time, BR_RF_inc, label= 'B&R Increasing', color = 'tab:blue', linestyle='dotted', linewidth=1)
plt.plot(df.Time, BR_RF_decr, label= 'B&R Decreasing', color = 'tab:orange', linestyle='dotted', linewidth=1)
plt.plot(df.Time, BR_RF_cnst, label= 'B&R Constant', color = 'tab:green', linestyle='dotted', linewidth=1)
plt.plot(df.Time, BR_RF_bc, label= 'B&R Base Case', color = 'tab:red', linestyle='dotted', linewidth=1)
plt.xlim(1940, 2120)
plt.ylim(0, max(SS_RF_inc)*1.05)
plt.title("Radiative Forcing Evolution)")
plt.xlabel('Time (years)')
plt.ylabel('Radiative Forcing (W/m^2)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'small')


## Question B ##
def RFcalccontrail(km):
    RF = (1.82E-9  * km)/1000
    return RF

SS_contrailRF_inc = RFcalccontrail(df.Flown_km_inc)
SS_contrailRF_decr = RFcalccontrail(df.Flown_km_decr)
SS_contrailRF_cnst = RFcalccontrail(df.Flown_km_cnst)
SS_contrailRF_bc = RFcalccontrail(df.Flown_km_bc)

plt.figure(3)
plt.plot(df.Time, SS_contrailRF_inc, label= 'S&S Increasing', color = 'tab:blue', linewidth=1)
plt.plot(df.Time, SS_contrailRF_decr, label= 'S&S Decreasing', color = 'tab:orange', linewidth=1)
plt.plot(df.Time, SS_contrailRF_cnst, label= 'S&S Constant', color = 'tab:green', linewidth=1)
plt.plot(df.Time, SS_contrailRF_bc, label= 'S&S Base Case', color = 'tab:red', linewidth=1)
plt.xlim(1940, 2120)
plt.ylim(0, max(SS_contrailRF_inc)*1.05)
plt.title("Radiative Forcing Evolution due to contrails)")
plt.xlabel('Time (years)')
plt.ylabel('Radiative Forcing (W/m^2)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'small')
plt.show()

### Part D ###

# S&S Temperature calculations
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
                y.append(f(df.Time[i]-df.Time[k])*RF[k])
            F = np.sum(y)
            dT.append(F)
    return dT

# SS_dT_inc = Tconv_SS(SS_RF_inc)
# SS_dT_decr = Tconv_SS(SS_RF_decr)
# SS_dT_cnst = Tconv_SS(SS_RF_cnst)
# SS_dT_bc = Tconv_SS(SS_RF_bc)
#
BR_dT_inc = Tconv_BR(BR_RF_inc)
BR_dT_decr = Tconv_BR(BR_RF_decr)
BR_dT_cnst = Tconv_BR(BR_RF_cnst)
BR_dT_bc = Tconv_BR(BR_RF_bc)

plt.figure(4)
plt.plot(df.Time, BR_dT_inc, label= 'B&R Increasing', color = 'tab:blue', linewidth=1)
plt.plot(df.Time, BR_dT_decr, label= 'B&R Decreasing', color = 'tab:orange', linewidth=1)
plt.plot(df.Time, BR_dT_cnst, label= 'B&R Constant', color = 'tab:green', linewidth=1)
plt.plot(df.Time, BR_dT_bc, label= 'B&R Base Case', color = 'tab:red', linewidth=1)
plt.xlim(1940, 2120)
plt.ylim(0, max(BR_dT_inc)*1.05)
plt.title("Temperature change over time")
plt.xlabel('Time (years)')
plt.ylabel('Delta T (C)')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc = 'best', fontsize = 'small')
plt.show()