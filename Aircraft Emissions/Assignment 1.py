from math import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv("NL_trajectories_2021.txt", sep=',', header=0)
#data.columns = ["time", "icao24","lat","lon","velocity","heading","vertrate","callsign","onground","alert","spi","squawk","baroaltitude","geoaltitude","lastposupdate","lastcontact"]
sounding = pd.read_csv("Sounding_DeBilt_20210405.txt", sep='\s+',skiprows=3,header=None,)
sounding.columns = "PRES","HGHT","TEMP","DWPT","RELH","MIXR","DRCT","SPD","THTA","THTE","THTV"
#data2.columns = "PRES (hPa)","HGHT (m)","TEMP (C)","DWPT (C)","RELH (%)","MIXR (g/kg)","DRCT (deg)","SPD (m/s)","THTA (K)","THTE (K)","THTV (K)"


# BB = ("lon:", min(data.lon),",", max(data.lon), "lat:", min(data.lat),",",max(data.lat))
# print(BB)

grouped = data.groupby(data.callsign)
flight1 = grouped.get_group("PGT27W")
flight2 = grouped.get_group("AFR96ZN")
flight3 = grouped.get_group("KLM862")
flight4 = grouped.get_group("AFR15AH")
flight5 = grouped.get_group("PGT53X")
flight6 = grouped.get_group("DLH908")
flight7 = grouped.get_group("GEC8201")
flight8 = grouped.get_group("KLM685")
flight9 = grouped.get_group("KLM867")

fig, ax = plt.subplots()
im = plt.imread("map2.png")
ax.imshow(im, extent=[min(data.lon), max(data.lon), min(data.lat), max(data.lat)])
ax.plot(flight1.lon, flight1.lat, label = flight1["callsign"].iloc[0])
ax.plot(flight2.lon, flight2.lat, label = flight2["callsign"].iloc[0])
ax.plot(flight3.lon, flight3.lat, label = flight3["callsign"].iloc[0])
ax.plot(flight4.lon, flight4.lat, label = flight4["callsign"].iloc[0])
ax.plot(flight5.lon, flight5.lat, label = flight5["callsign"].iloc[0])
ax.plot(flight6.lon, flight6.lat, label = flight6["callsign"].iloc[0])
ax.plot(flight7.lon, flight7.lat, label = flight7["callsign"].iloc[0])
ax.plot(flight8.lon, flight8.lat, label = flight8["callsign"].iloc[0])
ax.plot(flight9.lon, flight9.lat, label = flight9["callsign"].iloc[0])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title("Trajectories")
plt.legend(loc = 'upper right', fontsize = 'small')
plt.show()

plt.figure(1)
plt.plot(sounding.TEMP, sounding.HGHT)
plt.xlabel("Temperature (C)")
plt.ylabel("Altitude (m)")
plt.show()

def info(flight):
    duration = max(flight.time) - min(flight.time)
    icao = flight["icao24"].iloc[0]
    if str(flight["geoaltitude"].iloc[0]) == "nan":
        startalt = 0
    else:
        startalt = flight["geoaltitude"].iloc[0]

    if str(flight["geoaltitude"].iloc[-1]) == "nan":
        finalalt = 0
    else:
        finalalt = flight["geoaltitude"].iloc[-1]

    if np.sign(finalalt-startalt) == -1:
        phase = "Approach"
    else:
        phase = "Departure"
    print("Info flight",flight["callsign"].iloc[0],":   ICAO Code =", icao,", Start alt =", startalt,", Final alt =", finalalt,", Phase =", phase,", Duration =", duration)
    return

def contrailplt(flight):
    lat = []
    lon = []
    for i in range(0,len(flight.geoaltitude)):
        if flight["geoaltitude"].iloc[i] > 8037 and flight["geoaltitude"].iloc[i] < 8803:
            lat.append(flight["lat"].iloc[i])
            lon.append(flight["lon"].iloc[i])
        elif flight["geoaltitude"].iloc[i] > 8909 and flight["geoaltitude"].iloc[i] < 8952:
            lat.append(flight["lat"].iloc[i])
            lon.append(flight["lon"].iloc[i])
        elif flight["geoaltitude"].iloc[i] > 8984 and flight["geoaltitude"].iloc[i] < 9014:
            lat.append(flight["lat"].iloc[i])
            lon.append(flight["lon"].iloc[i])
        elif flight["geoaltitude"].iloc[i] > 9644 and flight["geoaltitude"].iloc[i] < 9855:
            lat.append(flight["lat"].iloc[i])
            lon.append(flight["lon"].iloc[i])
        elif flight["geoaltitude"].iloc[i] > 10625 and flight["geoaltitude"].iloc[i] < 10706:
            lat.append(flight["lat"].iloc[i])
            lon.append(flight["lon"].iloc[i])
    return ax.plot(lon,lat, label = flight["callsign"].iloc[0])

# info(flight1)
# info(flight2)
# info(flight3)
# info(flight4)
# info(flight5)
# info(flight6)
# info(flight7)
# info(flight8)
# info(flight9)


## Contrail calculations ##
c_p = 1003.5    #J/kg/K
EI = 0.14*18.016/(2*1.008)          #kg/kg
Q = 43E6        #J
eff = 0.8      #-

p_h2O = sounding.PRES*100 * sounding.MIXR/1000 * 28.97/18.01528
p_h2O_e = p_h2O + 28.97/18.01528 * sounding.PRES*100 * (EI - sounding.MIXR/1000) * c_p / ((1- eff)*Q)
TempK = sounding.TEMP + 273.15
e_w = np.exp(-6096.9385*TempK**-1 + 21.2409642 - 0.02711193*TempK + 0.00001673952*TempK**2 + 2.433502*np.log(TempK))
e_i = np.exp(-6024.5282*TempK**-1 + 29.32707 + 0.010613868*TempK - 0.000013198825*TempK**2 - 0.49382577*np.log(TempK))


# lowb = 1
# j = 1
# for i in range(1,len(sounding.TEMP)):
#     if p_h2O_e[i]>e_w[i] and p_h2O_e[i-1]<e_w[i-1]:
#         lowb = sounding.HGHT[i]
#     elif p_h2O_e[i-1]>e_w[i-1] and p_h2O_e[i]<e_w[i]:
#         upb = sounding.HGHT[i-1]
#         print("Criteria 1 is satisfied in range", j, "from", lowb,"m -", upb, "m")
#         j = j+1

lowb_i = 1
k = 1

for i in range(1, len(sounding.TEMP)):
    if p_h2O[i]>e_i[i] and p_h2O[i-1]<e_i[i-1]:
        lowb_i = sounding.HGHT[i]
    elif p_h2O[i-1]>e_i[i-1] and p_h2O[i]<e_i[i]:
        #if sounding.HGHT[i-1] > 8037:
        upb_i = sounding.HGHT[i-1]
        print("Criteria 2 is satisfied in range", k, "from", lowb_i,"m -", upb_i, "m")
        k = k+1

# for i in range(0,len(sounding.TEMP)):
#     if sounding.TEMP[i] < -38:
#         print("Criteria 3 is satisfied above", sounding.HGHT[i],"m")
#         break

fig, ax = plt.subplots()
im = plt.imread("map2.png")
ax.imshow(im, extent=[min(data.lon), max(data.lon), min(data.lat), max(data.lat)])
contrailplt(flight1)
contrailplt(flight2)
contrailplt(flight3)
contrailplt(flight4)
contrailplt(flight5)
contrailplt(flight6)
contrailplt(flight7)
contrailplt(flight8)
contrailplt(flight9)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title("Contrail formation along trajectories")
plt.legend(loc = 'upper right', fontsize = 'small')
plt.show()



