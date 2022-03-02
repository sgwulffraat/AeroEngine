from math import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv("NL_trajectories_2021.txt", sep=',', header=0)
#data.columns = ["time", "icao24","lat","lon","velocity","heading","vertrate","callsign","onground","alert","spi","squawk","baroaltitude","geoaltitude","lastposupdate","lastcontact"]
sounding = pd.read_csv("Sounding_DeBilt_20210405.txt", sep='\s+',skiprows=3,header=None,)
sounding.columns = "PRES","HGHT","TEMP","DWPT","RELH","MIXR","DRCT","SPD","THTA","THTE","THTV"
#data2.columns = "PRES (hPa)","HGHT (m)","TEMP (C)","DWPT (C)","RELH (%)","MIXR (g/kg)","DRCT (deg)","SPD (m/s)","THTA (K)","THTE (K)","THTV (K)"



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

#plt.figure(0)
# plt.plot(flight1.lon, flight1.lat)
# plt.plot(flight2.lon, flight2.lat)
# plt.plot(flight3.lon, flight3.lat)
# plt.plot(flight4.lon, flight4.lat)
# plt.plot(flight5.lon, flight5.lat)
# plt.plot(flight6.lon, flight6.lat)
# plt.plot(flight7.lon, flight7.lat)
# plt.plot(flight8.lon, flight8.lat)
# plt.plot(flight9.lon, flight9.lat)
#plt.axis('equal')
#plt.show()

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

# info(flight1)
# info(flight2)
# info(flight3)
# info(flight4)
# info(flight5)
# info(flight6)
# info(flight7)
# info(flight8)
# info(flight9)

# ph20 > e_w (voor normaal en efficient)
# pah20> e_i
# T<-38

for i in range(0,len(sounding.TEMP)):
    if sounding.TEMP[i] < -38:
        print("Criteria 3 is satisfied above", sounding.HGHT[i],"m")
        break

#for i in range(0,len(sounding.TEMP)):
p_ah2O = sounding.PRES * sounding.MIXR * 28.97/18.01528

e_w = exp(-6096.9385*sounding.TEMP**-1 + 21.2409642 - 0.02711193*sounding.TEMP + 0.00001673952*sounding.TEMP**2 + 2.433502*np.log(sounding.TEMP))
print(e_w)

