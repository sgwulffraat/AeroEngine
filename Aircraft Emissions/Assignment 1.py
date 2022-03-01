from math import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd


data = pd.read_csv("NL_trajectories_2021.txt", sep=',', header=0)
#data.columns = ["time", "icao24","lat","lon","velocity","heading","vertrate","callsign","onground","alert","spi","squawk","baroaltitude","geoaltitude","lastposupdate","lastcontact"]

map = gpd.read_file(gpd.datasets.get_path('nybb'))

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

plt.figure(0)
plt.plot(flight1.lon, flight1.lat)
plt.plot(flight2.lon, flight2.lat)
plt.plot(flight3.lon, flight3.lat)
plt.plot(flight4.lon, flight4.lat)
plt.plot(flight5.lon, flight5.lat)
plt.plot(flight6.lon, flight6.lat)
plt.plot(flight7.lon, flight7.lat)
plt.plot(flight8.lon, flight8.lat)
plt.plot(flight9.lon, flight9.lat)
plt.axis('equal')
#plt.show()


def info(flight):
    duration = max(flight.time) - min(flight.time)
    icao = flight["icao24"].iloc[0]
    startalt = flight["geoaltitude"].iloc[0]
    finalalt = flight["geoaltitude"].iloc[-1]
    if np.sign(finalalt-startalt) == -1:
        phase = "Approach"
    else:
        phase = "Departure"
    print("Info flight",flight["callsign"].iloc[0],":   ICAO Code =", icao,", Start alt =", startalt,", Final alt =", finalalt,", Phase =", phase,", Duration =", duration)
    return

info(flight1)
info(flight2)
info(flight3)
info(flight4)
info(flight5)
info(flight6)
info(flight7)
info(flight8)
info(flight9)
