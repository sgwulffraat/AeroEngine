import matplotlib.pyplot as plt
import numpy as np

Heat_coarse = np.genfromtxt("heatflux_coarse.csv", delimiter=', ', skip_header=5)
Cp_coarse = np.genfromtxt("cp_coarse.csv", delimiter=', ', skip_header=5)

Heat_coarse_upwind = np.genfromtxt("heatflux_coarse_upwind.csv", delimiter=', ', skip_header=5)
Cp_coarse_upwind = np.genfromtxt("cp_coarse_upwind.csv", delimiter=', ', skip_header=5)

Heat_fine = np.genfromtxt("heatflux_fine.csv", delimiter=', ', skip_header=5)
Cp_fine = np.genfromtxt("cp_fine.csv", delimiter=', ', skip_header=5)

Heat_fine_upwind = np.genfromtxt("heatflux_fine_upwind.csv", delimiter=', ', skip_header=5)
Cp_fine_upwind = np.genfromtxt("cp_fine_upwind.csv", delimiter=', ', skip_header=5)

Heat_highintensity = np.genfromtxt("heatflux_highintensity.csv", delimiter=', ', skip_header=5)
Cp_highintensity = np.genfromtxt("cp_highintensity.csv", delimiter=', ', skip_header=5)

Heat_kepsilon = np.genfromtxt("heatflux_kepsilon.csv", delimiter=', ', skip_header=5)
Cp_kepsilon = np.genfromtxt("cp_kepsilon.csv", delimiter=', ', skip_header=5)

yplus_fine = np.genfromtxt("yplus_fine.csv", delimiter=', ', skip_header=5)
yplus_fine_upwind = np.genfromtxt("yplus_fine_upwind.csv", delimiter=', ', skip_header=5)
yplus_highintensity = np.genfromtxt("yplus_highintensity.csv", delimiter=', ', skip_header=5)
yplus_kepsilon = np.genfromtxt("yplus_kepsilon.csv", delimiter=', ', skip_header=5)
yplus_coarse = np.genfromtxt("yplus_coarse.csv", delimiter=', ', skip_header=5)
yplus_coarse_upwind = np.genfromtxt("yplus_coarse_upwind.csv", delimiter=', ', skip_header=5)

Xy = yplus_fine[:, 0]
yf = yplus_fine[:, 1]
yfu = yplus_fine_upwind[:, 1]
yhi = yplus_highintensity[:, 1]
yke = yplus_kepsilon[:, 1]
yc = yplus_coarse[:, 1]
ycu = yplus_coarse_upwind[:, 1]

Xh = Heat_coarse[:, 0]
Hc = Heat_coarse[:, 1]
Hf = Heat_fine[:, 1]
Hcu = Heat_coarse_upwind[:,1]
Hfu = Heat_fine_upwind[:,1]
Hhi = Heat_highintensity[:,1]
Hke = Heat_kepsilon[:,1]

Xcp = Cp_coarse[:, 0]
Cpc = Cp_coarse[:, 1]
Cpf = Cp_fine[:, 1]
Cpcu = Cp_coarse_upwind[:,1]
Cpfu = Cp_fine_upwind[:,1]
Cphi = Cp_highintensity[:,1]
Cpke = Cp_kepsilon[:,1]


plt.figure(0)
# plt.plot(Xh, Hke, label= 'k-Epsilon (Fine Mesh)', color = 'tab:brown')
plt.plot(Xh, Hc, label = 'Coarse Mesh', color = 'tab:blue')
# plt.plot(Xh, Hhi, label= 'High Turbulence Intensity (Fine Mesh)', color = 'tab:purple')
# plt.plot(Xh, Hf, label= 'Fine Mesh', color = 'tab:orange')
# plt.plot(Xh, Hcu, label= 'Upwind scheme (Coarse Mesh)', color = 'tab:green')
# plt.plot(Xh, Hfu, label= 'Upwind scheme (Fine Mesh)', color = 'tab:red')
plt.xlabel('X [m]')
plt.ylabel('Heat Flux Coefficient [W m^-2 K^-1)]')
plt.grid()
plt.title('Heat Flux distribution')
plt.legend()

plt.figure(1)
# plt.plot(Xcp, Cphi, label= 'High Turbulence Intensity (Fine Mesh)', color = 'tab:purple')
# plt.plot(Xcp, Cpke, label= 'k-Epsilon (Fine Mesh)', color = 'tab:brown')
plt.plot(Xcp, Cpc, label = 'Coarse Mesh', color = 'tab:blue')
# plt.plot(Xcp, Cpf, label= 'Fine Mesh', color = 'tab:orange')
# plt.plot(Xcp, Cpcu, label= 'Upwind scheme (Coarse Mesh)', color = 'tab:green')
# plt.plot(Xcp, Cpfu, label= 'Upwind scheme (Fine Mesh)', color = 'tab:red')
plt.title('Cp distribution')
plt.xlabel('X [m]')
plt.ylabel('Cp [-]')
plt.grid()
plt.legend()

plt.figure(2)
#plt.plot(Xy, yhi, label= 'High Turbulence Intensity (Fine Mesh)', color = 'tab:purple')
#plt.plot(Xy, yc, label= 'Coarse Mesh', color = 'tab:blue')
# plt.plot(Xy, ycu, label= 'Upwind scheme (Coarse Mesh)', color = 'tab:red')
plt.plot(Xy, yf, label= 'Fine Mesh', color = 'tab:orange')
#plt.plot(Xy, yfu, label= 'Upwind scheme (Fine Mesh)', color = 'tab:red')
#plt.plot(Xy, yke, label= 'k-Epsilon (Fine Mesh)', color = 'tab:brown')
plt.title('y+ distribution')
plt.xlabel('X [m]')
plt.ylabel('y+ [-]')
plt.grid()
plt.legend()

plt.show()