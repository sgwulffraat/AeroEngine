import matplotlib.pyplot as plt
import numpy as np

Heat_coarse = np.genfromtxt("heatflux_coarse.csv", delimiter=', ', skip_header=5)
Heat_fine = np.genfromtxt("heatflux_fine.csv", delimiter=', ', skip_header=5)
Cp_coarse = np.genfromtxt("cp_coarse.csv", delimiter=', ', skip_header=5)
Cp_fine = np.genfromtxt("cp_fine.csv", delimiter=', ', skip_header=5)

Xh = Heat_coarse[:, 0]
Hc = Heat_coarse[:, 1]
Hf = Heat_fine[:, 1]

Xcp = Cp_coarse[:, 0]
Cpc = Cp_coarse[:, 1]
Cpf = Cp_fine[:, 1]

plt.figure(0)
plt.plot(Xh, Hc, label = 'Coarse Mesh')
plt.plot(Xh, Hf, label= 'Fine Mesh')
plt.xlabel('X [m]')
plt.ylabel('Heat Flux Coefficient [W m^-2 K^-1)]')
plt.grid()
plt.title('Heat Flux distribution')
plt.legend()


plt.figure(1)
plt.plot(Xcp, Cpc, label = 'Coarse Mesh')
plt.plot(Xcp, Cpf, label= 'Fine Mesh')
plt.title('Cp distribution')
plt.xlabel('X [m]')
plt.ylabel('Cp [-]')
plt.grid()
plt.legend()
plt.show()