import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('GAW_data')


data = np.loadtxt("GAW_data/GAW-2.txt", skiprows=1)
x = data[:, 0]
y = data[:, 1]

# data = np.loadtxt("GAW_data/GAW-2_withoutskirt.txt", skiprows=1)
# x = data[:, 0]
# y = data[:, 1]

data_flap = np.loadtxt("GAW_data/GAW-2_flap.txt", skiprows=1)
x_flap = data_flap[:, 0]
y_flap = data_flap[:, 1]

#Ref data
ref_0 = np.loadtxt("GAW_data/GAW2_ref_0.txt", skiprows=1)
ref_10 = np.loadtxt("GAW_data/GAW2_ref_10.txt", skiprows=1)
ref_20 = np.loadtxt("GAW_data/GAW2_ref_20.txt", skiprows=1)
ref_30 = np.loadtxt("GAW_data/GAW2_ref_30.txt", skiprows=1)
a_ref_0 = ref_0[:, 0]
cl_ref_0 = ref_0[:, 1]
a_ref_10 = ref_10[:, 0]
cl_ref_10 = ref_10[:, 1]
a_ref_20 = ref_20[:, 0]
cl_ref_20 = ref_20[:, 1]
a_ref_30 = ref_30[:, 0]
cl_ref_30 = ref_30[:, 1]

#JavaFoil data
jf_0_eppler = np.loadtxt("GAW_data/GAW2_0_eppler.txt", skiprows=5)
jf_10_eppler = np.loadtxt("GAW_data/GAW2_10_eppler.txt", skiprows=5)
jf_20_eppler = np.loadtxt("GAW_data/GAW2_20_eppler.txt", skiprows=5)
jf_30_eppler = np.loadtxt("GAW_data/GAW2_30_eppler.txt", skiprows=5)

jf_0_calcfoil = np.loadtxt("GAW_data/GAW2_0_calcfoil.txt", skiprows=5)
jf_10_calcfoil = np.loadtxt("GAW_data/GAW2_10_calcfoil.txt", skiprows=5)
jf_20_calcfoil = np.loadtxt("GAW_data/GAW2_20_calcfoil.txt", skiprows=5)
jf_30_calcfoil = np.loadtxt("GAW_data/GAW2_30_calcfoil.txt", skiprows=5)

a_jf_0_calc = jf_0_calcfoil[:, 0]
cl_jf_0_calc = jf_0_calcfoil[:, 1]
a_jf_10_calc = jf_10_calcfoil[:28, 0]
cl_jf_10_calc = jf_10_calcfoil[:28, 1]
a_jf_20_calc = jf_20_calcfoil[:27, 0]
cl_jf_20_calc = jf_20_calcfoil[:27, 1]
a_jf_30_calc = jf_30_calcfoil[:24, 0]
cl_jf_30_calc = jf_30_calcfoil[:24, 1]

a_jf_0_eppler = jf_0_eppler[:, 0]
cl_jf_0_eppler = jf_0_eppler[:, 1]
a_jf_10_eppler = jf_10_eppler[:28, 0]
cl_jf_10_eppler = jf_10_eppler[:28, 1]
a_jf_20_eppler = jf_20_eppler[:27, 0]
cl_jf_20_eppler = jf_20_eppler[:27, 1]
a_jf_30_eppler = jf_30_eppler[:24, 0]
cl_jf_30_eppler = jf_30_eppler[:24, 1]


fig1, axes = plt.subplots()
plt.plot(x, y, label='Airfoil', color='tab:blue', marker='o', markersize=2)
plt.plot(x_flap, y_flap, label='Flap', color='tab:orange', marker='o', markersize=2)
plt.xlabel('Chordwise position (x/c)')
plt.ylabel('(z/c)')
plt.axis('equal')
plt.grid()
plt.show()

fig2 = plt.figure()
plt.plot(a_ref_0, cl_ref_0, label='Reference data', color='dimgray', linestyle='dashed', linewidth='1', marker='v', markersize=4)
plt.plot(a_ref_10, cl_ref_10, color='dimgray', linestyle='dashed', linewidth='1', marker='v', markersize=4)
plt.plot(a_ref_20, cl_ref_20, color='dimgray', linestyle='dashed', linewidth='1', marker='v', markersize=4)
plt.plot(a_ref_30, cl_ref_30,  color='dimgray', linestyle='dashed', linewidth='1', marker='v', markersize=4)
plt.plot(a_jf_0_calc, cl_jf_0_calc, label='Calcfoil stall model', color='tab:blue', marker='o', markersize=2)
plt.plot(a_jf_0_eppler, cl_jf_0_eppler, label='Eppler stall model', color='tab:orange', marker='o', markersize=2)
plt.plot(a_jf_10_calc, cl_jf_10_calc, color='tab:blue', marker='o', markersize=2)
plt.plot(a_jf_10_eppler, cl_jf_10_eppler, color='tab:orange', marker='o', markersize=2)
plt.plot(a_jf_20_calc, cl_jf_20_calc, color='tab:blue', marker='o', markersize=2)
plt.plot(a_jf_20_eppler, cl_jf_20_eppler, color='tab:orange', marker='o', markersize=2)
plt.plot(a_jf_30_calc, cl_jf_30_calc, color='tab:blue', marker='o', markersize=2)
plt.plot(a_jf_30_eppler, cl_jf_30_eppler, color='tab:orange', marker='o', markersize=2)
plt.legend()
plt.xlabel(r'Angle of attack ($\alpha$)[]')
plt.ylabel(r'Lift ($C_L$) [-]')
plt.grid()
plt.show()