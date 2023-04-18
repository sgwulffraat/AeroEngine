import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('GAW_data')


data = np.loadtxt("GAW_data/GAW-2.txt", skiprows=1)
x = data[:, 0]
y = data[:, 1]

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
a_jf_10_calc = jf_10_calcfoil[:, 0]
cl_jf_10_calc = jf_10_calcfoil[:, 1]
a_jf_20_calc = jf_20_calcfoil[:, 0]
cl_jf_20_calc = jf_20_calcfoil[:, 1]
a_jf_30_calc = jf_30_calcfoil[:, 0]
cl_jf_30_calc = jf_30_calcfoil[:, 1]

a_jf_0_eppler = jf_0_eppler[:, 0]
cl_jf_0_eppler = jf_0_eppler[:, 1]
a_jf_10_eppler = jf_0_eppler[:, 0]
cl_jf_10_eppler = jf_0_eppler[:, 1]
a_jf_20_eppler = jf_20_eppler[:, 0]
cl_jf_20_eppler = jf_20_eppler[:, 1]
a_jf_30_eppler = jf_30_eppler[:, 0]
cl_jf_30_eppler = jf_30_eppler[:, 1]


fig1, axes = plt.subplots()
circle = plt.Circle((0.7119, -0.0071), 0.0117, fill= False)
plt.plot(x, y, label='Airfoil', color='tab:blue')
plt.plot(x_flap, y_flap, label='Flap', color='tab:orange')
plt.xlabel('Chordwise position (x/c)')
plt.ylabel('(z/c)')
plt.axis('equal')
axes.set_aspect( 1 )
axes.add_artist(circle)
plt.grid()
plt.show()

fig2 = plt.figure()
plt.plot(a_ref_0, cl_ref_0, label='alpha = 0', color='black')
plt.plot(a_ref_10, cl_ref_10, label='alpha = 10', color='black')
plt.plot(a_ref_20, cl_ref_20, label='alpha = 10', color='black')
plt.plot(a_ref_30, cl_ref_30, label='alpha = 30', color='black')
plt.plot(a_jf_0_calc, cl_jf_0_calc, label='alpha = 0 (calc)', color='tab:red')
plt.plot(a_jf_0_eppler, cl_jf_0_eppler, label='alpha = 0 (eppler)', color='tab:green')
plt.plot(a_jf_10_calc, cl_jf_10_calc, label='alpha = 10 (calc)', color='tab:red')
plt.plot(a_jf_10_eppler, cl_jf_10_eppler, label='alpha = 20 (eppler)', color='tab:green')
plt.plot(a_jf_20_calc, cl_jf_20_calc, label='alpha = 20 (calc)', color='tab:red')
plt.plot(a_jf_20_eppler, cl_jf_20_eppler, label='alpha = 20 (eppler)', color='tab:green')
plt.plot(a_jf_30_calc, cl_jf_30_calc, label='alpha = 30 (calc)', color='tab:red')
plt.plot(a_jf_30_eppler, cl_jf_30_eppler, label='alpha = 30 (eppler)', color='tab:green')
plt.legend()
plt.xlabel('Angle of attack ()[]')
plt.ylabel('Lift (c_l) [-]')
plt.grid()
plt.show()