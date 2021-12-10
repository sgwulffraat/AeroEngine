from math import *
import numpy as np

n_comb = 0.99   #[-]
c_p_g = 1150    #[-]
LHV = ((48*393.5+46*241.826)-211.46*4)/(167.31102*4)    #[MJ/kg]
T_3 = np.array([612, 723, 765, 820])    #[K]
T_4 = np.array([1234, 1423, 1523, 1631])    #[K]
m_dot_air = np.array([17.12, 29.17, 34.58, 39.58])  #[kg/s]
m_dot_f = (m_dot_air * c_p_g * (T_4 - T_3)) / (n_comb * LHV *10**6) #[kg/s]

AF_stoch = (71*3.76*28.01+71*32)/(4*167.31102)  #[-]
AF = m_dot_air/m_dot_f  #[-]
ER = AF_stoch/AF    #[-]

