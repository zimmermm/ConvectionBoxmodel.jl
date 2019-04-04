#===========================================
Emulates the dependence between u_w^3 on U10
by a linear interpolation of precalculated
values
===========================================#

u_w_discrete = [0, 1.35e-08, 3.28e-08, 5.49e-08, 7.93e-08, 1.05e-07, 1.33e-07, 1.62e-07, 2.13e-07, 3.20e-07, 4.62e-07, 6.46e-07, 8.78e-07, 1.17e-06, 1.52e-06, 1.95e-06, 2.46e-06, 3.07e-06, 3.78e-06, 4.61e-06, 5.56e-06, 6.66e-06, 7.91e-06, 9.32e-06, 1.09e-05, 1.27e-05, 1.47e-05, 1.69e-05, 1.94e-05, 2.21e-05, 2.50e-05, 2.83e-05, 3.19e-05, 3.58e-05, 4.00e-05, 4.46e-05, 4.95e-05, 5.49e-05, 6.06e-05, 6.68e-05, 7.34e-05, 8.05e-05, 8.81e-05, 9.62e-05, 0.000104755, 0.000113888, 0.000123579, 0.000133846, 0.000144709, 0.000156185, 0.000168293, 0.000181052, 0.000194479, 0.000208593, 0.00022341, 0.000238949, 0.000255226, 0.000272259, 0.000290063, 0.000308656, 0.000328052, 0.000348267, 0.000369315, 0.000391212, 0.00041397, 0.000437604, 0.000462126, 0.000487548, 0.000513882, 0.000541139, 0.000569329, 0.000598461, 0.000628545, 0.000659588, 0.000691599, 0.000724583, 0.000758546, 0.000793493, 0.000829428, 0.000866355, 0.000904275, 0.00094319, 0.0009831, 0.001024004, 0.001065901, 0.001108788, 0.001152661, 0.001197515, 0.001243345, 0.001290142, 0.0013379, 0.001386608, 0.001436257, 0.001486834, 0.001538325, 0.001590719, 0.001643997, 0.001698145, 0.001753143, 0.001808973]
U10_discrete = [0, 0.505050505, 1.01010101, 1.515151515, 2.02020202, 2.525252525, 3.03030303, 3.535353535, 4.04040404, 4.545454545, 5.050505051, 5.555555556, 6.060606061, 6.565656566, 7.070707071, 7.575757576, 8.080808081, 8.585858586, 9.090909091, 9.595959596, 10.1010101, 10.60606061, 11.11111111, 11.61616162, 12.12121212, 12.62626263, 13.13131313, 13.63636364, 14.14141414, 14.64646465, 15.15151515, 15.65656566, 16.16161616, 16.66666667, 17.17171717, 17.67676768, 18.18181818, 18.68686869, 19.19191919, 19.6969697, 20.2020202, 20.70707071, 21.21212121, 21.71717172, 22.22222222, 22.72727273, 23.23232323, 23.73737374, 24.24242424, 24.74747475, 25.25252525, 25.75757576, 26.26262626, 26.76767677, 27.27272727, 27.77777778, 28.28282828, 28.78787879, 29.29292929, 29.7979798, 30.3030303, 30.80808081, 31.31313131, 31.81818182, 32.32323232, 32.82828283, 33.33333333, 33.83838384, 34.34343434, 34.84848485, 35.35353535, 35.85858586, 36.36363636, 36.86868687, 37.37373737, 37.87878788, 38.38383838, 38.88888889, 39.39393939, 39.8989899, 40.4040404, 40.90909091, 41.41414141, 41.91919192, 42.42424242, 42.92929293, 43.43434343, 43.93939394, 44.44444444, 44.94949495, 45.45454545, 45.95959596, 46.46464646, 46.96969697, 47.47474747, 47.97979798, 48.48484848, 48.98989899, 49.49494949, 50]
eps_U10_emulator = LinearInterpolation(u_w_discrete, U10_discrete)