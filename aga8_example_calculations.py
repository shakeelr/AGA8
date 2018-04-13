#AGA8 Detail Method Example Code
#Shakeel Rajwani
#27/04/2017

#! Python

import aga8

c_id = [0 for x in range(21)]
x_i = [0.0 for x in range(21)]
x_ic = [0.0 for x in range(21)]

#set up temperature array
T_data = [65.00, 65.00, 65.05, 65.05, 65.05,
	65.05, 65.05, 65.05, 65.05, 65.10,
	65.15, 65.15, 65.15, 65.15, 65.15,
	65.15, 65.15, 65.15, 65.20, 65.25]

#set up pressure array
P_data = [750.0, 750.5, 751.0, 751.0, 751.0,
	751.0, 751.0, 751.5, 751.5, 752.0,
	752.5, 752.5, 752.5, 752.5, 752.5,
	753.0, 753.5, 753.5, 754.0, 754.0]

#set up composition array using compositions for Gulf Coast Gas and Amarillo Gas
x_k = [[0 for x in range(21)] for y in range(2)]
x_k[0] = [0.965222, 0.002595, 0.005956, 0.018186, 0.004596,
	0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
	0.000977, 0.001007, 0.000473, 0.000324, 0.000664, 
	0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
	0.000000]
x_k[1] = [0.906724, 0.031284, 0.004676, 0.045279, 0.008280,
	0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
	0.001037, 0.001563, 0.000321, 0.000443, 0.000393,
	0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
	0.000000]

i_comp = 0
ncc = 0
i_level = 0
for i in range(0,20):
	if (x_k[i_comp][i] != 0.0):
		c_id[ncc] = i
		x_i[ncc]=x_k[i_comp][i]
		ncc+=1
		
aga8.paramdl(ncc, c_id)
aga8.chardl(ncc, x_i)
aga8.temp(T_data[1])