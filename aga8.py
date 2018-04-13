#AGA 8 Calculation
#Ported from AGA Report 8 FORTRAN code listings, with some minor enhancements
#This code may seem convoluded in some places due to the direct port
#Significant room to optimize, especially in braket
#May build a native python version in the future
#revision 0
#Shakeel Rajwani
#02/05/2017

import math

#Dictionary of component IDs
components = {
	0: 'C1',
	1: 'N2',
	2: 'CO2',
	3: 'C2',
	4: 'C3',
	5: 'H2O',
	6: 'H2S',
	7: 'H2',
	8: 'CO',
	9: 'O2',
	10: 'iC4',
	11: 'nC4',
	12: 'iC5',
	13: 'nC5',
	14: 'C6',
	15: 'C7',
	16: 'C8',
	17: 'c9',
	18: 'C10',
	19: 'He',
	20: 'Ar'
	}

#Initialize Global Variables and Constants
rgas = 8.31451*(10**-3) #Gas constant
told = 0.0
tlow = 0.0
thigh = 100000
plow = 0.5*(10**-9)
phigh = 275.0
dhigh = 12.0

#Initialize fn array used in equation of state
fn = [0.0 for x in range(58)] 

#Initialize detla used in braket
delta = 0.0

"""
AGA8 Subroutine A.3.8
BLOCK DATA

Shak Notes:  In retrospect I should have redone the matricies in numpy, but I didn't want unnessesary dependencies.
"""

#Equation of State Parameters a[0:52]
a = [0.153832600, 1.341953000, -2.998583000, -0.048312280,
	0.375796500, -1.589575000, -0.053588470,					#a[0,6]
	0.886594630, -0.710237040, -1.471722000, 1.321850350,
	-0.786659250,												#a[8:11] from GRI1991 A[8,11]
	0.229129 * (10**-8),										#a[12]
	0.157672400, -0.436386400, -0.044081590, -0.003433888,
	0.032059050, 0.024873550, 0.073322790, -0.001600573,
	0.642470600, -0.416260100, -0.066899570, 0.279179500,
	-0.696605100, -0.002860589, -0.008098836, 3.150547000,
	0.007224479, -0.705752900, 0.534979200, -0.079314010,
	-1.418465000,												#a[13:33]
	-0.599905 * (10**-16),										#a[34]
	0.105840200, 0.034317290, -0.007022847, 0.024955870,
	0.042968180, 0.746545300, -0.291961300, 7.294616000,
	-9.936757000, -0.005399808, -0.243256700, 0.049870160,
	0.003733797, 1.874951000, 0.002168144, -0.658716400,
	0.000205518, 0.009776195, -0.020487080, 0.015573220,
	0.006862415, -0.001226752, 0.002850908]						#a[35:57]

#initialize first virial coefficients, bn, and second virial coefficient bmix
b1 = 0
b2 = 0
b3 = 0
b4 = 0
b5 = 0
b6 = 0
b7 = 0
b8 = 0
b9 = 0
b10 = 0
b11 = 0
b12 = 0
b13 = 0
b14 = 0
b15 = 0
b16 = 0
b17 = 0
b18 = 0
bmix = 0

#Initialize Individual Component Parameters
#qib = Qi, hib=Fi, rkib = Ki, eib=Ei, wib=Gi, cmwb=MWi, mib=Mi, dib=Delta
di = [0.0 for x in range(21)] 
dib = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] 

qi = [0.0 for x in range(21)] 
qib = [0.0, 0.0, 0.69, 0.0, 0.0, 1.06775, 0.633276, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
q2p0 = 0

hi = [0.0 for x in range(21)] 
hib = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
hh = 0

rki = [0.0 for x in range(21)] 
rkib = [0.4619255, 0.4479153, 0.4557489, 0.5279209,
	0.5837490, 0.3825868, 0.4618263, 0.3514916,
	0.4533894, 0.4186954, 0.6406937, 0.6341423,
	0.6738577, 0.6798307, 0.7175118, 0.7525189,
	0.7849550, 0.8152731, 0.8437826, 0.3589888,
	0.4216551]
rk3p0 = 0

ei = [0.0 for x in range(21)] 
eib = [151.318300, 99.737780, 241.960600, 244.166700,
	298.118300, 514.015600, 296.355000, 26.957940,
	105.534800, 122.766700, 324.068900, 337.638900,
	365.599900, 370.682300, 402.636293, 427.722630,
	450.325022, 470.840891, 489.558373, 2.610111,
	119.629900]
uu = 0
	
wi = [0.0 for x in range(21)] 
wib = [0.000000, 0.027815, 0.189065, 0.079300, 0.141239,
	0.332500, 0.088500, 0.034369, 0.038953, 0.021000,
	0.256692, 0.281835, 0.332267, 0.366911, 0.289731,
	0.337542, 0.383381, 0.427354, 0.469659, 0.000000,
	0.000000]
ww = 0
	
cmw = [0.0 for x in range(21)] 
cmwb = [16.0430, 28.0135, 44.0100, 30.0700, 44.0970,
	18.0153, 34.0820, 2.0159, 28.0100, 31.9988,
	58.1230, 58.1230, 72.1500, 72.1500, 86.1770, 100.2040, 114.2310,
	128.2580, 142.2850, 4.0026, 39.9480]

mi = [0.0 for x in range(21)] 
mib = [0.0, 0.0, 0.0, 0.0, 0.0, 1.5822, 0.390, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#Binary Interaction Parameters
#beijb = Eij, buijb = Uij, bkijb = Kij, bwijb = Gij
buij = [[0 for x in range(21)] for y in range(21)]
buijb = [[0 for x in range(21)] for y in range(21)]
buijb[0] = [1.000000, 0.886106, 0.963827, 1.000000, 0.990877,
	1.000000, 0.736833, 1.156390, 1.000000, 1.000000,
	1.000000, 0.992291, 1.000000, 1.003670, 1.302576,
	1.191904, 1.205769, 1.219634, 1.233498, 1.000000,
	1.000000]
buijb[1] = [1.000000, 0.835058, 0.816431, 0.915502, 1.000000,
	0.993476, 0.408838, 1.000000, 1.000000, 1.000000,
	0.993556, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
buijb[2] = [1.000000, 0.969870, 1.000000, 1.000000, 1.045290,
	1.000000, 0.900000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.066638, 1.077634, 1.088178,
	1.098291, 1.108021, 1.000000, 1.000000]
buijb[3] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
buijb[4] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
buijb[5] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
buijb[6] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.028973, 1.033754,
	1.038338, 1.042735, 1.046966, 1.000000, 1.000000]
buijb[7] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
buijb[8] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
buijb[9] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
buijb[10] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
	1.000000]
buijb[11] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
buijb[12] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
buijb[13] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
buijb[14] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
buijb[15] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
buijb[16] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
buijb[17] = [1.000000, 1.000000, 1.000000, 1.000000]
buijb[18] = [1.000000, 1.000000, 1.000000]
buijb[19] = [1.000000, 1.000000]
buijb[20] = [1.000000]

bkij = [[0 for x in range(21)] for y in range(21)] 
bkijb = [[0 for x in range(21)] for y in range(21)] 
bkijb[0] = [1.000000, 1.003630, 0.995933, 1.000000, 1.007619,
	1.000000, 1.000080, 1.023260, 1.000000, 1.000000,
	1.000000, 0.997596, 1.000000, 1.002529, 0.982962,
	0.983565, 0.982707, 0.981849, 0.980991, 1.000000,
	1.000000]
bkijb[1] = [1.000000, 0.982361, 1.007960, 1.000000, 1.000000, 
	0.942696, 1.032270, 1.000000, 1.000000, 1.000000, 
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
bkijb[2] = [1.000000, 1.008510, 1.000000, 1.000000, 1.007790, 
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
	1.000000, 1.000000, 0.910183, 0.895362, 0.881152, 
	0.867520, 0.854406, 1.000000, 1.000000]
bkijb[3] = [1.000000, 0.986893, 1.000000, 0.99969, 1.020340,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
bkijb[4] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
bkijb[5] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
bkijb[6] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 0.968130, 0.962870,
	0.957828, 0.952441, 0.948338, 1.000000, 1.000000]
bkijb[7] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
bkijb[8] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
bkijb[9] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
bkijb[10] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
bkijb[11] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
bkijb[12] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
bkijb[13] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
bkijb[14] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
bkijb[15] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
bkijb[16] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
bkijb[17] = [1.000000, 1.000000, 1.000000, 1.000000]
bkijb[18] = [1.000000, 1.000000, 1.000000]
bkijb[19] = [1.000000, 1.000000]
bkijb[20] = [1.000000]

beij = [[0 for x in range(21)] for y in range(21)] 
beijb = [[0 for x in range(21)] for y in range(21)] 
beijb[0] = [1.000000, 0.971640, 0.960644, 1.000000, 0.994635,
	0.708218, 0.931484, 1.170520, 0.990126, 1.000000,
	1.019530, 0.989844, 1.002350, 0.999268, 1.107274,
	0.880880, 0.880973, 0.881067, 0.881161, 1.000000,
	1.000000]
beijb[1] = [1.000000, 1.022740, 0.970120, 0.945939, 0.746954,
	0.902271, 1.086320, 1.005710, 1.021000, 0.946914,
	0.973384, 0.959340, 0.945520, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
beijb[2] = [1.000000, 0.925053, 0.960237, 0.849408, 0.955052,
	1.281790, 1.500000, 1.000000, 0.906849, 0.897362,
	0.726255, 0.859764, 0.855134, 0.831229, 0.808310,
	0.786323, 0.765171, 1.000000, 1.000000]
beijb[3] = [1.000000, 1.022560, 0.693168, 0.946871, 1.164460,
	1.000000, 1.000000, 1.000000, 1.013060, 1.000000,
	1.005320, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
beijb[4] = [1.000000, 1.000000, 1.000000, 1.034787, 1.000000,
	1.000000, 1.000000, 1.004900, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
beijb[5] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
beijb[6] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.008692, 1.010126,
	1.011501, 1.012821, 1.014089, 1.000000, 1.000000]
beijb[7] = [1.000000, 1.100000, 1.000000, 1.300000, 1.300000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
beijb[8] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
beijb[9] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
beijb[10] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
beijb[11] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
beijb[12] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
beijb[13] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
beijb[14] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
beijb[15] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
beijb[16] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
beijb[17] = [1.000000, 1.000000, 1.000000, 1.000000]
beijb[18] = [1.000000, 1.000000, 1.000000]
beijb[19] = [1.000000, 1.000000]
beijb[20] = [1.000000]

bwij = [[0 for x in range(21)] for y in range(21)] 
bwijb = [[0 for x in range(21)] for y in range(21)] 
bwijb[0] = [1.000000, 1.000000, 0.807653, 1.000000, 1.000000,
	1.000000, 1.000000, 1.957310, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
bwijb[1] = [1.000000, 0.982746, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
bwijb[2] = [1.000000, 0.370296, 1.000000, 1.673090, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
bwijb[3] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
bwijb[4] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
bwijb[5] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
bwijb[6] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
bwijb[7] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
bwijb[8] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
bwijb[9] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
bwijb[10] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
bwijb[11] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
bwijb[12] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000, 1.000000]
bwijb[13] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000, 1.000000]
bwijb[14] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000, 1.000000]
bwijb[15] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
	1.000000]
bwijb[16] = [1.000000, 1.000000, 1.000000, 1.000000, 1.000000]
bwijb[17] = [1.000000, 1.000000, 1.000000, 1.000000]
bwijb[18] = [1.000000, 1.000000, 1.000000]
bwijb[19] = [1.000000, 1.000000]
bwijb[20] = [1.000000]
	
def paramdl(ncc, cid):
	"""
	AGA8 A.3.8.2 Subroutine PARAMDL
	Sets up constants used by the DETAIL method
	
	Shak notes:  Is this code really nessesary?  Probably not, could rewrite the rest of 
	program to eliminate, however trying to be somewhat consistent with FORTRAN reference
	"""
	#ncc = number of components
	#cid = list containing component ID numbers.  The first ncc elements of the array
	#	define the composition of the mixture.  The remaining elements of the array
	#	are set to 0.  See components dictionary for ID numbers.
	global cmwb, eib, rkib, wib, qib, hib, beijb, bkijb, bwijb, buijb, mib, dib
	global cmw, ei, rki, wi, qi, hi, beij, bkij, bwij, buij, mi, di
	
	#resize working arrays for number of components (ncc)
	di = [0.0 for x in range(ncc)]
	qi = [0.0 for x in range(ncc)]
	hi = [0.0 for x in range(ncc)]
	rki = [0.0 for x in range(ncc)]
	ei = [0.0 for x in range(ncc)]
	wi = [0.0 for x in range(ncc)]
	cmw = [0.0 for x in range(ncc)]
	mi = [0.0 for x in range(ncc)]
	beij = [[0 for x in range(ncc)] for y in range(ncc)] 
	bkij = [[0 for x in range(ncc)] for y in range(ncc)] 
	bwij = [[0 for x in range(ncc)] for y in range(ncc)] 
	buij = [[0 for x in range(ncc)] for y in range(ncc)] 
	
	#transfer constants into working arrays
	for j in range(0, ncc):
		di[j] = dib[cid[j]]
		qi[j] = qib[cid[j]]
		hi[j] = hib[cid[j]]
		rki[j] = rkib[cid[j]] 
		ei[j] = eib[cid[j]]
		wi[j] = wib[cid[j]]
		cmw[j] = cmwb[cid[j]]
		mi[j] = mib[cid[j]]
		for k in range(j,ncc):
			beij[j][k] = beijb[cid[j]][cid[k]-cid[j]]
			bkij[j][k] = bkijb[cid[j]][cid[k]-cid[j]]
			bwij[j][k] = bwijb[cid[j]][cid[k]-cid[j]]
			buij[j][k] = buijb[cid[j]][cid[k]-cid[j]]		
	
def chardl(ncc, xi):
	"""
	AGA8 A.3.8.3 Subroutine CHARDL
	Sets up composition dependant terms.  Paramdl should
	be called before calling this function
	
	Shak Notes: Why not just use math functions for all exponentials?  Maybe a legacy reason in FORTRAN.
	"""
	#ncc = number of components
	#xi = array values for mole fraction composition
	#Z_b = compressibility at 60F and 14.73 psia
	#D_b = molar density at 60F and 14.73 psia
	global a
	global b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18
	global uu, rk3p0, ww, q2p0, hh
	
	#Normalize mole fractions and calculate molar mass
	tmfrac = 0
	for j in range(ncc):
		tmfrac = tmfrac + xi[j]
	for j in range(ncc):
		xi[j] = xi[j]/tmfrac
	
	#Re-initialize variables
	b1 = 0
	b2 = 0
	b3 = 0
	b4 = 0
	b5 = 0
	b6 = 0
	b7 = 0
	b8 = 0
	b9 = 0
	b10 = 0
	b11 = 0
	b12 = 0
	b13 = 0
	b14 = 0
	b15 = 0
	b16 = 0
	b17 = 0
	b18 = 0
	
	rk5p0 =0
	rk2p5 = 0
	u5p0 = 0
	u2p5 = 0
	ww = 0
	q1p0 = 0
	hh = 0
	
	mwx = 0
	for j in range(ncc):
		mwx = mwx + xi[j]*cmw[j]
	
	for i in range(ncc):
		rk2p5 = rk2p5 + xi[i]*rki[i]**(2.5)
		u2p5 = u2p5 + xi[i]*ei[i]**(2.5)
		ww = ww + xi[i]*wi[i]
		q1p0 = q1p0 + xi[i]*qi[i]
		hh = hh + xi[i]*xi[i]*hi[i]
		for j in range(i, ncc):
			if i != j:
				xij = 2.0*xi[i]*xi[j]
			else:
				xij = xi[i]*xi[j]
			if bkij[i][j] != 1.0:
				rk5p0 = rk5p0 + xij*(bkij[i][j]**5.0 - 1.0)*(math.sqrt((rki[i]**5.0)*(rki[j]**5.0)))
			if bkij[i][j] != 1.0:
				u5p0 = u5p0 + xij*(buij[i][j]**5.0 - 1.0)*(math.sqrt((ei[i]**5.0)*(ei[j]**5.0)))
			if bkij[i][j] != 1.0:
				ww = ww = xij*(bwij[i][j] - 1.0)*((wi[i]+wi[j])/2)
			#Terms in second virial coefficeints
			eij = beij[i][j]*math.sqrt(ei[i]*ei[j])
			wij = bwij[i][j]*(wi[i]+wi[j])/2.0
			e0p5 = math.sqrt(eij)
			e2p0 = eij*eij
			e3p0 = eij*e2p0
			e3p5 = e3p0*e0p5
			e4p5 = eij*e3p5
			e6p0 = e3p0*e3p0
			e11p0 = e4p5*e4p5*e2p0
			e7p5 = e4p5*eij*e2p0
			e9p5 = e7p5*e2p0
			e12p0 = e11p0*eij
			e12p5 = e12p0*e0p5
			s3 = xij*math.sqrt((rki[i]**3)*(rki[j]**3))
			b1 = b1 + s3
			b2 = b2 + s3*e0p5
			b3 = b3 + s3*eij
			b4 = b4 + s3*e3p5
			b5 = b5 + s3*wij/e0p5
			b6 = b6 + s3*wij*e4p5
			b7 = b7 + s3*qi[i]*qi[j]*e0p5
			b8 = b8 + s3*mi[i]*mi[j]*e7p5
			b9 = b9 + s3*mi[i]*mi[j]*e9p5
			b10 = b10 + s3*di[i]*di[j]*e6p0
			b11 = b11 + s3*di[i]*di[j]*e12p0
			b12 = b12 + s3*di[i]*di[j]*e12p5
			b13 = b13 + s3*hi[i]*hi[j]/e6p0
			b14 = b14 + s3*e2p0
			b15 = b15 + s3*e3p0
			b16 = b16 + s3*qi[i]*qi[j]*e2p0
			b17 = b17 + s3*e2p0
			b18 = b18 + s3*e11p0
	
	b1 = b1*a[1]
	b2 = b2*a[2]
	b3 = b3*a[3]
	b4 = b4*a[4]
	b5 = b5*a[5]
	b6 = b6*a[6]
	b7 = b7*a[7]
	b8 = b8*a[8]
	b9 = b9*a[9]
	b10 = b10*a[10]
	b11 = b11*a[11]
	b12 = b12*a[12]
	b13 = b13*a[13]
	b14 = b14*a[14]
	b15 = b15*a[15]
	b16 = b16*a[16]
	b17 = b17*a[17]
	b18 = b18*a[18]
	
	rk3p0 = (rk5p0 + rk2p5*rk2p5)**0.6
	uu = (u5p0 + u2p5*u2p5)**0.2
	q2p0 = q1p0*q1p0
	
	#Z_b = compressibility factor at base conditions, T_b and P_b
	#Base conditions = 60F, 14.73psia
	T_b = (60.0 + 459.67)/1.8
	P_b = 14.73*6894.757/100000.0
	D_b = ddetail(P_b, T_b)
	Z_b = zdetail(T_b, D_b)

def B(T):
	"""
	AGA8 A.3.8.4 Subroutine B
	Calculates the second viral coefficient using the AGA8 
	method given temperature.
	"""
	#T = Temperature in kelvins.
	#bmix = Second virial coefficient in mol/dm^3
	global b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18
	global bmix
	
	t0p5 = math.sqrt(T)
	t2p0 = T*T
	t3p0 = T*t2p0
	t3p5 = t3p0*t0p5
	t4p5 = T*t3p5
	t6p0 = t3p0*t3p0
	t11p0 = t4p5*t4p5*t2p0
	t7p5 = t6p0*T*t0p5
	t9p5 = t7p5*t2p0
	t12p0 = t9p5*t0p5*t2p0
	t12p5 = t12p0*t0p5
	bmix = (b1 + b2/t0p5 + b3/T + b4/t3p5 + b5*t0p5 + b6/t4p5 + b7/t0p5 + b8/t7p5 + b9/t9p5 + b10/t6p0 
		+ b11/t12p0 + b12/t12p5 + b13*t6p0 + b14/t2p0 + b15/t3p0 + b16/t2p0 + b17/t2p0 + b18/t11p0)
	
def temp(T):
	"""
	AGA8 A.3.8.5 Subroutine TEMP
	Sets up temperature dependent terms
	"""
	#T = Temperature in kelvins
	global a
	global uu, rk3p0, ww, q2p0, hh, bmix
	global fn
	
	B(T)

	tr = T/uu
	tr0p5 = math.sqrt(tr)
	tr1p5 = tr*tr0p5
	tr2p0 = tr*tr
	tr3p0 = tr*tr2p0
	tr4p0 = tr*tr3p0
	tr5p0 = tr*tr4p0
	tr6p0 = tr*tr5p0
	tr7p0 = tr*tr6p0
	tr8p0 = tr*tr7p0
	tr9p0 = tr*tr8p0
	tr11p0 = tr6p0*tr5p0
	tr13p0 = tr6p0*tr7p0
	tr21p0 = tr9p0*tr9p0*tr3p0
	tr22p0 = tr*tr21p0
	tr23p0 = tr*tr22p0
	
	fn[12] = a[12]*hh*tr6p0
	fn[13] = a[13]/tr2p0
	fn[14] = a[14]/tr3p0
	fn[15] = a[15]*q2p0/tr2p0
	fn[16] = a[16]*tr2p0
	fn[17] = a[17]*tr11p0
	fn[18] = a[18]*tr0p5
	fn[19] = a[19]/tr0p5
	fn[20] = a[20]
	fn[21] = a[21]/tr4p0
	fn[22] = a[22]/tr6p0
	fn[23] = a[23]/tr21p0
	fn[24] = a[24]*ww/tr23p0
	fn[25] = a[25]*q2p0/tr22p0
	fn[26] = a[26]*hh*tr
	fn[27] = a[27]*q2p0*tr0p5
	fn[28] = a[28]*ww/tr7p0
	fn[29] = a[29]*hh*tr
	fn[30] = a[30]/tr6p0
	fn[31] = a[31]*ww/tr4p0
	fn[32] = a[32]*ww/tr
	fn[33] = a[33]*ww/tr9p0
	fn[34] = a[34]*hh*tr13p0
	fn[35] = a[35]/tr21p0
	fn[36] = a[36]*q2p0/tr8p0
	fn[37] = a[37]*tr0p5
	fn[38] = a[38]
	fn[39] = a[39]/tr2p0
	fn[40] = a[40]/tr7p0
	fn[41] = a[41]*q2p0/tr9p0
	fn[42] = a[42]/tr22p0
	fn[43] = a[43]/tr23p0
	fn[44] = a[44]/tr
	fn[45] = a[45]/tr9p0
	fn[46] = a[46]*q2p0/tr3p0
	fn[47] = a[47]/tr8p0
	fn[48] = a[48]*q2p0/tr23p0
	fn[49] = a[49]/tr1p5
	fn[50] = a[50]*ww/tr5p0
	fn[51] = a[51]*q2p0*tr0p5
	fn[52] = a[52]/tr4p0
	fn[53] = a[53]*ww/tr7p0
	fn[54] = a[54]/tr3p0
	fn[55] = a[55]*ww
	fn[56] = a[56]/tr
	fn[57] = a[57]*q2p0
		
def ddetail(P, T):
	"""
	AGA8 A.3.8.6 Function Subprogram DDETAIL
	"""
	pass

def density_search(IT, rho1, rho2, P1, P2, T)
	"""
	Broke density search loop out of BRACKET subroutine to eliminate GOTOs
	"""
	global rgas
	global uu, rk3p0, ww, q2p0, hh, bmix
	global delta
	
	code = 0
	imax = 200	#max number of iterations
	rhomax = 1.0/rk3p0
	if T > 1.2593*uu:
		rhomax = 20.0*rhomax
	
	if IT > imax:	#Maximum number of iterations reached
		code = 3
		rho_l = rho1
		rho_h = rho2
		p_rho_l = P1
		p_rho_h = P2
		print("Code=3: Maximum number of iterations in BRAKET exceeded, default density used")
		return code, rho_l, rho_h, p_rho_l, p_rho_h
		
	if (code != 2) and rho2 > rhomax:	#Density in braket exceeds maximum allowable desnity
		code = 2
		delta = 0.1*(rhomax - rho1) + P/(rgas*t)/20.0
		print("Code=2: Density in BRAKET exceeds maximum")
		rho_l, rho_h, p_rho_l, p_rho_h = density_search(IT+1, rho1, rho1+delta, P1, P2, T)
	
	#calculate pressure P2 at density rho2
	P2 = pdetail(rho2, T)
	
	#test value of P2 relative to P and relative to P1
	if P2 > P:	#the density root is bracketed
		rho_l = rho1
		p_rho_l = P1
		rho_h = rho2
		p_rho_h = P2
		return code, rho_l, rho_h, p_rho_l, p_rho_h
	elif P2 > P1 and code == 2:
		rho1 = rho2
		P1 = P2
		rho_l, rho_h, p_rho_l, p_rho_h = density_search(IT+1, rho1, rho1+delta, P1, P2, T)
	elif P2 > P1 and code == 0:
		delta = 2.0*delta
		rho1 = rho2
		P1 = P2
		

def braket(T, P):
	"""
	AGA8 A.3.8.7 Subroutine BRAKET
	
	Shak Notes: Untangled spaghetti code by breaking out density search routine.  This was GOTO hell
	"""
	#returns code, rho, rho_l, rho_h, p_rho_l, p_rho_h
	global rgas
	global uu, rk3p0, ww, q2p0, hh, bmix
	global delta
	
	#set initial values for rho1, P1, rho2, P2
	rho1 = 0.0
	P1 = 0.0
	videal = rgas*T/P
	if math.abs(bmix) < 0.167*videal:
		rho2 = 0.95/(videal + bmix)
	else:
		rho2 = 1.15/videal
	P2 = 0.0
	
	#run iterative density search
	code, rho, rho_l, rho_h, p_rho_l, p_rho_h = density_search(0, rho1, rho2, P1, P2, T)
	
	return code, rho, rho_l, rho_h, p_rho_l, p_rho_h
	
def pdetail(D, T):
	"""
	AGA8 A.3.8.8 Function Subprogram PDETAIL
	"""
	global rgas
	
	pdetail = zdetail(D, T)*D*rgas*T
	
	return pdetail
	
def zdetail(D, T):
	"""
	AGA8 A.3.8.9 Function Subprogram ZDETAIL
	
	Calculates the compressibility factor from the AGA8
	as a function of density and temperature.
	"""
	#D = Molar density in mol/dm^3
	#T = Temperature in kelvins
	#returns Z = compressibility factor
	global fn
	global uu, rk3p0, ww, q2p0, hh, bmix
	
	if T != told:
		temp(T)
	told = T
	D1 = rk3p0*D1
	D2 = D1*D1
	D3 = D2*D1
	D4 = D3*D1
	D5 = D4*D1
	D6 = D5*D1
	D7 = D6*D1
	D8 = D7*D1
	D9 = D8*D1
	
	exp1 = math.exp(-1*D1)
	exp2 = math.exp(-1*D2)
	exp3 = math.exp(-1*D3)
	exp4 = math.exp(-1*D4)
	
	Z = (1.0 + bmix*D
		+ fn[12]*D1*(exp3 - 1.0 - 3.0*D3*exp3)
		+ (fn[13] + fn[14] + fn[15])*D1*(exp2 - 1.0 - 2.0*D2*exp2)
		+ (fn[16] + fn[17])*D1*(exp4 - 1.0 - 4.0*D4*exp4)
		+ (fn[18] + fn[19])*D2*2.0
		+ (fn[20] + fn[21] + fn[22]*D2*(2.0 - 2.0*D2)*exp2)
		+ (fn[23] + fn[24] + fn[25]*D2*(2.0 - 4.0*D4)*exp4)
		+ fn[26]*D2*(2.0 - 4.0*D4)*exp4
		+ fn[27]*D3*3.0
		+ (fn[28] + fn[29])*D3*(3.0 - D1)*exp1
		+ (fn[30] + fn[31])*D3*(3.0 - 2.0*D2)*exp2
		+ (fn[32] + fn[33])*D3*(3.0 - 3.0*D3)*exp3
		+ (fn[34] + fn[35] + fn[36])*D3*(3.0 - 4.0*D4)*exp4
		+ (fn[37] + fn[38])*D4*4.0
		+ (fn[39] + fn[40] + fn[41])*D4*(4.0 - 2.0*D2)*exp2
		+ (fn[42] + fn[43])*D4*(4.0 - 4.0*D4)*exp4
		+ fn[44]*D5*5.0
		+ (fn[45] + fn[46])*D5*(5.0 - 2.0*D2)*exp2
		+ (fn[47] + fn[48])*D5*(5.0 - 4.0*D4)*exp4
		+ fn[49]*D6*6.0
		+ fn[50]*D6*(6.0 - 2.0*D2)*exp2
		+ fn[51]*D7*7.0
		+ fn[52]*D7*(7.0 - 2.0*D2)*exp2
		+ fn[53]*D8*(8.0 - D1)*exp1
		+ (fn[54] + fn[55])*D8*(8.0 - 2.0*D2)*exp2
		+ (fn[56] + fn[57])*D8*(9.0 - 2.0*D2)*exp2)
	
	return Z
