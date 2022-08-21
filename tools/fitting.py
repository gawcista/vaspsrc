import numpy as np
#import matplotlib.pyplot as plt
from scipy import optimize
import re

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)
def ReadXY():

	fin = open('energy.dat','r')
	lines = fin.readlines()
	fin.close()
	x = []
	y = []
	for line in lines:
		x.append(float(getnum(line)[0]))
		y.append(float(getnum(line)[1]))
	return np.array(x),np.array(y)

def read_Area(k):
	fin = open('POSCAR','r')
	global POSCAR
	POSCAR = fin.readlines()
	fin.close()
	vector=[]
	vector.append(float(getnum(POSCAR[2])[0]))
	vector.append(float(getnum(POSCAR[3])[1]))
	vector.append(float(getnum(POSCAR[4])[2]))
	return vector[int(k)]*vector[2],vector[0],vector[1]

def Birch_Murnaghan(V,E0,V0,B0,B1):
	return E0+(9.*V0*B0/16.)*(((V0/V)**(2/3)-1)**3 *B1 +((V0/V)**(2/3)-1)**2 * (6-4* (V0/V)**(2/3)))
'''
def plotAll(X1,Y1,X2,Y2):
	plt.figure()
	plt.title('Birch Murnaghan Fitting')
	plt.xlabel('Volume (Angstrom^3)')  
	plt.ylabel('Energy (eV)') 
	plt.plot(X1,Y1)
	plt.plot(X2,Y2)
	plt.show()
'''
X,Y=ReadXY()
para,cov=optimize.curve_fit(Birch_Murnaghan,X,Y,[np.min(Y),X[len(X)>>2],1,1])
X1 = np.arange(X[0],X[len(X)-1],0.1)
#print('E0 =',para[0])
#k = input("Axis=\n")
print(para[1])
#print('B0 =',para[2])
#print( "B0'=",para[3])
#plotAll(X,Y,X1,Birch_Murnaghan(X1,para[0],para[1],para[2],para[3]))
