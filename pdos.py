import re
import os
import subprocess
from numpy import sqrt,array,zeros

from subprocess import Popen, PIPE, call
global ATOMS,NGRID,EF,E,DOS

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def correction():
	band=float(getnum(subprocess.getoutput('grep E-fermi '+'OUTCAR'))[0])
	return -band

def read_DOSCAR():
	global ATOMS,NGRID,EF,E,DOS
	fin = open('DOSCAR','r')
	lines = fin.readlines()
	fin.close()
	ATOMS = int(getnum(lines[0])[0])
	NGRID = int(getnum(lines[5])[2])
	EF = float(getnum(lines[5])[3]) 
	E = []
	tot = []
	for i in range(NGRID):
		E.append(float(getnum(lines[i+6])[0])-EF)
		tot.append(float(getnum(lines[i+6])[1])*10**int(getnum(lines[i+6])[2]))
	DOS = []
	DOS.append(tot)
	for i in range(ATOMS):
		s = []
		py = []
		pz = []
		px = []
		dxy = []
		dyz = []
		dz2 = []
		dxz = []
		dx2 = []
		t = []
		for j in range(NGRID):
			k = (i+1)*(NGRID+1)+6+j
			s.append(float(getnum(lines[k])[1])*10**int(getnum(lines[k])[2]))
			py.append(float(getnum(lines[k])[3])*10**int(getnum(lines[k])[4]))
			pz.append(float(getnum(lines[k])[5])*10**int(getnum(lines[k])[6]))
			px.append(float(getnum(lines[k])[7])*10**int(getnum(lines[k])[8]))
			dxy.append(float(getnum(lines[k])[9])*10**int(getnum(lines[k])[10]))
			dyz.append(float(getnum(lines[k])[11])*10**int(getnum(lines[k])[12]))
			dz2.append(float(getnum(lines[k])[13])*10**int(getnum(lines[k])[14]))
			dxz.append(float(getnum(lines[k])[15])*10**int(getnum(lines[k])[16]))
			dx2.append(float(getnum(lines[k])[17])*10**int(getnum(lines[k])[18]))
			t.append(s[-1]+py[-1]+pz[-1]+px[-1]+dxy[-1]+dyz[-1]+dz2[-1]+dxz[-1]+dx2[-1])
		DOS.append([s,py,pz,px,dxy,dyz,dz2,dxz,dx2,t])

def calc_Atom(a,b):
	temp = zeros(NGRID)
	for i in range(a,b+1):
		temp += array(DOS[i][-1])
	return temp

def calc_Orbit(a,b):
	temp = array(DOS[b])
	for i in range(a,b):
		temp+=array(DOS[i])
	return temp


def generate_A():
	Cu = calc_Atom(1,3)
	O = calc_Atom(4,9)
	C = calc_Atom(10,15)
	TOT = calc_Atom(1,15)
	for i in range(NGRID):
		print(E[i],DOS[0][i],Cu[i],O[i],C[i],TOT[i])

def generate_P():
	atom = calc_Orbit(1,15)
	for i in range(NGRID):
		print(E[i],end=' ')
		s = atom[0]
		p = atom[1]+atom[2]+atom[3]
		d = atom[4]+atom[5]+atom[6]+atom[7]+atom[8]
		#for val in atom:
		#	print(val[i],end=' ')
		print(s[i],p[i],d[i],atom[9][i])
read_DOSCAR()
generate_A()
#generate_P()