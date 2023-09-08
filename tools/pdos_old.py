import re
import os
import subprocess
from numpy import sqrt,array,zeros

from subprocess import Popen, PIPE, call
global ATOMS,NGRID,EF,E,DOS_up,DOS_down
def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def read_DOSCAR():
	global ATOMS,NGRID,EF,E,DOS_up,DOS_down
	fin = open('DOSCAR','r')
	lines = fin.readlines()
	fin.close()
	ATOMS = int(lines[0].split()[0])
	NGRID = int(lines[5].split()[2]) 
	EF = float(lines[5].split()[3])
	E = []
	tot = []
	for i in range(NGRID):
		E.append(float(getnum(lines[i+6])[0])-EF)
		tot.append(float(getnum(lines[i+6])[1])*10**int(getnum(lines[i+6])[2]))
	DOS_up = []
	DOS_down = []
	DOS_up.append(tot)
	DOS_down.append(tot)
	for i in range(ATOMS):
		s_up = [] 
		s_down = []
		py_up = [] 
		py_down = []
		pz_up = [] 
		pz_down = []
		px_up = [] 
		px_down = []
		dxy_up = [] 
		dxy_down = []
		dyz_up = [] 
		dyz_down = []
		dz2_up = [] 
		dz2_down = []
		dxz_up = [] 
		dxz_down = []
		dx2_up = [] 
		dx2_down = []
		t_up = []
		t_down = []
		for j in range(NGRID):
			k = (i+1)*(NGRID+1)+6+j
			s_up.append(float(getnum(lines[k])[1])*10**int(getnum(lines[k])[2]))
			s_down.append(float(getnum(lines[k])[3])*10**int(getnum(lines[k])[4]))
			py_up.append(float(getnum(lines[k])[5])*10**int(getnum(lines[k])[6]))
			py_down.append(float(getnum(lines[k])[7])*10**int(getnum(lines[k])[8]))
			pz_up.append(float(getnum(lines[k])[9])*10**int(getnum(lines[k])[10]))
			pz_down.append(float(getnum(lines[k])[11])*10**int(getnum(lines[k])[12]))
			px_up.append(float(getnum(lines[k])[13])*10**int(getnum(lines[k])[14]))
			px_down.append(float(getnum(lines[k])[15])*10**int(getnum(lines[k])[16]))
			dxy_up.append(float(getnum(lines[k])[17])*10**int(getnum(lines[k])[18]))
			dxy_down.append(float(getnum(lines[k])[19])*10**int(getnum(lines[k])[20]))
			dyz_up.append(float(getnum(lines[k])[21])*10**int(getnum(lines[k])[22]))
			dyz_down.append(float(getnum(lines[k])[23])*10**int(getnum(lines[k])[24]))
			dz2_up.append(float(getnum(lines[k])[25])*10**int(getnum(lines[k])[26]))
			dz2_down.append(float(getnum(lines[k])[27])*10**int(getnum(lines[k])[28]))
			dxz_up.append(float(getnum(lines[k])[29])*10**int(getnum(lines[k])[30]))
			dxz_down.append(float(getnum(lines[k])[31])*10**int(getnum(lines[k])[32]))
			dx2_up.append(float(getnum(lines[k])[33])*10**int(getnum(lines[k])[34]))
			dx2_down.append(float(getnum(lines[k])[35])*10**int(getnum(lines[k])[36]))
			#t_up.append(s_up[-1]+pz_up[-1]+dyz_up[-1]+dz2_up[-1]+dxz_up[-1])
			#t_down.append(s_down[-1]+pz_down[-1]+dyz_down[-1]+dz2_down[-1]+dxz_down[-1])
			t_up.append(s_up[-1]+py_up[-1]+pz_up[-1]+px_up[-1]+dxy_up[-1]+dyz_up[-1]+dz2_up[-1]+dxz_up[-1]+dx2_up[-1])
			t_down.append(s_down[-1]+py_down[-1]+pz_down[-1]+px_down[-1]+dxy_down[-1]+dyz_down[-1]+dz2_down[-1]+dxz_down[-1]+dx2_down[-1])
		DOS_up.append([s_up,py_up,pz_up,px_up,dxy_up,dyz_up,dz2_up,dxz_up,dx2_up,t_up])
		DOS_down.append([s_down,py_down,pz_down,px_down,dxy_down,dyz_down,dz2_down,dxz_down,dx2_down,t_down])

def calc_Atom(a,b):
	temp_up = zeros(NGRID)
	temp_down = zeros(NGRID)
	for i in range(a,b+1):
		temp_up += array(DOS_up[i][-1])
		temp_down += array(DOS_down[i][-1])
	return temp_up,temp_down

def calc_Orbit(a,b):
	temp_up = array(DOS_up[b])
	temp_down = array(DOS_down[b])
	for i in range(a,b):
		temp_up += array(DOS_up[i])
		temp_down += array(DOS_down[i])
	return temp_up,temp_down


def generate_A():
	#H = calc_Atom(1,2)
	#C = calc_Atom(3,8)
	#O = calc_Atom(9,12)
	#HITP_up,HITP_down = calc_Atom(1,72)
	H_up,H_down = calc_Atom(1,24)
	C_up,C_down = calc_Atom(25,60)
	N_up,N_down = calc_Atom(61,72)
	Fe_up,Fe_down = calc_Atom(73,75)
	for i in range(NGRID):
		#print(E[i],H_up[i],-H_down[i],C_up[i],-C_down[i],N_up[i],-N_down[i],Fe_up[i],-Fe_down[i])
		print(E[i],H_up[i]+H_down[i],C_up[i]+C_down[i],N_up[i]+N_down[i],Fe_up[i]+Fe_down[i])
		#print(E[i],DOS[0][i],M[i],Fe[i],TOT[i])

def generate_P():
	atom = calc_Orbit(1,13)
	for i in range(NGRID):
		print(E[i],end=' ')
		s = atom[0]
		p = atom[1]+atom[2]+atom[3]
		d = atom[4]+atom[5]+atom[6]+atom[7]+atom[8]
		#for val in atom:
		#	print(val[i],end=' ')
		print(s[i],p[i],d[i],atom[9][i])
def calc_d_orbital(ion1,ion2):
	atom_up,atom_down = calc_Orbit(ion1,ion2)
	print(E[0],0,0,0,0,0)
	for i in range(NGRID):
		if i>0:
			print(E[i],atom_up[4][i],atom_up[5][i],atom_up[6][i],atom_up[7][i],atom_up[8][i])
		#print(E[i],atom_up[4][i]+atom_down[4][i],atom_up[5][i]+atom_down[5][i],atom_up[6][i]+atom_down[6][i],atom_up[7][i]+atom_down[7][i],atom_up[8][i]+atom_down[8][i])
		#print(E[i],atom_up[4][i],-atom_down[4][i],atom_up[5][i],-atom_down[5][i],atom_up[6][i],-atom_down[6][i],atom_up[7][i],-atom_down[7][i],atom_up[8][i],-atom_down[8][i])
	print(E[-1],0,0,0,0,0)
	print(E[0],0,0,0,0,0)
	for i in range(NGRID):
		if i>0:
			print(E[i],-atom_down[4][i],-atom_down[5][i],-atom_down[6][i],-atom_down[7][i],-atom_down[8][i])
read_DOSCAR()
#generate_A()
#generate_P()
#calc_d_orbital(73,73)