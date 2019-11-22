import copy
import re
import os
import subprocess
import math
import numpy as np
from scipy import optimize
global POSCAR

global POSCAR
def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def read_POSCAR():
	fin = open('POSCAR','r')
	global POSCAR
	POSCAR = fin.readlines()
	fin.close()
	a=float(getnum(POSCAR[2])[0])
	b=float(getnum(POSCAR[3])[1])
	return a,b

def write_POSCAR():
	global POSCAR
	fout = open('POSCAR','w')
	for line in POSCAR:
		fout.write(line)
	fout.close()

def modify(name,k):
	global POSCAR
	dx = 0.01
	N=10
	os.chdir(name)
	axis = list(read_POSCAR())
	ab = copy.copy(axis)
	for i in range(2*N):
		print(ab[k])
		os.chdir(str(i))
		ab[k] = round(axis[k]+i*dx,3)
		POSCAR[2]=str(ab[0])+' 0.0 0.0\n'
		POSCAR[3]='0.0 '+str(ab[1])+' 0.0\n'
		write_POSCAR()
		os.chdir('..')
	os.system('bsub<series.sh')
	os.chdir('..')

modify('VOF',1)



