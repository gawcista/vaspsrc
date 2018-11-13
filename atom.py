import re
import os
import subprocess
from numpy import sqrt,array,zeros

from subprocess import Popen, PIPE, call

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def correction():
	band=float(getnum(subprocess.getoutput('grep E-fermi '+'OUTCAR'))[0])
	return -band

def read_PROCAR(file):
	fin = open(file,'r')
	lines = fin.readlines()
	fin.close()
	print(lines[1])
    kpoints = int(getnum(lines[1])[0])
	bands = int(getnum(lines[1])[1])
	natoms = int(getnum(lines[1])[2])
	K = []
	ENERGY = []
	DOS = []
	dE = correction()
	for i in range(KPOINTS):
		line_K = (bands*(natoms+5)+3)*i+3
		print(lines[line_K])

read_PROCAR('PROCAR')