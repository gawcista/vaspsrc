import copy
import re
import os
import subprocess
import math
import numpy as np
from scipy import optimize
global POSCAR

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def Birch_Murnaghan(V,E0,V0,B0,B1):
	return E0+(9.*V0*B0/16.)*(((V0/V)**(2/3)-1)**3 *B1 +((V0/V)**(2/3)-1)**2 * (6-4* (V0/V)**(2/3)))

def read_POSCAR():
	fin = open('POSCAR','r')
	global POSCAR
	POSCAR = fin.readlines()
	fin.close()
	a=float(getnum(POSCAR[2])[0])
	b=float(getnum(POSCAR[3])[1])
	return a,b

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
	return vector[1-k]*vector[2],vector[0],vector[1]
	

def write_POSCAR():
	global POSCAR
	fout = open('POSCAR','w')
	for line in POSCAR:
		fout.write(line)
	fout.close()

def write_series(name):
	fout = open('run.sh','w')
	fout.write('#!/bin/sh\n#BSUB -q all\n#BSUB -o %J.out\n#BSUB -R "span[ptile=24]"\n#BSUB -n 48\n')
	fout.write('#BSUB -J '+name+'_s\n')
	fout.write('cd $LS_SUBCWD\n')
	fout.write('source /etc/profile.d/intel.sh\n')
	fout.write('VASP="mpirun /share/soft/vasp_intelmpi/bin/vasp5.4.1_Wannier90/vasp_std"\n')
	fout.write('for ((x=0;x<=19;x=x+1))\ndo\n')
	fout.write('cd $x\n$VASP >> log\n')
	fout.write('e=`grep TOTEN OUTCAR|tail -1`\n')
	fout.write('g=`grep volume OUTCAR|tail -1`\n')
	fout.write('echo ${g:18:11} ${e:25:18}>> ../energy.dat\n')
	fout.write('cd ..\n done\n')
	fout.close()

def optimze(name,k):
	global POSCAR
	dx = 0.01
	N=10
	os.chdir(name)
	axis = list(read_POSCAR())
	ab = copy.copy(axis)
	for i in range(-N,N):
		ab[k] = round(axis[k]+i*dx,3)
		POSCAR[2]=str(ab[0])+' 0.0 0.0\n'
		POSCAR[3]='0.0 '+str(ab[1])+' 0.0\n'
		os.system('rm '+str(i+N))
		os.system('mkdir '+str(i+N))
		os.chdir(str(i+N))
		os.system('cp ../INCAR .')
		os.system('cp ../POTCAR .')
		os.system('cp ../KPOINTS KPOINTS')
		write_POSCAR()
		#os.system('mpirun /share/soft/vasp_intelmpi/bin/vasp5.4.1_Wannier90/vasp_std >>log')
		#energy = float(getnum(subprocess.getoutput('grep TOTEN OUTCAR|tail -1'))[0])
		#volume = float(getnum(subprocess.getoutput('grep volume OUTCAR|tail -1'))[0])
		os.chdir('..')
		#os.system('echo '+str(volume)+' '+str(energy)+' >>energy.dat')	
	os.system('bsub<series.sh')
	os.chdir('..')	

def fitting(name,k):
	global POSCAR
	dx = 0.01
	N=10
	os.chdir(name)
	X,Y=ReadXY()
	para,cov=optimize.curve_fit(Birch_Murnaghan,X,Y,[np.min(Y),X[len(X)>>2],1,1])
	X1 = np.arange(X[0],X[len(X)-1],0.1)
	vector=list(read_Area(k))
	vector[1+k]=para[1]/vector[0]
	#vector[1+k]=(para[1]-X[0])/(X[len(Y)-1]-X[0])
	#print(name,para[1],para[1]/vector[0],X[0],X[len(X)-1],vector[1+k])
	print(name,vector[1+k],(para[1]-X[0])/(X[len(Y)-1]-X[0]))
	
	print('Writing POSCAR',name,'with axis',k)
	POSCAR[2]=str(vector[1])+' 0.0 0.0\n'
	POSCAR[3]='0.0 '+str(vector[2])+' 0.0\n'
	os.system('cp POSCAR POSCAR2')
	write_POSCAR()
	os.system('bsub<run.sh')
	os.chdir('..')

def read_INCAR(name):
	os.chdir(name)
	fin = open('INCAR','r')
	INCAR = fin.readlines()
	fin.close()
	INCAR[3]='   ISTART = 0 ; ICHARG = 11\n'
	INCAR[14]='   LCHARG = .TRUE.\n'
	print(name,INCAR[3])
	fout = open('INCAR.band','w')
	for i in range(22):
		fout.write(INCAR[i])
	fout.close()
	#print(name,INCAR[7])
	#INCAR[7]='   EDIFF = 1E-8\n'
	#fout = open('INCAR','w')
	#for line in INCAR:
	#	fout.write(line)
	#fout.close()
	os.chdir('..')

def terminal():
	M = ['Co','Cr','Fe','Mn','Ni','Sc','Ti','V']
	O = ['N','O']
	X = ['Br','Cl','F']
	axis = 1
	for x in M:
		for y in O:
			for z in X:
				name = x+y+z
				read_INCAR(name)
				#if (name=='FeOF'):
			#		print('Optimizing '+name+' with axis '+str(axis) )
		#			optimze(name,axis)
			 	#print('Fitting '+name+' with axis '+str(axis))
				#fitting(name,axis) 
			#	if os.path.exists(name+'/energy.dat'):
			#		fitting(name,axis) 

terminal()
