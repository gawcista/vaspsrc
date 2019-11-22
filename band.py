import re
import os
import subprocess
import math

from subprocess import Popen, PIPE, call

global SPACE,NK,NBAND,ISPIN,SPIN_MODE,DIR,HEAD,FERMI
global E_fermi,Band,Energy,Occupancy,Kpoints,Special,Route

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def distance(a,b):
	dis = .0
	for i in range(3):
		dis += (a[i]-b[i])**2
	return math.sqrt(dis)

def initialize(): # initialize constants and global variables
	global SPACE,NK,NBAND,ISPIN,SPIN_MODE,DIR,HEAD,FERMI
	global E_fermi,Band,Energy,Occupancy,Kpoints,Special,Route
	NK = 0
	NBAND = 120
	SPIN_MODE = 0
	DIR = ''
	HEAD = 8
	FERMI = '0'
	SPACE = ' '
	ISPIN = 1
	E_fermi = 0.
	Band = []
	Energy = []
	Occupancy = []
	Kpoints = []
	Special = []
	Route = ['{/Symbol G}','M','K','{/Symbol G}']

def set_DIR():
	mode = ''
	global DIR
	while(mode!='1' and mode!='2'):
		mode = input("Please choose mode(1:filename;2:path)\n")
		if (mode=='1'):
			print('Filename mode chosen!') # developing
		if (mode=='2'):
			print('Path mode chosen!')
			global DIR
			DIR = input("Input the calculation directory:\n")
			if (len(DIR)!=0 and DIR[len(DIR)-1]!='/'):
				DIR += '/'
		if (mode!='1' and mode!='2'):
			print('Error input!')

def set_SPIN():
	global SPIN_MODE
	while(mode!='1' and mode!='2'):
		mode = input("Please choose output spin mode(1:Up; 2:Down)\n")
		if (mode=='1'):
			print('Filename mode chosen!')
		if (mode=='2'):
			print('Path mode chosen!')
		if (mode!='1' and mode!='2'):
			print('Error input!')
	SPIN_MODE = int(mode)-1

def correction():
	global FERMI,DIR,E_fermi
	band=float(getnum(subprocess.getoutput('grep E-fermi '+DIR+'band/OUTCAR'))[0])
	if FERMI == '0':
		E_fermi = 0.
		return -band
	else: 
		E_fermi = band
		scf=float(getnum(subprocess.getoutput('grep E-fermi '+DIR+'SCF/OUTCAR'))[0])
		return scf-band

def read_EIGENVAL():
	global NK,NBAND,HEAD,ISPIN,DIR,Energy,Occupancy,Kpoints,Special
	fin = open(DIR+'band/EIGENVAL','r')
	lines = fin.readlines()
	fin.close()
	NK = int(getnum(lines[5])[1])
	NBAND = int(getnum(lines[5])[2])
	ISPIN =int(getnum(subprocess.getoutput('grep ISPIN '+DIR+'band/OUTCAR'))[0])
	os.system('grep -A '+str(NK)+' 2pi/SCALE '+DIR+'band/OUTCAR >'+DIR+'kpoints.temp')
	fin = open(DIR+'kpoints.temp','r')
	klines = fin.readlines()
	fin.close()
	os.system('rm '+DIR+'kpoints.temp')
	del klines[0]
	dE = correction()
	route = 0.
	Special.append(0.)
	print(NK)
	for i in range(NK):
		Et = [[],[]]
		Ot = [[],[]]
		for j in range(NBAND):
			ls = i*(NBAND+2)+HEAD
			for spin in range(ISPIN):			
				Et[spin].append(float(getnum(lines[ls+j])[spin+1])+dE)
				#Ot[spin].append(float(getnum(lines[ls+j])[spin+3]))
		Energy.append(Et)
		Occupancy.append(Ot)

		point = getnum(klines[i])
		ktemp = [0,1,2]
		for j in ktemp:
			ktemp[j] = float(point[j])
		if i>0:
			route += distance(ktemp,Kpoints[i-1])
			ktemp.append(route)
			if distance(ktemp,Kpoints[i-1]) == 0:
				Special.append(route)
		else:
			ktemp.append(0)
		Kpoints.append(ktemp)
	Special.append(Kpoints[NK-1][3])

def print_SPIN(OUTPUT_SPIN=0):
	global NK,NBAND,Energy,Occupancy,Kpoints
	for i in range(NK):
		print(Kpoints[i][3],end=' ')
		for j in range(NBAND):
			print(Energy[i][OUTPUT_SPIN][j],end=' ')
		print()

def print2File(OUTPUT_SPIN,filename,dir=''):
	global NK,NBAND,Energy,SPACE,Occupancy,Kpoints
	fout = open(dir+filename,'w')
	for i in range(NK):
		fout.write(str(Kpoints[i][3])+SPACE)
		for j in range(NBAND):
			fout.write(str(Energy[i][OUTPUT_SPIN][j])+SPACE)
		fout.write('\n')
	fout.close()

def writePlot(y0,y1,color,inputfilename,outputfilename='output.eps',dir='',filename='plot_test'):
	global NK,Kpoints

	fout = open(dir+filename,'w')
	fout.write('set term post eps color enhanced "Helvetica" 20\n')
	fout.write('set output "'+outputfilename+'"\n')
	fout.write('unset key\n')
	fout.write('set ylabel "Energy (eV)"\n')
	fout.write('set xrange [0:'+str(Kpoints[NK-1][3])+']\n')
	fout.write('set yrange ['+str(y0)+':'+str(y1)+']\n')

	fout.write('set xtics (')
	for i in range(len(Special)):
		fout.write('"'+Route[i]+'" '+str(Special[i]))
		if i < len(Special)-1:
			fout.write(', ')
		else:
			fout.write(')\n')

	for i in range(len(Special)-2):
		fout.write('set arrow from '+str(Special[i+1])+','+str(y0)+' to '+str(Special[i+1])+','+str(y1)+' nohead lt 0 lw 2\n')
	
	fout.write('plot ')
	for i in range(NBAND-1):
		fout.write("'"+inputfilename+"'"+' u 1:'+str(i+2)+' w l lc '+str(color)+' lw 4,\\\n')
	fout.write("'"+inputfilename+"'"+' u 1:'+str(NBAND+1)+' w l lc '+str(color)+ ' lw 4\n')

	fout.close()

def plotSpin(y0,y1,color1,color2,spinf1,spinf2,outputfilename='output.eps',dir='',filename='plot_test'):
	global NK,Kpoints
	fout = open(dir+filename,'w')
	fout.write('set term post eps color enhanced "Helvetica" 20\n')
	fout.write('set output "'+outputfilename+'"\n')
	fout.write('unset key\n')
	fout.write('set ylabel "Energy (eV)"\n')
	fout.write('set xrange [0:'+str(Kpoints[NK-1][3])+']\n')
	fout.write('set yrange ['+str(y0)+':'+str(y1)+']\n')
	fout.write('set border lw 2 \n')

	fout.write('set xtics (')
	for i in range(len(Special)):
		fout.write('"'+Route[i]+'" '+str(Special[i]))
		if i < len(Special)-1:
			fout.write(', ')
		else:
			fout.write(')\n')

	for i in range(len(Special)-2):
		fout.write('set arrow from '+str(Special[i+1])+','+str(y0)+' to '+str(Special[i+1])+','+str(y1)+' nohead lt 0 lw 4\n')
	
	fout.write('set arrow from 0,0 to '+str(Special[len(Special)-1])+',0 nohead lt 0 lw 4\n')

	fout.write('plot ')
	for i in range(NBAND):
		fout.write("'"+dir+spinf1+"'"+' u 1:'+str(i+2)+' w l lw 2 lc '+str(color1)+',\\\n')
	for i in range(NBAND-1):
		fout.write("'"+dir+spinf2+"'"+' u 1:'+str(i+2)+' w l lw 2 lc '+str(color2)+',\\\n')
	fout.write("'"+dir+spinf2+"'"+' u 1:'+str(NBAND+1)+' w l lw 2 lc '+str(color2)+'\n')

	fout.close()

def plot(y0,y1,color,color2,Dir='',outputfilename='output.jpeg',filename='plot_test'):
	global ISPIN
	if ISPIN == 1:
		print2File(0,'plot.dat',Dir)
		writePlot(y0,y1,color,'plot.dat',outputfilename,Dir,filename)
		os.system('gnuplot<'+Dir+filename+'>'+outputfilename)
		os.system('rm '+Dir+filename)
		os.system('rm '+Dir+'plot.dat')
	if ISPIN == 2:
		print2File(0,'SPINUP.dat',Dir)
		print2File(1,'SPINDOWN.dat',Dir)
		plotSpin(y0,y1,color,color2,'SPINUP.dat','SPINDOWN.dat',outputfilename,Dir,filename)
		os.system('gnuplot<'+Dir+filename+'>'+outputfilename)
		os.system('rm '+Dir+filename)

def dirGenerator():
	global DIR
	os.system('mkdir pic')
	M = ['Co','Cr','Fe','Mn','Ni','Sc','Ti','V']
	O = ['N','O']
	X = ['Br','Cl','F']
	for x in M:
		for y in O:
			for z in X:
				name = x+y+z
				if os.path.exists(name+'/band/EIGENVAL'):
					initialize()
					DIR = name+'/'
					outputfilename = 'pic/'+name+'-band'
					read_EIGENVAL()
					plot(-23,12,8,7,DIR,outputfilename+'0.eps')
					plot(-4,4,8,7,DIR,outputfilename+'1.eps')

def run():
	initialize()
	read_EIGENVAL()
	#plot(3.02527,20.03929,8,7,'','band.eps')
	plot(-23,12,8,7,DIR,'band0.eps')
	plot(-4,4,8,7,DIR,'band1.eps')
	plot(0,0.75,8,7,DIR,'band2.eps')
	#dirGenerator()
run()
#print(Special)
#dirGenerator()
#print_SPIN()
#print_SPIN(1)
