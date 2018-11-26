import re
import os
import subprocess
from numpy import sqrt,array,zeros

from subprocess import Popen, PIPE, call

global KPOINTS,NBANDS,NATOMS,ENERGY,K,DOS,Special
ORBITAL=['s','py','pz','px','dxy','dyz','dz2','dxz','dx2-y2','tot']
HIGH_SYM = ['{/Symbol G}','M','K','{/Symbol G}']

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def correction():
	band=float(getnum(subprocess.getoutput('grep E-fermi '+'OUTCAR'))[0])
	return -band

def read_para_OUTCAR():
	global FERMI,SOC,ISPIN,NBANDS,NIONS
	FERMI = float(getnum(subprocess.getoutput('grep E-fermi '+'OUTCAR'))[0])
	SOC = bool(subprocess.getoutput('grep LSORBIT '+'OUTCAR')[18])
	ISPIN = int(getnum(subprocess.getoutput('grep ISPIN '+'OUTCAR'))[0])
	NBANDS = int(getnum(subprocess.getoutput('grep NBANDS '+'OUTCAR'))[-1])
	NIONS = int(getnum(subprocess.getoutput('grep NIONS '+'OUTCAR'))[1])

def read_kpath_OUTCAR():
	global NKPTS,KPOINTS,Kpath,Route,Special
	NKPTS = int(getnum(subprocess.getoutput('grep NKPTS '+'OUTCAR'))[0])
	K_lines = subprocess.getoutput('grep -A '+str(NKPTS)+' 2pi/SCALE '+'OUTCAR'+'|tail -'+str(NKPTS)).split('\n')
	KPOINTS = []
	Route = []
	k = 0.
	Special=[0.]
	Route.append(k)
	for i in range(NKPTS):
		kx = float(getnum(K_lines[i])[0])
		ky = float(getnum(K_lines[i])[1])
		kz = float(getnum(K_lines[i])[2])
		#occ = float(getnum(K_lines[i])[3])
		KPOINTS.append([kx,ky,kz])
		if i>0:
			k+=sqrt(sum((array(KPOINTS[-1])-array(KPOINTS[-2]))**2))
			Route.append(k)
			if Route[-1]==Route[-2]:
				Special.append(Route[-1])
	Special.append(Route[-1])

def read_dos(dos_line):
	s = float(getnum(dos_line)[-10])
	py = float(getnum(dos_line)[-9])
	pz = float(getnum(dos_line)[-8])
	px = float(getnum(dos_line)[-7])
	dxy = float(getnum(dos_line)[-6])
	dyz = float(getnum(dos_line)[-5])
	dz2 = float(getnum(dos_line)[-4])
	dxz = float(getnum(dos_line)[-3])
	dx2 = float(getnum(dos_line)[-2])
	tot = float(getnum(dos_line)[-1])
	return [s,py,pz,px,dxy,dyz,dz2,dxz,dx2,tot]

def read_PROCAR_spin():
	global NKPTS,NBANDS,NIONS,ENERGY_UP,ENERGY_DOWN,FERMI,DOS_UP,DOS_DOWN
	fin = open('PROCAR','r')
	lines = fin.readlines()
	fin.close()
	ENERGY_UP = []
	ENERGY_DOWN = []
	DOS_UP = []
	DOS_DOWN = []
	head = ((NIONS+5)*NBANDS+3)*NKPTS+1
	for i in range(NKPTS):
		print('loading kpoints # ',str(i),'...')
		line_K_UP = ((NIONS+5)*NBANDS+3)*i+1+1
		line_K_DOWN = line_K_UP + head
		E_k_UP = []
		E_k_DOWN = []
		DOS_k_UP = []
		DOS_k_DOWN = []
		for j in range(NBANDS):
			line_band_UP = line_K_UP + (NIONS+5)*j + 3
			line_band_DOWN = line_band_UP + head
			energy_band_UP = float(getnum(lines[line_band_UP])[1]) - FERMI
			energy_band_DOWN = float(getnum(lines[line_band_DOWN])[1]) - FERMI
			E_k_UP.append(energy_band_UP)
			E_k_DOWN.append(energy_band_DOWN)
			DOS_band_UP = []
			DOS_band_DOWN = []
			for ion in range(NIONS):
				line_ion_UP = line_band_UP + 3 + ion
				line_ion_DOWN = line_ion_UP + head
				DOS_band_UP.append(read_dos(lines[line_ion_UP]))
				DOS_band_DOWN.append(read_dos(lines[line_ion_DOWN]))
			DOS_k_UP.append(DOS_band_UP)
			DOS_k_DOWN.append(DOS_band_DOWN)
		ENERGY_UP.append(E_k_UP)
		ENERGY_DOWN.append(E_k_DOWN)
		DOS_UP.append(DOS_k_UP)
		DOS_DOWN.append(DOS_k_DOWN)

def sum_DOS(atom1,atom2):
	global NKPTS,DOS_UP,DOS_DOWN
	TDOS_UP = []
	TDOS_DOWN = []
	for i in range(NKPTS):
		TDOS_K_UP = []
		TDOS_K_DOWN = []
		for j in range(NBANDS):
			TDOS_band_UP = zeros(10)
			TDOS_band_DOWN = zeros(10)
			for ion in range(atom1-1,atom2):
				TDOS_band_UP += array(DOS_UP[i][j][ion])
				TDOS_band_DOWN += array(DOS_DOWN[i][j][ion])
			TDOS_K_UP.append(TDOS_band_UP)
			TDOS_K_DOWN.append(TDOS_band_DOWN)
		TDOS_UP.append(TDOS_K_UP)
		TDOS_DOWN.append(TDOS_K_DOWN)
	return TDOS_UP,TDOS_DOWN

def print2File(x,y,z,orbit,file):
	fout = open(file,'w')
	for i in range(NKPTS):
		for j in range(NBANDS):
			fout.write(str(x[i])+' '+str(y[i][j])+' '+str(z[i][j][orbit])+' ')
		fout.write('\n')
	fout.close()

def writePlot_spin(x,y0,y1,spinfile1,spinfile2,color1=8,color2=7):
	global Special,HIGH_SYM,NBANDS
	filename = 'pro_plot.gnu'
	fout = open(filename,'w')
	fout.write('set term post eps color enhanced "Helvetica" 20\n')
	fout.write('unset key\n')
	fout.write('set ylabel "Energy (eV)"\n')
	fout.write('set xrange [0:'+str(x[-1])+']\n')
	fout.write('set yrange ['+str(y0)+':'+str(y1)+']\n')
	fout.write('set arrow from 0,0 to '+str(x[-1])+',0 nohead lt 0 lw 2\n')#set fermi level
	#Plot High symmetry points (lines and xtics)
	fout.write('set xtics (')
	for i in range(len(Special)-1):
		fout.write('"'+HIGH_SYM[i]+'" '+str(Special[i])+', ')
	fout.write('"'+HIGH_SYM[-1]+'" '+str(Special[-1])+')\n')
	for i in range(1,len(Special)):
		fout.write('set arrow from '+str(Special[i])+','+str(y0)+' to '+str(Special[i])+','+str(y1)+' nohead lt 0 lw 2\n')
	#Plot data
	fout.write("plot ")
	for i in range(NBANDS):
		fout.write("'"+spinfile1+"' u "+str(i*3+1)+":"+str(i*3+2)+":(0.005*sqrt($"+str(i*3+3)+")) w circles lc 6 fs transparent solid 0.8 noborder,\\\n")
		fout.write("'"+spinfile2+"' u "+str(i*3+1)+":"+str(i*3+2)+":(0.005*sqrt($"+str(i*3+3)+")) w circles lc 6 fs transparent solid 0.8 noborder,\\\n")
	for i in range(NBANDS-1):
		fout.write("'"+spinfile1+"' u "+str(i*3+1)+":"+str(i*3+2)+" w l lc "+str(color1)+" lw 4,\\\n")
		fout.write("'"+spinfile2+"' u "+str(i*3+1)+":"+str(i*3+2)+" w l lc "+str(color2)+" lw 4,\\\n")
	i+=1
	fout.write("'"+spinfile1+"' u "+str(i*3+1)+":"+str(i*3+2)+" w l lc "+str(color1)+" lw 4,\\\n")
	fout.write("'"+spinfile2+"' u "+str(i*3+1)+":"+str(i*3+2)+" w l lc "+str(color2)+" lw 4\n")
	fout.close()

def PlotElement_spin(atom1,atom2,ATOM):
	TDOS1,TDOS2=sum_DOS(atom1,atom2)
	for i in range(10):
		print2File(Route,ENERGY_UP,TDOS1,i,'SPINUP.dat')
		print2File(Route,ENERGY_DOWN,TDOS2,i,'SPINDOWN.dat')
		writePlot_spin(Route,-1.5,1.5,'SPINUP.dat','SPINDOWN.dat')
		os.system('gnuplot<pro_plot.gnu>'+ATOM+'_'+ORBITAL[i]+'.eps')
		os.system('convert -density 600 '+ATOM+'_'+ORBITAL[i]+'.eps '+ATOM+'_'+ORBITAL[i]+'.png')
		os.system('rm '+ATOM+'_'+ORBITAL[i]+'.eps')
		os.system('rm SPINUP.dat SPINDOWN.dat')
		os.system('rm pro_plot.gnu')

def Auto_PlotElement(): #Auto Plot (POSCAR required)
	atom = re.findall(r"[A-Z]+[a-z]*",subprocess.getoutput('sed -n "6,6p" '+'POSCAR'))
	natom = re.findall(r"[0-9]+",subprocess.getoutput('sed -n "7,7p" '+'POSCAR'))
	print(atom,natom)
	atom_start = 0
	for i in range(len(atom)):
		atom_start += 1
		PlotElement_spin(atom_start,atom_start+int(natom[i])-1,atom[i])
		atom_start += int(natom[i])-1

read_para_OUTCAR()
read_kpath_OUTCAR()
read_PROCAR_spin()
Auto_PlotElement()
#PlotElement_spin(1,2,'Ag')
#PlotElement_spin(3,8,'Cl')
#PlotElement_spin(9,9,'Re')

