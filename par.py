import re
import os
import subprocess
from numpy import sqrt,array,zeros

from subprocess import Popen, PIPE, call

global KPOINTS,NBANDS,NATOMS,ENERGY,K,DOS
ORBITAL=['s','py','pz','px','dxy','dyz','dz2','dxz','dx2','tot']

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def correction():
	band=float(getnum(subprocess.getoutput('grep E-fermi '+'OUTCAR'))[0])
	return -band

def read_PROCAR():
	global KPOINTS,NBANDS,NATOMS,ENERGY,K,DOS
	fin = open('PROCAR','r')
	lines = fin.readlines()
	fin.close()
	KPOINTS = int(getnum(lines[1])[0])
	NBANDS = int(getnum(lines[1])[1])
	NATOMS = int(getnum(lines[1])[2])
	NATOMS = 76*4-1
	K = []
	ENERGY = []
	DOS = []
	dE = correction()
	#print(NATOMS)
	for i in range(KPOINTS):
		line_K = (NBANDS*(NATOMS+5)+3)*i+3
		kx = float(getnum(lines[line_K])[1])
		ky = float(getnum(lines[line_K])[2])
		kz = float(getnum(lines[line_K])[3])
		K.append([kx*0.121506683,kx*0.070151916+ky*0.140303832,kz*0.050000000])
		Ek = []
		DOS_k = []
		for j in range(NBANDS):
			line_E = line_K+(NATOMS+5)*j+2
			energy = float(getnum(lines[line_E])[1])+dE
			Ek.append(energy)
			DOS_B = []
			for ion in range(NATOMS+1):
				line_I = line_E+ion+3
				print(line_I,line_E,line_K)
				s = float(getnum(lines[line_I])[-10])
				py = float(getnum(lines[line_I])[-9])
				pz = float(getnum(lines[line_I])[-8])
				px = float(getnum(lines[line_I])[-7])
				dxy = float(getnum(lines[line_I])[-6])
				dyz = float(getnum(lines[line_I])[-5])
				dz2 = float(getnum(lines[line_I])[-4])
				dxz = float(getnum(lines[line_I])[-3])
				dx2 = float(getnum(lines[line_I])[-2])
				tot = float(getnum(lines[line_I])[-1])
				DOS_B.append([s,py,pz,px,dxy,dyz,dz2,dxz,dx2,tot])
			DOS_k.append(DOS_B)
		ENERGY.append(Ek)
		DOS.append(DOS_k)

def generate(atom1,atom2):
	x = [0.]
	y = ENERGY
	z = []
	for i in range(len(K)):
		if i>0:
			x.append(x[i-1]+sqrt(sum((array(K[i])-array(K[i-1]))**2)))
		zj=[]
		for j in range(NBANDS):
			t = zeros(10)
			for n in range(atom1,atom2):
				t+=array(DOS[i][j][n])
			zj.append(t)
		z.append(zj)
	return x,y,z

def print2File(x,y,z,orbit):
	fout = open('plot.dat','w')
	for i in range(KPOINTS):
		for j in range(NBANDS):
			fout.write(str(x[i])+' '+str(y[i][j])+' '+str(z[i][j][orbit])+' ')
		fout.write('\n')
	fout.close()

def writePlot(x,y0,y1):
	filename = 'pro_plot'
	outputfile = 'plot.dat'
	fout = open(filename,'w')
	fout.write('set term post eps color enhanced "Helvetica" 20\n')
	fout.write('unset key\n')
	fout.write('set ylabel "Energy (eV)"\n')
	fout.write('set xrange [0:'+str(x[-1])+']\n')
	fout.write('set yrange ['+str(y0)+':'+str(y1)+']\n')
	fout.write('set xtics ("{/Symbol G}" '+str(x[0])+', "M" '+str(x[int(KPOINTS/3)-1])+', "K" '+str(x[int(KPOINTS*2/3)-1])+', "{/Symbol G}" '+str(x[-1])+')\n')
	fout.write('set arrow from '+str(x[int(KPOINTS/3)-1])+','+str(y0)+' to '+str(x[int(KPOINTS/3)-1])+','+str(y1)+' nohead lt 0 lw 2\n')
	fout.write('set arrow from '+str(x[int(KPOINTS*2/3)-1])+','+str(y0)+' to '+str(x[int(KPOINTS*2/3)-1])+','+str(y1)+' nohead lt 0 lw 2\n')
	fout.write('set arrow from 0,0 to '+str(x[-1])+',0 nohead lt 0 lw 2\n')
	fout.write("plot ")
	for i in range(NBANDS):
		fout.write("'"+outputfile+"' u "+str(i*3+1)+":"+str(i*3+2)+":(0.005*sqrt($"+str(i*3+3)+")) w circles lc 6 fs transparent solid 0.8 noborder,\\\n")
	for i in range(NBANDS-1):
		fout.write("'"+outputfile+"' u "+str(i*3+1)+":"+str(i*3+2)+" w l lc 8 lw 4,\\\n")
	i+=1
	fout.write("'"+outputfile+"' u "+str(i*3+1)+":"+str(i*3+2)+" w l lc 8 lw 4\n")
	fout.close()

def PlotElement(atom1,atom2,ATOM):
	read_PROCAR()
	x,y,z=generate(atom1,atom2)
	for i in range(10):
		print2File(x,y,z,i)
		writePlot(x,-15,10)
		os.system('gnuplot<pro_plot>'+ATOM+'_'+ORBITAL[i]+'.eps')
		os.system('convert -density 600 '+ATOM+'_'+ORBITAL[i]+'.eps '+ATOM+'_'+ORBITAL[i]+'.png')
		#os.system('rm '+ATOM+'_'+ORBITAL[i]+'.eps')

#PlotElement(3,9,'O')
#PlotElement(0,3,'Cu')
PlotElement(0,24,'H')
PlotElement(24,60,'C')
PlotElement(60,72,'N')
PlotElement(72,75,'Ta')
