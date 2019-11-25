#
# Band structure plot script for VASP calculations 
# gnuplot, ImageMagick is required to process images
# code in Python3
#
# Coder: Gao (Gawcista) Yifan (SUSTech & HKUST)
#

from subprocess import getoutput
from os import system,path,remove

# calculate distance between two sites
def distance(a,b):
	dis = .0
	for i in range(3):
		dis += (a[i]-b[i])**2
	return dis**0.5

# initialize constants and global variables
def initialize():
	global Kpoints,Special,Kpath,Energy,Breakpoints
	# initialization for list, don't make change here 
	Special = []
	Kpoints = []
	Energy = [[],[]]
	Kpath = []
	Breakpoints = []
	# Global constants you can make change here
	global DIR,Xtics
	DIR = 'band/' #where your OUTCAR file locates
	Xtics = ['{/Symbol G}','L','B_1','B','Z','{/Symbol G}','X','Q','F','P_1','Z','L','P'] # the high-symmetry point tics


# read data from OUTCAR
def read_OUTCAR():
	global DIR,NKPTS,NBANDS,E_fermi,ISPIN
	global Kpints,Special,Kpath,Energy,Breakpoints
	# read constant parameters
	NKPTS = int(getoutput("grep NKPTS %sOUTCAR"%DIR).split()[3])
	NBANDS = int(getoutput("grep NBANDS %sOUTCAR"%DIR).split()[-1])
	ISPIN = int(getoutput("grep ISPIN %sOUTCAR"%DIR).split()[2])
	E_fermi = float((getoutput("grep E-fermi %sOUTCAR"%DIR)).split()[2])
	# read k-path and find high-symmetry points
	kpoint_list = getoutput("grep 2pi/SCALE %sOUTCAR -A %d|tail -n %d"%(DIR,NKPTS,NKPTS)).split()
	Special.append(0.)
	route = 0.
	for i in range(NKPTS):
		Kpoints.append([float(kpoint_list[4*i]),float(kpoint_list[4*i+1]),float(kpoint_list[4*i+2])])
		if i > 0:
			dr = distance(Kpoints[i],Kpoints[i-1])
			if i > 2 and dr > 5*(Kpath[-1]-Kpath[-2]) and Kpath[-1]-Kpath[-2] != 0:
				dr = 0
				Breakpoints.append(route)
				Xtics[len(Special)] = Xtics[len(Special)]+'|'+Xtics[len(Special)+1]
				del Xtics[len(Special)+1]
			route += dr
			if dr == 0:
				Special.append(route)
		Kpath.append(route)
	Special.append(route)
	# read eigenvalues
	global VBM,CBM
	VBM = [-100,0,0]
	CBM = [100.0,0]
	for k in range(1,NKPTS+1):
		klines = getoutput("grep 'k-point%6d' %sOUTCAR -A %d|egrep -v 'k-point|band'|awk '{print $2}'"%(k,DIR,NBANDS+1)).split()
		Ek = [[],[]] # eigenvalues for kth k-point
		for band in range(NBANDS):
			for spin in range(ISPIN):
				Ek[spin].append(float(klines[band+NBANDS*spin])-E_fermi)
				if Ek[spin][-1]>E_fermi and Ek[spin][-1]<CBM[0]:
					CBM = [Ek[spin][-1],band,k]
				if Ek[spin][-1]<E_fermi and Ek[spin][-1]>VBM[0]:
					VBM = [Ek[spin][-1],band,k]
		for spin in range(ISPIN):
			Energy[spin].append(Ek[spin])

# write plot data
def write_plot(type=1,head='',dir=''):
	if type == 1:
		for spin in range(ISPIN):
			with open('%s%d.dat'%(head,spin),'w') as fout:
				for band in range(NBANDS):
					for k in range(NKPTS):
						print('%f\t%f'%(Kpath[k],Energy[spin][k][band]),file=fout)
						if k<NKPTS-1 and Kpath[k] in Breakpoints and Kpath[k+1] in Breakpoints:
							print(file=fout)
					print(file=fout)
	if type == 2:
		for spin in range(ISPIN):
			with open('%s%d.dat'%(head,spin),'w') as fout:
				for k in range(NKPTS):
					print(Kpath[k],end='\t',file=fout)
					for band in range(NBANDS):
						print(Energy[spin][k][band],end='\t',file=fout)
					print(file=fout)

#write gnuplot command
def write_gnuplot(y0,y1,color=['blue','red'],outputfilename='band'):
	font = 'Times'
	border_width = 4
	arrow_width = 2
	line_width = 4
	xtics_size = 16
	ytics_size = 20
	ylabel_size = 28
	with open('plot.gnu','w') as fout:
		print('set term post eps color enhanced "%s"'%font,file=fout)
		print('set output "%s.eps"'%outputfilename,file=fout)
		print('unset key',file=fout)
		print('set size 0.6,1',file=fout)
		print('set border lw %d'%border_width,file=fout)
		print('set ylabel "Energy (eV)" font "%s,%d"'%(font,ylabel_size),file=fout)
		print('set xrange [%f:%f]'%(Kpath[0],Kpath[-1]),file=fout)
		print('set yrange [%f:%f]'%(y0,y1),file=fout)
		print('set ytics font "%s,%d"'%(font,ytics_size),file=fout)
		print('set xtics () font "%s,%d"'%(font,xtics_size),file=fout)
		for i in range(len(Special)):
			if '_' not in Xtics[i]:
				Xtics[i]+='_'
			print('set xtics add ("%s" %f)'%(Xtics[i],Special[i]),file=fout)
		print('set arrow from 0,0 to %f,0 nohead lt 0 lw %d'%(Special[-1],arrow_width),file=fout)
		for x in Special[1:-1]:
			if x not in Breakpoints:
				print('set arrow from %f,%f to %f,%f nohead front lt 0 lw %d'%(x,y0,x,y1,arrow_width),file=fout)
			else:
				print('set arrow from %f,%f to %f,%f nohead front lt -1 lw %d'%(x,y0,x,y1,arrow_width),file=fout)
		print('plot ',end=' ',file=fout)
		for spin in range(ISPIN):
			print('"%d.dat" w l lw %d lc rgb "%s"'%(spin,line_width,color[spin]),file=fout,end='')
			if spin+1 < ISPIN:	
				print(',\\',file=fout)

# convert eps to png
def eps_to_png(name):
	system("convert -density 600 %s.eps %s.png"%(name,name))

# clean all unnecessary files
def clean_all(name):
	file_list = ['0.dat','1.dat','plot.gnu','%s.eps'%name]
	for file in file_list:
		if path.exists(file):
			remove(file)

# quick plot 
def plot(y0,y1,color=['blue','red'],file='band'): 
	initialize()
	read_OUTCAR()
	write_plot()
	write_gnuplot(y0,y1,color,file)
	system('gnuplot<plot.gnu')
	eps_to_png(file)
	clean_all(file)

if __name__ == '__main__':
	plot(-2,2)
	print(CBM[0]-VBM[0])