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
	global Kpoints,Special,Kpath,Energy
	# initialization for list, don't make change here 
	Special = []
	Kpoints = []
	Energy = [[],[]]
	Kpath = []
	# Global constants you can make change here
	global DIR,Xtics
	DIR = 'band/' #where your OUTCAR file locates
	Xtics = ['{/Symbol G}','M','K','{/Symbol G}'] # the high-symmetry point tics


# read data from OUTCAR
def read_OUTCAR():
	global DIR,NKPTS,NBANDS,E_fermi,ISPIN
	global Kpints,Special,Kpath,Energy
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
			route += dr
			if dr == 0:
				Special.append(route)
		Kpath.append(route)
	Special.append(route)
	# read eigenvalues
	for k in range(1,NKPTS+1):
		klines = getoutput("grep 'k-point%6d' %sOUTCAR -A %d|egrep -v 'k-point|band'|awk '{print $2}'"%(k,DIR,NBANDS+1)).split()
		Ek = [[],[]] # eigenvalues for kth k-point
		for band in range(NBANDS):
			for spin in range(ISPIN):
				Ek[spin].append(float(klines[band+NBANDS*spin])-E_fermi)
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
	border_width = 4
	arrow_width = 4
	line_width = 4
	with open('plot.gnu','w') as fout:
		#print('set term png font "times.ttf,14"')
		print('set term post eps color enhanced "Times,20"',file=fout)
		print('set output "%s.eps"'%outputfilename,file=fout)
		print('unset key',file=fout)
		print('set size 0.6,1',file=fout)
		print('set border lw %d'%border_width,file=fout)
		print('set ylabel "Energy (eV)"',file=fout)
		print('set xrange [%f:%f]'%(Kpath[0],Kpath[-1]),file=fout)
		print('set yrange [%f:%f]'%(y0,y1),file=fout)
		print('set xtics (',end='',file=fout)
		for i in range(len(Special)):
			if i < len(Special)-1:
				print('"%s" %f'%(Xtics[i],Special[i]),file=fout,end=',')
			else:
				print('"%s" %f'%(Xtics[i],Special[i]),file=fout,end=')\n')
		print('set arrow from 0,0 to %f,0 nohead lt 0 lw %d'%(Special[-1],arrow_width),file=fout)
		for x in Special[1:-1]:
			print('set arrow from %f,%f to %f,%f nohead lt 0 lw %d'%(x,y0,x,y1,arrow_width),file=fout)
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