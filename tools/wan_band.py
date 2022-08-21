import re
import os
import subprocess
import math

global Route,Special,x0,x1

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def initialize(): 
	global Special,Route
	Route = ['{/Symbol G}','M','K','{/Symbol G}']
	Special=[]

def read_GNU():
	global Special,Xrange,x0,x1
	fin = open('wan/wannier90_band.gnu','r')
	lines = fin.readlines()
	fin.close()
	Special = getnum(lines[7])
	x0 = float(getnum(lines[2])[0])
	x1 = float(getnum(lines[2])[1])

def writePlot(y0,y1,color,datafile='wannier90_band.dat',outputfilename='output.eps',filename='plot_wan.gnu'):
	fout = open(filename,'w')
	fout.write('set term post eps color enhanced "Helvetica" 20\n')
	fout.write('set output "'+outputfilename+'"\n')
	fout.write('unset key\n')
	fout.write('set ylabel "Energy (eV)"\n')
	fout.write('set yrange ['+str(x0)+':'+str(x1)+']\n')
	fout.write('set yrange ['+str(y0)+':'+str(y1)+']\n')
	fout.write('set xtics (')
	for i in range(len(Special)):
		fout.write('"'+Route[i]+'" '+Special[i])
		if i < len(Special)-1:
			fout.write(', ')
		else:
			fout.write(')\n')
	for i in range(len(Special)-2):
		fout.write('set arrow from '+Special[i+1]+','+str(y0)+' to '+Special[i+1]+','+str(y1)+' nohead lt 0 lw 2\n')
	fout.write('set arrow from '+str(x0)+',0 to '+str(x1)+',0 nohead lt 0 lw 2\n')
	fout.write('plot "'+datafile+'" w l lw 4 lc '+str(color))
	fout.close()

def plot(y0,y1,color,outputfilename='output.eps',datafile='plot_wan.dat',filename='plot_wan.gnu'):
	writePlot(y0,y1,color,datafile,outputfilename,filename)
	os.system('gnuplot<'+filename+'>'+outputfilename)
	os.system('convert -density 800 '+outputfilename+' '+outputfilename+'.png')
	os.system('rm '+outputfilename)


def edit_plot(raw='wan/wannier90_band.dat',output='plot_wan.dat'):
	E_fermi = float(getnum(subprocess.getoutput('grep E-fermi bandsoc/OUTCAR'))[0])
	fin = open(raw,'r')
	lines = fin.readlines()
	fin.close()
	#print(E_fermi)
	fout = open(output,'w')
	for line in lines:
		xy = getnum(line)
		if len(xy)>3:
			fout.write(str(float(xy[0])*10**float(xy[1]))+' '+str(float(xy[2])*10**float(xy[3])-E_fermi))
		fout.write('\n')
	fout.close()

initialize()
read_GNU()
edit_plot()
print(Special)
line_color = 6
plot(-15,20,line_color,'band_wan')
plot(-1,1,line_color,'band_wan1')
plot(-0.005,0.005,line_color,'band_wan2')
