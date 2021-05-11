#
# Band structure plot module for VASP calculations 
# gnuplot, ImageMagick is required to process images
# code in Python3
#
# Author: Gao (Gawcista) Yifan (SUSTech & HKUST)
#
from subprocess import getoutput
from os import system,path,remove
from numpy import array,linspace,cross,linalg
from copy import deepcopy

global NUM_SUB_PLOTS

ORBITAL=['s','py','pz','px','dxy','dyz','dz2','dxz','dx2-y2','tot']

class band_point:
	energy = 0.
	weight = 0.
	kx = 0.
	ky = 0.
	kz = 0.
	kpath = 0.
	def __init__(self,energy=0.,weight=0.):
		self.energy = energy
		self.weight = weight

class projected_band:
	atom1 = 0
	atom2 = 0
	orbital = []
	color = 'red'
	outputfilename = ''
	proj = [[],[]]
	def __init__(self,atom1=0,atom2=0,orbital=[],color='red',outputfilename=''):
		self.atom1 = atom1
		self.atom2 = atom2
		self.orbital = orbital
		self.color = color
		self.outputfilename = outputfilename
		self.proj = [[],[]]

	def calculate(self,ISPIN,NKPTS,NBANDS,Projection):
		for spin in range(ISPIN):
			for k in range(NKPTS):
				proj_k=[]
				for band in range(NBANDS):
					proj_band = 0.
					for ion in range(self.atom1-1,self.atom2):
						for orbit in self.orbital:
							proj_band += Projection[spin][k][band][ion][orbit]
					proj_k.append(proj_band)
				self.proj[spin].append(proj_k)

class band_structure:
	color = []
	color_circle = []
	DIR = ''
	Special = []
	Kpoints = []
	Energy = [[],[]]
	Projection = [[],[]]
	Kpath = []
	NKPTS = 0
	NBANDS = 0
	ISPIN = 0
	NIONS = 0
	NSPECIAL = 0
	E_fermi = 0
	Xtics = []
	VBM = [-100,0,0]
	CBM = [100,0,0]
	plot_gap = False
	set_plot_projection = False
	projection_list = []
	gap = 0
	METAGGA = 'F'
	LHFCALC = False
	def __init__(self,file='band/OUTCAR',selected_bands=[],xtics=['{/Symbol G}','M','K','{/Symbol G}','A','L','H','A','L','M','H','K'],color=['black','red'],dashtype=[1,1]):
		#initializing
		self.color = color
		self.dashtype = dashtype
		self.file = file
		self.Xtics = xtics
		self.Special = []
		self.Kpoints = []
		self.Kpoints_rec = []
		self.Energy = [[],[]]
		self.Projection = [[],[]]
		self.Kpath = []
		self.NSPECIAL = 0
		self.VBM = [0,-100,0,0]
		self.CBM = [0,100,0,0]
		self.gap = 0
		self.plot_gap = False
		self.plot_projection = False
		self.projection_list = []
		self.bandindex = []
		# read constant parameters
		self.NKPTS = 0#int(getoutput("grep NKPTS %s"%self.file).split()[3])
		if selected_bands!=[]:
			self.NBANDS=len(selected_bands)
			self.bandindex=selected_bands
		else:
			self.NBANDS = int(getoutput("grep 'number of bands    NBANDS=' %s"%self.file).split()[-1])
		self.NIONS = int(getoutput("grep NIONS %s"%self.file).split()[-1])
		self.ISPIN = int(getoutput("grep ISPIN %s"%self.file).split()[2])
		self.E_fermi = float((getoutput("grep E-fermi %s"%self.file)).split()[2])
		self.METAGGA = str(getoutput("grep 'METAGGA=' %s"%self.file).split()[1])
		self.LHFCALC = str(getoutput("grep 'LHFCALC =' %s"%self.file).split()[2])=='T'
		self.flag_read = 'EIGENVAL'
	
	def set_E_fermi(self,E):
		self.E_fermi = E

	def set_plot_gap(self,flag):
		self.plot_gap = flag

	def load(self,OUTCAR='band/OUTCAR',PROCAR='band/PROCAR',EIGENVAL='band/EIGENVAL'): # add band in OUTCAR, along k-axis
		NKPTS_new = int(getoutput("grep NKPTS %s"%OUTCAR).split()[3])
		NBANDS_new = int(getoutput("grep 'number of bands    NBANDS=' %s"%OUTCAR).split()[-1])
		if self.bandindex!=[]:
			band_list=self.bandindex
		else:
			band_list=[i for i in range(0,NBANDS_new,1)]
		NK_skip = 0
		if self.METAGGA != 'F' or self.LHFCALC:
			NK_skip = 1
			while float((getoutput("grep 2pi/SCALE %s -A %d|tail -n %d|tail -1"%(OUTCAR,NK_skip,NK_skip)).split()[-1])) > 0:
				NK_skip += 1
			NK_skip = NK_skip-1
		# read k-path and find high-symmetry points
		kpoint_list = getoutput("grep '2pi/SCALE' %s -A %d|tail -n %d"%(OUTCAR,NKPTS_new,NKPTS_new-NK_skip)).split()
		kpoint_list_rec = getoutput("grep 'k-points in reciprocal lattice and weights' %s -A %d|tail -n %d"%(OUTCAR,NKPTS_new,NKPTS_new-NK_skip)).split()
		if len(self.Special)==0:
			route = 0
			self.Special.append([route,0,self.NSPECIAL])
		else:
			route = self.Kpath[-1]
		flag_special = False
		for i in range(NKPTS_new-NK_skip):
			self.Kpoints.append([float(kpoint_list[4*i]),float(kpoint_list[4*i+1]),float(kpoint_list[4*i+2])])
			self.Kpoints_rec.append([float(kpoint_list_rec[4*i]),float(kpoint_list_rec[4*i+1]),float(kpoint_list_rec[4*i+2])])
			if len(self.Kpoints) > 1:
				dr = distance(self.Kpoints[-1],self.Kpoints[-2])
				if len(self.Kpath)>1 and not is_parallel(array(self.Kpoints[-1])-array(self.Kpoints[-2]),array(self.Kpoints[-2])-array(self.Kpoints[-3])):
					dr = 0
				if flag_special:
					flag_special = False
					if distance(self.Kpoints[-2],self.Kpoints[-3])>1e-8:
						self.NSPECIAL +=1
						self.Special.append([route,i+self.NKPTS-1,self.NSPECIAL])
					elif not is_parallel(array(self.Kpoints[-1])-array(self.Kpoints[-2]),array(self.Kpoints[-3])-array(self.Kpoints[-4])):
						self.Special.append([route,i+self.NKPTS-1,self.NSPECIAL])
					else:
						self.NSPECIAL -=1
						del(self.Special[-1])
				if dr == 0:
					flag_special = True
				route += dr
			self.Kpath.append(route)
		if route != self.Special[-1][0] :
			self.NSPECIAL += 1
		self.Special.append([route,i+self.NKPTS,self.NSPECIAL])

		if self.flag_read=='OUTCAR':# read eigenvalues from OUTCAR
			for spin in range(self.ISPIN):
				for k in range(1+NK_skip,NKPTS_new+1):
					klines = getoutput("grep 'k-point%6d' %s -A %d|egrep -v 'k-point|band'|awk '{print $2}'"%(k,OUTCAR,NBANDS_new+1)).split()
					band_k = [[],[]] # eigenvalues for kth k-point
					for band in band_list:
						bandpoint=band_point(energy=float(klines[band+NBANDS_new*spin])-self.E_fermi)
						band_k[spin].append(bandpoint)
						if bandpoint.energy>0 and bandpoint.energy<self.CBM[1]:
							self.CBM = [self.Kpath[k+self.NKPTS-1],bandpoint.energy,k+self.NKPTS-1,band]
						if bandpoint.energy<0 and bandpoint.energy>self.VBM[1]:
							self.VBM = [self.Kpath[k+self.NKPTS-1],bandpoint.energy,k+self.NKPTS-1,band]
					self.Energy[spin].append(band_k[spin])

		elif self.flag_read=='EIGENVAL':# read eigenvalues from EIGENVAL
			for spin in range(self.ISPIN):
				temp_k = []
				for band in band_list:
					bandlines = getoutput("grep '%5d      ' %s|awk '{print $%d}'"%(band+1,EIGENVAL,spin+2)).split()
					for k in range(1+NK_skip,NKPTS_new+1):
						bandpoint=band_point(energy=float(bandlines[k-NK_skip-1])-self.E_fermi)
						if k-NK_skip>len(temp_k):
							temp_k.append([bandpoint])
						else:
							temp_k[k-NK_skip-1].append(bandpoint)
						if bandpoint.energy>0 and bandpoint.energy<self.CBM[1]:
							self.CBM = [self.Kpath[k+self.NKPTS-1],bandpoint.energy,k+self.NKPTS-1,band]
						if bandpoint.energy<0 and bandpoint.energy>self.VBM[1]:
							self.VBM = [self.Kpath[k+self.NKPTS-1],bandpoint.energy,k+self.NKPTS-1,band]
				for k_list in temp_k:
					self.Energy[spin].append(k_list)

		# read projections from PROCAR
		if self.plot_projection:
			with open(PROCAR,'r') as fin:
				lines = fin.readlines()
			band_lines = self.NIONS + 5
			k_lines = band_lines*self.NBANDS + 3
			spin_lines = k_lines*NKPTS_new + 1
			for spin in range(self.ISPIN):
				for k in range(NK_skip,NKPTS_new):
					line_kpoints = k_lines*k+ 3 + spin_lines*spin
					energy_kpoints = []
					projection_kpoint = []
					for band in band_list:
						line_bands = line_kpoints+band_lines*band+2
						energy = float(lines[line_bands].split()[4])-self.E_fermi
						energy_kpoints.append(energy)
						projection_band = []
						tot_i = float(lines[line_bands+2+self.NIONS+1].split()[-1])
						for ion in range(self.NIONS+1):
							line_ion = line_bands+ion+3
							s = float(lines[line_ion].split()[-10])/tot_i
							py = float(lines[line_ion].split()[-9])/tot_i
							pz = float(lines[line_ion].split()[-8])/tot_i
							px = float(lines[line_ion].split()[-7])/tot_i
							dxy = float(lines[line_ion].split()[-6])/tot_i
							dyz = float(lines[line_ion].split()[-5])/tot_i
							dz2 = float(lines[line_ion].split()[-4])/tot_i
							dxz = float(lines[line_ion].split()[-3])/tot_i
							dx2 = float(lines[line_ion].split()[-2])/tot_i
							tot = float(lines[line_ion].split()[-1])/tot_i
							projection_band.append([s,py,pz,px,dxy,dyz,dz2,dxz,dx2,tot])
							#correction
							#for orbital in range(len(projection_band[-1])):
								#if projection_band[-1][orbital]<0.6:
								#	projection_band[-1][orbital] = 0
						projection_kpoint.append(projection_band)
					self.Projection[spin].append(projection_kpoint)
		
		self.gap = self.CBM[0]-self.VBM[0]
		self.NKPTS += NKPTS_new-NK_skip
	def rearrange(self,newlist):
		Special_new = []
		Kpoints_new = []
		Kpoints_rec_new = []
		Energy_new = [[],[]]
		Projection_new = [[],[]]
		Kpath_new = []
		VBM_new = [0,-100,0,0]
		CBM_new = [0,100,0,0]
		route_new = 0
		k_new = 0
		#print(self.Special)
		for i in range(len(newlist)):
			j = abs(newlist[i])
			if newlist[i]>0:
				Special_new.append(deepcopy(self.Special)[(j-1)*2])
				Special_new.append(deepcopy(self.Special)[(j-1)*2+1])
				Kpoints_new += self.Kpoints[Special_new[-2][1]:Special_new[-1][1]+1]
				Kpoints_rec_new += self.Kpoints_rec[Special_new[-2][1]:Special_new[-1][1]+1]
				Kpath_temp = self.Kpath[Special_new[-2][1]:Special_new[-1][1]+1]
				for spin in range(self.ISPIN):
					Energy_new[spin] += self.Energy[spin][Special_new[-2][1]:Special_new[-1][1]+1]
					Projection_new[spin] += self.Projection[spin][Special_new[-2][1]:Special_new[-1][1]+1]
			else:	
				Special_new.append(deepcopy(self.Special)[(j-1)*2+1])
				Special_new.append(deepcopy(self.Special)[(j-1)*2])
				Kpoints_new += list(reversed(self.Kpoints[Special_new[-1][1]:Special_new[-2][1]+1]))
				Kpoints_rec_new += list(reversed(self.Kpoints_rec[Special_new[-1][1]:Special_new[-2][1]+1]))
				Kpath_temp = self.Kpath[Special_new[-1][1]:Special_new[-2][1]+1]
				for spin in range(self.ISPIN):
					Energy_new[spin] += list(reversed(self.Energy[spin][Special_new[-1][1]:Special_new[-2][1]+1]))
					Projection_new[spin] += list(reversed(self.Projection[spin][Special_new[-1][1]:Special_new[-2][1]+1]))
			#print(abs(self.Special[(j-1)*2][0]-self.Special[(j-1)*2+1][0]))
			Special_new[-2][0] = route_new
			Special_new[-2][1] = k_new
			for i in range(len(Kpath_temp)):
				if i==0:
					Kpath_new.append(route_new)
				else:
					Kpath_new.append(Kpath_new[-1]+abs(Kpath_temp[i]-Kpath_temp[i-1]))
			route_new += abs(self.Special[(j-1)*2][0]-self.Special[(j-1)*2+1][0])
			k_new += abs(self.Special[(j-1)*2][1]-self.Special[(j-1)*2+1][1])
			Special_new[-1][0] = route_new
			Special_new[-1][1] = k_new
			k_new += 1

		self.Special = Special_new
		self.Kpoints = Kpoints_new
		self.Kpoints_rec = Kpoints_rec_new
		self.Energy = Energy_new
		self.Projection = Projection_new
		self.Kpath = Kpath_new
		#print(Special_new)
	
	def exchange(self,ibanda,ibandb,ik,ispin=0):
		self.Energy[ispin][ik][ibanda],self.Energy[ispin][ik][ibandb]=self.Energy[ispin][ik][ibandb],self.Energy[ispin][ik][ibanda]
		if self.plot_projection:
			self.Projection[ispin][ik][ibanda],self.Projection[ispin][ik][ibandb]=self.Projection[ispin][ik][ibandb],self.Projection[ispin][ik][ibanda]


	def read_bandplot(self,OUTCAR='OUTCAR',PROCAR='PROCAR',EIGENVAL='EIGENVAL'):
		self.load(OUTCAR=OUTCAR,PROCAR=PROCAR,EIGENVAL=EIGENVAL)



# calculate distance between two sites
def distance(a,b):
	dis = .0
	for i in range(3):
		dis += (a[i]-b[i])**2
	return dis**0.5

def is_parallel(a,b,error=1e-8):
	tmp = linalg.norm(cross(a,b))
	if tmp<error:
		return True
	else:
		return False

def write_plot(bands=[],head='',dir=''):
	filenum = 0
	for bandstruct in bands:
		#print(bandstruct.NBANDS,bandstruct.NKPTS,bandstruct.ISPIN)
		for spin in range(bandstruct.ISPIN):
			filenum += 1
			with open('%s%d.dat'%(head,filenum),'w') as fout:
				for band in range(bandstruct.NBANDS):
					cc=0
					for k in range(bandstruct.NKPTS):
						print('%f\t%f'%(bandstruct.Kpath[k],bandstruct.Energy[spin][k][band].energy+cc),file=fout,end=' ')
						if bandstruct.plot_projection:
							for proj in bandstruct.projection_list:
								print('\t%f'%proj.proj[spin][k][band],end=' ',file=fout)
						print(file=fout)
						if (k>0 and k in array(bandstruct.Special).T[1] and k+1 in array(bandstruct.Special).T[1]) or k+1 == bandstruct.NKPTS:
							print(file=fout)
						#if k<bandstruct.NKPTS-1 and bandstruct.Kpath[k] in bandstruct.Breakpoints and bandstruct.Kpath[k+1] in bandstruct.Breakpoints:
						#	print(file=fout)
					#print(file=fout)
	global NUM_SUB_PLOTS
	NUM_SUB_PLOTS = filenum

def print_plot(bands=[]):
	filenum = 0
	for bandstruct in bands:
		for spin in range(bandstruct.ISPIN):
			filenum += 1
			for k in range(bandstruct.NKPTS):
				print("%f\t%f\t%f\t%f"%(bandstruct.Kpoints_rec[k][0],bandstruct.Kpoints_rec[k][1],bandstruct.Kpoints_rec[k][2],bandstruct.Kpath[k]),end='\t')
				for band in range(bandstruct.NBANDS):
					print(bandstruct.Energy[spin][k][band].energy,end="\t")
				print()

def write_gnuplot(y0,y1,band_list=[],filename='',outputfilename=''):
	font = 'Times'
	x0 = band_list[0].Special[0][0]
	x1 = band_list[0].Special[-1][0]
	border_width = 4
	arrow_width = 2
	line_width = 16
	xtics_size = 22
	ytics_size = 22
	ylabel_size = 24
	with open('%s.gnu'%filename,'w') as fout:
		print('set term post eps dl 0.4 color enhanced "%s"'%font,file=fout)
		print('set output "%s.eps"'%outputfilename,file=fout)
		print('unset key',file=fout)
		print('set size ratio 0.8',file=fout)
		print('set border lw %d'%border_width,file=fout)
		print('set ylabel "Energy (eV)" font "%s,%d"'%(font,ylabel_size),file=fout)
		print('emin=%f'%y0,file=fout)
		print('emax=%f'%y1,file=fout)
		print('set xrange [%f:%f]'%(x0,x1),file=fout)
		print('set yrange [emin:emax]',file=fout)
		print('set ytics font "%s,%d"'%(font,ytics_size),file=fout)
		print('set xtics () font "%s,%d"'%(font,xtics_size),file=fout)
		print(band_list[0].Special)
		for bandstruct in band_list:
			if bandstruct.plot_gap:
				print('set object 1 circle center %f,%f size 0.0004 lc 1 fs transparent solid 0.8 noborder fc "red"'%(bandstruct.VBM[0],bandstruct.VBM[1]),file=fout)
				print('set object 2 circle center %f,%f size 0.0004 lc 1 fs transparent solid 0.8 noborder fc "blue"'%(bandstruct.CBM[0],bandstruct.CBM[1]),file=fout)		
		
		for i in range(len(band_list[0].Special)):
			if i==0 or i==len(band_list[0].Special)-1:
				print('set xtics add ("%s" %f)'%(band_list[0].Xtics[band_list[0].Special[i][-1]],band_list[0].Special[i][0]),file=fout)
			elif band_list[0].Special[i][0] == band_list[0].Special[i-1][0]:
				if band_list[0].Xtics[band_list[0].Special[i][2]] == band_list[0].Xtics[band_list[0].Special[i-1][2]]:
					print('set xtics add ("%s" %f)'%(band_list[0].Xtics[band_list[0].Special[i][-1]],band_list[0].Special[i][0]),file=fout)
					print('set arrow from %f,emin to %f,emax nohead front lt 0 lw %d'%(band_list[0].Special[i][0],band_list[0].Special[i][0],arrow_width),file=fout)
				else:
					print('set xtics add ("%s|%s" %f)'%(band_list[0].Xtics[band_list[0].Special[i-1][-1]],band_list[0].Xtics[band_list[0].Special[i][-1]],band_list[0].Special[i][0]),file=fout)	
					print('set arrow from %f,emin to %f,emax nohead front lt -1 lw %d'%(band_list[0].Special[i][0],band_list[0].Special[i][0],arrow_width),file=fout)

		filenum = 0
		print('plot ',end=' ',file=fout)
		for bandstruct in band_list:
			for spin in range(bandstruct.ISPIN):
				filenum += 1
				for p in range(len(bandstruct.projection_list)):
					print('"%s%d.dat" u 1:2:(0.001*sqrt($%d)) w circles lc rgb "%s" fs transparent solid 0.8 noborder,\\'%(filename,filenum,p+3,bandstruct.projection_list[p].color),file=fout)
			#	print(filenum,bandstruct.dashtype[spin],line_width,bandstruct.color[spin])
				print('"%s%d.dat" u 1:2 w l dt %d  lw %d lc rgb "%s"'%(filename,filenum,bandstruct.dashtype[spin],line_width,bandstruct.color[spin]),file=fout,end='')
				if filenum < NUM_SUB_PLOTS:
						print(',\\',file=fout)

# convert eps to png
def eps_to_png(name):
	system("convert -density 600 %s.eps %s.png"%(name,name))

# clean all unnecessary files
def clean_all(name):
	file_list = ['plot.gnu','%s.eps'%name]
	for i in range(NUM_SUB_PLOTS):
		file_list.append("%s.dat"%(i+1))
	for file in file_list:
		if path.exists(file):
			remove(file)