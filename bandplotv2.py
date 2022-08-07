#
# Band structure plot script for VASP calculations 
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
	def __init__(self,file='band/OUTCAR',xtics=['{/Symbol G}','M','K','{/Symbol G}','A','L','H','A','L','M','H','K'],color=['black','red'],dashtype=[1,1]):
		#initializing
		self.color = color
		self.dashtype = dashtype
		self.file = file
		self.Xtics = xtics
		self.Special = []
		self.Kpoints = []
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
		# read constant parameters
		self.NKPTS = 0#int(getoutput("grep NKPTS %s"%self.file).split()[3])
		self.NBANDS = int(getoutput("grep 'number of bands    NBANDS=' %s"%self.file).split()[-1])
		self.NIONS = int(getoutput("grep NIONS %s"%self.file).split()[-1])
		self.ISPIN = int(getoutput("grep ISPIN %s"%self.file).split()[2])
		self.E_fermi = float((getoutput("grep E-fermi %s"%self.file)).split()[2])
		#self.METAGGA = str(getoutput("grep 'METAGGA=' %s"%self.file).split()[1])
		self.LHFCALC = str(getoutput("grep 'LHFCALC =' %s"%self.file).split()[2])=='T'
		self.flag_read = 'EIGENVAL'
		self.autoKLABEL = False

	def set_E_fermi(self,E):
		self.E_fermi = E

	def set_plot_gap(self,flag):
		self.plot_gap = flag

	def set_plot_projection(self,flag):
		self.set_plot_projection = flag

	def set_projection(self,atom1=0,atom2=0,orbital=[],color='red',outputfilename=''):
		new_proj = projected_band(atom1=atom1,atom2=atom2,orbital=orbital,color=color,outputfilename=outputfilename)
		new_proj.calculate(self.ISPIN,self.NKPTS,self.NBANDS,self.Projection)
		self.projection_list.append(new_proj)
	
	def generate_newobj(self,obj,indexlist):
		# reorder certain object list according to given index list
		obj_new = []
		for i in indexlist:
			obj_new.append(obj[i])
		return obj_new

	def generate_k_path(self,Kpoints):
		# generate Kpath and Special for given Kpoints list
		route = 0.
		Kpath = [0.]
		Special = []
		tmp = 0.
		NSPECIAL = 0
		for i in range(1,len(Kpoints)):
			dr = distance(Kpoints[i],Kpoints[i-1])
			if abs(dr-tmp)>1e-6:
				Special.append([Kpath[i-1],i-1,NSPECIAL])
				if abs(dr)>1e-6:
					NSPECIAL+=1
			route += dr
			Kpath.append(route)
			tmp = dr
		Special.append([Kpath[-1],len(Kpoints)-1,NSPECIAL])
		return Kpath,Special
	
		

	def load(self,OUTCAR='band/OUTCAR',PROCAR='band/PROCAR',EIGENVAL='band/EIGENVAL',iband=[]): # add band in OUTCAR, along k-axis
		NKPTS_new = int(getoutput("grep NKPTS %s"%OUTCAR).split()[3])
		NBANDS_new = int(getoutput("grep 'number of bands    NBANDS=' %s"%OUTCAR).split()[-1])
		NK_skip = 0
		if iband == []:
			iband = [0,NBANDS_new-1]
		if self.METAGGA != 'F' or self.LHFCALC:
			NK_skip = 1
			while float((getoutput("grep 2pi/SCALE %s -A %d|tail -n %d|tail -1"%(OUTCAR,NK_skip,NK_skip)).split()[-1])) > 0:
				NK_skip += 1
			NK_skip = NK_skip-1
			
		# read k-path and find high-symmetry points
		kpoint_list = getoutput("grep '2pi/SCALE' %s -A %d|tail -n %d"%(OUTCAR,NKPTS_new,NKPTS_new-NK_skip)).split()
		for i in range(NKPTS_new-NK_skip):
			self.Kpoints.append([float(kpoint_list[4*i]),float(kpoint_list[4*i+1]),float(kpoint_list[4*i+2])])
		self.Kpath, self.Special = self.generate_k_path(self.Kpoints)
		##########################################################################################
		# read eigenvalues from EIGENVAL
		if self.flag_read=='EIGENVAL':
			for spin in range(self.ISPIN):
				temp_k = []
				for band in range(iband[0],iband[-1]+1):
					bandlines = getoutput("grep '%5d      ' %s|awk '{print $%d}'"%(band+1,EIGENVAL,spin+2)).split()[-NKPTS_new:]
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
					for band in range(NBANDS_new):
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
		# Reorder the bandstructure according to given k-points index list
		self.Kpoints = self.generate_newobj(self.Kpoints,newlist)
		self.Kpath,self.Special = self.generate_k_path(self.Kpoints)
		for i in range(self.ISPIN):
			self.Energy[i] = self.generate_newobj(self.Energy[i],newlist)
		if self.plot_projection:
			for i in range(self.ISPIN):
				self.Projection[i] = self.generate_newobj(self.Projection[i],newlist)	
	
	def exchange(self,ibanda,ibandb,ik,ispin=0):
		self.Energy[ispin][ik][ibanda],self.Energy[ispin][ik][ibandb]=self.Energy[ispin][ik][ibandb],self.Energy[ispin][ik][ibanda]
		if self.plot_projection:
			self.Projection[ispin][ik][ibanda],self.Projection[ispin][ik][ibandb]=self.Projection[ispin][ik][ibandb],self.Projection[ispin][ik][ibanda]

	def read_bandplot(self,OUTCAR='OUTCAR',PROCAR='PROCAR',EIGENVAL='EIGENVAL'):
		self.load(OUTCAR=OUTCAR,PROCAR=PROCAR,EIGENVAL=EIGENVAL)

def write_plot(bands=[],head='',dir='',iband=[]):
	filenum = 0
	if iband==[]:
		iband=[0,bands.NBANDS-1]
	for bandstruct in bands:
		#print(bandstruct.NBANDS,bandstruct.NKPTS,bandstruct.ISPIN)
		for spin in range(bandstruct.ISPIN):
			filenum += 1
			with open('%s%d.dat'%(head,filenum),'w') as fout:
				for band in range(iband[-1]-iband[0]+1):
					for k in range(bandstruct.NKPTS):
						print('%f\t%f'%(bandstruct.Kpath[k],bandstruct.Energy[spin][k][band].energy),file=fout,end=' ')
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

def write_gnuplot(y0,y1,band_list=[],filename='',outputfilename=''):
	font = 'Times'
	x0 = band_list[0].Special[0][0]
	x1 = band_list[0].Special[-1][0]
	border_width = 4
	arrow_width = 2
	line_width = 4
	xtics_size = 22
	ytics_size = 22
	ylabel_size = 32
	with open('%s.gnu'%filename,'w') as fout:
		print('set term post eps dl 0.4 color enhanced "%s"'%font,file=fout)
		print('set output "%s.eps"'%outputfilename,file=fout)
		print('unset key',file=fout)
		print('set size ratio 1.2',file=fout)
		print('set border lw %d'%border_width,file=fout)
		print('set ylabel "Energy (eV)" font "%s,%d"'%(font,ylabel_size),file=fout)
		print('emin=%f'%y0,file=fout)
		print('emax=%f'%y1,file=fout)
		print('set xrange [%f:%f]'%(x0,x1),file=fout)
		print('set yrange [emin:emax]',file=fout)
		print('set ytics font "%s,%d"'%(font,ytics_size),file=fout)
		print('set xtics () font "%s,%d"'%(font,xtics_size),file=fout)
		for bandstruct in band_list:
			if bandstruct.plot_gap:
				print('set object 1 circle center %f,%f size 0.0004 lc 1 fs transparent solid 0.8 noborder fc "red"'%(bandstruct.VBM[0],bandstruct.VBM[1]),file=fout)
				print('set object 2 circle center %f,%f size 0.0004 lc 1 fs transparent solid 0.8 noborder fc "blue"'%(bandstruct.CBM[0],bandstruct.CBM[1]),file=fout)		
		
		for i in range(len(band_list[0].Special)):
			if band_list[0].Xtics[band_list[0].Special[i][-1]] != "":
				if i==0 or i==len(band_list[0].Special)-1:
					print('set xtics add ("%s" %f)'%(band_list[0].Xtics[band_list[0].Special[i][-1]],band_list[0].Special[i][0]),file=fout)
				elif band_list[0].Special[i][0] != band_list[0].Special[i-1][0]:
					print('set xtics add ("%s" %f)'%(band_list[0].Xtics[band_list[0].Special[i][-1]],band_list[0].Special[i][0]),file=fout)
					print('set arrow from %f,emin to %f,emax nohead front lt 0 lw %d'%(band_list[0].Special[i][0],band_list[0].Special[i][0],arrow_width),file=fout)
		
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

def fat_band_plot(y0,y1,atom1,atom2,orbit,filename,outputfilename):
	band = band_structure(file='OUTCAR.%s'%filename,color=['blue','red'],xtics=['X','{/Symbol G}','Y'])
	band.plot_projection=True
	band.read_bandplot(OUTCAR='OUTCAR.%s'%filename,PROCAR='PROCAR.%s'%filename)
	band.set_projection(atom1=atom1,atom2=atom2,orbital=orbit,color='green')
	band_list = [band]
	write_plot(band_list)
	write_gnuplot(y0=y0,y1=y1,band_list=band_list,outputfilename=outputfilename)
	system('gnuplot<plot.gnu')
	eps_to_png(outputfilename)
	#clean_all(outputfilename)


def plot(y0,y1,E_fermi=-0.76917987,xtics=['{/Symbol G}','M','K','{/Symbol G}',"K\'"],dir='./',filename='',outputfilename='band'):
	band = band_structure(file=dir+'OUTCAR.0',color=['red','red'],xtics=xtics)
	iband=[8000,8075]
	band.set_E_fermi(E_fermi)

	for i in range(9):
		print('Reading file set%d...'%i)
		band.load(OUTCAR=dir+'OUTCAR.%d'%i,EIGENVAL=dir+'EIGENVAL.%d'%i,iband=iband)

	new_klist = [x for x in range(18)]+[x for x in range(26,17,-1)]
	band.rearrange(new_klist)
	band_list = [band]
	write_plot(band_list,filename,iband=iband)
	write_gnuplot(y0=y0,y1=y1,band_list=band_list,filename=filename,outputfilename=outputfilename)
	system('gnuplot<%s.gnu'%filename)
	eps_to_png(outputfilename)
	clean_all(outputfilename)

plot(y0=-0.5,y1=0.5,E_fermi=0.44883190,dir='band/',filename='band',outputfilename='bandall',xtics=['{/Symbol G}',"M","K",'{/Symbol G}'])
