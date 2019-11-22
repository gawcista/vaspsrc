import re,os,subprocess

def getnum(s):
	return re.findall(r"[-+]?[0-9]*\.?[0-9]+",s)

def read_para_OUTCAR():
	global FERMI,SOC,ISPIN,NKPTS,NBANDS,NIONS
	FERMI = float(getnum(subprocess.getoutput('grep E-fermi '+'OUTCAR'))[0])
	SOC = bool(subprocess.getoutput('grep LSORBIT '+'OUTCAR')[18])
	ISPIN = int(getnum(subprocess.getoutput('grep ISPIN '+'OUTCAR'))[0])
	NKPTS = int(getnum(subprocess.getoutput('grep NKPTS '+'OUTCAR'))[0])
	NBANDS = int(getnum(subprocess.getoutput('grep NBANDS '+'OUTCAR'))[0])
	NIONS = int(getnum(subprocess.getoutput('grep NIONS '+'OUTCAR'))[1])

def read_kpath_OUTCAR():
	global NKPTS,KPOINTS,Kpath,Route
	NKPTS = int(getnum(subprocess.getoutput('grep NKPTS '+'OUTCAR'))[0])
	K_lines = subprocess.getoutput('grep -A '+str(NKPTS)+' 2pi/SCALE '+'OUTCAR'+'|tail -'+str(NKPTS)).split('\n')
	KPOINTS = []
	Route = []
	k = 0.
	Route.append(k)
	for i in range(NKPTS):
		kx = float(getnum(K_lines[i])[0])
		ky = float(getnum(K_lines[i])[1])
		kz = float(getnum(K_lines[i])[2])
		#occ = float(getnum(K_lines[i])[3])
		KPOINTS.append([kx,ky,kz])
		if i>0:
			k+=sqrt(sum((array(K[-1])-array(K[-2]))**2))
			Route.append(k)