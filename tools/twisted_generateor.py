from numpy import matrix,array,arccos,linalg,degrees,floor,gcd
from copy import deepcopy
class atom_unit:
	element=''
	dynamics=''
	position = []
	def __init__(self,position=[],element='',dynamics=['T','T','T']):
		self.position = position
		self.element = element
		self.dynamics = dynamics

	def read_poslines(self,pos_line):
		pos_list = pos_line.split()
		self.position = [float(pos_list[0]),float(pos_list[1]),float(pos_list[2])]
		if len(pos_list) >= 6:
			self.dynamics =[float(pos_list[3]),float(pos_list[4]),float(pos_list[5])]

	def set_position(self,position):
		self.position=position

	def set_position_c(self,c):
		self.position[2]=c

	def set_element(self,element):
		self.element=element

	def set_dynamics(self,dynamics):
		self.dynamics=dynamics

class atomic_structure:
	file = ''
	title = ''
	frac = 1.0
	basis = []
	elements = []
	natoms = []
	selective_dynamics = 0
	tag = 'Direct'
	atom = []
	nions = 0
	thickness = 0
	def __init__(self,file='POSCAR'):
		with open(file,'r') as fin:
			POSCAR_lines = fin.readlines()
		#print(POSCAR_lines)
		self.title = POSCAR_lines[0][:-2]
		self.frac = float(POSCAR_lines[1].split()[0])
		self.basis = []
		for i in range(3):
			x = float(POSCAR_lines[i+2].split()[0])
			y = float(POSCAR_lines[i+2].split()[1])
			z = float(POSCAR_lines[i+2].split()[2])
			self.basis.append([x,y,z])
		self.elements = POSCAR_lines[5].split()
		self.natoms = []
		for s in POSCAR_lines[6].split():
			self.natoms.append(int(s))
		if POSCAR_lines[7].split()[0][0] in ['S','s']:
			self.selective_dynamics = 1
		if POSCAR_lines[7+self.selective_dynamics].split()[0][0] in ['D','d']:
			self.tag = 'Direct'
		elif POSCAR_lines[7+self.selective_dynamics].split()[0][0] in ['C','c']:
			self.tag = 'Cartesian'
		self.nions = 0
		self.atom = []
		for i in range(len(self.elements)):
			for j in range(self.natoms[i]):
				self.nions += 1
				atom_ij = atom_unit(element=self.elements[i])
				atom_ij.read_poslines(POSCAR_lines[8+self.selective_dynamics+self.nions-1])
				self.atom.append(atom_ij)
		self.standardize()

	def reset(self):
		self.frac = 1.0
		self.basis = []
		self.elements = []
		self.natoms = []
		self.selective_dynamics = 0
		self.tag = 'Direct'
		self.atom = []
		self.nions = 0

	def set_basis(self,basis):
		self.basis=basis

	def atom_in_cell(self,position):
		tol=1e-4
		flag = True
		trans=[0,0,0]
		for i in range(3):
			if abs(position[i])<tol:
				position[i]=0.
			if abs(position[i]-1)<tol:
				position[i]=1.0	
			if not(position[i]>=0 and position[i]<1.):
				flag = False
				trans[i] = int(floor(position[i]))
		return flag,trans

	def add_atom(self,atom):
		if atom.element not in self.elements:
			self.elements.append(atom.element)
			self.natoms.append(0)
		i = self.elements.index(atom.element)
		pos = 0
		for j in range(i+1):
			pos += self.natoms[j]
		self.nions += 1
		self.natoms[i] += 1
		self.atom.insert(pos,atom)

	def del_atom(self,i):
		self.nions -= 1
		index = self.elements.index(self.atom[i].element)
		self.natoms[index] -= 1
		del(self.atom[i])

	def cell_angle(self):
		alpha = round(calc_Angle(self.basis[0],self.basis[2]),6)
		beta = round(calc_Angle(self.basis[1],self.basis[2]),6)
		gamma = round(calc_Angle(self.basis[0],self.basis[1]),6)
		return [alpha,beta,gamma]

	def cell_length(self):
		a = linalg.norm(self.basis[0])
		b = linalg.norm(self.basis[1])
		c = linalg.norm(self.basis[2])
		return [a,b,c]

	def strain(self,strain_list):
		for i in range(3):
			self.basis[i]=(array(self.basis[i])*strain_list[i]).tolist()

	def add_vaccum(self,d):
		old_c = self.basis[-1][-1]
		new_c = old_c + d
		if self.tag=='Direct':
			for i in range(self.nions):
				self.atom[i].set_position_c(self.atom[i].position[-1]*old_c/new_c)
			self.basis[-1][-1]=new_c
		self.standardize()

	def standardize(self):
		c_min = 1
		c_max = -1
		for atom in self.atom:   
			if atom.position[-1]>c_max:
				c_max=deepcopy(atom.position[2])
			if atom.position[-1]<c_min:
				c_min=deepcopy(atom.position[2])
		if c_max-c_min>0.5:
			self.thickness = 1-(c_max-c_min) 
		else:
			self.thickness = c_max-c_min
		center = c_min+self.thickness/2
		#print(c_min,c_max,center,self.thickness)
		for atom in self.atom:
			atom.set_position_c((atom.position[2]+0.5-center)%1)

	def checkcenter(self):
		c_min = 1
		c_max = -1
		for atom in self.atom:   
			if atom.position[-1]>c_max:
				c_max=deepcopy(atom.position[2])
			if atom.position[-1]<c_min:
				c_min=deepcopy(atom.position[2])
		if c_max-c_min>0.5:
			thickness = 1-(c_max-c_min) 
		else:
			thickness = c_max-c_min
		center = c_min+thickness/2
		print(c_min,c_max,center,thickness)

	def transform(self,P,p=[0,0,0]):
		newPOSCAR=deepcopy(self)
		newPOSCAR.reset()
		if linalg.det(array(P)) == 0:
			print("-Error: The determination of the transform matrix is 0!")
		else:
			if linalg.det(array(P)) < 0:
				print("-Warning: The transform matrix changes the coordinate system from right- to left-handed.")
			if linalg.det(array(P)) != 1:
				print("-The transform matrix changes the cell volume.")
			#print(array(P))
			newPOSCAR.set_basis((array(self.basis).T.dot(array(P))).T)
			atom_new = []
			inv_p = linalg.pinv(array(P)) #the atom position vector should dot this
			#################search atoms in the new basis#####################
			stack = [[0,0,0],[1,0,0],[0,1,0],[0,-1,0],[-1,0,0],[0,0,1],[0,0,-1]]
			search_map = []
			while len(stack)!=0:
				n_new  = 0
				for i in range(self.nions):
					newatom = deepcopy(self.atom[i])
					position = ((array(self.atom[i].position)+array(stack[0])).dot(array(inv_p).T)+array(p)).tolist()
					flag,trans = self.atom_in_cell(position)

					if flag:
						newatom.set_position(position)
						newPOSCAR.add_atom(newatom)
						n_new += 1

				if n_new>0:
					for nextcell in [[-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,1],[0,0,-1]]:
						cell = array(nextcell)+array(stack[0])
						if cell.tolist() not in search_map and cell.tolist() not in stack:
							stack.append(cell.tolist())

				search_map.append(stack.pop(0))
		newPOSCAR.standardize()
		return newPOSCAR


	def print_POSCAR(self,file='CONTCAR'):
		with open(file,'wt') as fout:
			print(self.title,file=fout)
			print("%19.14f"%(self.frac),file=fout)
			for i in range(3):
				print(" %22.16f%22.16f%22.16f"%(self.basis[i][0],self.basis[i][1],self.basis[i][2]),file=fout)
			for val in self.elements:
				print("   %-2s"%(val),end='',file=fout)
			print(file=fout)
			for val in self.natoms:
				print("   %-2d"%(val),end='',file=fout)  #vasp POSCAR standard "%6d"
			print(file=fout)
			if self.selective_dynamics == 1:
				print('Selective Dynamics',file=fout)
			print(self.tag,file=fout)
			for i in range(self.nions):
				print("%20.16f%20.16f%20.16f"%(self.atom[i].position[0],self.atom[i].position[1],self.atom[i].position[2]),end='',file=fout)
				if self.selective_dynamics == 1:
					print(" %s %s %s"%(self.atom[i].dynamics[0],self.atom[i].dynamics[1],self.atom[i].dynamics[2]),file=fout)
				else:
					print(file=fout)


def calc_Angle(a,b): #Calculate the angle between vectors <a,b>
	return degrees(arccos(array(a).dot(array(b))/(linalg.norm(a)*linalg.norm(b))))

def calc_Mismatch(lattA,lattB,printflag=True):
	#consider strain was applied on lattB
	cellA = array(lattA.cell_length())
	cellB = array(lattB.cell_length())
	mismatch = (cellA-cellB)/cellB
	if printflag:
		if mismatch[0] - mismatch[1]<1e-8:
			print("%.4f%% biaxial strain should be applied!"%(100*mismatch[0]))
		else:
			print("%.4f%% strain should be applied on a!"%(100*mismatch[0]))
			print("%.4f%% strain should be applied on b!"%(100*mismatch[1]))
		if mismatch[2]!=0:
			print("%.4f%% mismatch was found on c!"%(100*mismatch[2]))
	return mismatch

def choost(a,b):
	if a>b:
		return(b)
	else:
		return(a)

def calculate_strain(mp,nq,a,b,printflag=True):
	a1=(-nq*b**2/(mp))**0.5
	b1=(-mp*a**2/(nq))**0.5
	#print(b**2/a**2)
	delta_a = (a1-a)/a
	delta_b = (b1-b)/b
	if abs(delta_a)<abs(delta_b):
		if printflag:
			print("%.4f%% strain was applied on a"%(100*delta_a))
			return [a1/a,1,1]
		else:
			return 'a',a1/a-1
	else:
		if printflag:
			print("%.4f%% strain was applied on b"%(100*delta_b))
			return [1,b1/b,1]
		else:
			return 'b',b1/b-1


def equals_vec(a,b,tol=1e-8):
	flag=True
	for i in range(len(a)):
		if abs(a[i]-b[i])>tol:
			flag=False
			break
	return flag


def combine_POSCAR(POSCARa,POSCARb,translation=[0,0,0]):
	newPOSCAR=deepcopy(POSCARa)
	#print(POSCARa.cell_angle(),POSCARb.cell_angle(),POSCARa.cell_length(),POSCARb.cell_length())
	if equals_vec(POSCARa.cell_angle(),POSCARb.cell_angle()) and equals_vec(POSCARa.cell_length(),POSCARb.cell_length()):
		for atom in POSCARb.atom:
			newatom=deepcopy(atom)
			position=array(newatom.position)+array(translation)
			newatom.set_position(position.tolist())
			newPOSCAR.add_atom(newatom)
	else:
		print('WARNING: cell angles or sized does not match!! The CONTCAR may be incorrect!!' )
	newPOSCAR.standardize()
	return newPOSCAR


def find_rectangle(MPcut,tolerance,lattice,writeflag=False):
	a = lattice.cell_length()[0]
	b = lattice.cell_length()[1]
	mp = -1
	print("[mp,nq]\t[m,n,p,q]\tAngle(deg)\tNIONS\tstrain")
	while -mp <= MPcut:
		nq = -round(a**2 * mp/b**2)
		axis,strain = calculate_strain(mp,nq,a,b,printflag=False)
		if abs(100*strain) < tolerance:
			#print('%d\t%d\t%f'%(mp,nq,100*strain))
			# find m,n,p,q:
			for m in range(1,-mp+1):
				for q in range(1,nq+1):
					if mp % m == 0 and nq % q == 0:
						p = mp // m
						n = nq // q
						if gcd(p,q) == 1 and gcd(m,n)==1:
							if axis == 'a':
								a_0 = (array(lattice.basis[0])*m*(strain+1) + array(lattice.basis[1])*n).tolist()
								a_1 = (array(lattice.basis[0])*m*(strain+1) - array(lattice.basis[1])*n).tolist()
								angle = calc_Angle(a_0,a_1)
							if axis == 'b':
								b_0 = (array(lattice.basis[0])*p + array(lattice.basis[1])*q*(1+strain)).tolist()
								b_1 = (-array(lattice.basis[0])*p + array(lattice.basis[1])*q*(1+strain)).tolist()
								angle = calc_Angle(b_0,b_1)
							#print(b_0,b_1,lattice.basis[0],lattice.basis[1],strain+1)
							print("[%d,%d]\t[%d,%d,%d,%d]\t%.3f\t%d\t%.4f%%"%(mp,nq,m,n,p,q,angle,(m*q-n*p)*2*lattice.nions,100*strain))
							if writeflag:
								generate_rectangle(m,n,p,q,lattice,filename='%.3f_%d.vasp'%(angle,(m*q-n*p)*2*lattice.nions))
							#print('\t',m,n,p,q,angle,(m*q-n*p)*2*lattice.nions)

		mp -= 1

def find_hexagonal(Mcut,lattice,writeflag=False):
	a = lattice.cell_length()[0]
	for m in range(1,Mcut+1):
		for n in range(1,m):
			if gcd(m,n) == 1:
				v1 = (array(lattice.basis[0])*m + array(lattice.basis[1])*n).tolist()
				v2 = (array(lattice.basis[0])*(m-n) - array(lattice.basis[1])*n).tolist()
				angle = calc_Angle(v1,v2)
				print("%d\t%d\t%.3f\t%d"%(m,n,angle,(m**2+n**2-m*n)*lattice.nions))
				if writeflag:
					generate_hexagonal(m,n,lattice,filename='%.3f_%d.vasp'%(angle,(m**2+n**2-m*n)*3*lattice.nions))

def find_hexagonal_210(Mcut,lattice,writeflag=False):
	a = lattice.cell_length()[0]
	for m in range(1,Mcut+1):
		for n in range(-(m-1),2*m):
			if gcd(m,n) == 1 and m != 2*n:
				v1 = (array(lattice.basis[0])*m + array(lattice.basis[1])*n).tolist()
				v2 = (array(lattice.basis[0])*m - array(lattice.basis[1])*(n-m)).tolist()
				angle = calc_Angle(v1,v2)
				print("%d\t%d\t%.3f\t%d"%(m,n,angle,(m**2+n**2-m*n)*lattice.nions))
				if writeflag:
					generate_hexagonal(m,n,lattice,filename='%.3f_%d.vasp'%(angle,(m**2+n**2-m*n)*3*lattice.nions))
	

def generate_square(m,n,lattice):
	P1=[[n,m,0],[-m,n,0],[0,0,1]]
	P2=[[n,-m,0],[m,n,0],[0,0,1]]
	d=0.3
	build_bilayer(P1,P2,d,lattice)

def build_bilayer(P1,P2,d,lattice,translation=[0.,0.,0.]):
	latticeA=lattice.transform(P1)
	latticeB=lattice.transform(P2).transform([[1,0,0],[0,1,0],[0,0,1]],translation)
	print("The rotation angle is %.3f deg"%(calc_Angle(latticeA.basis[0],latticeB.basis[0])))
	bilayer=combine_POSCAR(latticeA,latticeB,[0,0,d/latticeA.cell_length()[-1]+(latticeA.thickness+latticeB.thickness)/2])
	bilayer.add_vaccum(d+latticeA.thickness*latticeA.cell_length()[-1])
	return bilayer

def build_trilayer_ABA_heter(P1,P2,d,lattA,lattB,translation=[0.,0.,0.]):
	latticeA=lattA.transform(P1)
	latticeB=lattB.transform(P2).transform([[1,0,0],[0,1,0],[0,0,1]],translation)
	latticeB.add_vaccum(latticeA.cell_length()[-1]-latticeB.cell_length()[-1])
	mismatch=calc_Mismatch(latticeA,latticeB)
	latticeB=latticeB.transform([[1+mismatch[0],0,0],[0,1+mismatch[1],0],[0,0,1]],translation)
	bilayer=combine_POSCAR(latticeA,latticeB,[0,0,(d/latticeA.cell_length()[-1]+(latticeA.thickness+latticeB.thickness)/2)])
	bilayer.add_vaccum(d+latticeB.thickness*latticeB.cell_length()[-1])
	latticeB.add_vaccum(d+latticeB.thickness*latticeB.cell_length()[-1])
	trilayer = combine_POSCAR(latticeB,bilayer,[0,0,(d/bilayer.cell_length()[-1]+(bilayer.thickness+latticeB.thickness)/2)])
	trilayer.add_vaccum(d+latticeB.thickness*latticeB.cell_length()[-1])
	return trilayer

def build_bilayer_heter(P1,P2,d,lattA,lattB,translation=[0.,0.,0.],frac=1.0):
	latticeA=lattA.transform(P1)
	latticeB=lattB.transform(P2).transform([[1,0,0],[0,1,0],[0,0,1]],translation)
	latticeB.add_vaccum(latticeA.cell_length()[-1]-latticeB.cell_length()[-1])
	mismatch=calc_Mismatch(latticeA,latticeB)
	latticeB=latticeB.transform([[1+mismatch[0],0,0],[0,1+mismatch[1],0],[0,0,1]],translation)
	bilayer=combine_POSCAR(latticeA,latticeB,[0,0,frac*(d/latticeA.cell_length()[-1]+(latticeA.thickness+latticeB.thickness)/2)])
	bilayer.add_vaccum(d+latticeB.thickness*latticeB.cell_length()[-1])
	print("The rotation angle is %.3f deg"%(calc_Angle(latticeA.basis[0],latticeB.basis[0])))
	return bilayer

def build_trilayer_ABA(P1,P2,d,lattice):
	latticeA=lattice.transform(P1)
	latticeB=lattice.transform(P2)#.transform([[1,0,0],[0,1,0],[0,0,1]],[0,0.5,0.])
	print("The rotation angle is %.3f deg"%(calc_Angle(latticeA.basis[0],latticeB.basis[0])))
	bilayer=combine_POSCAR(latticeA,latticeB,[0,0,d/latticeA.cell_length()[-1]+latticeA.thickness])
	bilayer.add_vaccum(d+latticeA.thickness*latticeA.cell_length()[-1])
	latticeA.add_vaccum(d+latticeA.thickness*latticeA.cell_length()[-1])
	trilayer=combine_POSCAR(bilayer,latticeA,[0,0,d/bilayer.cell_length()[-1]+(bilayer.thickness+latticeA.thickness)/2])
	trilayer.add_vaccum(d+latticeA.thickness*latticeA.cell_length()[-1])
	return trilayer

def build_trilayer_AAB(P1,P2,d,lattice):
	latticeA=lattice.transform(P1)
	latticeB=lattice.transform(P2)#.transform([[1,0,0],[0,1,0],[0,0,1]],[1/3,1/3,0.])
	print("The rotation angle is %.3f deg"%(calc_Angle(latticeA.basis[0],latticeB.basis[0])))
	bilayer=combine_POSCAR(latticeA,latticeA,[0,0,d/latticeA.cell_length()[-1]+latticeA.thickness])
	bilayer.add_vaccum(d+latticeA.thickness*latticeA.cell_length()[-1])
	latticeB.add_vaccum(d+latticeB.thickness*latticeB.cell_length()[-1])
	trilayer=combine_POSCAR(bilayer,latticeB,[0,0,d/bilayer.cell_length()[-1]+(bilayer.thickness+latticeB.thickness)/2])
	trilayer.add_vaccum(d+latticeB.thickness*latticeB.cell_length()[-1])
	return trilayer



def generate_rectangle(m,n,p,q,lattice,d=4,filename='CONTCAR'):
	#############################################################
	# d: distance between two layers with d in unit of angstrom #
	#############################################################
	latt = deepcopy(lattice)
	correction=calculate_strain(m*p,n*q,latt.cell_length()[0],latt.cell_length()[1])
	latt.strain(correction)
	P1=[[m,p,0],[n,q,0],[0,0,1]]
	P2=[[m,-p,0],[-n,q,0],[0,0,1]]
	bilayer=build_bilayer(P1,P2,d,latt)
	bilayer.print_POSCAR(filename)
	#print((m*q-n*p)*4*2)


def generate_hexagonal_ABA(m,n,lattice,d=4,filename='CONTCAR'):
	#############################################################
	# d: distance between two layers with d in unit of angstrom #
	#############################################################
	P1=[[m,-n,0],[n,m-n,0],[0,0,1]]
	P2=[[m-n,n,0],[-n,m,0],[0,0,1]]
	trilayer=build_trilayer_ABA(P1,P2,d,lattice)
	trilayer.print_POSCAR(filename)
	#print((m*q-n*p)*4*2)

def generate_hexagonal_ABA_heter(m,n,latticeA,latticeB,d=4,filename='CONTCAR'):
	#############################################################
	# d: distance between two layers with d in unit of angstrom #
	#############################################################
	P1=[[1,0,0],[0,1,0],[0,0,1]]
	P2=[[m,-n,0],[n,m-n,0],[0,0,1]]
	P3=[[m,n-m,0],[m-n,n,0],[0,0,1]]
	for atom in latticeB.atom:
		atom.set_dynamics(["F","F","F"])
	bilayer=build_bilayer_heter(P1,P2,d,latticeA,latticeB)
	trilayer=build_bilayer_heter(P1,P3,d,bilayer,latticeB,frac=-1.0)
	trilayer.selective_dynamics=1
	trilayer.print_POSCAR(filename)
	print(trilayer.nions)
	#print((m*q-n*p)*4*2)

def generate_hexagonal_ABA_210(m,n,lattice,d=4,filename='CONTCAR'):
	#############################################################
	# d: distance between two layers with d in unit of angstrom #
	#############################################################
	P1=[[m,-n,0],[n,m-n,0],[0,0,1]]
	P2=[[m,n-m,0],[m-n,n,0],[0,0,1]]
	trilayer=build_trilayer_ABA(P1,P2,d,lattice)
	trilayer.print_POSCAR(filename)

def generate_hexagonal_AAB(m,n,lattice,d=4,filename='CONTCAR'):
	#############################################################
	# d: distance between two layers with d in unit of angstrom #
	#############################################################
	P1=[[m,-n,0],[n,m-n,0],[0,0,1]]
	P2=[[m-n,n,0],[-n,m,0],[0,0,1]]
	trilayer=build_trilayer_AAB(P1,P2,d,lattice)
	trilayer.print_POSCAR(filename)
	#print((m*q-n*p)*4*2)

def generate_hexagonal_AAB_210(m,n,lattice,d=4,filename='CONTCAR'):
	#############################################################
	# d: distance between two layers with d in unit of angstrom #
	#############################################################
	P1=[[m,-n,0],[n,m-n,0],[0,0,1]]
	P2=[[m,n-m,0],[m-n,n,0],[0,0,1]]
	trilayer=build_trilayer_AAB(P1,P2,d,lattice)
	trilayer.print_POSCAR(filename)

def generate_hexagonal(m,n,lattice,d=4,translation=[0.,0.,0.],filename='CONTCAR'):
	#############################################################
	# d: distance between two layers with d in unit of angstrom #
	#############################################################
	P1=[[m,-n,0],[n,m-n,0],[0,0,1]]
	P2=[[m-n,n,0],[-n,m,0],[0,0,1]]
	bilayer=build_bilayer(P1,P2,d,lattice,translation)
	bilayer.print_POSCAR(filename)

def generate_hexagonal_210(m,n,lattice,d=4,filename='CONTCAR'):
	#############################################################
	# d: distance between two layers with d in unit of angstrom #
	#############################################################
	P1=[[m,-n,0],[n,m-n,0],[0,0,1]]
	P2=[[m,n-m,0],[m-n,n,0],[0,0,1]]
	bilayer=build_bilayer(P1,P2,d,lattice)
	bilayer.print_POSCAR(filename)

def find_mn(Mcut,structA,structB):
	a = structA.cell_length()[0]
	min_trans = [0.001,0,0]
	for m in range(1,Mcut+1):
		for n in range(1,m):
			v = array(structB.basis[0])*m + array(structB.basis[1])*n
			b = linalg.norm(v)
			if abs(b-a)/a < min_trans[0]:
				min_trans[0] = abs(b-a)/a
				min_trans[1],min_trans[2] = m,n
				print(*min_trans)
	if min_trans[1:] != [0,0]:
		print(*min_trans)
	else:
		print("Not found!")

structure=atomic_structure('monolayer.vasp')
#generate_hexagonal(13,6,structure,filename='bilayer_13_6.vasp')
#generate_hexagonal(11,1,structure,filename='bilayer_11_1.vasp')
generate_hexagonal_210(7,3,structure,d=4,filename='bilayer_210_7_3.vasp')
#generate_hexagonal_210(15,7,structure,d=4,filename='bilayer_210_15_7.vasp')
#generate_hexagonal_210(17,8,structure,d=4,filename='bilayer_210_17_8.vasp')
#generate_hexagonal_210(21,10,structure,d=4,filename='bilayer_210_21_10.vasp')
#generate_hexagonal_ABA_210(13,6,structure,d=4,filename='trilayer_ABA_210_13_6.vasp')



#tWSe2 = atomic_structure('1.00.vasp')
#hBN = atomic_structure('hBN.vasp')
#find_mn(20,tWSe2,hBN)
#m=17
#n=10
#generate_hexagonal_ABA_heter(17,7,tWSe2,hBN,d=3,filename="test.vasp")