def print_KPOINTS_line(k_i,k_f,nk,title):
	with open(title,'w') as fout:
		print("Automesh",file=fout)
		print("%d ! number of intersections"%nk,file=fout)
		print("Line-mode",file=fout)
		print("rec",file=fout)
		print("%f %f %f"%(k_i[0],k_i[1],k_i[2]),file=fout)
		print("%f %f %f"%(k_f[0],k_f[1],k_f[2]),file=fout)

def print_KPOINTS(k_list,title):
	weight = 1./len(k_list)
	with open(title,'w') as fout:
		print("Automesh",file=fout)
		print("%d ! number of kpt"%len(k_list),file=fout)
		print("Reciprocal lattice",file=fout)
		for kpt in k_list:
			print("%f %f %f %f"%(kpt[0],kpt[1],kpt[2],weight),file=fout)	

def generate_mesh(nka,nkb,nkc,k_i,k_f,nk,title):
	ktemp0 = k_i
	k_list=[]
	nf = 0
	for i in range(nka):
		if nka == 1:
			ka = k_i[0]+(k_f[0]-k_i[0])*i/nka
		else:
			ka = k_i[0]+(k_f[0]-k_i[0])*i/(nka-1)
		for j in range(nkb):
			if nkb == 1:
				kb = k_i[1]+(k_f[1]-k_i[1])*j/nkb
			else:
				kb = k_i[1]+(k_f[1]-k_i[1])*j/(nkb-1)
			for k in range(nkc):
				if nkc == 1:
					kc = k_i[2]+(k_f[2]-k_i[2])*k/nkc
				else:
					kc = k_i[2]+(k_f[2]-k_i[2])*k/(nkc-1)
				k_list.append([ka,kb,kc])
	while k_list!=[]:
		list_temp=[]
		for i in range(nk):
			list_temp.append(k_list.pop(0))
		print_KPOINTS(list_temp,title+str(nf))
		nf +=1
		pass


def generate_line(nklines,k_i,k_f,nk,title):
	ktemp0 = k_i
	nf = -1
	if nklines == 1:
		print_KPOINTS_line(k_i,k_f,nk,title)
	else:
		for i in range(nklines):
			ka = k_i[0]+(k_f[0]-k_i[0])*i/(nklines-1)
			kb = k_i[1]+(k_f[1]-k_i[1])*i/(nklines-1)
			kc = k_i[2]+(k_f[2]-k_i[2])*i/(nklines-1)
			ktemp1 = [ka,kb,kc]
			if nf>=0:
				print_KPOINTS_line(ktemp0,ktemp1,nk,title+str(nf))
			nf += 1
			ktemp0 = ktemp1
			

generate_mesh(11,11,1,[-0.5,0,0],[0.5,0.5,0],11,'KPOINTS.mesh')
