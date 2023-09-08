from numpy import array
def read_DOSCAR(file="DOSCAR"):
    with open(file,"r") as fin:
        dat = fin.readlines()
    NIONS = int(dat[0].split()[0])
    NGRID = int(dat[5].split()[2]) 
    E_fermi = float(dat[5].split()[3])
    dat_tot = array([x.split() for x in dat[6:NGRID+6]]).astype(float)
    E = dat_tot[:,0]-E_fermi
    DOS = []
    DOS.append(dat_tot[:,1])
    for i in range(NIONS):
        pass


read_DOSCAR("/mnt/d/Working/CuHAT/dos/DOSCAR.Cu4")
