from subprocess import getoutput
from numpy import array,zeros

ORBITAL=['s','py','pz','px','dxy','dyz','dz2','dxz','dx2-y2','tot']

def distance(a,b):
# calculate distance between two sites
    dis = .0
    for i in range(3):
        dis += (a[i]-b[i])**2
    return dis**0.5

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

class band_projection:
    color = 'red'
    outputfilename = ''
    label = ''
    proj = []
    def __init__(self,Projection,ion_list=[],orb_list=[],color='red',label='',channel=0):
        self.color = color
        self.label = label
        self.proj = self.generateProjection(Projection,ion_list,orb_list,channel=channel)

    def generateProjection(self,Projection,ion_list=[],orb_list=[],channel=0):
        proj = zeros((Projection.shape[0],Projection.shape[1],Projection.shape[2]))
        if ion_list == []:
            ion_list = [-1]
        if orb_list == [-1]:
            orb_list = list(range(Projection.shape[5]))
        for ion in ion_list:
            for orb in orb_list:
                proj += Projection[:,:,:,channel,ion,orb]
        return proj

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
    LMAXMIX = 2
    E_fermi = 0
    Xtics = []
    VBM = [-100,0,0]
    CBM = [100,0,0]
    plot_gap = False
    plot_projection = False
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
        self.Projection = []
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
        self.ISPIN = int(getoutput("grep ISPIN %s|tail -1"%self.file).split()[2])
        self.E_fermi = float((getoutput("grep E-fermi %s"%self.file)).split()[2])
        if int(getoutput("grep 'METAGGA =' %s|wc -l"%self.file)) > 0:
            self.METAGGA = str(getoutput("grep 'METAGGA =' %s|tail -1"%self.file).split()[2])
        self.LHFCALC = str(getoutput("grep 'LHFCALC =' %s|tail -1"%self.file).split()[2])=='T'
        self.LSORBIT = str(getoutput("grep 'LSORBIT =' %s|tail -1"%self.file).split()[2])=='T'
        self.LMAXMIX = int(getoutput("grep LMAXMIX %s|tail -1"%self.file).split()[2])
        self.flag_read = 'EIGENVAL'
        #self.autoKLABEL = False
        self.plotmode = "line"
        self.label = ""
        self.bandcut = 100.

    def set_E_fermi(self,E):
        self.E_fermi = E
    
    def set_plotmode(self,mode):
        self.plotmode = mode

    def set_plot_gap(self,flag):
        self.plot_gap = flag

    def set_xtics(self,xtics):
        self.Xtics = xtics
    
    def set_label(self,label):
        self.label = label

    def set_plot_projection(self,flag):
        self.plot_projection = flag
    
    def set_projection(self,ion_list=[],orb_list=[],channel=0,color='red',label=''):
        self.projection_list.append(band_projection(self.Projection,ion_list,orb_list,channel=channel,color=color,label=label))
    
    def generate_newobj(self,obj,indexlist):
        # reorder certain object list according to given index list
        obj_new = []
        for i in indexlist:
            obj_new.append(obj[i])
        return array(obj_new)

    def generate_k_path(self,Kpoints,tol=1e-6):
        # generate Kpath and Special for given Kpoints list
        route = 0.
        Kpath = [0.]
        Special = []
        tmp = 0.
        NSPECIAL = 0
        for i in range(1,len(Kpoints)):
            dr = distance(Kpoints[i],Kpoints[i-1])
            if abs(dr-tmp) > tol:
                Special.append([Kpath[i-1],i-1,NSPECIAL])
                if abs(dr) > self.bandcut*abs(tmp) and abs(tmp) > tol:
                    NSPECIAL += 1
                    dr = 0
                if abs(dr) > tol:
                    NSPECIAL += 1
            route += dr
            tmp = dr
            Kpath.append(route)
        Special.append([Kpath[-1],len(Kpoints)-1,NSPECIAL])
        return Kpath,Special
          

    def load(self,OUTCAR='band/OUTCAR',PROCAR='band/PROCAR',EIGENVAL='band/EIGENVAL',iband=[],bandcut=10): # add band in OUTCAR, along k-axis
        self.bandcut = bandcut
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
        NKPTS_new = NKPTS_new - NK_skip 
        # read k-path and find high-symmetry points
        kpoint_list = getoutput("grep '2pi/SCALE' %s -A %d|tail -n %d"%(OUTCAR,NKPTS_new+NK_skip,NKPTS_new)).split()
        for i in range(NKPTS_new):
            self.Kpoints.append([float(kpoint_list[4*i]),float(kpoint_list[4*i+1]),float(kpoint_list[4*i+2])])
        self.Kpath, self.Special = self.generate_k_path(self.Kpoints)
        ##########################################################################################
        # read eigenvalues from EIGENVAL
        if self.flag_read=='EIGENVAL':
            for spin in range(self.ISPIN):
                temp_k = []
                for band in range(iband[0],iband[-1]+1):
                    bandlines = getoutput("grep '%5d      ' %s|awk '{print $%d}'"%(band+1,EIGENVAL,spin+2)).split()[-NKPTS_new:]
                    for k in range(NKPTS_new):
                        bandpoint=band_point(energy=float(bandlines[k])-self.E_fermi)
                        if k>=len(temp_k):
                            temp_k.append([bandpoint])
                        else:
                            temp_k[k].append(bandpoint)
                        if bandpoint.energy>0 and bandpoint.energy<self.CBM[1]:
                            self.CBM = [self.Kpath[k],bandpoint.energy,k,band]
                        if bandpoint.energy<0 and bandpoint.energy>self.VBM[1]:
                            self.VBM = [self.Kpath[k],bandpoint.energy,k,band]
                for k_list in temp_k:
                    self.Energy[spin].append(k_list)

        # read projections from PROCAR
        if self.plot_projection:
            if self.LSORBIT:
                channels = 4
            else:
                channels = 1
            projlines=[]
            for band in range(iband[0],iband[-1]+1):
                projlines += getoutput("grep 'band%6d' %s -A %d|awk '!/ion|--|band/ {$1=\"\";print $0}'"%(band+1,PROCAR,2+channels*(self.NIONS+1))).split()
            self.Projection = array(projlines).astype(float).reshape((self.ISPIN,iband[-1]-iband[0]+1,NKPTS_new-NK_skip,channels,self.NIONS+1,10)).transpose((0,2,1,3,4,5))      
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
    
    def read_bandplot(self,OUTCAR='OUTCAR',PROCAR='PROCAR',EIGENVAL='EIGENVAL'):
        self.load(OUTCAR=OUTCAR,PROCAR=PROCAR,EIGENVAL=EIGENVAL)

def write_plot(bands=[],head='band',dir=''):
    ###### to be reconstructed
    filenum = 0
    for bandstruct in bands:
        for spin in range(bandstruct.ISPIN):
            filenum += 1
            with open('%s%d.dat'%(head,filenum),'w') as fout:
                for band in range(array(bandstruct.Energy[spin]).shape[1]):
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

def plot_matplotlib(band_list,erange=[-1.0,1.0],outputfile="band.png",title="",plot_legend=False,**kargs):
    def set_params(kargs):
        default_arg = {}
        default_arg['border_width'] = 3
        default_arg['fontfamily'] = ['Arial']
        default_arg['fontsize'] = 18
        default_arg['ylabel_fontsize'] = 24
        default_arg['figure_figsize'] = (4,4)
        default_arg['figure_autolayout'] = True
        default_arg['band_linewidth'] = 1.5
        default_arg['band_markersize'] = 5
        default_arg['projection_weight'] = 50
        default_arg['projection_alpha'] = 0.3
        for key,value in kargs.items():
            default_arg[key] = value
        return default_arg
    import matplotlib.pyplot as plt
    kpath = band_list[0].Kpath
    param = set_params(kargs)
    plt.rcParams['font.family'] = param['fontfamily']
    plt.rcParams['font.size'] = param['fontsize']
    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 8
    plt.rcParams['figure.autolayout'] = param['figure_autolayout']
    plt.rcParams['figure.figsize'] = param['figure_figsize']
    fig, ax = plt.subplots(1, 1)
    ax.spines['bottom'].set_linewidth(param['border_width'])
    ax.spines['left'].set_linewidth(param['border_width'])
    ax.spines['top'].set_linewidth(param['border_width'])
    ax.spines['right'].set_linewidth(param['border_width'])
    ax.tick_params(direction='in', width=param['border_width'])
    ax.set_xlim(kpath[0],kpath[-1])
    ax.set_ylim(erange[0],erange[1])
    ax.set_ylabel('Energy (eV)',fontsize=param['ylabel_fontsize'])
    xticks_position=[]
    xticks_label=[]
    xticks_position.append(band_list[0].Special[0][0])
    xticks_label.append(band_list[0].Xtics[band_list[0].Special[0][2]])
    for n in range(1,len(band_list[0].Special)):
        if band_list[0].Xtics[band_list[0].Special[n][2]] !="":
            if band_list[0].Special[n-1][0] == band_list[0].Special[n][0] and  band_list[0].Special[n-1][2] != band_list[0].Special[n][2]:
                xticks_position.append(band_list[0].Special[n][0])
                ax.vlines(band_list[0].Special[n][0],erange[0],erange[1],linewidth=1.0,edgecolor="black")
                xticks_label.append("%s|%s"%(band_list[0].Xtics[band_list[0].Special[n-1][2]],band_list[0].Xtics[band_list[0].Special[n][2]]))
            elif band_list[0].Special[n-1][0] != band_list[0].Special[n][0]:
                xticks_position.append(band_list[0].Special[n][0])
                ax.vlines(band_list[0].Special[n][0],erange[0],erange[1],linewidth=1.0,edgecolor="gray",linestyles="--")
                xticks_label.append(band_list[0].Xtics[band_list[0].Special[n][2]])



    ax.set_xticks(xticks_position)
    ax.set_xticklabels(xticks_label)
    for bandstruct in band_list:
        kpath = bandstruct.Kpath
        for spin in range(bandstruct.ISPIN):
            energy = array(bandstruct.Energy[spin]).T
            for iband in range(len(energy)):
                band = [e.energy for e in energy[iband]]
                if bandstruct.plotmode == "line":
                    for n in range(1,len(band_list[0].Special)):
                        if band_list[0].Special[n][0] != band_list[0].Special[n-1][0]:
                            ki = band_list[0].Special[n-1][1]
                            kf = band_list[0].Special[n][1]
                            line,=ax.plot(kpath[ki:kf+1],band[ki:kf+1],color=bandstruct.color[spin],linewidth=param['band_linewidth'])
                elif bandstruct.plotmode == "scatter":
                    line,=ax.plot(kpath,band,marker='o',ms=param['band_markersize'],mec=bandstruct.color[spin],color='none')

                if bandstruct.plot_projection and len(bandstruct.projection_list)>0:
                    for projection in bandstruct.projection_list:
                        proj = param['projection_weight']*projection.proj[spin].T[iband]
                        #print(proj)
                        ax.scatter(kpath,band,proj,alpha=param['projection_alpha'],marker='o',color=projection.color)
        if bandstruct.label != "":
            line.set_label(bandstruct.label)

    if title!= "":
        ax.set_title(title,pad=15)

    if "lines_sup" in param:
        for line in param['lines_sup']:
            ax.plot(line["x"],line["y"],color=line["color"],linewidth=line["linewidth"])

    if plot_legend:
        ax.legend()
    ax.hlines(0.,kpath[0],kpath[-1],linewidth=1.0,edgecolor="gray",linestyles="--")
    plt.savefig(outputfile,transparent=True,dpi=600)
    plt.close()


def quickplot(erange, xtics=[" "," "," "," "], E_fermi=0., dir='./', outputfilename='band', title="",**kargs):
    band = band_structure(file=dir+'OUTCAR',color=['tab:blue','tab:red'],xtics=xtics)
    band.set_E_fermi(E_fermi)
    band.load(OUTCAR=dir+'OUTCAR',EIGENVAL=dir+'EIGENVAL')
    plot_matplotlib(band_list=[band],erange=erange,outputfile=outputfilename+".png",title=title,**kargs)


def quickplot_mult(erange, xtics=[" ", " ", " ", " "], E_fermi=0., dir='./', outputfilename='band', title="", nfile=1, iband=[]):
    band = band_structure(
        file=dir+'OUTCAR.0', color=['tab:blue', 'tab:red'], xtics=xtics)
    band.set_E_fermi(E_fermi)
    for i in range(nfile):
        band.load(OUTCAR=dir+'OUTCAR.%d' %
                  i, EIGENVAL=dir+'EIGENVAL.%d' % i, iband=iband)
    plot_matplotlib(band_list=[band], erange=erange,
                    outputfile=outputfilename+".png", title=title)

def quickplot_proj(erange, ion_list, orb_list, xtics=[" ", " ", " ", " "], E_fermi=0., dir='./', outputfilename='band', title="",newlist=[],label=""):
    band = band_structure(
        file=dir+'OUTCAR', color=['tab:grey', 'tab:grey'], xtics=xtics)
    band.set_E_fermi(E_fermi)
    band.set_plot_projection(True)
    band.load(OUTCAR=dir+'OUTCAR',EIGENVAL=dir+'EIGENVAL',PROCAR=dir+"PROCAR")
    if newlist !=[]:
        band.rearrange(newlist)
    band.set_projection(ion_list=ion_list,orb_list=orb_list,color="tab:red",label=label)
    plot_matplotlib(band_list=[band],erange=erange,outputfile=outputfilename+".png",title=title,band_linewidth=0.5)

def quickplot_proj_mult(erange, ion_list, orb_list, xtics=[" ", " ", " ", " "], E_fermi=0., dir='./', outputfilename='band', title="",newlist=[],label="",nfile=1):
    band = band_structure(file=dir+'OUTCAR.0', color=['tab:grey', 'tab:grey'], xtics=xtics)
    band.set_E_fermi(E_fermi)
    band.set_plot_projection(True)
    for i in range(nfile):
        band.load(OUTCAR=dir+'OUTCAR.%d'%i,EIGENVAL=dir+'EIGENVAL.%d'%i,PROCAR=dir+"PROCAR.%d"%i)
    if newlist !=[]:
        band.rearrange(newlist)
    band.set_projection(ion_list=ion_list,orb_list=orb_list,color="tab:red",label=label)
    plot_matplotlib(band_list=[band],erange=erange,outputfile=outputfilename+".png",title=title,band_linewidth=0.5)

def quickplot_proj_auto(erange, xtics=[" ", " ", " ", " "], E_fermi=0., dir='./', outputfilename='band', POSCAR='POSCAR',newlist=[]):
    elements = getoutput("awk 'NR==6' %s" % (dir+POSCAR)).split()
    ions = getoutput("awk 'NR==7' %s" % (dir+POSCAR)).split()
    flag = 0
    print(elements,ions)
    for i in range(len(elements)):
        ion_list = range(flag,flag+int(ions[i]))
        for orb in range(10):
            print("Plotting projection on %s %s"%(elements[i],ORBITAL[orb]))
            quickplot_proj(erange=erange, ion_list=ion_list, orb_list=[orb], xtics=xtics, E_fermi=E_fermi, dir=dir, outputfilename=outputfilename+'.%s.%s'%(elements[i],ORBITAL[orb]), title="%s %s"%(elements[i],ORBITAL[orb]),newlist=newlist)
        flag += int(ions[i])
