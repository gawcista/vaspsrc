from numpy import array,zeros,sum

def read_DOSCAR(file="DOSCAR"):
    with open(file,"r") as fin:
        dat = fin.readlines()
    NIONS = int(dat[0].split()[0])
    NGRID = int(dat[5].split()[2]) 
    E_fermi = float(dat[5].split()[3])
    dat_tot = array([x.split() for x in dat[6:NGRID+6]]).astype(float)
    E = dat_tot[1:,0]-E_fermi
    DOS = []
    DOS.append(dat_tot[:,1])
    for i in range(NIONS):
        k = (i+1)*(NGRID+1)+6
        DOS.append(array([x.split()[1:] for x in dat[k+1:NGRID+k]]).astype(float))
    return E,DOS

def gen_atom(atom_list=[],file="DOSCAR"):
    E, dat = read_DOSCAR(file)
    DOS = []
    for l in atom_list:
        tmp = zeros((dat[1].shape[0],dat[1].shape[1]))
        for n in l:
            tmp += dat[n]
        DOS.append(sum(tmp,axis=1))
    return E,DOS


def plot(E, dos_list, erange, color=["tab:blue", "tab:red"], outputfile="dos.png", title="", plot_legend=False, **kargs):
    def set_params(kargs):
        default_arg = {}
        default_arg['border_width'] = 3
        default_arg['fontfamily'] = ['Arial']
        default_arg['fontsize'] = 18
        default_arg['label_fontsize'] = 24
        default_arg['figure_figsize'] = (4, 4)
        default_arg['figure_autolayout'] = True
        default_arg['band_linewidth'] = 1.5
        default_arg['band_markersize'] = 5
        default_arg['projection_weight'] = 50
        default_arg['projection_alpha'] = 0.3
        for key, value in kargs.items():
            default_arg[key] = value
        return default_arg
    import matplotlib.pyplot as plt
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
    ax.set_xlim(erange[0], erange[-1])
    ax.set_xlabel('Energy (eV)',fontsize=param['label_fontsize'])
    if title != "":
        ax.set_title(title, pad=15)
    if plot_legend:
        ax.legend()
    plt.savefig(outputfile, transparent=True, dpi=600)

gen_atom(atom_list=[[1]], file="/mnt/d/Working/CuHAT/dos/DOSCAR.Cu4")
