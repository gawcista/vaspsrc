from numpy import array,loadtxt
from subprocess import getoutput

def read_band(file="wannier90_band.dat",fermi=0.):
    dat = loadtxt(file, dtype=float).T
    return dat[0],dat[1]-fermi

def read_nkpts(file="wannier90_band.kpt"):
    return int(getoutput("head wannier90_band.kpt -n 1"))

def read_label(file="wannier90_band.labelinfo.dat"):
    label = loadtxt(file, dtype=str, usecols=0)
    index = loadtxt(file, dtype=int, usecols=1)
    k = loadtxt(file, dtype=float, usecols=2)
    return {"label":label,"index":index,"k":k}

def plot(erange=[-1.0,1.0],fermi=0.,input="wannier90",outputfile="band_wan.png",title="",**kargs):
    def set_params(kargs):
        default_arg = {}
        default_arg['border_width'] = 3
        default_arg['fontfamily'] = ['Arial']
        default_arg['fontsize'] = 18
        default_arg['ylabel_fontsize'] = 24
        default_arg['figure_figsize'] = (4, 4)
        default_arg['figure_autolayout'] = True
        default_arg['band_linewidth'] = 1.5
        default_arg['band_color'] = "tab:blue"
        default_arg['band_markersize'] = 5
        default_arg['projection_weight'] = 50
        default_arg['projection_alpha'] = 0.3
        default_arg['background_transparent'] = True
        for key, value in kargs.items():
            default_arg[key] = value
        return default_arg
    # load data
    kpath_raw, band_raw = read_band(file=input+"_band.dat", fermi=fermi)
    nkpts = read_nkpts(file=input+"_band.kpt")
    special = read_label(file=input+"_band.labelinfo.dat")
    kpath = kpath_raw.reshape(-1, nkpts)
    band = band_raw.reshape(-1, nkpts)
    nbands = band.shape[0]
    # set plot sytle
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
    ax.set_xlim(kpath[0][0], kpath[0][-1])
    ax.set_ylim(erange[0], erange[1])
    ax.set_ylabel('Energy (eV)', fontsize=param['ylabel_fontsize'])
    ax.set_xticks(special["k"])
    ax.set_xticklabels(special["label"])
    if title != "":
        ax.set_title(title, pad=15)
    ax.hlines(0., kpath[0][0], kpath[0][-1], linewidth=1.0,
              edgecolor="gray", linestyles="--")
    for k in special["k"]:
        if k > kpath[0][0] and k < kpath[0][-1]:
            ax.vlines(k, erange[0], erange[1],
                      linewidth=1.0, edgecolor="gray", linestyles="--")
    for i in range(nbands):
        ax.plot(kpath[i], band[i],
                color=param['band_color'], linewidth=param['band_linewidth'])
    plt.savefig(
        outputfile, transparent=param['background_transparent'], dpi=600)

