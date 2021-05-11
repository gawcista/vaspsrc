import band as h
def plot(y0,y1,E_fermi=-0.76917987,xtics=['{/Symbol G}','X','S','Y','{/Symbol G}','S'],dir='./',filename='',outputfilename='band'):
	#a=[i for i in range(1816,1820,1)]
	a=[i for i in range(3632,3640,1)]
	band = h.band_structure(file=dir+'OUTCAR.0',selected_bands=a,color=['blue','red'],xtics=xtics)
	band.set_E_fermi(E_fermi)
	#for i in range(25):
	#	print('Reading file set%d...'%i)
	#	band.read_bandplot(OUTCAR=dir+'OUTCAR.%d'%i,EIGENVAL=dir+'EIGENVAL.%d'%i)
	for i in range(5):
		band.read_bandplot(OUTCAR=dir+'OUTCAR.%d'%i,EIGENVAL=dir+'EIGENVAL.%d'%i)
	#band.rearrange([-4,-3,-5,1,2])
	band_list = [band]
	h.write_plot(band_list,filename)
	#h.print_plot(band_list)
	#print(len(band.Special))
	h.write_gnuplot(y0=y0,y1=y1,band_list=band_list,filename=filename,outputfilename=outputfilename)
	h.system('gnuplot<%s.gnu'%filename)
	h.eps_to_png(outputfilename)
	#clean_all(outputfilename)
def band_plot_print(banda,bandb,nfiles,E_fermi=-0.76917987,xtics=['{/Symbol G}','X','S','Y','{/Symbol G}','S'],dir='./'):
	a=[i for i in range(banda,bandb,1)]
	band = h.band_structure(file=dir+'OUTCAR.0',selected_bands=a,color=['blue','red'],xtics=xtics)
	band.set_E_fermi(E_fermi)
	for i in range(nfiles):
		band.read_bandplot(OUTCAR=dir+'OUTCAR.%d'%i,EIGENVAL=dir+'EIGENVAL.%d'%i)
	band.rearrange([-4,-3,-5,1,2])
	band_list = [band]
	h.print_plot(band_list)
band_plot_print(1456,1460,5,E_fermi=-0.70343092,dir='SnS.13.443_new/')
#plot(y0=-0.262,y1=-0.25,E_fermi=-0.77080878,dir='SnS.6.03.new/scp/',filename='SnS.6.03.new',outputfilename='SnS.6.03.new')
#plot(y0=-0.31,y1=-0.23,E_fermi=-0.79738994,xtics=['{/Symbol G}','X','S','Y','{/Symbol G}','S'],dir='SnS.12.03.p2/',filename='SnS.12.03.p2',outputfilename='SnS.12.03.p2.1')
#plot(y0=-0.335,y1=-0.325,E_fermi=-.76723717,xtics=['S','G'],dir='SnS.12.03.SG/',filename='SnS.12.03.SG',outputfilename='SnS.12.03.SG.1')
#plot(y0=-0.285,y1=-0.265,E_fermi=-.76723717,xtics=['S','G'],dir='SnS.12.03.SG/',filename='SnS.12.03.SG',outputfilename='SnS.12.03.SG.2')
#plot(y0=-0.335,y1=-0.265,E_fermi=-.76723717,xtics=['C','D'],dir='SnS.12.03.line1/',filename='SnS.12.03.line1',outputfilename='SnS.12.03.line1.1')
#plot(y0=-0.335,y1=-0.265,E_fermi=-.76723717,xtics=['A','B'],dir='SnS.12.03.line2/',filename='SnS.12.03.line2',outputfilename='SnS.12.03.line2.1')
plot(y0=-0.278,y1=-0.274,E_fermi=-.76723717,xtics=['S','G'],dir='SnS.12.03.0.4XS/',filename='SnS.12.03.0.4XS',outputfilename='SnS.12.03.0.4XS.1')
#plot(y0=-0.332,y1=-0.329,E_fermi=-.76723717,xtics=['S','G'],dir='SnS.12.03.0.4XS/',filename='SnS.12.03.0.4XS',outputfilename='SnS.12.03.0.4XS.2')
#plot(y0=-0.285,y1=-0.265,E_fermi=-.76723717,xtics=['S','G'],dir='SnS.12.03.SG/',filename='SnS.12.03.SG',outputfilename='SnS.12.03.SG.2')