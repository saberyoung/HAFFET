from sdapy import snerun, constants, functions
import numpy as np
import matplotlib.pyplot as plt
import os, sys

    
if False:
    a = snerun.snelist()
    a.run(syntax='type in ["Ib"]')
    a.add_subset(syntax='z>.01')
    a.add_parameter(parname='mag', filt1='r', phase1=0)
    a.add_parameter(parname='color', filt1='g', filt2='r', phase1=10)
    a.add_parameter(parname='mni', model='Arnett', quant=[.05,.5,.95])
    print (a.table(syntax='z>.01'))
if False:
    a = snerun.snobject(objid='ZTF20aajcdad', aliasid='SN2020bcq',
        z=0.0186, dm=34.6, mkwebv=0.01387, hostebv=0, sntype='Ib', 
        ra='13:26:29.65', dec='+36:00:31.1', jdpeak=2458888.02)
    #a.run()
    a.get_fp_ztf()
    print (a.lc)

#from scipy.optimize import curve_fit
#def func1(x, a, b, sntype='Ib'):
#    if sntype=='Ib': c=10
#    else: c=1
#    print (c)
#    return a * np.exp(-b * x) + c
#xdata = np.linspace(0, 4, 50)
#y = func1(xdata, 2.5, 1.3, 0.5)
#ydata = y + 0.2 * np.random.normal(size=len(xdata))
#func = lambda x, a, b:func1(x, a, b, sntype='Ib')
#popt, pcov = curve_fit(func, xdata, ydata)
#print ( popt, pcov )
#sys.exit()

fig,ax = plt.subplots(1,1,figsize=(6,4))
colors=dict()
colors['r'] = 'r'
colors['g'] = 'g'
colors['i'] = 'k'


if False:
    fig, ax = plt.subplots(1,1)
    sncls = snerun.snelist(ax=ax)
    sncls.parse_meta(withnew='skip', source='BTS', debug=True,
                     syntax='type in ["Ib"]', metafile=None)
    sncls.run(lctype= ['ztffp'], spectype=['marshal', 'fritz'], 
              fit_methods=['bazin1', 'arnett_fit_taum', 'gauss', 'exp'])
    #for objid in sncls.meta.index:
    #    sncls.load_data(objid, force=False, datafile='%s_data.clf')
    sncls.add_subset('z<.05',color='y',ls='',label='dataset 1',)
    sncls.add_subset('z>.05',color='k',ls='',label='dataset 2',nbin=20)
    sncls.add_subset(color='r',ls='',label='dataset full')    
    sncls.add_parameter(parname='ra')
    sncls.add_parameter(parname='mag', filt1='r')
    sncls.add_parameter(parname='color', phase1=0)
    sncls.add_parameter(parname='mni')
    #tab = sncls.table()
    #print (tab)
    
    #sncls.show2d(1,2)
    sncls.shownd()
    plt.savefig('test.png', dpi=400)
    os.system('open test.png')

    sys.exit()

ra=127.3615417
dec=63.2618889
ztfp = snerun.snelist()
ztfp.parse_meta(syntax='type in ["Ib", "Ic"]')
ztfp.parse_params()
#ztfp.run(clobber=False)

for objid in [
        #'ZTF21aciuhqw',
        #'ZTF22aaezyos',
        #'ZTF21acjxzsa',
        #'ZTF21achgfbk',
        #'ZTF22aaolhhl',
        #'ZTF22aakpyqw',
        #'ZTF22aarbard',
        #'ZTF22aaebkwd',
        #'ZTF22aadqlhy',
        #'ZTF21acdxhuv',
        #'ZTF21acdqgns',
        #'ZTF22aalfmec',
        'ZTF22aadnjsh',
        #'ZTF22aagzmcz',
        #'ZTF21aceehxt',
        #'ZTF22aaoolua',
        #'ZTF21aciwixf',
        #'ZTF22aadoihx'
]:
    ztfp.parse_meta_all(ztfp.kwargs, objid)
    ztfp.data[objid] = snerun.snobject(objid, z=ztfp.z, ra=ztfp.ra, dec=ztfp.dec, mkwebv=ztfp.mkwebv, hostebv=ztfp.hostebv, sntype=ztfp.sntype, dm=ztfp.dm, jdpeak=ztfp.jdpeak)
    #ztfp.data[objid].query_fp_ztf(dend=200)
    ra0, dec0 = ztfp.data[objid].parse_coo(deg=True)[0], ztfp.data[objid].parse_coo(deg=True)[1]
    dist = np.sqrt((ra-ra0)**2+(dec-dec0)**2)
    print (objid, ra0,dec0,dist)
sys.exit()



objid = 'ZTF20aajcdad'
#objid = 'ZTF20abswdbg'

ztfp.parse_meta_all(ztfp.kwargs, objid)
ztfp.data[objid] = snerun.snobject(objid, z=ztfp.z, ra=ztfp.ra, dec=ztfp.dec,
    mkwebv=ztfp.mkwebv, hostebv=ztfp.hostebv, sntype=ztfp.sntype,
    dm=ztfp.dm, jdpeak=ztfp.jdpeak)

#ztfp.data[objid].get_alert_ztf(source='marshal')
ztfp.data[objid].get_fp_ztf()
#ztfp.data[objid].get_local_spectra()
#ztfp.data[objid].get_external_phot('notebook/ZTF20abswdbg_nir.txt', source='myfavour')

if False:
    ztfp.data[objid].run_gp(gp_bands=['r','g','i'])
    ztfp.data[objid].set_peak_gp('r')

if False:
    t0, z = ztfp.data[objid].t0, ztfp.data[objid].z
    ll = open('Shreya_IR_2007gr.txt', 'w')
    for jd in [2459089.98, 2459121.99, 2459151.93]:
        phase = (jd-t0)/(1+z)
        ll.write('%.2f ' % phase)
        for filt in ['r','g','H','K','J']:
            m, em = ztfp.data[objid]._mag_at(filt, phase, interpolation='bin',tdbin=3, corr_mkw=True, corr_host=True)
            if m is None and filt in ['r','g']:
                m, em = ztfp.data[objid]._mag_at(filt, phase, interpolation=None, corr_mkw=True, corr_host=True)
            if m is None: m=99
            if em is None: em=99
            ll.write ('%.2f %.2f ' % (m, em))
        ll.write('\n')
    ll.close()

if True:
    ztfp.data[objid].bin_lc(bin_lc=True, tdbin=1, returnv=False)
    
    #ztfp.data[objid].clip_lc(clipsigma=5)    
    ztfp.data[objid].fig,(ztfp.data[objid].ax, ztfp.data[objid].ax2)=plt.subplots(2,1)
    ztfp.data[objid].fig,ztfp.data[objid].ax=plt.subplots(1,1)
    ztfp.data[objid]._ax(plot_bands=['r', 'g', 'i', 'c', 'o'], show_title=False, show_legend=False)
    #ztfp.data[objid]._ax2(plot_bands=['r', 'g', 'i', 'c', 'o'])
    
if False:
    ztfp.data[objid].run_fit('multiband_main', fit_methods=['bz'],
        multiband_main_bands=['g','r','i'], multiband_main_routine='minimize')

    #ztfp.data[objid].run_fit('multiband_early', fit_methods=['pl'],
    #        multiband_early_bands=['g','r','i'], multiband_early_routine='mcmc',
    #        nsteps_burnin=5000, nsteps=30000)    
    #par = ztfp.data[objid].get_par('tfall', index=0, engine=None, model=None,
    #        source=None, filt='r', quant=[.05,.5,.95])
    #print (par)
    #ztfp.data[objid].fig,ztfp.data[objid].ax2=plt.subplots(1,1)
    #ztfp.data[objid]._ax2(plot_bands=None, xstyle='rp', ax_ystyle='app', show_legend=True, show_fit_error=False)
    
if False:
    ztfp.data[objid].ax = ax
    ztfp.data[objid]._ax(show_title=False, show_gp=True, show_texp=False)
    
if False:
    #ztfp.data[objid].calc_colors(color_bands=['g','r'], tdbin=1, color_interp=['bin','gp','fit'])
    ztfp.data[objid].lyman_bol()        
    
    ztfp.data[objid].set_texp_midway()
    ztfp.data[objid].run_fit('bol_main', fit_methods=['arnett_fit_taum'], fit_redo=True,
                             bol_main_xrangep=[-18, 60], bol_main_routine='minimize')

    # spectra
    #ztfp.data[objid].run_fit('specline', fit_methods=['gauss'])
    #ztfp.data[objid].run_fit('specv_evolution', fit_methods=['exp'])
        
    #par = ztfp.data[objid].get_par('mej', index=0, engine=None, model=None,
    #                            source=None, filt='r', quant=[.05,.5,.95])
    #print (par)    

    ztfp.data[objid].fig,ztfp.data[objid].ax4=plt.subplots(1,1)    
    ztfp.data[objid].ax4.errorbar(-2, 4.82e+42, marker='o', color='r', yerr=4.29e+41, label='BB')
    ztfp.data[objid].ax4.errorbar(27, 2.07e+42, marker='o', color='r', yerr=1.91e+41)
    ztfp.data[objid].ax4.errorbar(56, 1.28e+42, marker='o', color='r', yerr=5.81e+41)
    ztfp.data[objid]._ax4(show_title=True, show_legend=True)
    print ( ztfp.data[objid].get_par('mni'), ztfp.data[objid].get_par('taum') )

if False:

    ztfp.data[objid].bb_colors(bb_bands=['g','r','i'], tdbin=1, color_interp=['bin','gp','fit'],)

    ztfp.data[objid].run_fit('sed', fit_methods=['bb'], sed_routine='minimize', fit_redo=True)
    #ztfp.data[objid].all_fittings()
    ztfp.data[objid].bb_bol()
    #print (ztfp.data[objid].mbolbb)
    
    #ztfp.data[objid].fig,ztfp.data[objid].ax4=plt.subplots(1,1)
    fig, ax = plt.subplots(1,1)
    ztfp.data[objid]._ax5(ax,0,show_title=True, show_legend=True)
    #input('...here')
    
if False:
    ztfp.data[objid].ax1 = ax
    ztfp.data[objid]._ax1(show_text=False)
if False:    
    ztfp.data[objid].ax4 = ax
    ztfp.data[objid]._ax4(show_texp=False)

#for ax in plt.gcf().get_axes():
#    ax.grid(False)
#    ax.axis('off')
plt.savefig('test.png', dpi=400)
os.system('open test.png')
