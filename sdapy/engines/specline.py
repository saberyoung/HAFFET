from sdapy.model_fitters import fit_model
from sdapy.constants import line_location, line_forsn, \
    line_velocity_p0, line_velocity_bounds, gh_p0, gh_bounds, \
    ga_p0, ga_bounds, gs_p0, gs_bounds
from sdapy.functions import *
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")
            
def engine(self, model_name, engine_name, **kwargs):
    ''' engine used to fit spectral lines
    '''
    if 'fitcls' not in self.__dict__: self.fitcls = dict()
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])

    assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
    source = kwargs['%s_type'%engine_name]
    assert source in self.__dict__, 'Error: having %s first before fitting'%source
    source = eval('self.%s'%source)
    
    if len(kwargs['sn_line']) == 0:        
        if self.sntype in line_forsn:
            _sline = line_forsn[self.sntype]                          
        else:
            print ('!!! No spectral fitting line defined for %s, use Halpha'%self.sntype)
            _sline = 'H~$\alpha$'
    else:
        _sline = kwargs['sn_line']
    assert _sline in line_location, 'define line location for %s, otherwise use %s'%\
        (_sline, line_location.keys())
    
    cw = line_location[_sline]
    # blueside: twice of the allowed maximum velocity
    # redside : line intrinstic wavelength
    blueside, redside = calc_wav(max(line_velocity_bounds)*2,cw), cw    
    region = [blueside, cw + abs(blueside-redside)/10.]
    guessw = [calc_wav(line_velocity_p0, cw),
              [calc_wav(max(line_velocity_bounds),cw),
               calc_wav(min(line_velocity_bounds),cw)]]
    for _ in source.data:
        spec   = source.data[_]['data']
        phase  = float(source.data[_]['phase']) 
        if phase < kwargs['specfit_phase'][0] or phase > kwargs['specfit_phase'][1]:
            continue        
        _source = '_'.join(_.split())        
        if _source not in self.fitcls[engine_name]:
            self.fitcls[engine_name][_source] = dict() 
        _sengine(self, _source, spec, model_name, engine_name, cw, region, guessw, **kwargs)
    self.cw, self.region = cw, region
    
def _sengine(self, source, spec, model_name, engine_name, cw, region, guessw, **kwargs):
    '''
    df is a dictionary, instead of pandas dataframe
    '''
    bestw = guessw[0]
    bestr = guessw[1]
    
    # prepare data, with flux-continuum, which is zoomed and normalized
    wave, flux = spec._norm_spectrum(region=region, stype='flat', **kwargs)
    
    # smooth data locally
    fsmooth = spec._calc_continuum(
        wave, flux, kwargs['continuum_method'], kwargs['continuum_degree'],
    )
    
    # find peaks on smoothed spectra
    findp = spec._calc_peaks(wave, fsmooth, kwargs['pfactor'])    
    
    if 'specpeaks' not in self.__dict__: self.specpeaks = dict()
    self.specpeaks[source] = findp
    
    # check peaks
    if not -1 in findp: # no valley
        print ('no absorption feature detected for %s %s, check spectrum, pfactor, velocity range and intrinstic wavelength'%(self.objid, source))
        self.fitcls[engine_name][source][model_name] = None
        return
    
    # bestv is the closest valley to the guess
    _ = np.argmin(abs(findp[-1]-bestw))
    bestw = findp[-1][_]
    
    if not 1 in findp: # no peaks
        print ('no absorption boundries detected for %s %s, check spectrum, pfactor, velocity range and intrinstic wavelength'%(self.objid, source))
        self.fitcls[engine_name][source][model_name] = None
        return    
    _ = np.where(findp[1]<bestw)        
    bestr_left = findp[1][_]
    if len(bestr_left) == 0:
        print ('no blueend of absorption detected for %s %s, check spectrum, pfactor, velocity range and intrinstic wavelength'%(self.objid, source))
        self.fitcls[engine_name][source][model_name] = None
        return
    bestr[0] = max(bestr_left)
    
    _ = np.where(findp[1]>bestw)        
    bestr_right = findp[1][_]
    if len(bestr_right) > 0: bestr[1] = max(bestr_right)
    
    # cut spectrum again based on peak/valley information
    _ = np.logical_and(wave>bestr[0], wave<bestr[1])
    w, f = wave[_], flux[_]    
    
    # relative to intrinstic wavelength
    w -= cw    
    
    # generate flux errors
    f_unc = spec._add_noise(f, kwargs['spec_snr'])    
    
    if False:
        import matplotlib.pyplot as plt
        import sys
        fig,ax = plt.subplots(1,1)
        ax.plot(wave-cw,flux,'k-', alpha=.6)
        ax.plot(wave-cw,fsmooth,'g-', alpha=.6)
        ax.plot(w,f,'r-', alpha=1)
        #ax.fill_between(w,f-f_unc,f+f_unc,alpha=.2,color='g')        
        #for _ in findp[1]: ax.axvline(_-cw, color='orange', ls='--')
        #for _ in findp[-1]: ax.axvline(_-cw, color='cyan', ls='--')        
        #ax.axvline(bestw-cw, color='y')
        #ax.axvline(bestr[0]-cw, color='b')
        #ax.axvline(bestr[1]-cw, color='b')
        plt.savefig('test.png')
        os.system('open test.png')        
        input(source)
        
    # start fit        
    self.fitcls[engine_name][source][model_name] = fit_model(
        w, f, f_unc, filters=None
    )
    self.fitcls[engine_name][source][model_name].train(
        opt_routine=kwargs['%s_routine'%engine_name],
        fit_mean=model_name, nwalkers=kwargs['nwalkers'],
        nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'],
        ncores=kwargs['ncores'], thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
        mcmc_h5_file='%s_%s_%s_%s'%\
        (self.objid, model_name, kwargs['%s_routine'%engine_name], source),
        emcee_burnin=kwargs['emcee_burnin'], datadir='%s/cache/'%LOCALSOURCE,
        use_emcee_backend=kwargs['use_emcee_backend'], scipysamples=kwargs['scipysamples'],
        clobber=kwargs['fit_redo'], verbose=kwargs['verbose'],sigma=kwargs['fsigma'],
        t0=cw, timedilation=1, xpredict=None,
        bestv=[gh_p0, ga_p0, guessw[0]-cw, gs_p0],
        bounds=((min(gh_bounds), min(ga_bounds), min(guessw[1])-cw, min(gs_bounds)),
                (max(gh_bounds), max(ga_bounds), max(guessw[1])-cw, max(gs_bounds))),
    )
    self.fitcls[engine_name][source][model_name].predict(quant=kwargs['quantile'])
    self.fitcls[engine_name][source][model_name].get_random_samples()
