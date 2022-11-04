from sdapy.model_fitters import fit_model, get_pars
from sdapy import constants
from sdapy.functions import *
from sdapy.specline_handler import handle_spectrum
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")
            
def engine(self, model_name, engine_name, sourcename=None, **kwargs):
    ''' engine used to fit spectral lines
    '''    
    if 'fitcls' not in self.__dict__: self.fitcls = dict()
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()    
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
    
    assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
    source = kwargs['%s_type'%engine_name]
    assert source in self.__dict__, 'Error: having %s first before fitting'%source
    source = eval('self.%s'%source)
    if sourcename is not None: sourcename = '_'.join(sourcename.split())

    _sline=None
    if kwargs['sn_line']=='full':
        print ('!!! Warning: define a line first')
        return
    elif kwargs['sn_line']=='sntype' and self.sntype in constants.line_forsn:
        _sline = constants.line_forsn[self.sntype]
    elif len(kwargs['sn_line']) == 0 or kwargs['sn_line'] is None:
        if self.sntype in constants.line_forsn:
            _sline = constants.line_forsn[self.sntype]
        else:
            print ('!!! Warning: define a line first')
            return
    else:
        _sline = kwargs['sn_line']
        assert _sline in constants.line_location, 'define line location for %s, otherwise use %s'%\
            (_sline, constants.line_location.keys())
        
    for _ in source.data:
        spec   = source.data[_]['data']
        phase  = float(source.data[_]['phase'])        
        if phase < kwargs['specfit_phase'][0] or phase > kwargs['specfit_phase'][1]:
            #print ('Skipped _ (phase=%.2f) out of phase range: %s'%(phase, kwargs['specfit_phase']))
            continue
        _source = '_'.join(_.split())                
        if sourcename is not None and _source != sourcename: continue
        __source = '%s_%s' % (_source, _sline)
        if _source not in self.fitcls[engine_name]:
            self.fitcls[engine_name][__source] = dict()
        _sengine(self, _source, spec, model_name, engine_name, _sline, **kwargs)
        
def _sengine(self, source, spec, model_name, engine_name, sline, **kwargs):
    '''
    df is a dictionary, instead of pandas dataframe
    '''
    # find peaks on smoothed spectra
    cw, region = handle_spectrum.parse_element(sline, kwargs['v_bounds'])
    
    # prepare data, with flux-continuum, which is zoomed and normalized
    wave, flux = spec._norm_spectrum(region=region, stype='flat')

    # fit range
    if kwargs['force_range'] is None:
        findp = spec.find_peaks(region=region, stype='flat', smooth=True,
                pfactor=kwargs['pfactor'], continuum_method=kwargs['continuum_method'],
                continuum_degree=kwargs['continuum_degree'], line_velocity_bounds=kwargs['v_bounds'])
        if findp is None: return
    
        # check peak informations for fitting range
        if not -1 in findp: # no valley
            print ('no absorption feature detected for %s %s, check spectrum, pfactor, velocity range and intrinstic wavelength'%(self.objid, source))
            self.fitcls[engine_name]['%s_%s' % (source, sline)][model_name] = None
            return

        # mean (the valley) cannot be too close to the intrinstci wavelength
        _ = np.where( abs(findp[-1] - cw) > min(kwargs['v_bounds'])*1e8/constants.c*cw )
        if len(findp[-1][_]) == 0:
            print ('no absorption feature at reasonable location for %s %s, check spectrum, pfactor, velocity range and intrinstic wavelength'%(self.objid, source))
            self.fitcls[engine_name]['%s_%s' % (source, sline)][model_name] = None
            return
    
        # bestv is the closest valley to the guess
        guess_v_loc = calc_wav(kwargs['v_p'], cw)
        __ = np.argmin( abs(findp[-1][_] - guess_v_loc) )
        bestw = findp[-1][_][__]
    
        # bestr_left/right defines the fitting range
        if not 1 in findp: # no peaks
            print ('no red endabsorption boundries detected for %s %s, check spectrum, pfactor, velocity range and intrinstic wavelength'%(self.objid, source))
            self.fitcls[engine_name]['%s_%s' % (source, sline)][model_name] = None
            return    
        _ = np.where(findp[1]<bestw)        
        bestr_left = findp[1][_]
        if len(bestr_left) == 0:
            print ('no blue end of absorption detected for %s %s, check spectrum, pfactor, velocity range and intrinstic wavelength'%(self.objid, source))
            self.fitcls[engine_name]['%s_%s' % (source, sline)][model_name] = None
            return    
        _ = np.where(findp[1]>bestw)        
        bestr_right = findp[1][_]
        if len(bestr_right) > 0:
            # cut spectrum again based on peak/valley information    
            _ = np.logical_and(wave>bestr_left[-1], wave<bestr_right[0]) 
            w, f = wave[_], flux[_]
        else:
            _ = np.where(wave>bestr_left[-1])
            w, f = wave[_], flux[_]        

    else:
        assert type(kwargs['force_range']) is list        
        # cut spectrum based on provided range
        _ = np.logical_and(wave>min(kwargs['force_range']), wave<max(kwargs['force_range']))
        w, f = wave[_], flux[_]
    
    # relative to intrinstic wavelength
    w -= cw    
    
    # generate flux errors
    f_unc = spec._add_noise(f, kwargs['spec_snr'])    
    
    # start fit        
    self.fitcls[engine_name]['%s_%s' % (source, sline)][model_name] = fit_model(
        w, f, f_unc, filters=None
    )
    self.fitcls[engine_name]['%s_%s' % (source, sline)][model_name].train(
        opt_routine=kwargs['%s_routine'%engine_name],
        fit_mean=model_name, nwalkers=kwargs['nwalkers'],
        nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'],
        ncores=kwargs['ncores'], thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
        mcmc_h5_file='%s_%s_%s_%s_%s'%\
        (self.objid, model_name, kwargs['%s_routine'%engine_name], source, sline),
        emcee_burnin=kwargs['emcee_burnin'], datadir='%s/cache/'%LOCALSOURCE,
        use_emcee_backend=kwargs['use_emcee_backend'], scipysamples=kwargs['scipysamples'],
        clobber=kwargs['fit_redo'], verbose=kwargs['verbose'],sigma=kwargs['fsigma'],
        t0=cw, timedilation=1, xpredict=None,
    )
    self.fitcls[engine_name]['%s_%s' % (source, sline)][model_name].predict(quant=kwargs['quantile'])
    self.fitcls[engine_name]['%s_%s' % (source, sline)][model_name].get_random_samples()
