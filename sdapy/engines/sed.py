from sdapy.model_fitters import fit_model, get_pars
from sdapy.functions import *
from sdapy.filters import *
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def engine(self, model_name, engine_name, sourcename=None, **kwargs):
    ''' engine used to fit SED on multiband photometry or spectra for bolometric LCs
    '''    
    if 'fitcls' not in self.__dict__: self.fitcls = dict()
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])    
    assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
    
    if sourcename is not None: sources = [sourcename]
    else: sources = kwargs['%s_type'%engine_name].split(',')
    
    for source in sources:
        source = source.strip()        
        
        if source == 'bb' and 'bb' in kwargs['make_bol']:  # photometry sed
            if 'cbb' in self.__dict__:
                cbb = self.cbb
            else:
                cbb = self.bb_colors(xpred=None, index=0, returnv=True, **kwargs)                        
            filters = kwargs['sed_bands']
            if filters is None: filters = np.unique(self.lc['filter'])            
            if kwargs['sed_color_interp_mix']:
                _fluxes, _efluxes, _ws, _bs = dict(), dict(), dict(), dict()
            for interpolation in kwargs['sed_color_interp']: 
                jds = cbb[interpolation]['jd']
                for _ in range(len(jds)):
                    jd = jds[_]
                    phase = (jd - self.t0) / (1+self.z)
                    if kwargs['%s_range' % interpolation] is not None:
                        if phase<min(kwargs['%s_range' % interpolation]) or phase>max(kwargs['%s_range' % interpolation]):
                            continue
                    if not kwargs['sed_color_interp_mix']:
                        fluxes, efluxes, ws, bs = [], [], [], []
                    for n,f in enumerate(filters):
                        if f not in cbb[interpolation]: continue
                        m, em = cbb[interpolation][f]
                        m, em = m[_], em[_]                                        
                        if m is None or em is None: continue
                        _flux, _dflux, wlref, bandwidths = self._to_luminosity(m, f, em=em, **kwargs)
                        if not kwargs['sed_color_interp_mix']:
                            fluxes.append(_flux)
                            efluxes.append(_dflux)
                            ws.append(wlref)
                            bs.append(bandwidths)
                        else:
                            if not jd in _fluxes:
                                _fluxes[jd], _efluxes[jd], _ws[jd], _bs[jd] = [], [], [], []
                            _fluxes[jd].append(_flux)
                            _efluxes[jd].append(_dflux)
                            _ws[jd].append(wlref)
                            _bs[jd].append(bandwidths)
                    if not kwargs['sed_color_interp_mix']:
                        fluxes, efluxes, ws, bs = np.array(fluxes), np.array(efluxes), np.array(ws), np.array(bs)
                        
                        # if any band cannot provide fluxes, just skip them                
                        __ = np.isnan(fluxes)
                        if len(fluxes[__]) > 0 or len(fluxes) < 2: continue
                    
                        # Set flux to zero at red and blue extrema matching wlref1
                        #flux1 = np.insert(fluxes,0,0)
                        #flux1 = np.append(flux1,0)
                        
                        _source = '%s_%s_%.1f' % (source, interpolation, jd)
                        if _source not in self.fitcls[engine_name]:
                            self.fitcls[engine_name][_source] = dict()                                            
                            _sengine(self, _source, ws, fluxes, efluxes, model_name, engine_name, **kwargs)
                            
            if kwargs['sed_color_interp_mix']:
                for jd in _fluxes:
                    fluxes, efluxes, ws, bs = np.array(_fluxes[jd]), np.array(_efluxes[jd]), np.array(_ws[jd]), np.array(_bs[jd])
                    
                    # if any band cannot provide fluxes, just skip them  
                    __ = np.isnan(fluxes)
                    if len(fluxes[__]) > 0 or len(fluxes) < 2: continue
                    
                    _source = '%s_mix_%.1f' % (source, jd)
                    if _source not in self.fitcls[engine_name]:
                        self.fitcls[engine_name][_source] = dict()                                            
                        _sengine(self, _source, ws, fluxes, efluxes, model_name, engine_name, **kwargs)
            
        if source == 'spec' and 'spec' in self.__dict__ and 'spec' in kwargs['make_bol']:
            for _ in self.spec.data: 
                spec   = self.spec.data[_]['data']
                phase  = float(self.spec.data[_]['phase'])
                jd = phase * (1+self.z) + self.t0                
                if kwargs['specfit_phase'] is not None:
                    if phase<min(kwargs['specfit_phase']) or phase>max(kwargs['specfit_phase']):
                        continue                    
                # get original data
                ws, fluxes = spec.wave, spec.flux
                efluxes = spec._add_noise(fluxes, kwargs['spec_snr'])
                if len(ws) == 0: continue
                
                fluxes, efluxes = self._cal_spectrum(ws, fluxes, efluxes, phase, **kwargs)                    
                _source = '%s_%.1f' % (source, jd)
                if _source not in self.fitcls[engine_name]:
                    self.fitcls[engine_name][_source] = dict()                                            
                    _sengine(self, _source, ws, fluxes, efluxes, model_name, engine_name, **kwargs)                        
    
def _sengine(self, source, ws, fluxes, efluxes, model_name, engine_name, **kwargs):
    '''
    df is a dictionary, instead of pandas dataframe
    '''       
    # prepare data
    xx, yy, yye = ws, fluxes, efluxes
    if xx is None or yy is None or yye is None: return
    
    # select data
    if kwargs['%s_xrange'%engine_name] is not None:
        pmin,pmax = min(kwargs['%s_xrange'%engine_name]), max(kwargs['%s_xrange'%engine_name])
        __ = np.logical_and(xx>=pmin, xx<=pmax)
        xx,yy,yye = xx[__],yy[__],yye[__]

    # clear data
    xx  = xx[~np.isnan(yy)]    
    yye = yye[~np.isnan(yy)]
    yy  = yy[~np.isnan(yy)]
    
    xx  = xx[~np.isinf(yy)]    
    yye = yye[~np.isinf(yy)]
    yy  = yy[~np.isinf(yy)]

    if len(xx) < 2: return
    # start fit        
    self.fitcls[engine_name][source][model_name] = fit_model(
        xx, yy, yye, filters=None,
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
        t0=0, timedilation=1, xpredict=kwargs['%s_xrangep'%engine_name],
    )
    self.fitcls[engine_name][source][model_name].predict(quant=kwargs['quantile'])
    self.fitcls[engine_name][source][model_name].get_random_samples()
