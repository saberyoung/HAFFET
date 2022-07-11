from sdapy.model_fitters import fit_model, get_pars
from sdapy.functions import *
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def engine(self, model_name, engine_name, sourcename=None, **kwargs):
    ''' engine used to fit bolometric lc around the main peak
    '''
    if 'fitcls' not in self.__dict__: self.fitcls = dict()
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])

    assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
    if sourcename is not None: sources = [sourcename]
    else: sources = kwargs['%s_type'%engine_name].split(',')
    #if len( sources ) == 1:
    #    source = sources[0]
    #    assert source in self.__dict__, 'Error: having %s first before fitting'%source
    #    if source not in self.fitcls[engine_name]:
    #        self.fitcls[engine_name][source] = dict()        
    #    _sengine(self, source, model_name, engine_name, **kwargs)
    #else:        
    for source in sources:
        if not source in self.__dict__:
            print ('Error: having %s first before fitting'%source)
            continue            
        if source not in self.fitcls[engine_name]:
            self.fitcls[engine_name][source] = dict()            
        _sengine(self, source, model_name, engine_name, **kwargs)
    
def _sengine(self, source, model_name, engine_name, **kwargs):
    '''
    df is a dictionary, instead of pandas dataframe
    '''
    df = eval('self.%s'%source)
    
    # prepare data
    xx, yy, yye = [], [], []
    for _ in kwargs['%s_color_interp'%engine_name]:        
        if not _ in df: continue        
        xx = np.append(xx, get_numpy(df[_][0]))
        yy = np.append(yy, get_numpy(df[_][1]))
        yye = np.append(yye, get_numpy(df[_][2]))
    # time to rest frame phase relaitve to peak
    xx = (xx-self.t0)/(1+self.z)
        
    # select data
    if kwargs['%s_xrange'%engine_name] is not None:
        pmin,pmax = min(kwargs['%s_xrange'%engine_name]), max(kwargs['%s_xrange'%engine_name])
        __ = np.logical_and(xx>=pmin, xx<=pmax)       
        xx,yy,yye = xx[__],yy[__],yye[__] 
    if kwargs['%s_yrange'%engine_name] is not None:
        fmin,fmax = min(kwargs['%s_xrange'%engine_name]), max(kwargs['%s_xrange'%engine_name])
        __ = np.logical_and(yy>=fmin, yy<=fmax)
        xx,yy,yye = xx[__],yy[__],yye[__]
        
    # if fit texp,
    fitpars = get_pars(model_name)['par']
    if not 'texp' in fitpars:
        # time to explosion epoch for fits
        assert self.texp is not None, 'set texp before run fits for %s'%model_name
        xx -= self.texp[1]
        t0 = self.t0 + (1+self.z) * self.texp[1]
    else:  t0 = self.t0
    
    # start fit        
    self.fitcls[engine_name][source][model_name] = fit_model(
        xx, yy, yye, filters=None,
    )   
    self.fitcls[engine_name][source][model_name].train(
        opt_routine=kwargs['%s_routine'%engine_name],
        fit_mean=model_name, nwalkers=kwargs['nwalkers'],
        nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'],
        ncores=kwargs['ncores'], thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
        mcmc_h5_file='%s_%s_%s'%\
           (self.objid, model_name, kwargs['%s_routine'%engine_name]),
        emcee_burnin=kwargs['emcee_burnin'], datadir='%s/cache/'%LOCALSOURCE,
        use_emcee_backend=kwargs['use_emcee_backend'], scipysamples=kwargs['scipysamples'],
        clobber=kwargs['fit_redo'], verbose=kwargs['verbose'],sigma=kwargs['fsigma'],
        t0=t0, timedilation=1+self.z, xpredict=kwargs['%s_xrangep'%engine_name],
    )
    self.fitcls[engine_name][source][model_name].predict(quant=kwargs['quantile'])
    self.fitcls[engine_name][source][model_name].get_random_samples()
    self.fitcls[engine_name][source][model_name].set_peak()
    # update t0 and fpeak
    if len(kwargs['set_tpeak_method'])==0 and (self.t0 ==0 or len(self.fpeak)==0):
        self.set_peak_bol_main(model_name=model_name)
    elif kwargs['set_tpeak_method'] == 'bol':
        self.set_peak_bol_main(model_name=model_name)
    if 'texp' in fitpars:  # update texp
        if len(kwargs['set_texp_method'])==0 and self.texp is None:
            self.set_texp_bol_main(model_name=model_name)  
        elif kwargs['set_texp_method'] == 'bolmain':
            self.set_texp_bol_main(model_name=model_name)
