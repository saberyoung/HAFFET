from sdapy.model_fitters import fit_model, get_pars
from sdapy.functions import *
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def engine(self, model_name, engine_name, sourcename=None, **kwargs):
    ''' engine used to fit bolometric lc at the tail
    '''
    if 'fitcls' not in self.__dict__: self.fitcls = dict()
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])

    assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'    
    if sourcename is not None:
        sources = [sourcename]
    else:
        sources = kwargs['%s_type'%engine_name].split(',')
    for source in sources:
        source = source.strip()
        if not source in self.__dict__: continue            
        if source not in self.fitcls[engine_name]:
            self.fitcls[engine_name][source] = dict()            
        _sengine(self, source, model_name, engine_name, **kwargs)
    
def _sengine(self, source, model_name, engine_name, **kwargs):
    '''
    df is a dictionary, instead of pandas dataframe
    '''
    df = eval('self.%s'%source)
    print (source)
    
    # prepare data
    xx, yy, yye = [], [], []
    for _ in kwargs['%s_color_interp'%engine_name]:
        if not _ in df: continue
        xx = np.append(xx, get_numpy(df[_][0]))
        yy = np.append(yy, get_numpy(df[_][1]))
        yye = np.append(yye, get_numpy(df[_][2]))

    if len(xx) == 0: return
    # time to rest frame phase relaitve to peak
    xx = (xx-self.t0)/(1+self.z)
        
    # select data
    if kwargs['%s_xrange'%engine_name] is not None:
        pmin,pmax = min(kwargs['%s_xrange'%engine_name]), max(kwargs['%s_xrange'%engine_name])
        __ = np.logical_and(xx>=pmin, xx<=pmax)

        if len(xx[__]) == 0: return
        xx,yy,yye = xx[__],yy[__],yye[__]
    if kwargs['%s_yrange'%engine_name] is not None:
        fmin,fmax = min(kwargs['%s_xrange'%engine_name]), max(kwargs['%s_xrange'%engine_name])
        __ = np.logical_and(yy>=fmin, yy<=fmax)

        if len(xx[__]) == 0: return
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
    self.fitcls[engine_name][source][model_name].get_random_samples(
        limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples']
    )
    self.fitcls[engine_name][source][model_name].set_peak()    
