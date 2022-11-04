from sdapy.model_fitters import fit_model
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def engine(self, model_name, engine_name, sourcename=None, **kwargs):
    ''' engine used to fit multiband lcs independently around the main peak
    '''
    if 'fitcls' not in self.__dict__: self.fitcls = dict()
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])            
    source = kwargs['%s_type'%engine_name]
    assert source in self.__dict__, 'Error: having %s first before fitting'%source
    source = eval('self.%s'%source)        
    assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'    
    
    df = source
    if sourcename is not None: df = df.query('source==@sourcename')
    if kwargs['%s_bands'%engine_name] is not None:
        df = df.query('filter in %s' % kwargs['%s_bands'%engine_name])    
    if kwargs['%s_xrange'%engine_name] is not None:
        pmin,pmax = min(kwargs['%s_xrange'%engine_name]), max(kwargs['%s_xrange'%engine_name])
        df = df.query('jdobs<=@self.t0+@pmax and jdobs>=@self.t0+@pmin')
    if kwargs['%s_yrange'%engine_name] is not None:
        fmin,fmax = min(kwargs['%s_xrange'%engine_name]), max(kwargs['%s_xrange'%engine_name])
        df = df.query('flux<=@fmax and jdobs>=@fmin')
    # t: rest frame phase
    t = (df['jdobs'].to_numpy()-self.t0)/(1+self.z)
    # f: original fluxes
    f = df['flux'].to_numpy()
    f_unc = df['eflux'].to_numpy()
    bands = df['filter'].to_numpy()
    
    # start fit        
    self.fitcls[engine_name][model_name] = fit_model(
        t, f, f_unc, filters=bands
    )
    self.fitcls[engine_name][model_name].train(
        opt_routine=kwargs['%s_routine'%engine_name],
        fit_mean=model_name, nwalkers=kwargs['nwalkers'],
        nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'],
        ncores=kwargs['ncores'], thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
        mcmc_h5_file='%s_%s_%s'%\
           (self.objid, model_name, kwargs['%s_routine'%engine_name]),
        emcee_burnin=kwargs['emcee_burnin'], datadir='%s/cache/'%LOCALSOURCE,
        use_emcee_backend=kwargs['use_emcee_backend'], scipysamples=kwargs['scipysamples'],
        clobber=kwargs['fit_redo'], verbose=kwargs['verbose'],sigma=kwargs['fsigma'],
        t0=self.t0, timedilation=1+self.z, xpredict=kwargs['%s_xrangep'%engine_name],
    )
    self.fitcls[engine_name][model_name].predict(quant=kwargs['quantile'])
    self.fitcls[engine_name][model_name].get_random_samples(
        limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples']
    )
    self.fitcls[engine_name][model_name].set_peak()    
    # update t0 and fpeak
    if len(kwargs['set_tpeak_method'])==0 and (self.t0 ==0 or len(self.fpeak)==0):
        self.set_peak_multiband_main(model_name=model_name, **kwargs)
    elif kwargs['set_tpeak_method'] == 'fit':
        self.set_peak_multiband_main(model_name=model_name, **kwargs) 
