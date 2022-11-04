from sdapy.model_fitters import fit_model
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")
    
def engine(self, model_name, engine_name, sourcename=None, **kwargs):
    ''' engine used to fit multiband lcs independently before peak
    '''
    if 'fitcls' not in self.__dict__: self.fitcls = dict()
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])            
    source = kwargs['%s_type'%engine_name]
    assert source in self.__dict__, 'Error: having %s first before fitting'%source
    source = eval('self.%s'%source)        
    assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'

    ''' I need peak flux here
    '''
    if len(self.fpeak) == 0:
        print ('Error: set peak flux with GP or main LC fits before run early LC fits')
        return        
    df = source
    if sourcename is not None: df = df.query('source==@sourcename')
    if kwargs['%s_bands'%engine_name] is not None:
        df = df.query('filter in %s' % kwargs['%s_bands'%engine_name])    
    assert kwargs['%s_xrange'%engine_name] is not None
    pmin,pmax = min(kwargs['%s_xrange'%engine_name]), max(kwargs['%s_xrange'%engine_name])
    df = df.query('jdobs<=@self.t0+@pmax and jdobs>=@self.t0+@pmin')
        
    # t: rest frame phase
    t = (df['jdobs'].to_numpy()-self.t0)/(1+self.z)
    # f: original fluxes
    f = df['flux'].to_numpy()
    f_unc = df['eflux'].to_numpy()
    bands = df['filter'].to_numpy()
    
    assert kwargs['%s_yrange'%engine_name] is not None
    fmin,fmax = min(kwargs['%s_yrange'%engine_name]), max(kwargs['%s_yrange'%engine_name])
    
    # normalize flux scale for all bands
    fscale = kwargs['flux_scale']
    f_zp = np.zeros_like(f)
    f_zp_unc = np.zeros_like(f_unc)
    bands = np.array(['']*len(t))
    for filt in np.unique(df['filter']):            
        this_chip = np.where(df['filter'] == filt)
        assert filt in self.fpeak, 'Error: run GP or multiband main fits before to guess the peak for %s, so that it can be afterwards normalized'%filt
        tm = self.fpeak[filt][0]
        f_zp[this_chip] = f[this_chip]/tm*fscale
        f_zp_unc[this_chip] = f_unc[this_chip]/tm*fscale
        bands[this_chip] = filt        
    __ = np.logical_and(f_zp > fmin * fscale, f_zp < fmax * fscale)    
    t = t[__]
    f_zp = f_zp[__]
    f_zp_unc = f_zp_unc[__]
    bands = bands[__]
    
    # start fit
    self.fitcls[engine_name][model_name] = fit_model(
        t, f_zp, f_zp_unc, filters=bands
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
    self.fitcls[engine_name][model_name].get_random_samples(
        limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples']
    )
    self.fitcls[engine_name][model_name].predict(quant=kwargs['quantile'])
    # update texp
    if len(kwargs['set_texp_method'])==0 and self.texp is None:
        self.set_texp_pl(model_name=model_name, **kwargs)
    elif kwargs['set_texp_method'] == 'fit':
        self.set_texp_pl(model_name=model_name, **kwargs)
