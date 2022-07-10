from sdapy.model_fitters import fit_model, get_pars
from sdapy.functions import *
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def engine(self, model_name, engine_name, **kwargs):
    ''' engine used to fit bolometric lc at the tail
    '''
    if 'fitcls' not in self.__dict__: self.fitcls = dict()
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])

    assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'    
    sources = kwargs['%s_type'%engine_name].split(',')
    if len( sources ) == 1:
        source = sources[0]
        assert source in self.__dict__, 'Error: having %s first before fitting'%source
        if source not in self.fitcls[engine_name]:
            self.fitcls[engine_name][source] = dict()        
        _sengine(self, source, model_name, engine_name, **kwargs)
    else:        
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
    
    

    
'''
def engine(self,**kwargs):           
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
    if kwargs['Tail_type'] == 0: return
    elif kwargs['Tail_type'] == 1: clobber = False
    elif kwargs['Tail_type'] == 2: clobber = True
    else: return
        
    if kwargs['Tail_bolopt'] == 1:
        assert 'mbol' in self.__dict__, print ('constrcut bolometric lc first')
        mbol = self.mbol
    elif kwargs['Tail_bolopt'] == 2:
        assert 'mbolbb' in self.__dict__, print ('constrcut bolometric lc first')
        mbol = self.mbolbb
    else:
        if kwargs['verbose']: print ('skip tail fit')
        return
    if not 'texp' in self.__dict__:
        if kwargs['verbose']: print ('estimate explosion epoch first')
        return
    assert self.t0 > 2400000
        
    ncores = kwargs['ncores']
    if ncores is None: ncores = cpu_count() - 1
        
    xx, yy, yye = [], [], []
    for _ in kwargs['Tail_copt']:
        if _ in mbol:
            xx = np.append(xx, get_numpy(mbol[_][0]))
            yy = np.append(yy, get_numpy(mbol[_][1]))
            yye = np.append(yye, get_numpy(mbol[_][2]))
                
    # Tail fits
    xx = (xx-self.t0)/(1+self.z) - self.texp[1]
    p1, p2 = min(kwargs['Tail_fitr']), max(kwargs['Tail_fitr'])
    __ = np.logical_and(xx>=p1, xx<=p2)
    
    if len( xx[__] ) <=3:
        if kwargs['verbose']: print ('too few tail points for fit')
        return
    
    # style=1 fit Mni and t0, no v
    # style=2 fit Mni, Ek and Mej, v can be decided
    # style=3 same fit as style=1, with ts as free parameter as well
    # style=4 same fit as style=2, with ts as free parameter as well
    if kwargs['Tail_style'] == 1:
        fit_mean='tail_fit_t0'        
    elif kwargs['Tail_style'] == 2:
        fit_mean='tail_fit_Mej_Ek'
    elif kwargs['Tail_style'] == 3:
        fit_mean='tail_fit_t0_ts'
    elif kwargs['Tail_style'] == 4:
        fit_mean='tail_fit_Mej_Ek_ts'
    else: return
        
    self.tailcls = fit_model(xx[__],yy[__],yye[__],filters=None)        
    self.tailcls.train(opt_routine=kwargs['Tail_routine'],
                    fit_mean=fit_mean, nwalkers=kwargs['nwalkers'],
                    nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'],ncores=ncores,
                    thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
                    mcmc_h5_file=kwargs['Tail_cache']%
                           (self.objid, kwargs['Tail_routine'], kwargs['Tail_style']),
                    emcee_burnin=kwargs['emcee_burnin'], p0=None, datadir=LOCALCACHE,
                    use_emcee_backend=kwargs['use_emcee_backend'], scipysamples=kwargs['scipysamples'],
                    clobber=clobber, verbose=kwargs['verbose'], sguess=kwargs['sguess'])
    #self.tailcls.predict(quant=kwargs['quantile'])
    self.tail_theta, self.tail_rtheta = self.get_random_samples(self.tailcls.samples,
                self.tailcls.lnprob, cachefile=kwargs['Tail_cache1']%
                            (self.objid, kwargs['Tail_routine'], kwargs['Tail_style']),                    
                    clobber=clobber, verbose=kwargs['verbose'], datadir=LOCALCACHE,)
'''
