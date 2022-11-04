from sdapy.model_fitters import fit_model
from sdapy.constants import line_location, line_forsn
from sdapy.functions import *
from sdapy.models import *
import numpy as np
import os
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def engine(self, model_name, engine_name,
           sourcename=None, modelname=None, quant=[.05,.5,.95], **kwargs):
    ''' engine used to fit spectral line velocity evolution
    '''
    if 'fitcls' not in self.__dict__: return
    if not 'specline' in self.fitcls.keys(): return
    if engine_name not in self.fitcls: self.fitcls[engine_name] = dict()
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key]) 

    assert self.t0 > 2400000, 'set t0 first'
    if sourcename is None: # sourcename is spectral line name
        sntype = self.sntype
        if sntype in line_forsn: sourcename = line_forsn[sntype]
    if sourcename is None: return 
    intrinstic_wave = line_location[sourcename]
    
    xx, yy, yye = [], [], []
    def func(a): return -a/intrinstic_wave*299792.458/1000.
    _dict = self._get_par('x0', func, None, 'specline', None, None, None, quant)    
    
    for sourcename in _dict:        
        epoch = sourcename.split()[1].split('_')[0]
        vexp = _dict[sourcename]
        t = self.spec.define_phase(epoch)
        xx.append(float(t))
        yy.append(vexp[1])
        yerror = (abs(vexp[0]-vexp[1])+abs(vexp[2]-vexp[1]))/2.
        yye.append( max(.1, yerror) )
    if len(xx) == 0: return
    xx, yy, yye = np.array(xx), np.array(yy), np.array(yye)
    
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
    
    # start fit
    self.fitcls[engine_name][model_name] = fit_model(
        xx, yy, yye, filters=None
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
        t0=0, timedilation=1+self.z, xpredict=kwargs['%s_xrangep'%engine_name],
    )
    self.fitcls[engine_name][model_name].predict(quant=kwargs['quantile'])
    self.fitcls[engine_name][model_name].get_random_samples()
    
    # set vexp
    self.set_vexp()
