#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : gaussian_process.py
# License           : BSD-3-Clause
# Author            : syang <sheng.yang@astro.su.se>
# Date              : 01.01.2020
# Last Modified Date: 15.05.2020
# Last Modified By  : syang <sheng.yang@astro.su.se>

import numpy as np
import george, os, emcee
from scipy.optimize import least_squares, minimize
import matplotlib.pyplot as plt
from multiprocessing import Pool
from sdapy.filters import central_wavelengths
from sdapy.functions import *
from sdapy.corner_hack import corner_hack

class fit_gp:
    """Fits data with gaussian process.

    The package 'george' is used for the gaussian process fit.

    Parameters
    ----------
    x_data : array
        Independent values.
    y_data : array
        Dependent values.
    yerr_data : array, int
        Dependent value errors.
    filters : array
        If filters available, will convolve wavelengths to x_data to train gaussian process.
    kernel : str, default 'matern52'
        Kernel to be used with the gaussian process. 
        Possible choices are: 'matern52', 'matern32', 'squaredexp'.
    fix_scale : bool
        If fix default gaussian process param
    gp_mean: str, default 'mean'
        Mean y_data function.
        Possible choices are: 'mean', 'gaussian', 'bazin', 'villar'.
    opt_routine : str, 
        Which technic to be used to realize optimization.
        Possible choices are: 'minimize', 'mcmc', 'leastsq'.
    nwalkers : int
        if mcmc adopted, set walker number
    nsteps:  int
        if mcmc adopted, set step
    nsteps_burnin: int
        if mcmc adopted, set burnin step
    clobber: bool
        if gp already done, if redo it or not

    Returns
    -------
    Returns the interpolated independent and dependent values with the 1-sigma standard deviation.

    Examples
    -------
    gp = fit_gp(jd, flux, fluxerr, central_wavelength)
    gp.train(gp_mean='bazin', opt_routine = 'mcmc')
    gp_jd, gp_flux, gp_flux_errors, gp_ws = gp.predict()
    gp.save_corner(saveplotas='tmp.png')

    """
    def __init__(self,  x_data, y_data, yerr_data=1e-8, filters=None):
        try: # if pandas.dataframe
            self.x = x_data.to_numpy()
            self.y = y_data.to_numpy()
            self.yerr = yerr_data.to_numpy()
        except: # otherwise, numpy array
            self.x = x_data
            self.y = y_data
            self.yerr = yerr_data    
        self.filters = filters
        self.parse_filters()        
        self.y_norm = self.y.max()
        self.gp = None
        self.meanp = None

    def parse_filters(self):
        if self.filters is None:
            self.wavelengths = None
        else:
            self.wavelengths = []
            for _band in self.filters:
                self.wavelengths.append( central_wavelengths[_band] )
    
    def train(self,kernel='matern52', fix_scale=False, gp_mean='mean', 
              opt_routine = 'minimize', nwalkers=30, nsteps=1000,
              nsteps_burnin=50, mcmc_h5_file='tmp', clobber=False,
              verbose=False, datadir='./', thin_by=1,
              t0=0, timedilation=1, xpredict=None,
              emcee_burnin=True, use_emcee_backend=True):
        
        assert opt_routine in ['mcmc', 'minimize', 'leastsq']
        assert gp_mean in ['mean', 'gaussian', 'bazin', 'villar']
        assert kernel in ['matern52', 'matern32', 'squaredexp']
        self.datadir = datadir
        self.t0 = t0
        self.timedilation = timedilation
        self.xpredict = xpredict
        
        _s = ''
        for f in central_wavelengths:            
            if f in self.filters: _s += f
        mcmc_h5_file = '%s/%s_%s'%(datadir, mcmc_h5_file, _s)
        if not '.h5' in mcmc_h5_file:  mcmc_h5_file = '%s.h5'%mcmc_h5_file
        
        # lc mean model
        if gp_mean=='gaussian':
            class lcMeanModel(george.modeling.Model):
                """Gaussian light curve model."""
                parameter_names = ("A", "mu", "log_sigma2")

                def get_value(self, t):
                    return self.A * np.exp(-0.5*(t-self.mu)**2 * np.exp(-self.log_sigma2))

                # This method is to compute the gradient of the objective function below.
                def compute_gradient(self, t):
                    e = 0.5*(t-self.mu)**2 * np.exp(-self.log_sigma2)
                    dA = np.exp(-e)
                    dmu = self.A * dA * (t-self.mu) * np.exp(-self.log_sigma2)
                    dlog_s2 = self.A * dA * e
                    return np.array([dA, dmu, dlog_s2])
            
        elif gp_mean=='bazin':
            class lcMeanModel(george.modeling.Model):
                """Bazin et al. (2011) light curve model."""
                parameter_names = ("A", "t0", "tf", "tr")

                def get_value(self, t):
                    np.seterr(all='ignore')
                    bazin_model = self.A * np.exp(-(t-self.t0)/self.tf) / (1 + np.exp(-(t-self.t0)/self.tr))
                    np.seterr(all='warn')

                    return bazin_model

                def compute_gradient(self, t):
                    np.seterr(all='ignore')
                    Tf = (t-self.t0)/self.tf
                    Tr = (t-self.t0)/self.tr
                    B = np.exp(-Tf)
                    C = 1 + np.exp(-Tr)
                    
                    dA = B/C
                    dtf = A*B/C * Tf/self.tf
                    dtr = A*B/C**2 * Tf/self.tr * np.exp(-Tr)
                    dt0 = (np.exp(-t/self.tf + self.t0/self.tf + t/self.tr) *
                           (self.tr*(np.exp(t/self.tr) + np.exp(self.t0/self.tr)) - self.tr*np.exp(self.t0/self.tr)) /
                           (self.tr*self.tf*(np.exp(t/self.tr) + np.exp(self.t0/self.tr))**2)
                    )
                    np.seterr(all='warn')

                    return np.array([dA, dt0, dtf, dtr])

        elif gp_mean=='zheng':
            class lcMeanModel(george.modeling.Model):
                """Zheng et al. (2018) light curve model."""
                parameter_names = ("A", "t0", "tb", "ar", "ad", "s")
            
                def get_value(self, t):
                    np.seterr(all='ignore')
                    Tb = (t-self.t0)/self.tb
                    zheng_model = self.A * Tb**self.ar * (1 + Tb**(self.s*self.ad))**(-2/self.s)
                    np.seterr(all='warn')

                    return zheng_model

                def compute_gradient(self, t):
                    np.seterr(all='ignore')
                    Tb = (t-self.t0)/self.tb
                    Tb_ar = Tb**self.ar
                    Tb_sad = Tb**(self.s*self.ad)
                    
                    dA = Tb_ar * (1 + Tb_sad)**(-2/self.s)
                    dt0 = ((2*self.A*self.ad * (Tb_sad + 1)**(-2/self.s - 1) * Tb_sad*Tb_ar/Tb)/self.tb -
                           (self.A*self.ar*Tb_ar/Tb * (Tb_sad + 1)**(-2/self.s))/self.tb
                    )
                    dtb = dt0*Tb
                    dar = self.A * Tb_ar * np.log(Tb) * (Tb_sad + 1)**(-2/self.s)
                    dad = -2*self.A * np.log(Tb) * (Tb_sad + 1)**(-2/self.s - 1) * Tb_sad*Tb_ar
                    ds = self.A * Tb_ar * (Tb_sad + 1)**(-2/self.s) * (2*np.log(Tb_sad + 1)/self.s**2
                                                        - 2*self.ad*np.log(Tb)*Tb_sad/(self.s*(Tb_sad + 1)))
                    np.seterr(all='warn')

                    return np.array([dA, dt0, dtb, dar, dad, ds])
            
        elif gp_mean=='villar':
            class lcMeanModel(george.modeling.Model):
                """Villar et al. (2019) light curve model."""
                parameter_names = ("A", "t0", "t1", "tf", "tr")
                
                def get_value(self, t):
                    np.seterr(all='ignore')
                    for _t in t:
                        if _t < self.t1:                        
                            villar_model = self.A*(_t-self.t0) / (1 + np.exp(-(_t - self.t0) / self.tr))
                        else:
                            villar_model = self.A*(self.t1-self.t0) * np.exp(-(_t - self.t1) / self.tf) / (1 + np.exp(-(_t - self.t0) / self.tr))
                    np.seterr(all='warn')

                    return villar_model

                def compute_gradient(self, t):
                    np.seterr(all='ignore')
                    Tf = (t-self.t0)/self.tf
                    Tr = (t-self.t0)/self.tr
                    B = np.exp(-Tf)
                    C = 1 + np.exp(-Tr)
                    
                    dA = B/C
                    dtf = A*B/C * Tf/self.tf
                    dtr = A*B/C**2 * Tf/self.tr * np.exp(-Tr)
                    dt0 = (np.exp(-t/self.tf + self.t0/self.tf + t/self.tr) *
                           (self.tr*(np.exp(t/self.tr) + np.exp(self.t0/self.tr)) - self.tr*np.exp(self.t0/self.tr)) /
                           (self.tr*self.tf*(np.exp(t/self.tr) + np.exp(self.t0/self.tr))**2)
                    )
                    dt1 = dt0
                    np.seterr(all='warn')
                
                    return np.array([dA, dt0, dt1, dtf, dtr])
            
        # define the objective function (negative log-likelihood in this case)
        def neg_log_like(params):
            """Negative log-likelihood."""
            gp.set_parameter_vector(params)
            log_like = gp.log_likelihood(y, quiet=True)        
            if np.isfinite(log_like):  return -log_like
            else:  return np.inf

        # and the gradient of the objective function
        def grad_neg_log_like(params):
            """Gradient of the negative log-likelihood."""
            gp.set_parameter_vector(params)
            return -gp.grad_log_likelihood(y, quiet=True)

        # for mcmc
        def lnprob(p):
            gp.set_parameter_vector(p)
            return gp.log_likelihood(y, quiet=True) + gp.log_prior()
        
        x, ndim = np.copy(self.x), 1
        y, yerr = np.copy(self.y), np.copy(self.yerr)
        
        if not self.wavelengths is None:
            assert len(x) == len(self.wavelengths)
            x = np.vstack([x, self.wavelengths]).T
            ndim = 2    

        # normalise the data for better results
        y /= self.y_norm
        yerr /= self.y_norm
        
        # choose GP mean model with the respective bounds
        if gp_mean=='mean':
            mean_model = y.mean()

        elif gp_mean=='gaussian':
            A, mu, log_sigma2 = y.max(), x[y==y.max()][0], np.log(10)
            mean_bounds = {'A':(1e-3, 1e2),
                           'mu':(mu-50, mu+50),
                           'log_sigma2':(np.log(10), np.log(60)),
            }
            mean_model = lcMeanModel(A=A, mu=mu, log_sigma2=log_sigma2, bounds=mean_bounds)

        elif gp_mean=='bazin':
            A, tf, tr = y.max(), 40, 20
            t0 = 20 + tr*np.log(tf/tr-1)
            mean_bounds = {'A':(1e-3, 1e2),
                           't0':(t0-50, t0+50),
                           'tf':(tf-10, tf+10),
                           'tr':(tr-10, tr+10),
            }
            mean_model = lcMeanModel(A=A, t0=t0, tf=tf, tr=tr, bounds=mean_bounds)

        elif gp_mean=='zheng':
            A, t0, tb, ar, ad, s = y.max(), x[y==y.max()][0]-20, 20, 2, 2.5, 1.5
            mean_bounds = {'A':(1e-3, 1e2),
                           't0':(t0-50, t0+50),
                           'tb':(tb-10, tb+10),
                           'ar':(ar-1.0, ar+1.0),
                           'ad':(ad-1.5, ad+1.5),
                           's':(s-1.0, s+1.0),
            }
            mean_model = lcMeanModel(A=A, t0=t0, tb=tb, ar=ar, ad=ad, s=s, bounds=mean_bounds)

        elif gp_mean=='villar':
            A, tf, tr = y.max(), 40, 20
            t0 = 20 + tr*np.log(tf/tr-1)
            mean_bounds = {'A':(1e-3, 1e2),
                           't0':(t0-50, t0+50),
                           't1':(t0-50, t0+50),
                           'tf':(tf-10, tf+10),
                           'tr':(tr-10, tr+10),
            }
            mean_model = lcMeanModel(A=A, t0=t0, t1=t0, tf=tf, tr=tr, bounds=mean_bounds)
        
        var, length_scale = np.var(y), np.diff(self.x).max()
        bounds_var, bounds_length = [(np.log(1e-6), np.log(10))], [(np.log(1), np.log(1e6))]

        if not self.wavelengths is None:
            prio = [length_scale ** 2, 6000 ** 2]
            bounds_length += [(0, np.log(10000 ** 2))]
        else:
            prio = length_scale**2
        
        # a constant kernel is used to allow adding bounds
        k1 = george.kernels.ConstantKernel(np.log(var), bounds=bounds_var, ndim=ndim)
    
        if kernel == 'matern52':
            k2 = george.kernels.Matern52Kernel(prio, metric_bounds=bounds_length, ndim=ndim)
        elif kernel == 'matern32':
            k2 = george.kernels.Matern32Kernel(prio, metric_bounds=bounds_length, ndim=ndim)
        elif kernel == 'squaredexp':
            k2 = george.kernels.ExpSquaredKernel(prio, metric_bounds=bounds_length, ndim=ndim)
        else:
            raise ValueError(f'"{kernel}" is not a valid kernel.')

        ker = k1*k2
        if fix_scale: ker.freeze_parameter("k2:metric:log_M_0_0") 
        #ker.freeze_parameter("k1:log_constant")

        if not self.wavelengths is None:  gp = george.GP(kernel=ker)
        else: gp = george.GP(kernel=ker, mean=mean_model, fit_mean=True)
        
        # initial guess
        gp.compute(x, yerr)
        
        # check if cache exists or not        
        if not os.path.exists(mcmc_h5_file):
            pass
        elif clobber:
            print ('remove %s' % mcmc_h5_file)
            os.remove( mcmc_h5_file )
        else:
            if opt_routine=='mcmc':
                self.samples, self.lnprob = get_samples_mc(mcmc_h5_file)               
                p_final = np.mean(self.samples, axis=0)
                gp.set_parameter_vector(p_final)
                self.meanp = gp.get_parameter_dict()        
                self.gp = gp
                return
        
        # optimization routine for hyperparameters
        self.samples=None
        if opt_routine == 'minimize':                
            p0 = gp.get_parameter_vector()
            bounds = gp.get_parameter_bounds()       
            results = minimize(neg_log_like, p0, jac=grad_neg_log_like,
                               method="L-BFGS-B", options={'maxiter':30},
                               bounds=bounds)        
            gp.set_parameter_vector(results.x)
            self.meanp = gp.get_parameter_dict()
            self.gp = gp
            
        elif opt_routine == 'leastsq':                
            p0 = gp.get_parameter_vector()        
            results = least_squares(neg_log_like, p0, jac=grad_neg_log_like,
                                    method="trf")   
            gp.set_parameter_vector(results.x)
            self.meanp = gp.get_parameter_dict()
            self.gp = gp
            
        elif opt_routine == 'mcmc':                    
            initial = gp.get_parameter_vector()
            ndim = len(initial)
            p0 = initial + 1e-8 * np.random.randn(nwalkers, ndim)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
            
            if emcee_burnin:                    
                # Running burn-in...                    
                p0, _, _ = sampler.run_mcmc(p0, nsteps_burnin,
                                            thin_by=thin_by, progress=verbose)
                
            if use_emcee_backend:            
                # file to save samples
                filename = mcmc_h5_file
                backend = emcee.backends.HDFBackend(filename)
                backend.reset(nwalkers, ndim)
                    
                sampler = emcee.EnsembleSampler(nwalkers, ndim, 
                                                lnprob, backend=backend)

            else:
                sampler = emcee.EnsembleSampler(nwalkers, ndim, lbprob)
            
            sampler.reset()
            # Running production...
            sampler.run_mcmc(p0, nsteps, thin_by=thin_by, progress=verbose)
            
            #self.samples = sampler.flatchain                
            self.samples, self.lnprob = get_samples_mc(mcmc_h5_file)
            
            p_final = np.mean(self.samples, axis=0)
            gp.set_parameter_vector(p_final)
            self.meanp = gp.get_parameter_dict()        
            self.gp = gp
        
    def save_corner(self, figpath, datadir=None, quantiles=[.05,.95], clobber=False):
        ''' generate corner plots '''
        if self.samples is None:
            print ('contours for mcmc mode')
            return
        if datadir is None: datadir = self.datadir        
        _fign = '{}/{}.png'.format(datadir, figpath)
        if os.path.exists(_fign) and not clobber: pass
        else:
            if os.path.exists(_fign): os.remove(_fign)
            cfig = corner_hack(self.samples, labels=list(self.meanp.keys()),
                               label_kwargs={'fontsize':16}, ticklabelsize=13,
                               show_titles=True, quantiles=quantiles,
                               title_fmt=".2f", title_kwargs={'fontsize':16},
                               plot_datapoints=True, plot_contours=True)                          
            cfig.savefig(_fign, bbox_inches='tight')        
        return [_fign]
    
    def predict(self, x_pred=None, step = 1, clobber=False, returnv=False):
        ''' output GP products '''
        if 'x_pred' in self.__dict__ and not clobber and not returnv: return
        if x_pred is None:
            if self.xpredict is not None:
                x_pred = np.arange(min(self.xpredict), max(self.xpredict), step) * self.timedilation+self.t0
            else:
                x_pred = np.arange(self.x.min(), self.x.max(), step)                
        x_predl, y_predl, y_pred_el, _ws = [], [], [], []                
        if not self.wavelengths is None:            
            for _wavelength in np.unique(self.wavelengths):
                _wavelengths = np.ones(len(x_pred)) * [_wavelength]            
                _x_pred = np.vstack([x_pred, _wavelengths]).T
                mu, var = self.gp.predict(self.y/self.y_norm, _x_pred, return_var=True)
                std = np.sqrt(var)
                x_predl.append( x_pred )
                y_predl.append( mu * self.y_norm )
                y_pred_el.append( std * self.y_norm )                
                _ws.append( np.array( len(x_pred) * [self.parse_wavelength( _wavelength )]) )
        else:
            mu, var = self.gp.predict(self.y/self.y_norm, x_pred, return_var=True)
            std = np.sqrt(var)
            x_predl, y_predl, y_pred_el = x_pred, mu*self.y_norm, std*self.y_norm          

        if returnv:            
            return np.array(x_predl), np.array(y_predl), np.array(y_pred_el), np.array(_ws)        
        self.x_pred = np.array(x_predl)
        self.y_pred = np.array(y_predl)
        self.y_prede = np.array(y_pred_el)
        self.f_pred = np.array(_ws)     

    @staticmethod
    def parse_wavelength(w):        
        for f in central_wavelengths:            
            _w = central_wavelengths[f]            
            if _w == w: return  f
        return

    def set_peak(self, clobber=False):        
        if 'tpeak' in self.__dict__ and not clobber: return
        self.tpeak = dict() 
        self.fpeak = dict()        
        for filt in np.unique(self.f_pred):
            __ = np.where(self.f_pred==filt)
            xx = self.x_pred[__]
            yy = self.y_pred[__]
            yye = self.y_prede[__]
            __max = np.argmax(yy)            
            self.tpeak[filt] = xx[__max]
            self.fpeak[filt] = ( yy[__max], yye[__max] )
