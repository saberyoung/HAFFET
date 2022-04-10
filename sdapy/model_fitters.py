#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : model_fitters.py
# License           : BSD-3-Clause
# Author            : syang <sheng.yang@astro.su.se>
# Date              : 01.01.2020
# Last Modified Date: 14.02.2022
# Last Modified By  : syang <sheng.yang@astro.su.se>

import os
import numpy as np
from scipy.optimize import minimize, curve_fit
import matplotlib.pyplot as plt
import emcee, time
from multiprocessing import Pool
from .functions import *
from .corner_hack import quantile
from . import models
from .filters import central_wavelengths

class fit_model:
    """Fits data with power law.
    power law parts were from:
       Miller et al, https://ui.adsabs.harvard.edu/abs/2020ApJ...902...47M/abstract
    
    Parameters
    ----------
    x_data : array
        Independent values, e.g. (rest frame) phase relative to peak
    y_data : array
        Dependent values, e.g. fluxes
    yerr_data : array, int
        Dependent value errors, e.g. flux errors
    filters : array
        If filters available, will fit for each band simultaneously.   
    opt_routine : str, 
        Which technic to be used to realize optimization.
        Possible choices are: 'mcmc', lm’, ‘trf’, ‘dogbox’
    nwalkers : int
        if mcmc adopted, set walker number
    ncores  : int
        core numbers to be used for multi processing
    nsteps:  int
        if mcmc adopted, set step
    clobber: bool
        if power law already done, if redo it or not
    verbose: bool
        show progress or not

    Returns
    -------
    Returns the interpolated independent and dependent values with the 1-sigma standard deviation.
    """
    def __init__(self,  x_data, y_data, yerr_data=1e-8, filters=None):
        assert len(x_data) >= 3 and len(x_data) == len(y_data), 'Error: check input data...'
        try: # if pandas.dataframe
            self.x = x_data.to_numpy()
            self.y = y_data.to_numpy()
            self.yerr = yerr_data.to_numpy()
        except: # otherwise, numpy array
            self.x = x_data
            self.y = y_data
            self.yerr = yerr_data
        self.filters = filters
        self.y_norm = y_data.max()
        self.x_norm = x_data[np.argmax(y_data)]        
        
    def train(self, t_fl=18, opt_routine='curvefit', fit_mean='powerlaw',
              nwalkers=30, nsteps=20000, nsteps_burnin=50, ncores=27, thin_by=1, 
              maxfev=20000, mcmc_h5_file='tmp', emcee_burnin=True,
              use_emcee_backend=True, p0=None, bounds=None, clobber=False,
              verbose=True, datadir='./', sigma=3, scipysamples=100, sguess=True):        
        assert opt_routine in ['mcmc', 'lm', 'trf', 'dogbox']
        assert fit_mean in ['powerlaw', 'bazin', 'villar', 'arnett', 'tail_fitt',
                            'tail', 'arnett_fitv', 'arnett_fitt', 'arnett_tail']
        mcmc_h5_file = '%s/%s'%(datadir, mcmc_h5_file)
        if not '.h5' in mcmc_h5_file:  mcmc_h5_file = '%s.h5'%mcmc_h5_file        
        self.opt_routine = opt_routine
        self.fit_mean = fit_mean
        
        # check if cache exists or not
        if 'samples' in self.__dict__:
            if clobber: pass
            else:  return
        elif not os.path.exists(mcmc_h5_file):  pass
        elif clobber:  os.remove( mcmc_h5_file )
        else:
            if opt_routine=='mcmc':
                self.samples, self.lnprob = get_samples_mc(mcmc_h5_file)
            else:
                self.samples, self.lnprob = get_samples_scipy(mcmc_h5_file)
            return
        
        # optimization routine for hyperparameters
        def run_mcmc(nllfunc, lnposteriorfunc, t_data, f_data,
                     f_unc_data, guess_0, filters=None, mcmc_h5_file='tmp.h5'):
            t_mcmc_start = time.time()

            
            if filters is None: args = (f_data, t_data, f_unc_data)
            else: args = (f_data, t_data, f_unc_data, filters)

            if sguess:
                ml_res = minimize(nllfunc, guess_0, method='Powell', args=args)
                ml_guess = ml_res.x            
                ndim = len(ml_guess)                 
                pos = ml_guess + 1e-4 * np.random.randn(nwalkers, ndim)
            else:
                ndim = len(guess_0)            
                #initial position of walkers                          
                pos = guess_0 + 1e-4 * np.random.randn(nwalkers, ndim)
            
            with Pool(ncores) as pool:
                if emcee_burnin:                   
                    burn_sampler = emcee.EnsembleSampler(nwalkers, ndim, 
                                                         lnposteriorfunc, 
                                                         args=args,
                                                         pool=pool)
                    #burn_sampler.run_mcmc(pos, nsteps=nsteps_burnin, 
                    #                      thin_by=thin_by, progress=False)
                    #flat_burn_chain = burn_sampler.get_chain(flat=True)
                    #flat_burn_prob = np.argmax(burn_sampler.get_log_prob(flat=True))
                    #max_prob = flat_burn_chain[flat_burn_prob]
                    pos, _, _ = burn_sampler.run_mcmc(pos, nsteps=nsteps_burnin, thin_by=thin_by, progress=verbose)
                    

                if use_emcee_backend:            
                    # file to save samples
                    filename = mcmc_h5_file
                    backend = emcee.backends.HDFBackend(filename)
                    backend.reset(nwalkers, ndim)
                    
                    sampler = emcee.EnsembleSampler(nwalkers, ndim, 
                                                    lnposteriorfunc, 
                                                    args=args,
                                                    backend=backend,
                                                    pool=pool)

                else:
                    sampler = emcee.EnsembleSampler(nwalkers, ndim, 
                                                    lnposteriorfunc, 
                                                    args=args,
                                                    pool=pool)
            
                max_samples = nsteps
        
                old_tau = np.inf
                steps_so_far, tau = 0, None                
                
                for sample in sampler.sample(pos, iterations=max_samples,
                                             thin_by=thin_by, progress=False):            
                    if sampler.iteration <= int(1e3/thin_by):
                        continue
                    elif ((int(1e3/thin_by) < sampler.iteration <= int(1e4/thin_by)) 
                          and sampler.iteration % int(1e3/thin_by)):
                        continue
                    elif ((int(1e4/thin_by) < sampler.iteration <= int(1e5/thin_by)) 
                          and sampler.iteration % int(1e4/thin_by)):
                        continue
                    elif ((int(1e5/thin_by) < sampler.iteration) and 
                          sampler.iteration % int(2e4/thin_by)):
                        continue

                    tstart = time.time()
                    tau = sampler.get_autocorr_time(tol=0)
                    tend = time.time()            
                    steps_so_far = sampler.iteration

                    if verbose:
                        print('''After {:d} steps, 
                        autocorrelation takes {:.3f} s ({} total FFTs)                
                        acceptance fraction = {:.4f}, and
                        tau = {}'''.format(steps_so_far, 
                                           tend-tstart, nwalkers*ndim,
                                           np.mean(sampler.acceptance_fraction), 
                                           tau))
            
                    # Check convergence
                    converged = np.all(tau * 100 < sampler.iteration)
                    converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                    if converged:
                        break
                    old_tau = tau

                if verbose:
                    print("Ran {} steps; final tau= {}".format(steps_so_far*thin_by, tau))
                    t_mcmc_end = time.time()
                    print("All in = {:.2f} s on {} cores".format(t_mcmc_end - t_mcmc_start, ncores))
                        
        if opt_routine in ['lm', 'trf', 'dogbox']:   # fit successively
            assert self.filters is None, 'for multi-band data, fit different filters successively'           
            if fit_mean == 'powerlaw':
                func = models.powerlaw_full
                if p0 is None: p0 = np.array([self.y_norm, -t_fl, 2, 6e1])
                assert len(p0) == 4                
                samples = np.ndarray(shape=(scipysamples,4), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)
                if bounds is None: bounds = ([0,-100,0,0], [np.inf,0,10,np.inf])
            elif fit_mean == 'bazin':
                func = models.bazin
                if p0 is None: p0 = np.array([self.y_norm, self.x_norm, 10, 10, 5])
                assert len(p0) == 5               
                samples = np.ndarray(shape=(scipysamples,5), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)
                if bounds is None: bounds = ([0,-100,0,0,0], [np.inf,100,120,60,np.inf])
            elif fit_mean == 'villar':
                func = models.villar
                if p0 is None: p0 = np.array([self.y_norm, self.x_norm, 10, 10, 5])
                assert len(p0) == 5               
                samples = np.ndarray(shape=(scipysamples,5), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)
                if bounds is None: bounds = ([0,-100,0,0,0], [np.inf,100,120,60,np.inf])
            elif fit_mean == 'arnett':
                func = models.Arnett_fit
                if p0 is None: p0 = np.array([.2, 10])
                assert len(p0) == 2                
                samples = np.ndarray(shape=(scipysamples,2), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)
                if bounds is None: bounds = ([0,0], [10,1000])
            elif fit_mean == 'arnett_fitt':
                func = models.Arnett_fit_withoutt
                if p0 is None: p0 = np.array([.2, 10, -t_fl])
                assert len(p0) == 3                
                samples = np.ndarray(shape=(scipysamples,3), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)
                if bounds is None: bounds = ([0,0,-100], [10,1000,0])
            elif fit_mean == 'arnett_fitv':
                func = models.Arnett_fit_withoutv
                if p0 is None: p0 = np.array([.2, 1, 1])
                assert len(p0) == 3                
                samples = np.ndarray(shape=(scipysamples,3), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)
                if bounds is None: bounds = ([0,0,0], [10,1000,50])
                print ("Warning: you'd better do mc for arnett fits without v, since it's quite hard to get coverged via scipy")
            elif fit_mean == 'tail':
                func = models.tail_nickel
                if p0 is None: p0 = np.array([.2, 100])
                assert len(p0) == 2                
                samples = np.ndarray(shape=(scipysamples,2), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)  
                if bounds is None: bounds = ([0,0], [10,1000])
            elif fit_mean == 'tail_fitt':
                func = models.tail_nickel_fitt
                if p0 is None: p0 = np.array([.2, 100, 60])
                assert len(p0) == 3             
                samples = np.ndarray(shape=(scipysamples,3), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)  
                if bounds is None: bounds = ([0,0,20], [10,1000,150])
            elif fit_mean == 'arnett_tail':
                func = models.arnett_tail
                if p0 is None: p0 = np.array([.2, 10, 100, -t_fl, 60])
                assert len(p0) == 5
                samples = np.ndarray(shape=(scipysamples,5), dtype=float)
                lnprob = np.ndarray(shape=(scipysamples,1), dtype=float)  
                if bounds is None: bounds = ([0,0,0,-100, 20], [10,1000,1000,0, 150])
            
            BBparams, covar = curve_fit(func,self.x,self.y, method=opt_routine,
                        sigma=self.yerr, p0=p0, absolute_sigma=True, maxfev=maxfev, bounds=bounds)
            perr = np.sqrt(np.diag(covar))
            
            for _n, (_1,_2) in enumerate(zip(BBparams + perr*sigma, BBparams - perr*sigma)):                
                __ = np.random.uniform(low=_1, high=_2, size=(scipysamples))
                for _nn in range(scipysamples):  samples[_nn, _n] = __[_nn]
            for _theta in samples:                
                model = func(self.x, *_theta)
                lnprob[_nn, 0] = -0.5*np.sum((self.y - model)**2 / ((self.yerr)**2)) - np.sum(np.log(self.yerr)) - 0.5*len(model)*np.log(2*np.pi)

            hf = h5py.File(mcmc_h5_file, 'w')
            hf.create_dataset('samples', data=samples)
            hf.create_dataset('lnprob', data=lnprob)
            hf.close()

            self.samples, self.lnprob = get_samples_scipy(mcmc_h5_file)
            
        elif opt_routine == 'mcmc':    # fit simultaneously                            
            if self.filters is None:
                if fit_mean == 'powerlaw':
                    nf, pf = models.nll_pls, models.lnposterior_pls
                    if p0 is None: p0 = np.array([-t_fl, 6e1, self.y_norm, 2, 1])
                    else: assert len(p0) == 5
                elif fit_mean == 'bazin':
                    nf, pf = models.nll_bzs, models.lnposterior_bzs
                    if p0 is None: p0 = np.array([self.y_norm, self.x_norm, 10, 10, 5, 1])
                    else: assert len(p0) == 6
                elif fit_mean == 'arnett':
                    nf, pf = models.nll_arnett, models.lnposterior_arnett
                    if p0 is None: p0 = np.array([.2, 10, 1])
                    else: assert len(p0) == 3
                elif fit_mean == 'arnett_fitv':
                    nf, pf = models.nll_arnett_fitv, models.lnposterior_arnett_fitv
                    if p0 is None: p0 = np.array([.2, 10, 1, 1])
                    else: assert len(p0) == 4
                elif fit_mean == 'arnett_fitt':
                    nf, pf = models.nll_arnett_fitt, models.lnposterior_arnett_fitt
                    if p0 is None: p0 = np.array([.2, 10, -t_fl, 1])
                    else: assert len(p0) == 4
                elif fit_mean == 'tail':
                    nf, pf = models.nll_tail, models.lnposterior_tail
                    if p0 is None: p0 = np.array([.2, 100, 1])
                    else: assert len(p0) == 3
                elif fit_mean == 'tail_fitt':
                    nf, pf = models.nll_tail_fitt, models.lnposterior_tail_fitt
                    if p0 is None: p0 = np.array([.2, 100, 60, 1])
                    else: assert len(p0) == 4
                elif fit_mean == 'arnett_tail':
                    nf, pf = models.nll_arnett_tail, models.lnposterior_arnett_tail
                    if p0 is None: p0 = np.array([.2, 10, 100, -t_fl, 60])
                    else: assert len(p0) == 5
            else:
                assert fit_mean not in ['arnett', 'arnett_fitv',
                    'arnett_fitt', 'tail', 'tail_fitt', 'arnett_tail'], 'fit only for bolometric LCs with out filters'
                if fit_mean == 'powerlaw':
                    nf, pf = models.nll_pl, models.lnposterior_pl
                    if p0 is None:
                        p0 = np.array([-t_fl] + [6e1, 1, 2, 1]*len(np.unique(self.filters)))
                    else: assert len(p0) == 1 + 4*len(np.unique(self.filters))
                elif fit_mean == 'bazin':
                    nf, pf = models.nll_bz, models.lnposterior_bz
                    if p0 is None:
                        p0 = np.array([self.y_norm, self.x_norm, 10, 10, 5, 1]*len(np.unique(self.filters)))
                    else: assert len(p0) == 6*len(np.unique(self.filters))             
                    
            run_mcmc(nf, pf, self.x, self.y, self.yerr, p0,
                     filters=self.filters, mcmc_h5_file=mcmc_h5_file)
            
            self.samples, self.lnprob = get_samples_mc(mcmc_h5_file)           

    def get_par(self, filt=None, quant=[.05,.5,.95]):
        assert 'samples' in self.__dict__
        nf = 0 
        if self.filters is not None and filt is not None:
            _forder = dict()               
            for __filt in central_wavelengths:
                if __filt in self.filters:
                    _forder[__filt] = nf
                    nf += 1
            nf = _forder[filt]
            
        mean, p1, p2 = [], [], []
        if self.fit_mean == 'powerlaw':
            for _ in [0, 1+4*nf, 2+4*nf, 3+4*nf]:                
                if self.opt_routine == 'mcmc': _ += nf
                xs = self.samples[:,_]
                ql,q,qu = quantile(np.atleast_1d(xs), quant, weights=None)
                mean.append(q)
                p1.append(ql)
                p2.append(qu)
        elif self.fit_mean == 'bazin':
            for _ in [0+5*nf, 1+5*nf, 2+5*nf, 3+5*nf, 4+5*nf]:
                if self.opt_routine == 'mcmc': _ += nf
                xs = self.samples[:,_]                
                ql,q,qu = quantile(np.atleast_1d(xs), quant, weights=None)
                mean.append(q)
                p1.append(ql)
                p2.append(qu)
        elif self.fit_mean in ['arnett','tail']:
            for _ in [0, 1]:
                xs = self.samples[:,_]
                ql,q,qu = quantile(np.atleast_1d(xs), quant, weights=None)
                mean.append(q)
                p1.append(ql)
                p2.append(qu)
        elif self.fit_mean in ['arnett_fitv', 'arnett_fitt', 'tail_fitt']:
            for _ in [0, 1, 2]:
                xs = self.samples[:,_]
                ql,q,qu = quantile(np.atleast_1d(xs), quant, weights=None)
                mean.append(q)
                p1.append(ql)
                p2.append(qu)
        elif self.fit_mean in ['arnett_tail',]:
            for _ in [0, 1, 2, 3, 4]:
                xs = self.samples[:,_]
                ql,q,qu = quantile(np.atleast_1d(xs), quant, weights=None)
                mean.append(q)
                p1.append(ql)
                p2.append(qu)
        return p1, mean, p2

    def predict(self, x_pred=None, step = 1e-3, clobber=False,
                returnv=False, quant=[.05,.5,.95]):
        ''' output fitting products '''
        if 'x_pred' in self.__dict__ and not clobber and not returnv: return
        if x_pred is None:
            x_pred = np.arange(self.x.min(), self.x.max(), step)        
        if self.fit_mean == 'powerlaw':  func = models.powerlaw_full        
        elif self.fit_mean == 'bazin':   func = models.bazin
        elif self.fit_mean == 'arnett':   func = models.Arnett_fit
        elif self.fit_mean == 'arnett_fitv': func = models.Arnett_fit_withoutv
        elif self.fit_mean == 'arnett_fitt': func = models.Arnett_fit_withoutt
        elif self.fit_mean == 'tail':   func = models.tail_nickel
        elif self.fit_mean == 'tail_fitt':   func = models.tail_nickel_fitt
        elif self.fit_mean == 'arnett_tail':   func = models.arnett_tail
        
        x_predl, y_predl, y_pred_1, y_pred_2, _ws = [], [], [], [], []        
        if not self.filters is None:            
            for _f in np.unique(self.filters):               
                q1,q,q2 = self.get_par(filt=_f, quant=quant)
                y = func( x_pred, *q )
                y1 = func( x_pred, *q1 )
                y2 = func( x_pred, *q2 )                  
                x_predl.append( x_pred )
                y_predl.append( y )
                y_pred_1.append( y1 )
                y_pred_2.append( y2 )                
                _ws.append( np.array( len(x_pred) * [_f]) )            
        else:
            q1,q,q2 = self.get_par(filt=None, quant=quant)            
            y = func( x_pred, *q )
            y1 = func( x_pred, *q1 )
            y2 = func( x_pred, *q2 )            
            x_predl, y_predl, y_pred_1, y_pred_2 = x_pred, y, y1, y2

        if returnv:
            return np.array(x_predl), np.array(y_predl), np.array(y_pred_1), np.array(y_pred_2), np.array(_ws)       
        self.x_pred = np.array(x_predl)
        self.y_pred = np.array(y_predl)
        self.y_pred1 = np.array(y_pred_1)
        self.y_pred2 = np.array(y_pred_2)
        self.f_pred = np.array(_ws)  
