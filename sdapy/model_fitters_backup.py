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
from multiprocessing import Pool, cpu_count
from sdapy.functions import *
from sdapy.corner_hack import quantile, corner_hack
from sdapy.filters import central_wavelengths
from sdapy import constants

def get_engine():
    mod = __import__('sdapy', fromlist=['engines'])
    enginelist = dict()
    for _engine in mod.engines.__all__:
        if not _engine in mod.engines.__dict__:
            print ('Warning: skip engine %s since not registered in engines'%\
                   _engine)
            continue
        fit_engine = mod.engines.__dict__[_engine].engine
        enginelist[_engine] = fit_engine
    return enginelist

def get_pars(which):
    enginelist = get_engine()
    which = which.lower()    
    mod = __import__('sdapy', fromlist=['models'])
    # accumulate all alias and engine
    allwhichs = dict()
    for _which in mod.models.__all__:
        if not _which in mod.models.__dict__:
            print ('Warning: skip model %s since not registered in models'%\
                   _which)
            continue        
        fit_pars = mod.models.__dict__[_which].parameters.modelmeta        
        for _subwhich in fit_pars:            
            assert 'engine' in fit_pars[_subwhich],\
                'Error: define an engine for %s.%s'%(_which, _subwhich)
            fit_engine = fit_pars[_subwhich]['engine']
            assert fit_engine in enginelist.keys(), \
                'Error: wrong engine found for %s, select from %s'%\
                (fit_engine, enginelist)        
            if 'alias' in fit_pars[_subwhich]:
                for __ in fit_pars[_subwhich]['alias']:
                    allwhichs[__.lower()] = (_which, _subwhich, fit_engine)
            else:
                allwhichs[_subwhich.lower()] = (_which, _subwhich, fit_engine)
    # make sure which is in alias
    if not which in allwhichs.keys():
        print('Error: model %s not exists, chosen from %s'%\
              (which, allwhichs.keys()))
        return    
    lv1, lv2, engine = allwhichs[which]
    _mod = mod.models.__dict__[lv1].parameters.modelmeta[lv2]    
    # func
    func = None
    for k in _mod['func'].split('.'):
        if func is None: func = mod.models.__dict__[lv1].__dict__[k]
        else: func = func.__dict__[k]              
    # output dict
    output = dict()
    output['func'] = func    
    output['engine'] = enginelist[engine]
    output['enginename'] = engine      
    output['parname'] = _mod['parname']
    output['par'] = _mod['par']     
    # if fit together?
    output['same_parameter'] = []
    if not 'same_parameter' in _mod: _mod['same_parameter']=[]
    for p in _mod['same_parameter']:
        assert p in _mod['par'], \
            'common parameter %s not in parlist %s?'%(p, _mod['par'])
        output['same_parameter'].append(p)    
    # if fit error?
    output['fit_error'] = _mod['fit_error']
    #output['fit_error'] = False
    #output['mcparname'] = _mod['parname']
    #output['mcpar'] = _mod['par'] 
    #if 'fit_error' in _mod:
    #    output['fit_error'] = _mod['fit_error']
    #    if output['fit_error'] and not 'fsigma' in output['mcpar']:
    #        output['mcpar'] = output['par'] + ['fsigma']
    #        output['mcparname'] = output['parname'] + [r'$f_\mathrm{\sigma}$']    
    # par bestv and bounds
    bestv, bounds1, bounds2 = [], [], []
    for _ in output['par']:
        bestv.append( eval('constants.%s_p0'%_) )
        bounds1.append( eval('constants.%s_bounds[0]'%_) )
        bounds2.append( eval('constants.%s_bounds[1]'%_) )
    output['bestv'] = bestv
    output['bounds'] = (bounds1, bounds2)
    #bestv, bounds1, bounds2 = [], [], []
    #for _ in output['mcpar']:
    #    bestv.append( eval('constants.%s_p0'%_) )
    #    bounds1.append( eval('constants.%s_bounds[0]'%_) )
    #    bounds2.append( eval('constants.%s_bounds[1]'%_) )
    #output['mcbestv'] = bestv
    #output['mcbounds'] = (bounds1, bounds2)
    return output

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
        Possible choices are: 'mcmc', 'minimize', 'leastsq'.
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
    def __init__(self,  x_data, y_data, yerr_data, filters=None):
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
        
    def train(self, opt_routine='curvefit', fit_mean='powerlaw', nwalkers=30,
              nsteps=20000, nsteps_burnin=50, ncores=27, thin_by=1, maxfev=20000,
              mcmc_h5_file='tmp', emcee_burnin=True, use_emcee_backend=True,
              clobber=False, verbose=True, datadir='./', sigma=3, scipysamples=100,
              t0=0, timedilation=1, xpredict=None):
        assert opt_routine in ['mcmc', 'minimize', 'leastsq']
        self.opt_routine = opt_routine
        self.nwalkers = nwalkers
        self.nsteps = nsteps
        self.nsteps_burnin = nsteps_burnin
        if ncores is None: self.ncores = cpu_count() - 1
        else:  self.ncores = ncores
        self.thin_by = thin_by
        self.maxfev = maxfev
        self.emcee_burnin = emcee_burnin
        self.use_emcee_backend = use_emcee_backend
        self.clobber = clobber
        self.verbose = verbose
        self.sigma = sigma
        self.scipysamples = scipysamples
        #self.sguess = sguess
        self.datadir = datadir
        self.t0 = t0
        self.timedilation = timedilation
        self.xpredict = xpredict
        
        _which = get_pars(fit_mean)            
        self.func = _which['func']
        self.bestv = _which['bestv']
        self.bounds = _which['bounds']        
        self.par = _which['par']
        self.parname = _which['parname']
        #self.mcbestv = _which['mcbestv']
        #self.mcbounds = _which['mcbounds']        
        #self.mcpar = _which['mcpar']
        #self.mcparname = _which['mcparname']
        self.fit_error =  _which['fit_error']       
        assert self.func is not None, 'define func correctly'
        assert self.bestv is not None, 'define bestv correctly'
        assert self.bounds is not None, 'define bounds correctly'        
        self.cachefile = '%s/%s'%(self.datadir, mcmc_h5_file)        
        if self.filters is not None:
            self.suffix = '_'
            for f in central_wavelengths:            
                if f in self.filters: self.suffix += f
        else:
            self.suffix = ''
        self.samples, self.lnprob, self.cl = dict(), dict(), []
        if self.filters is None:
            # bolometric lc
            #self.nll, self.lnposterior = self.nll1, self.lnposterior1
            self.lnposterior = self.lnposterior1
            self.samples, self.lnprob = self.run()
        elif len(_which['same_parameter'])>0 and opt_routine=='mcmc':
            # fit mcmc with all bands simultaneously, with some free parameters fitted in all bands
            #self.nll, self.lnposterior = self.nll2, self.lnposterior2
            self.lnposterior = self.lnposterior2
            # rearange parameters            
            for p in _which['same_parameter']:
                assert p in self.par                
                for _n, _p in enumerate(self.par):
                    if _p == p: self.cl.append(_n)
            assert len(self.cl) > 0
            #self.mcbestv *= len(np.unique(self.filters))
            #self.mcbounds = [self.mcbounds[0] * len(np.unique(self.filters)),
            #                 self.mcbounds[1] * len(np.unique(self.filters))]
            # run
            samples, lnprob = self.run()
            
            # assign samples to each bands
            nf = 0
            for f in central_wavelengths:
                if f in self.filters:                    
                    samp = np.ndarray(shape=(len(samples),len(self.mcpar)), dtype=float)
                    for nrow, row in enumerate(samples):                        
                        for col in range(nf*len(self.mcpar), (nf+1)*len(self.mcpar)):
                            if col%len(self.mcpar) in self.cl:
                                samp[nrow, col%len(self.mcpar)] = row[col%len(self.mcpar)]
                            else:
                                samp[nrow, col%len(self.mcpar)] = row[col]
                    self.samples[f] = samp
                    self.lnprob[f]  = lnprob
                    nf += 1
        else:
            # fit mcmc or scipy in each band successively
            #self.cl, self.nll, self.lnposterior = [], self.nll1, self.lnposterior1
            self.cl, self.lnposterior = [], self.lnposterior1
            for f in np.unique(self.filters):
                self.samples[f], self.lnprob[f] = self.run(f)

    def run_scipy(self, x, y, yerr):
        samples = np.ndarray(shape=(self.scipysamples,len(self.par)), dtype=float)
        lnprob = np.ndarray(shape=(self.scipysamples,1), dtype=float)
        BBparams, covar = curve_fit(self.func,x, y, method='trf',
                                    sigma=yerr, p0=self.bestv, absolute_sigma=True,
                                    maxfev=self.maxfev, bounds=self.bounds)
        perr = np.sqrt(np.diag(covar))
        return BBparams, perr
    
    def run(self, filt=None):        
        if self.filters is None:
            cachefile = '%s%s'%(self.cachefile, self.suffix)
            x, y, yerr = self.x, self.y, self.yerr
            filters=None
        elif len(self.cl)>0:
            cachefile =  '%s%s'%(self.cachefile, self.suffix)
            x, y, yerr = self.x, self.y, self.yerr
            filters=self.filters
        else:            
            assert filt is not None
            cachefile = '%s_%s'%(self.cachefile, filt)
            _ = np.where(self.filters == filt)
            x = self.x[_]
            y = self.y[_]
            yerr = self.yerr[_]
            filters = self.filters[_]
        if not '.h5' in cachefile: cachefile = '%s.h5'%cachefile
        
        # check if cache exists or not       
        if os.path.exists(cachefile) and self.clobber:
            os.remove( cachefile )
        elif os.path.exists(cachefile) and not self.clobber:
            if self.opt_routine=='mcmc':
                # check steps, if sample steps less than demanded, continue the chain
                sample_nsteps = get_samples_nstep(cachefile, self.thin_by)
                if self.nsteps <= sample_nsteps:
                    print ('Warning: cached file exists includes %i steps, larger than %i, reload'%
                           (sample_nsteps,self.nsteps))
                    samples, lnprob = get_samples_mc(cachefile)                    
                    return samples, lnprob
                else:
                    print ('Warning: cached file exists includes %i steps, less than %i, go on for the further %i steps'% (sample_nsteps,self.nsteps,self.nsteps-sample_nsteps))
                    self.nsteps -= sample_nsteps                    
            else:                
                samples, lnprob = get_samples_scipy(cachefile)               
                return samples, lnprob
        
        if self.opt_routine in ['minimize', 'leastsq']: # run scipy fit
            samples = np.ndarray(shape=(self.scipysamples,len(self.par)), dtype=float)
            lnprob = np.ndarray(shape=(self.scipysamples,1), dtype=float)
            BBparams, covar = curve_fit(self.func,x, y, method='trf',
                        sigma=yerr, p0=self.bestv, absolute_sigma=True,
                        maxfev=self.maxfev, bounds=self.bounds)
            perr = np.sqrt(np.diag(covar))            
            
            for _n, (_1,_2) in enumerate(zip(BBparams + perr*self.sigma, BBparams - perr*self.sigma)):
                __ = np.random.uniform(low=_1, high=_2, size=(self.scipysamples))
                for _nn in range(self.scipysamples):  samples[_nn, _n] = __[_nn]
            for _theta in samples:                
                model = self.func(self.x, *_theta)
                lnprob[_nn, 0] = -0.5*np.sum((self.y - model)**2 / ((self.yerr)**2)) - \
                    np.sum(np.log(self.yerr)) - 0.5*len(model)*np.log(2*np.pi)
            # store samplles
            hf = h5py.File(cachefile, 'w')
            hf.create_dataset('samples', data=samples)
            hf.create_dataset('lnprob', data=lnprob)
            hf.close()
            return samples, lnprob
        
        elif self.opt_routine == 'mcmc':   # emcee fit
            args = (x, y, yerr, filters, cachefile)            
            if not os.path.exists(cachefile):  self.run_mcmc(*args)
            else:  self.continue_chains(*args)
            samples, lnprob = get_samples_mc(cachefile)         
            return samples, lnprob
        
    # optimization routine for hyperparameters
    def run_mcmc(self, t_data, f_data, f_unc_data, filters=None, cachefile='tmp.h5'):
        ''' initial fit '''
        t_mcmc_start = time.time()        
        if len(self.cl) > 0:
            args = (f_data, t_data, f_unc_data, filters, self.cl, self.func, self.mcbounds, self.fit_error)
        else:
            args = (f_data, t_data, f_unc_data, self.func, self.mcbounds, self.fit_error)
        #if self.sguess:            
        #    ml_res = minimize(self.nll, self.mcbestv, method='Powell', args=args)
        #    ml_guess = ml_res.x
        #    ndim = len(ml_guess)                 
        #    pos = ml_guess + 1e-4 * np.random.randn(self.nwalkers, ndim)            
        #    
        #else:
        #    ndim = len(self.mcbestv)            
        #    #initial position of walkers                          
        #    pos = self.mcbestv + 1e-4 * np.random.randn(self.nwalkers, ndim)

        with Pool(self.ncores) as pool:
            if self.emcee_burnin:
                burn_sampler = emcee.EnsembleSampler(self.nwalkers, ndim, 
                                        self.lnposterior, args=args, pool=pool)               
                pos, _, _ = burn_sampler.run_mcmc(pos, nsteps=self.nsteps_burnin,
                                                thin_by=self.thin_by, progress=self.verbose)
            
            if self.use_emcee_backend:            
                # file to save samples
                filename = cachefile
                backend = emcee.backends.HDFBackend(filename)
                backend.reset(self.nwalkers, ndim)
                
                sampler = emcee.EnsembleSampler(self.nwalkers, ndim, 
                                self.lnposterior, args=args, backend=backend, pool=pool)

            else:
                sampler = emcee.EnsembleSampler(self.nwalkers, ndim, 
                                self.lnposterior, args=args, pool=pool)
                
            max_samples = self.nsteps        
            old_tau = np.inf
            steps_so_far, tau = 0, None                
            
            for sample in sampler.sample(pos, iterations=max_samples,
                                         thin_by=self.thin_by, progress=False):            
                if sampler.iteration <= int(1e3/self.thin_by):
                    continue
                elif ((int(1e3/self.thin_by) < sampler.iteration <= int(1e4/self.thin_by)) 
                      and sampler.iteration % int(1e3/self.thin_by)):
                    continue
                elif ((int(1e4/self.thin_by) < sampler.iteration <= int(1e5/self.thin_by)) 
                      and sampler.iteration % int(1e4/self.thin_by)):
                    continue
                elif ((int(1e5/self.thin_by) < sampler.iteration) and 
                      sampler.iteration % int(2e4/self.thin_by)):
                    continue

                tstart = time.time()
                tau = sampler.get_autocorr_time(tol=0)
                tend = time.time()            
                steps_so_far = sampler.iteration

                if self.verbose:
                    print('''After {:d} steps, 
                    autocorrelation takes {:.3f} s ({} total FFTs)                
                    acceptance fraction = {:.4f}, and
                    tau = {}'''.format(steps_so_far, 
                                       tend-tstart, self.nwalkers*ndim,
                                       np.mean(sampler.acceptance_fraction), 
                                       tau))
                    
                # Check convergence
                converged = np.all(tau * 100 < sampler.iteration)
                converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                if converged:
                    break
                old_tau = tau

            if self.verbose:
                print("Ran {} steps; final tau= {}".format(steps_so_far*self.thin_by, tau))
                t_mcmc_end = time.time()
                print("All in = {:.2f} s on {} cores".format(t_mcmc_end - t_mcmc_start, self.ncores))

    def continue_chains(self, t_data, f_data, f_unc_data, filters=None, cachefile='tmp.h5'):
        '''Run MCMC for longer than initial fit'''
        t_mcmc_start = time.time()                

        if len(self.cl) > 0:            
            args = (f_data, t_data, f_unc_data, filters, self.cl, self.func, self.mcbounds, self.fit_error)
        else:
            args = (f_data, t_data, f_unc_data, self.func, self.mcbounds, self.fit_error)       
        with Pool(self.ncores) as pool:
            # file to save samples
            filename = cachefile
            new_backend = emcee.backends.HDFBackend(filename) 
            _, nwalkers, ndim = np.shape(new_backend.get_chain())
            new_sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnposterior, 
                                                args=args, pool=pool, backend=new_backend)
            max_samples = self.nsteps
            if max_samples <= 2e3:
                print ('Warning: nsteps should be larger than 2000')
                return
            steps_so_far = new_sampler.iteration
            old_tau = new_sampler.get_autocorr_time(tol=0)
            for i in range(int(max_samples/(2e3/self.thin_by))):
                new_sampler.run_mcmc(None, int(2e3/self.thin_by), 
                                     thin_by=self.thin_by, progress=False)
                tstart = time.time()
                tau = new_sampler.get_autocorr_time(tol=0)
                tend = time.time()
                steps_so_far = new_sampler.iteration
                print('''After {:d} steps, 
                autocorrelation takes {:.3f} s ({} total FFTs)                
                acceptance fraction = {:.4f}, and
                tau = {}'''.format(steps_so_far, 
                        tend-tstart, nwalkers*ndim,
                        np.mean(new_sampler.acceptance_fraction), 
                        tau))
                # Check convergence
                converged = np.all(tau * 100 < new_sampler.iteration)
                converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                if converged:
                    break
                old_tau = tau
                
        print("Ran {} steps; final tau= {}".format(steps_so_far*self.thin_by, tau))
        t_mcmc_end = time.time()
        print("All in = {:.2f} s on {} cores".format(t_mcmc_end - t_mcmc_start, self.ncores))
        
    def save_corner(self, show=False, figpath=None, limit=0, quantiles=[.05,.95]):
        ''' generate corner plots '''
        assert 'samples' in self.__dict__
        assert 'lnprob' in self.__dict__
        if self.opt_routine == 'mcmc': parname = self.mcparname
        else: parname = self.parname
        
        if self.filters is None:            
            samples = self.filter_samples(self.samples, self.lnprob, limit=limit)            
            cfig = corner_hack(samples, labels=parname,
                        label_kwargs={'fontsize':16}, ticklabelsize=13,
                        show_titles=True, quantiles=quantiles,
                        title_fmt=".2f", title_kwargs={'fontsize':16},
                        plot_datapoints=True, plot_contours=True)
            if figpath is not None:                
                cfig.savefig('{}/{}.png'.format(self.datadir, figpath), bbox_inches='tight')  
        else:
            if len(self.cl) == 0:
                # make corner plots for different bands seperately
                for f in np.unique(self.filters):
                    samples = self.filter_samples(self.samples[f], self.lnprob[f], limit=limit)
                    cfig = corner_hack(samples, labels=parname,
                        label_kwargs={'fontsize':16}, ticklabelsize=13,
                        show_titles=True, quantiles=quantiles,
                        title_fmt=".2f", title_kwargs={'fontsize':16},
                        plot_datapoints=True, plot_contours=True)                                        
                    if figpath is not None:
                        cfig.savefig('{}/{}_{}.png'.format(self.datadir, figpath, f), bbox_inches='tight')  
            else:
                # make corner plots for different bands together
                npar = (len(parname)-len(self.cl)) * len(np.unique(self.filters)) + len(self.cl)
                for f in central_wavelengths:
                    if f in self.filters:
                        fc = f
                        break
                samples = np.zeros((np.shape(self.samples[fc])[0],npar))
                paramsNames, nc = [], 0                                     
                for n in self.cl:                   
                    samples[:,nc] = self.samples[fc][:,n]
                    paramsNames += [parname[n]]
                    nc += 1
                for f in np.unique(self.filters):
                    for nn, n in enumerate(parname):                        
                        if nn % len(parname) in self.cl: continue                        
                        samples[:,nc] = self.samples[f][:,nn]
                        paramsNames += ['%s_%s'%(n, f)]
                        nc += 1
                cfig = corner_hack(samples, labels=paramsNames,
                        label_kwargs={'fontsize':16}, ticklabelsize=13,
                        show_titles=True, quantiles=quantiles,
                        title_fmt=".2f", title_kwargs={'fontsize':16},
                        plot_datapoints=True, plot_contours=True)                                        
                if figpath is not None:
                    cfig.savefig('{}/{}.png'.format(self.datadir, figpath), bbox_inches='tight')  
        if show: plt.show()
            
    def get_par(self, filt=None, quant=[.05,.5,.95], parname=None):
        assert 'samples' in self.__dict__
        mean, p1, p2 = [], [], []        
        par = self.par
        for _ in range(len(par)):
            if self.filters is None:
                xs = self.samples[:,_]
            else:
                assert filt is not None
                xs = self.samples[filt][:,_]                
            ql,q,qu = quantile(np.atleast_1d(xs), quant, weights=None)
            mean.append(q)
            p1.append(ql)
            p2.append(qu)
        if parname is None:  return p1, mean, p2
        
        assert parname in self.par
        __ = np.where(np.array(self.par) == parname)        
        return np.array(p1)[__][0], np.array(mean)[__][0], np.array(p2)[__][0]
        
    def predict(self, x_pred=None, step = 1, returnv=False, quant=[.05,.5,.95]):
        ''' output fitting products '''        
        if x_pred is None:
            if self.xpredict is not None:
                x_pred = self.xpredict
            else:
                x_pred = np.arange(self.x.min(), self.x.max(), step) 
        x_predl, y_predl, y_pred_1, y_pred_2, _ws = [], [], [], [], []      
        if self.filters is not None:            
            for _f in np.unique(self.filters):               
                q1,q,q2 = self.get_par(filt=_f, quant=quant)                
                y = self.func( x_pred, *q )
                y1 = self.func( x_pred, *q1 )
                y2 = self.func( x_pred, *q2 )
                x_predl.append( x_pred )
                y_predl.append( y )
                y_pred_1.append( y1 )
                y_pred_2.append( y2 )
                _ws.append( np.array( len(x_pred) * [_f]) )            
        else:
            q1,q,q2 = self.get_par(filt=None, quant=quant)            
            y = self.func( x_pred, *q )
            y1 = self.func( x_pred, *q1 )
            y2 = self.func( x_pred, *q2 )            
            x_predl, y_predl, y_pred_1, y_pred_2 = x_pred, y, y1, y2

        if returnv:
            return np.array(x_predl)*self.timedilation+self.t0, np.array(y_predl), \
                np.array(y_pred_1), np.array(y_pred_2), np.array(_ws)       
        self.x_pred = np.array(x_predl)*self.timedilation+self.t0
        self.y_pred = np.array(y_predl)
        self.y_pred1 = np.array(y_pred_1)
        self.y_pred2 = np.array(y_pred_2)
        self.f_pred = np.array(_ws)

    def predict_random(self, limit=0., plotnsamples=100, x_pred=None, step = 1e-3):
        ''' output fitting products '''        
        if x_pred is None:
            if self.xpredict is not None:
                x_pred = self.xpredict
            else:
                x_pred = np.arange(self.x.min(), self.x.max(), step)
        if not 'theta_samp' in self.__dict__:
            self.get_random_samples(limit=limit, plotnsamples=plotnsamples)
        x_predl, y_predl, _ws = [], [], []
        if self.filters is not None:
            for _f in np.unique(self.filters):                
                nsamples = min(len(self.theta_samp[_f]), plotnsamples)
                for i in range(nsamples):                    
                    par = self.theta_samp[_f][i]
                    #y = self.func( x_pred, *par[0:len(self.bestv)] )
                    y = self.func( x_pred, *par )
                    x_predl.append( x_pred )
                    y_predl.append( y )                               
                    _ws.append( _f )            
        else:
            nsamples = min(len(self.theta_samp), plotnsamples)
            for i in range(nsamples):
                par = self.theta_samp[i]
                y = self.func( x_pred, *par )                               
                x_predl.append( x_pred )
                y_predl.append( y )    
        return np.array(x_predl)*self.timedilation+self.t0, np.array(y_predl), np.array(_ws)

    def get_random_samples(self, limit=0., plotnsamples=100):
        assert 'samples' in self.__dict__
        assert 'lnprob' in self.__dict__
        
        if self.filters is None: 
            # best sample
            self.theta = self.samples[np.argmax(self.lnprob)]
        
            # other samples
            _samples = self.filter_samples(self.samples, self.lnprob, limit=limit)
            nsamples = min(len(_samples), plotnsamples)
            
            self.theta_samp = []
            for samp_num in np.random.choice(range(len(_samples)),nsamples, replace=False):
                self.theta_samp.append( _samples[samp_num] )
                
        else:
            self.theta, self.theta_samp = dict(), dict()
            for _f in np.unique(self.filters):
                # best sample
                self.theta[_f] = self.samples[_f][np.argmax(self.lnprob[_f])]
                
                # other samples
                _samples = self.filter_samples(self.samples[_f], self.lnprob[_f], limit=limit)
                nsamples = min(len(_samples), plotnsamples)
        
                self.theta_samp[_f] = []        
                for samp_num in np.random.choice(range(len(_samples)),nsamples, replace=False):
                    self.theta_samp[_f].append( _samples[samp_num] )
                    
    @staticmethod
    def filter_samples(samples, lnprob, limit=0.):        
        thre = min(lnprob) + (max(lnprob) - min(lnprob)) * limit        
        theta_samp = []
        for nn, samp_num in enumerate(samples):
            if lnprob[nn] >= thre: theta_samp.append( samp_num )        
        return np.array(theta_samp)

    def set_peak(self):
        if self.filters is None:
            xx = self.x_pred
            yy = self.y_pred
            y_pred1 = self.y_pred1
            y_pred2 = self.y_pred2
            __ = np.argmax(yy)
            yye = penc_to_errors(yy[__], y_pred1[__], y_pred2[__])                               
            self.tpeak = xx[np.argmax(yy)]        
            self.fpeak = ( max(yy), yye )
        else:
            self.tpeak, self.fpeak = dict(), dict()
            for _f in np.unique(self.filters):
                _ = np.where(self.f_pred==_f)
                xx = self.x_pred[_]
                yy = self.y_pred[_]
                y_pred1 = self.y_pred1[_]
                y_pred2 = self.y_pred2[_]
                __ = np.argmax(yy)
                yye = penc_to_errors(yy[__], y_pred1[__], y_pred2[__])                               
                self.tpeak[_f] = xx[np.argmax(yy)]        
                self.fpeak[_f] = ( max(yy), yye )
    
    #@staticmethod
    #def nll1(theta, f, t, f_err, func, bounds, fit_error):
    #    return -1*fit_model.lnlikelihood1(theta, f, t, f_err, func, fit_error)
    
    @staticmethod
    def lnlikelihood1(theta, f, t, f_err, func, fit_error):
        if fit_error:
            theta, f_sigma = theta[:-1], theta[-1]
            model = func(t, *theta)
            assert np.all(model > -np.inf)
            ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - \
                np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
        else:
            model = func(t, *theta)
            assert np.all(model > -np.inf)
            ln_l = -0.5*np.sum((f - model)**2 / f_err**2) + np.sum(np.log(1/np.sqrt(2*np.pi*f_err**2)))
        return ln_l
    
    @staticmethod
    def lnprior1(theta, bounds):        
        for n, t in enumerate(theta):
            if t < bounds[0][n] or t > bounds[1][n]:
                return -np.inf
        return 0
        
    @staticmethod
    def lnposterior1(theta, f, t, f_err, func, bounds, fit_error):
        lnp = fit_model.lnprior1(theta, bounds)
        lnl = fit_model.lnlikelihood1(theta, f, t, f_err, func, fit_error)
        if not np.isfinite(lnl):  return -np.inf
        if not np.isfinite(lnp):  return -np.inf
        return lnl + lnp

    #@staticmethod
    #def nll2(theta, f, t, f_err, filters, cl, func, bounds, fit_error):
    #    return -1*fit_model.lnlikelihood2(theta, f, t, f_err, filters, cl, func, fit_error)

    @staticmethod
    def lnlikelihood2(theta, f, t, f_err, filters, cl, func, fit_error):
        n_filt = len(np.unique(filters))
        n_theta = int(len(theta) / n_filt)        
        ln_l = 0
        f_num = 0        
        for filt in central_wavelengths:
            # !!! arange filter order depending on central_wavelengths
            if not filt in filters: continue            
            __theta = theta[f_num*n_theta : (f_num+1)*n_theta]
            # common parameters
            for _ in cl:  __theta[_] = theta[_]
            f_obs = np.where(filters == filt)
            ln_l += fit_model.lnlikelihood1(__theta, f[f_obs], t[f_obs], f_err[f_obs], func, fit_error)
            f_num += 1
        return ln_l
    
    @staticmethod
    def lnprior2(theta, filters, cl, bounds):
        n_filt = len(np.unique(filters))
        n_theta = int(len(theta) / n_filt)
        ln_p = 0
        f_num = 0        
        for filt in central_wavelengths:
            if not filt in filters: continue                     
            __theta = theta[f_num*n_theta : (f_num+1)*n_theta]
            # common parameters
            for _ in cl:  __theta[_] = theta[_]
            __bounds = [bounds[0][f_num*n_theta : (f_num+1)*n_theta],
                        bounds[1][f_num*n_theta : (f_num+1)*n_theta]]
            ln_p += fit_model.lnprior1(__theta, __bounds)
            f_num += 1
        return ln_p
    
    @staticmethod
    def lnposterior2(theta, f, t, f_err, filters, cl, func, bounds, fit_error):
        lnp = fit_model.lnprior2(theta, filters, cl, bounds)    
        lnl = fit_model.lnlikelihood2(theta, f, t, f_err, filters, cl, func, fit_error)
        if not np.isfinite(lnl):  return -np.inf
        if not np.isfinite(lnp):  return -np.inf
        return lnl + lnp