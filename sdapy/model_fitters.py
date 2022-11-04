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
    ''' get all possible engines of HAFFET
    
    Returns
    ---------- 
    engines:         `list`
        return a list of engines
    '''
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

def get_model(which_engine=None, with_alias=False):
    ''' get all possible models of a given engine
    

    Parameters
    ----------
    which_engine:         `str`
        if None, will return all models of HAFFET    
    with_alias:         `bool`
        if return only model names, or names + alias names
    
    Returns
    ---------- 
    models:         `list`
        return a list of models
    '''
    enginelist = get_engine()
    mod = __import__('sdapy', fromlist=['models'])
    # accumulate all alias and engine
    allwhichs = dict()    
    for _which in mod.models.__all__:            
        if not _which in mod.models.__dict__:
            print ('Warning: skip model %s since not registered in models'%_which)
            continue        
        fit_pars = mod.models.__dict__[_which].parameters.modelmeta        
        for _subwhich in fit_pars:            
            assert 'engine' in fit_pars[_subwhich],\
                'Error: define an engine for %s.%s'%(_which, _subwhich)
            fit_engine = fit_pars[_subwhich]['engine']
            assert fit_engine in enginelist.keys(), \
                'Error: wrong engine found for %s, select from %s'%\
                (fit_engine, enginelist)
            if which_engine is not None and which_engine != fit_engine: continue
            if fit_engine not in allwhichs:  allwhichs[fit_engine] = []
            allwhichs[fit_engine].append(_subwhich.lower())
            if 'alias' in fit_pars[_subwhich] and with_alias:
                for __ in fit_pars[_subwhich]['alias']: allwhichs[fit_engine].append(__.lower())
    return allwhichs

def get_pars(which, with_alias=True):
    ''' get parameter of a given model
    

    Parameters
    ----------
    which:         `str`
        model name
    with_alias:         `bool`
        if check only model names, or check names as well as alias names
    
    Returns
    ---------- 
    pars:         `dict`
        return a list of parameters, e.g. model function, parameter name, range, etc
    '''
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
            allwhichs[_subwhich.lower()] = (_which, _subwhich, fit_engine)
            if 'alias' in fit_pars[_subwhich] and with_alias:
                for __ in fit_pars[_subwhich]['alias']:
                    allwhichs[__.lower()] = (_which, _subwhich, fit_engine)
    
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
    output['fit_error'] = False
    if 'fit_error' in _mod:
        output['fit_error'] = _mod['fit_error']       
    # guess and bounds
    bestv, bounds1, bounds2 = [], [], []
    for _ in output['par']:        
        bestv.append( _mod['bestv'][_] )
        bounds1.append( _mod['bounds'][_][0] )
        bounds2.append( _mod['bounds'][_][1] )
    output['bestv'] = bestv
    output['bounds'] = (bounds1, bounds2)    
    return output

class fit_model:
    """Fits data with analytic models.
    
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
    """
    def __init__(self,  x_data, y_data, yerr_data, filters=None):
        assert len(x_data) >= 2 and len(x_data) == len(y_data), 'Error: check input data...'
        try: # if pandas.dataframe
            self.x = x_data.to_numpy()
            self.y = y_data.to_numpy()            
        except: # otherwise, numpy array
            self.x = x_data
            self.y = y_data            
        if yerr_data is not None:                                    
            assert len(y_data) == len(yerr_data), 'Error: check input data...'
            try: # if pandas.dataframe            
                self.yerr = yerr_data.to_numpy()
            except: # otherwise, numpy array            
                self.yerr = yerr_data
        else:
            self.yerr = self.y/10.            
        self.filters = filters
        
    def train(self, opt_routine='curvefit', fit_mean='powerlaw', nwalkers=30,
              nsteps=20000, nsteps_burnin=50, ncores=27, thin_by=1, maxfev=20000,
              mcmc_h5_file='tmp', emcee_burnin=True, use_emcee_backend=True,
              clobber=False, verbose=True, datadir='./', sigma=3, scipysamples=100,
              t0=0, timedilation=1, xpredict=None, bestv=None, bounds=None):
        """Fits data with analytic models.
    
        Parameters
        ----------             
        opt_routine : str
               Which technic to be used to realize optimization.
               Possible choices are: 'mcmc', 'minimize', 'leastsq'.
        fit_mean : str
               model mean function
        nwalkers : int
               if mcmc adopted, set walker number
        nsteps : int
               steps for mcmc run
        nsteps_burnin : int
               steps for mcmc burnin process
        ncores  : int
               core numbers to be used for multi processing
        thin_by : int
               steps devided by a factor for running reports
        maxfev : int
               parameter for scipy minimize
        mcmc_h5_file : str
               cached sampling file name
        emcee_burnin : bool
               if run burnin process or not for mcmc
        use_emcee_backend : bool
               if use emcee backend file or not
        datadir : str
               data directory of cached files
        sigma  : int
               for minimize routine, how to get samplings with mean and std
        scipysamples : int
               sampling numbers of minimize routine
        t0    : int
               zeropoint in x axis
        timedilation : int
               scaling factor in x axis
        xpredict : list
               prediction range
        bestv :  list
               bestfit for parameters
        bounds : list
               boundaries for parameters
        clobber: bool
               if power law already done, if redo it or not
        verbose: bool
               show progress or not
        """
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
        self.datadir = datadir
        self.t0 = t0
        self.timedilation = timedilation
        self.xpredict = xpredict
        
        _which = get_pars(fit_mean)
        self.func = _which['func']
        if bestv is None:  self.bestv = _which['bestv']
        else: self.bestv = bestv
        if bounds is None: self.bounds = _which['bounds']
        else: self.bounds = bounds
        
        self.par = _which['par']
        self.parname = _which['parname']        
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
            self.lnposterior = self.lnposterior1
            # guess best fit
            try:
                p0, _ = self.run_scipy(self.func, self.bestv, self.bounds,
                    self.x, self.y, self.yerr, sigma=self.sigma, maxfev=self.maxfev)
            except:
                p0 = self.bestv
            self.samples, self.lnprob = self.run(p0, self.bounds)
        elif len(_which['same_parameter'])>0 and opt_routine=='mcmc':
            # fit mcmc with all bands simultaneously, with some free parameters fitted in all bands
            self.lnposterior = self.lnposterior2
            # rearange parameters            
            for p in _which['same_parameter']:
                assert p in self.par                
                for _n, _p in enumerate(self.par):
                    if _p == p: self.cl.append(_n)
            assert len(self.cl) > 0
            
            # run
            pos = []
            for f in np.unique(self.filters):
                _ = np.where(self.filters == f)
                try:
                    p0, _ = self.run_scipy(self.func, self.bestv, self.bounds,
                        self.x[_], self.y[_], self.yerr[_], sigma=self.sigma, maxfev=self.maxfev)
                except:
                    p0 = self.bestv
                _bounds1, _bounds2 = self.bounds[0], self.bounds[1]
                if self.fit_error: p0 = np.append( p0, constants.fsigma_p0 )                    
                pos = np.append(pos, p0)
            if not self.fit_error:
                bounds = [self.bounds[0] * len(np.unique(self.filters)),
                          self.bounds[1] * len(np.unique(self.filters))]
                nn = len(self.par)
            else:
                bounds = [(self.bounds[0]+[constants.fsigma_bounds[0]]) * len(np.unique(self.filters)),
                          (self.bounds[1]+[constants.fsigma_bounds[1]]) * len(np.unique(self.filters))]
                nn = len(self.par) + 1
            samples, lnprob = self.run(pos, bounds)            
            
            # assign samples to each bands            
            nf = 0
            for f in central_wavelengths:
                if f in self.filters:                    
                    samp = np.ndarray(shape=(len(samples),nn), dtype=float)
                    for nrow, row in enumerate(samples):                        
                        for col in range(nf*nn, (nf+1)*nn):
                            if col%nn in self.cl: samp[nrow, col%nn] = row[col%nn]
                            else:  samp[nrow, col%nn] = row[col]
                    self.samples[f] = samp
                    self.lnprob[f]  = lnprob
                    nf += 1
        else:
            # fit mcmc or scipy in each band successively
            #self.cl, self.nll, self.lnposterior = [], self.nll1, self.lnposterior1
            self.cl, self.lnposterior = [], self.lnposterior1
            for f in np.unique(self.filters):
                _ = np.where(self.filters == f)
                try:
                    p0, _ = self.run_scipy(self.func, self.bestv, self.bounds,
                        self.x[_], self.y[_], self.yerr[_], sigma=self.sigma, maxfev=self.maxfev)
                except:
                    p0 = self.bestv
                self.samples[f], self.lnprob[f] = self.run(p0, self.bounds, filt=f)
    
    def run(self, bestv, bounds, filt=None):
        ''' run fitting
        '''
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
        if os.path.exists(cachefile) and self.clobber: os.remove( cachefile )
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
                    print ('Warning: cached file exists includes %i steps, less than %i, go on for the further %i steps'%
                           (sample_nsteps,self.nsteps,self.nsteps-sample_nsteps))
                    self.nsteps -= sample_nsteps                    
            else:                
                samples, lnprob = get_samples_scipy(cachefile)               
                return samples, lnprob
        
        if self.opt_routine in ['minimize', 'leastsq']: # run scipy fit            
            samples, lnprob = self.run_scipy(self.func, bestv, bounds, x, y, yerr,
                                sigma=self.sigma, maxfev=self.maxfev, nsamp=self.scipysamples)
            # store samplles
            hf = h5py.File(cachefile, 'w')
            hf.create_dataset('samples', data=samples)
            hf.create_dataset('lnprob', data=lnprob)
            hf.close()
            return samples, lnprob
        
        elif self.opt_routine == 'mcmc':   # emcee fit                        
            args = (x, y, yerr, bestv, bounds, filters, cachefile)
            if not os.path.exists(cachefile):  self.run_mcmc(*args)
            else:  self.continue_chains(*args)
            samples, lnprob = get_samples_mc(cachefile)         
            return samples, lnprob
        
    # optimization routine for hyperparameters
    def run_mcmc(self, t_data, f_data, f_unc_data, bestv, bounds, filters=None, cachefile='tmp.h5'):
        ''' Initial MCMC fit '''
        t_mcmc_start = time.time() 
        if len(self.cl) > 0:            
            args = (f_data, t_data, f_unc_data, filters, self.cl, self.func, bounds, self.fit_error)
        else:
            args = (f_data, t_data, f_unc_data, self.func, bounds, self.fit_error)            
        ndim = len(bestv)
        pos = bestv + 1e-4 * np.random.randn(self.nwalkers, ndim)        
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

    def continue_chains(self, t_data, f_data, f_unc_data, bestv, bounds, filters=None, cachefile='tmp.h5'):
        '''Run MCMC for longer than initial fit'''
        t_mcmc_start = time.time()   
        
        if len(self.cl) > 0:            
            args = (f_data, t_data, f_unc_data, filters, self.cl, self.func, bounds, self.fit_error)
        else:
            args = (f_data, t_data, f_unc_data, self.func, bounds, self.fit_error)       
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
        
    def save_corner(self, figpath, datadir=None, filts=None, limit=0, quantiles=[.05,.95], clobber=False):
        ''' generate corner plots '''
        assert 'samples' in self.__dict__
        assert 'lnprob' in self.__dict__
        
        parname = self.parname
        if self.fit_error and not r'$f_\mathrm{\sigma}$' in parname:
            parname += [r'$f_\mathrm{\sigma}$']
        if datadir is None: datadir = self.datadir
        
        fignames = []
        if self.filters is None:
            _fign = '{}/{}.png'.format(datadir, figpath)            
            if os.path.exists(_fign) and not clobber: pass
            else:
                if os.path.exists(_fign): os.remove(_fign)
                samples = self.filter_samples(self.samples, self.lnprob, limit=limit)            
                cfig = corner_hack(samples, labels=parname,
                        label_kwargs={'fontsize':16}, ticklabelsize=13,
                        show_titles=True, quantiles=quantiles,
                        title_fmt=".2f", title_kwargs={'fontsize':16},
                        plot_datapoints=True, plot_contours=True)                          
                cfig.savefig(_fign, bbox_inches='tight')
            fignames.append(_fign)
        else:
            if len(self.cl) == 0:
                # make corner plots for different bands seperately                
                for f in np.unique(self.filters):
                    if filts is not None and f not in filts: continue
                    _fign = '{}/{}_{}.png'.format(datadir, figpath, f)
                    if os.path.exists(_fign) and not clobber: pass
                    else:
                        if os.path.exists(_fign): os.remove(_fign)
                        samples = self.filter_samples(self.samples[f],
                                                      self.lnprob[f], limit=limit)
                        cfig = corner_hack(samples, labels=parname,
                                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                                           show_titles=True, quantiles=quantiles,
                                           title_fmt=".2f", title_kwargs={'fontsize':16},
                                           plot_datapoints=True, plot_contours=True)
                        cfig.savefig(_fign, bbox_inches='tight')
                    fignames.append(_fign)
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
                _fign = '{}/{}.png'.format(datadir, figpath)
                if os.path.exists(_fign) and not clobber: pass
                else:
                    if os.path.exists(_fign): os.remove(_fign)
                    cfig = corner_hack(samples, labels=paramsNames,
                                       label_kwargs={'fontsize':16}, ticklabelsize=13,
                                       show_titles=True, quantiles=quantiles,
                                       title_fmt=".2f", title_kwargs={'fontsize':16},
                                       plot_datapoints=True, plot_contours=True)
                    cfig.savefig(_fign, bbox_inches='tight')
                fignames.append(_fign)
        return fignames
    
    def get_par(self, filt=None, quant=[.05,.5,.95], parname=None):
        ''' get bestfit and quantiles of a fitted parameter '''
        assert 'samples' in self.__dict__
        if self.samples is None: return None, None, None
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
        ''' re-produce the data with the fitted models '''        
        if x_pred is None:
            if self.xpredict is not None:
                x_pred = np.arange(min(self.xpredict), max(self.xpredict), step)
            else:
                x_pred = np.arange(self.x.min(), self.x.max(), step) 
        x_predl, y_predl, y_pred_1, y_pred_2, _ws = [], [], [], [], []      
        if self.filters is not None:            
            for _f in np.unique(self.filters):               
                q1,q,q2 = self.get_par(filt=_f, quant=quant)
                if q1 is None: continue
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
            if q1 is not None:
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

    def predict_random(self, limit=0., plotnsamples=100, x_pred=None, step = 1):
        ''' re-produce the data randomly with the fitted models '''        
        if x_pred is None:
            if self.xpredict is not None:
                x_pred = np.arange(min(self.xpredict), max(self.xpredict), step)
            else:
                x_pred = np.arange(self.x.min(), self.x.max(), step)        
        if not 'theta_samp' in self.__dict__:
            self.get_random_samples(limit=limit, plotnsamples=plotnsamples)        
        x_predl, y_predl, _ws = [], [], []
        if 'theta_samp' in self.__dict__:            
            if self.filters is not None:
                for _f in np.unique(self.filters):
                    if _f not in self.theta_samp: continue
                    nsamples = min(len(self.theta_samp[_f]), plotnsamples)
                    for i in range(nsamples):                    
                        par = self.theta_samp[_f][i]                                        
                        y = self.func( x_pred, *par[0:len(self.bestv)] )                    
                        x_predl.append( x_pred )
                        y_predl.append( y )                               
                        _ws.append( _f ) 
            else:
                nsamples = min(len(self.theta_samp), plotnsamples)
                for i in range(nsamples):                
                    par = self.theta_samp[i]                
                    y = self.func( x_pred, *par[0:len(self.bestv)] )                
                    x_predl.append( x_pred )
                    y_predl.append( y )        
        return np.array(x_predl)*self.timedilation+self.t0, np.array(y_predl), np.array(_ws)

    def get_random_samples(self, limit=0.2, plotnsamples=100):
        ''' get parameters for the random samplings '''
        assert 'samples' in self.__dict__
        assert 'lnprob' in self.__dict__
        if self.samples is None: return
        
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
                if _f not in self.samples:continue
                
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
        ''' get samplings with good chi square '''
        thre = min(lnprob) + (max(lnprob) - min(lnprob)) * limit        
        theta_samp = []
        for nn, samp_num in enumerate(samples):
            if lnprob[nn] >= thre: theta_samp.append( samp_num )        
        return np.array(theta_samp)

    def set_peak(self):
        ''' get peak infos of the fitting '''
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

    @staticmethod
    def run_scipy(func, p0, bounds, x, y, yerr, sigma=3, maxfev=2000, nsamp=None):                
        ''' run model fitting with scipy packages '''
        try:
            params, covar = curve_fit(func, x, y, method='trf',
                sigma=yerr, p0=p0, absolute_sigma=True, maxfev=maxfev, bounds=bounds)
        except:
            #import matplotlib.pyplot as plt
            #plt.plot(x,y,'ko')
            #plt.savefig('tmp.png')
            params, covar = curve_fit(func, x, y, method='trf', p0=p0, maxfev=maxfev, bounds=bounds)        
        perr = np.sqrt(np.diag(covar))                
        if nsamp is None: return params, perr
        
        try:
            # https://stackoverflow.com/questions/70263573/generate-200-data-points-drawn-from-a-multivariate-normal-distribution-with-mean
            w, v = np.linalg.eig(covar)    
            sigmas = np.sqrt(w) * v
            A = sigma @ np.random.randn(len(params), nsamp) + np.array([params]).T
            for _n in range(len(params)):
                for _nn in range(nsamp):
                    samples[_nn, _n] = A.T[_nn][_n]        
        except:        
            # ignore covariance that are not diagonal
            samples = np.ndarray(shape=(nsamp,len(p0)), dtype=float)
            lnprob = np.ndarray(shape=(nsamp,1), dtype=float)
            for _n, (_1,_2) in enumerate(zip(p0 + perr*sigma, p0 - perr*sigma)):
                __ = np.random.uniform(low=_1, high=_2, size=(nsamp))
                for _nn in range(nsamp):  samples[_nn, _n] = __[_nn]
        for _theta in samples:
            model = func(x, *_theta)
            lnprob[_nn, 0] = -0.5*np.sum((y - model)**2 / ((yerr)**2)) - \
                np.sum(np.log(yerr)) - 0.5*len(model)*np.log(2*np.pi)
        return samples, lnprob
    
    @staticmethod
    def lnlikelihood1(theta, f, t, f_err, func, fit_error):
        ''' define likelhood for MCMC '''
        if fit_error:
            theta, f_sigma = theta[:-1], theta[-1]
            model = func(t, *theta)
            #assert np.all(model > -np.inf)
            ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - \
                np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
        else:
            model = func(t, *theta)
            #assert np.all(model > -np.inf)
            ln_l = -0.5*np.sum((f - model)**2 / f_err**2) + np.sum(np.log(1/np.sqrt(2*np.pi*f_err**2)))
        return ln_l
    
    @staticmethod
    def lnprior1(theta, bounds):
        ''' define prior for MCMC '''
        for n, t in enumerate(theta):
            if t < bounds[0][n] or t > bounds[1][n]:
                return -np.inf
        return 0
        
    @staticmethod
    def lnposterior1(theta, f, t, f_err, func, bounds, fit_error):
        ''' define posterior for MCMC '''
        lnp = fit_model.lnprior1(theta, bounds)
        lnl = fit_model.lnlikelihood1(theta, f, t, f_err, func, fit_error)
        if not np.isfinite(lnl):  return -np.inf
        if not np.isfinite(lnp):  return -np.inf
        return lnl + lnp

    @staticmethod
    def lnlikelihood2(theta, f, t, f_err, filters, cl, func, fit_error):
        ''' define likelhood for MCMC, with filters '''
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
        ''' define prior for MCMC, with filters '''
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
        ''' define posterior for MCMC, with filters '''
        lnp = fit_model.lnprior2(theta, filters, cl, bounds)    
        lnl = fit_model.lnlikelihood2(theta, f, t, f_err, filters, cl, func, fit_error)
        if not np.isfinite(lnl):  return -np.inf
        if not np.isfinite(lnp):  return -np.inf
        return lnl + lnp
