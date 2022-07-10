def engine(self,**kwargs):           
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
                
    if kwargs['shock_type'] == 0: return
    elif kwargs['shock_type'] == 1: clobber = False
    elif kwargs['shock_type'] == 2: clobber = True
    else: return

    if kwargs['shock_bolopt'] == 1:
        assert 'mbol' in self.__dict__, print ('constrcut bolometric lc first')
        mbol = self.mbol
    elif kwargs['shock_bolopt'] == 2:
        assert 'mbolbb' in self.__dict__, print ('constrcut bolometric lc first')
        mbol = self.mbolbb
    else:
        if kwargs['verbose']: print ('skip tail fit')
        return 
    assert self.t0 > 2400000, 'set t0 first'
        
    ncores = kwargs['ncores']
    if ncores is None: ncores = cpu_count() - 1
        
    xx, yy, yye = [], [], []
    for _ in kwargs['shock_copt']:
        if _ in mbol:
            xx = np.append(xx, get_numpy(mbol[_][0]))
            yy = np.append(yy, get_numpy(mbol[_][1]))
            yye = np.append(yye, get_numpy(mbol[_][2]))
        
    # shock fits
    # style=1 fit Me, Re, Ee
    # style=2 fit Me, Re, Ee, and tc
    xx = (xx-self.t0)/(1+self.z)
    if kwargs['shock_style'] == 1:
        fit_mean='shock_fit'                    
    elif kwargs['shock_style'] == 2:
        fit_mean='shoch_fit_tc'           
    else: return
    p1, p2 = min(kwargs['shock_fitr']), max(kwargs['shock_fitr'])
    __ = np.logical_and(xx>=p1, xx<=p2)
    
    self.shockcls = fit_model(xx[__],yy[__],yye[__],filters=None)         
    self.shockcls.train(opt_routine=kwargs['shock_routine'],
                    fit_mean=fit_mean, nwalkers=kwargs['nwalkers'],
                    nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'], ncores=ncores,
                    thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
                    mcmc_h5_file=kwargs['shock_cache']%
                             (self.objid, kwargs['shock_routine'], kwargs['shock_style']),
                    emcee_burnin=kwargs['emcee_burnin'], datadir=LOCALCACHE,
                    use_emcee_backend=kwargs['use_emcee_backend'], clobber=clobber,
                    verbose=kwargs['verbose'], p0=None, bounds=None, scipysamples=kwargs['scipysamples'],
                    sguess=kwargs['sguess'])
    self.shock_theta, self.shock_rtheta = self.get_random_samples(self.shockcls.samples,
                    self.shockcls.lnprob, cachefile=kwargs['shock_cache1']%
                            (self.objid, kwargs['shock_routine'], kwargs['shock_style']),                    
                    clobber=clobber, verbose=kwargs['verbose'],datadir=LOCALCACHE,) 
