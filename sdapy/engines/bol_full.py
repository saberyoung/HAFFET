def engine(self,**kwargs):        
    for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
    if kwargs['joint_type'] == 0: return
    elif kwargs['joint_type'] == 1: clobber = False
    elif kwargs['joint_type'] == 2: clobber = True
    else: return
    
    if kwargs['joint_bolopt'] == 1:
        assert 'mbol' in self.__dict__, print ('constrcut bolometric lc first')
        mbol = self.mbol
    elif kwargs['joint_bolopt'] == 2:
        assert 'mbolbb' in self.__dict__, print ('constrcut bolometric lc first')
        mbol = self.mbolbb
    else:
        if kwargs['verbose']: print ('skip tail fit')
        return 
    assert self.t0 > 2400000, 'set t0 first'
        
    ncores = kwargs['ncores']
    if ncores is None: ncores = cpu_count() - 1
        
    xx, yy, yye = [], [], []
    for _ in kwargs['joint_copt']:
        if _ in mbol:
            xx = np.append(xx, get_numpy(mbol[_][0]))
            yy = np.append(yy, get_numpy(mbol[_][1]))
            yye = np.append(yye, get_numpy(mbol[_][2]))
            
    # joint fits
    # style=1 fit Arnett (Mni, taum) + tail (Mni, t0) + ts
    # style=2 fit Arnett (Mni, Mej, Ekin, vm) + tail (Mni, Mej, Ekin, vm) + ts
    # style=3 same fit as style=1, with texp as free parameter as well
    # style=4 same fit as style=2, with texp as free parameter as well        
    xx = (xx-self.t0)/(1+self.z)
    if kwargs['joint_style'] == 1:
        fit_mean='joint_fit_taum_t0'
        assert 'texp' in self.__dict__
        xx -= self.texp[1]
    elif kwargs['joint_style'] == 2:
        fit_mean='joint_fit_Mej_Ek'
        assert 'texp' in self.__dict__
        xx -= self.texp[1]
    elif kwargs['joint_style'] == 3:
        fit_mean='joint_fit_taum_t0_texp'
    elif kwargs['joint_style'] == 4:
        fit_mean='joint_fit_Mej_Ek_texp'                   
    else: return
        
    p1, p2 = min(kwargs['joint_fitr']), max(kwargs['joint_fitr'])
    __ = np.logical_and(xx>=p1, xx<=p2)
        
    self.jointcls = fit_model(xx[__],yy[__],yye[__],filters=None)         
    self.jointcls.train(opt_routine=kwargs['joint_routine'],
                    fit_mean=fit_mean, nwalkers=kwargs['nwalkers'],
                    nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'], ncores=ncores,
                    thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
                    mcmc_h5_file=kwargs['joint_cache']%
                             (self.objid, kwargs['joint_routine'], kwargs['joint_style']),
                    emcee_burnin=kwargs['emcee_burnin'], datadir=LOCALCACHE,
                    use_emcee_backend=kwargs['use_emcee_backend'], clobber=clobber,
                    verbose=kwargs['verbose'], p0=None, bounds=None, scipysamples=kwargs['scipysamples'],
                    sguess=kwargs['sguess'])
    self.joint_theta, self.joint_rtheta = self.get_random_samples(self.jointcls.samples,
                    self.jointcls.lnprob, cachefile=kwargs['joint_cache1']%
                            (self.objid, kwargs['joint_routine'], kwargs['joint_style']),                    
                    clobber=clobber, verbose=kwargs['verbose'],datadir=LOCALCACHE,)
    if kwargs['joint_style'] in [3,4]:
        # define texp here
        _, _1, _2 = self.jointcls.get_par(filt=None, quant=kwargs['quantile'])               
        self.texp = [_[2], _1[2], _2[2]]
