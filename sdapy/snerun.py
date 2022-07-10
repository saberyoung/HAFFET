#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/sne.py
# Author            : syang <sheng.yang@astro.su.se>
# Date              : 12.11.2019
# Last Modified Date: 11.02.2022
# Last Modified By  : syang <sheng.yang@astro.su.se>

from __future__ import print_function
import os, sys, requests, emcee, corner, random, urllib,\
    shlex, subprocess, time, json, math, re, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from scipy import integrate
from scipy.integrate import simps
from astropy.cosmology import Planck13 as cosmo
from astropy.time import Time
import astropy.units as u
from astropy import coordinates
from astropy.stats.sigma_clipping import sigma_clip
from joblib import dump, load
from io import StringIO
from pathlib import Path
from collections import OrderedDict
from astrorapid.classify import Classify
from bs4 import BeautifulSoup
import ztfquery
from ztfquery import marshal, fritz

# self import
from sdapy.gaussian_process import fit_gp
from sdapy.model_fitters import get_pars
from sdapy.filters import *
from sdapy.functions import *
from sdapy.specline_handler import handle_spectra
from sdapy.pbar import get_progress_bar
from sdapy import __path__, corner_hack, read_default

# all
__all__ = ('snelist', 'snobject')

# sdapy source path
srcpath = __path__[0]

# ztfquery source path
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def check_dir():
    answ = None
    if not os.path.isdir( LOCALSOURCE ):
        if answ is None:
            answ = input('as set, %s was the folder used to deal with data however not exists, create directory?(Y/N)'%LOCALSOURCE)
        if answ.strip() in ['y', 'Y']:
            try:
                os.mkdir( LOCALSOURCE )
            except:
                raise OSError("Can't create destination directory (%s)!" % (LOCALSOURCE))
        else:
            sys.exit()
    for filestype in ['fritz', 'marshal', 'ForcePhot', 'ForcePhot_atlas', 'cache', 'plots', 'images']:
        destdir = '%s/%s'%(LOCALSOURCE, filestype)
        if not os.path.isdir( destdir ):
            if answ is None:
                answ = input('as set, %s was the folder used to deal with data however some sub folder not exists, create them?(Y/N)'%(LOCALSOURCE))
            if answ.strip() in ['y', 'Y']:
                try:
                    os.mkdir( destdir )
                except:
                    raise OSError("Can't create destination directory (%s)!" % (destdir))
            else:
                sys.exit()
                
# parameters
snelist_pars = read_default.get_parameters(keylist='snelist')['snelist']
snobject_pars = read_default.get_parameters(keylist='snobject')['snobject']

# auth
keypairs = read_default.get_keypairs()
keypairs_type = read_default.get_keypairs(return_type=True)

def read_kwargs(kwargs, clsname):
    assert clsname in ['snelist', 'snobject']
    partype = read_default.get_parameters(return_type=True)[clsname]
    _kwargs = dict()
    for _key in kwargs:
        if not _key in partype: continue
        if partype[_key] == 'eval':
            _kwargs[_key] = eval( str(kwargs[_key]) )
        elif partype[_key] == 'str':
            _kwargs[_key] = str(kwargs[_key])
    return _kwargs

class snelist(object):

    # Static version info
    version = 1.0
    
    def __init__(self, ax=None, **kwargs):                
        """ initialize """
        # check dir
        check_dir()
        
        # ----- read default parameters ----- #        
        defkwargs = read_kwargs(snelist_pars, 'snelist')
        
        # ----- read meta ----- #
        self.kwargs = read_kwargs(kwargs, 'snelist')
        for _key in defkwargs: self.kwargs.setdefault(_key, defkwargs[_key])
        
        # ----- define products ----- #
        if not 'data' in self.__dict__: self.data = {}
        if not 'meta' in self.__dict__: self.meta = {}
        if not 'params' in self.__dict__: self.params = {}

        self.ax = ax
        self.dset = dict()
        
    def parse_meta(self, force=False, **kwargs):
        if len(self.meta) > 0 and not force: return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        metafile = '%s/data/meta.txt'%srcpath
        if not os.path.exists(metafile):
            print ('!!! Warning: no meta parsed, check if %s exsist' % metafile)
            return
        # read meta as pandas.dataframe
        skiprows = []
        for n,ll in enumerate(open(metafile).readlines()):
            if ll[0] == '#' or len(ll.strip()) == 0: skiprows.append(n)
        meta = pd.read_csv(metafile,sep=',', skiprows=skiprows).drop_duplicates()
        if len(kwargs['syntax']) != 0: meta = meta.query(kwargs['syntax'])
        self.meta = meta
        if kwargs['verbose']: print ( 'meta %i objs'%( len(meta)) )
        if len(meta) == 0:
            print ('!!! Warning: no meta data parsed')
            return        
        # format meta
        assert kwargs['idkey'] in self.meta.keys()
        self.meta = self.meta.set_index(kwargs['idkey'])
        if len(kwargs['sortkey']) != 0:
            assert kwargs['sortkey'] in self.meta.keys()
            self.meta = self.meta.sort_values(kwargs['sortkey'], ascending=False)

    def parse_meta_one(self, idkey, objid, key):        
        _meta = self.meta.query('%s==@objid'%idkey)        
        if len(_meta) == 0: return       
        if key in _meta:  return _meta[key][0]
        return

    def parse_meta_all(self, kwargs, objid):        
        # parse meta infos        
        self.ra = self.parse_meta_one(kwargs['idkey'], objid, kwargs['rakey'])
        self.dec = self.parse_meta_one(kwargs['idkey'], objid, kwargs['deckey'])
        self.z = self.parse_meta_one(kwargs['idkey'], objid, kwargs['zkey'])
        self.dist = self.parse_meta_one(kwargs['idkey'], objid, kwargs['distkey'])
        self.dm = self.parse_meta_one(kwargs['idkey'], objid, kwargs['dmkey'])
        self.mkwebv = self.parse_meta_one(kwargs['idkey'], objid, kwargs['mkwebvkey'])
        self.hostebv = self.parse_meta_one(kwargs['idkey'], objid, kwargs['hostebvkey'])
        self.mkwav = self.parse_meta_one(kwargs['idkey'], objid, kwargs['mkwavkey'])
        self.hostav = self.parse_meta_one(kwargs['idkey'], objid, kwargs['hostavkey'])
        self.sntype = self.parse_meta_one(kwargs['idkey'], objid, kwargs['typekey'])
        self.jdpeak = self.parse_meta_one(kwargs['idkey'], objid, kwargs['peaktkey'])        
        try:
            self.z = float(self.z)
        except:
            self.z = float(input('!!!  Error: %s redshift of %s should be a number\n\nmanully set z=' % (self.z, objid)))
        try:
            self.dist = float(self.dist)
        except:
            print ('Warning: %s distance set with standard cosmology' % (objid))
            self.dist = cosmo.luminosity_distance( self.z ).value
        try:
            self.dm = float(self.dm)
        except:
            self.dm = 5*np.log10(self.dist) + 25        
        try:
            self.mkwebv = float(self.mkwebv)
        except:
            try:
                self.mkwebv = float(self.mkwav) / kwargs['rv']
            except:
                print ('Warning: %s no mkw ebv and av was found, set mkw ebv=0' % (objid))
                self.mkwebv = 0        
        try:
            self.hostebv = float(self.hostebv)
        except:
            try:               
                self.hostebv = float(self.hostav) / kwargs['rv']
            except:
                print ('Warning: %s no host ebv and av was found, set host ebv=0' % (objid))
                self.hostebv = 0            
        try:
            self.jdpeak += kwargs['jdpeak_shift']
        except:
            pass
        
    def parse_params(self, force=False, **kwargs):
        ''' read individual_par.txt for particular SNe
        '''
        if len(self.params) > 0 and not force: return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        parfile = '%s/data/individual_par.txt'%srcpath
        if not os.path.exists(parfile):
            print ('!!! Warning: no prior information found, check if %s exsist' % parfile)
            return
        partype = read_default.get_parameters(return_type=True)['snobject']
        for line in open(parfile).readlines():
            if line[0]=='#' or len(line.strip())==0: continue
            obj = line.split()[0]                    
            if obj not in self.params: self.params[obj] = dict()
            for _ in line.split()[1:]:
                _key, _val = _.split('=')[0].strip(), _.split('=')[1].strip()
                if not _key in partype: continue
                if partype[_key] == 'eval':
                    self.params[obj][_key]=eval(_val)
                elif partype[_key] == 'str':                                
                    self.params[obj][_key]=str(_val)
                    
    def load_data(self, objid, force=False, datafile='%s_data.clf', **kwargs):
        if objid in self.data and not force: return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        datafile = '%s/cache/%s' % (LOCALSOURCE, datafile%objid ) 
        if os.path.exists(datafile):
            if kwargs['verbose']: print (  'load:', datafile )
            self.data[objid] = load(datafile)
            return True
        return False
    
    def save_data(self, objid, force=False, datafile='%s_data.clf', **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        datafile = '%s/cache/%s' % (LOCALSOURCE, datafile%objid)
        if kwargs['verbose']: print ('saved cache')
        if not os.path.exists(datafile) or force:
            if kwargs['verbose']: print (  'save:', datafile )
            dump(self.data[objid], datafile)
            return True
        return False
    
    def run(self, axes=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        self.parse_meta()
        self.parse_params()               
        with get_progress_bar(kwargs['verbose'], len(self.meta.index)) as pbar:            
            for i, objid in enumerate(self.meta.index):
                # parse meta infos
                self.parse_meta_all(kwargs, objid)
                
                # load data
                loaded = self.load_data(objid)
                
                if not loaded or kwargs['clobber']:
                    par = dict()
                    if objid in self.params: par = self.params[objid]
                    
                    # for each object                         
                    self.data[objid] = snobject(objid, z=self.z, ra=self.ra, dec=self.dec,
                            mkwebv=self.mkwebv, hostebv=self.hostebv, sntype=self.sntype,
                            dm=self.dm, jdpeak=self.jdpeak, axes=axes, **par)
                    
                    # run
                    self.data[objid].test_run()
                    
                    # save data
                    #saved = self.save_data(objid,force=True)
                    
                    sys.exit(objid)
                pbar.update(1)

    '''
    run/restore objects first, then distribute them in plots
    '''    
    def add_subset(self, syntax=None):
        _meta = self.meta
        if syntax is not None: _meta = _meta.query(syntax)
        self.dset[syntax] = {objid:self.data[objid] if objid in list(_meta.index) else None for objid in self.data}

    def get_par(self, objid, parname):        
        if parname == 'z': return self.data[objid].z
        if parname == 'ra': return self.data[objid].parse_coo(deg=True)[0]
        if parname == 'dec': return self.data[objid].parse_coo(deg=True)[1]
        if parname == 'mkwebv': return self.data[objid].mkwebv
        if parname == 'hostebv': return self.data[objid].hostebv
        if parname == 'dm': return self.data[objid].dm
        if parname == 't0': return self.data[objid].t0
        if parname == 'texp': return self.data[objid].texp[1]        
        
    def show1d(self, k, syntax=None, nbin=10, fontsize=12, labelpad=12, **kwargs):
        assert self.ax is not None
        self.add_subset(syntax=syntax)
        x = []
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue
            _x = self.get_par(objid, k)
            if _x is None:
                print ('can not get %s for %s'%(k, objid))
                continue
            x.append( float(_x) )            
        if len(x) == 0:
            print ('no data parsed')
            return
        binwidthx = (np.max(np.abs(x))-np.min(np.abs(x)))/nbin
        binx = np.arange(np.min(np.abs(x)), np.max(np.abs(x))+binwidthx, binwidthx)        
        self.ax.hist(x, bins=binx, histtype='step', **kwargs)
        self.ax.set_xlabel(k, fontsize=fontsize, labelpad=labelpad )
        self.ax.set_ylabel('counts', fontsize=fontsize, labelpad=labelpad )        
        
    def show2d(self, k1, k2, syntax=None, fontsize=12, labelpad=12, **kwargs):
        assert self.ax is not None
        self.add_subset(syntax=syntax)
        x, y = [], []
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue
            _x = self.get_par(objid, k1)
            if _x is None:
                print ('can not get %s for %s'%(k1, objid))
                continue
            _y = self.get_par(objid, k2)
            if _y is None:
                print ('can not get %s for %s'%(k2, objid))
                continue
            x.append( float(_x) )
            y.append( float(_y) )
        if len(x) == 0:
            print ('no data parsed')
            return
        self.ax.plot( x, y, **kwargs)
        self.ax.set_xlabel(k1, fontsize=fontsize, labelpad=labelpad )
        self.ax.set_ylabel(k2, fontsize=fontsize, labelpad=labelpad )
        self.ax.set_xlim([min(x)/1.1, max(x)*1.1])
        self.ax.set_ylim([min(y)/1.1, max(y)*1.1])

        # histograms        
        self.init_hist_axes(pad=0.1, labelbottom=False, labelleft=False)
        self.add_hist(x, y, syntax=None, nbinx = 10, nbiny = 10, xticks=None, yticks=None, **kwargs)
        
    def init_hist_axes(self, pad=0.1, labelbottom=False, labelleft=False):             
        assert self.ax is not None
        # create new axes on the right and on the top of the current axes
        divider = make_axes_locatable(self.ax)
        # below height and pad are in inches
        self.ax_histx = divider.append_axes("top", 1.2, pad=pad, sharex=self.ax)
        self.ax_histy = divider.append_axes("right", 1.2, pad=pad, sharey=self.ax)
        # make some labels invisible
        self.ax_histx.xaxis.set_tick_params(labelbottom=labelbottom)
        self.ax_histy.yaxis.set_tick_params(labelleft=labelleft)
        
    def add_hist(self, x, y, syntax=None, nbinx = 10, nbiny = 10, 
                 xticks=None, yticks=None, **kwargs):        
        assert 'ax_histx' in self.__dict__ and 'ax_histy' in self.__dict__, 'init_hist_axes() first'               
        assert len(x) > 0 and len(x) == len(y)
        
        # now determine nice limits by hand:
        binwidthx = (np.max(np.abs(x))-np.min(np.abs(x)))/nbinx
        binwidthy = (np.max(np.abs(y))-np.min(np.abs(y)))/nbiny
        binx = np.arange(np.min(np.abs(x)), np.max(np.abs(x))+binwidthx, binwidthx)
        biny = np.arange(np.min(np.abs(y)), np.max(np.abs(y))+binwidthy, binwidthy)
        self.ax_histx.hist(x, bins=binx, histtype='step', **kwargs)
        self.ax_histy.hist(y, bins=biny, orientation='horizontal', histtype='step', **kwargs)
        if yticks is not None: self.ax_histx.set_yticks(yticks)
        if xticks is not None: self.ax_histy.set_xticks(xticks)
        self.ax_histx.set_ylabel('Counts')
        self.ax_histy.set_xlabel('Counts')

    def showall(self, datatype, syntax=None, fontsize=12, labelpad=12,
                lc_cond='mag<99 and filter=="r"', color_interp=[1],
                **kwargs):
        assert self.ax is not None
        self.add_subset(syntax=syntax)
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue                        
            t0 = self.dset[syntax][objid].t0
            z = self.dset[syntax][objid].z
            dm = self.dset[syntax][objid].dm

            if datatype == 'lc':
                assert 'lc' in self.dset[syntax][objid].__dict__
                lc = self.dset[syntax][objid].lc
                if lc_cond is not None: lc = lc.query(lc_cond)
                self.ax.plot( (lc['jdobs']-t0)/(1+z), lc['mag']-dm, **kwargs )
                self.ax.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$', 
                                   fontsize=fontsize, labelpad=labelpad )
                self.ax.set_ylabel('M$_{abs}$ (mag)', fontsize=fontsize, labelpad=labelpad )
                #self.dset[syntax][objid]._ax(cond=lc_cond)

            if datatype == 'colour':
                assert 'colors' in self.dset[syntax][objid].__dict__
                colors = self.dset[syntax][objid].colors
                for _copt in color_interp:
                    if _copt in colors:
                        jd,g,r,ge,re = colors[_copt]
                        self.ax.plot( (jd-t0)/(1+z), g-r, **kwargs )
                self.ax.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$', 
                                   fontsize=fontsize, labelpad=labelpad )
                self.ax.set_ylabel('Colour (mag)', fontsize=fontsize, labelpad=labelpad )
                
            if datatype == 'mbol':
                assert 'mbol' in self.dset[syntax][objid].__dict__
                mbol = self.dset[syntax][objid].mbol            
                for _copt in color_interp:
                    if _copt in mbol:
                        self.ax.plot( (mbol[_copt][0]-t0)/(1+z), mbol[_copt][1], **kwargs )
                self.ax.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',
                                   fontsize=fontsize, labelpad=labelpad )
                self.ax.set_ylabel('L$_{bol}$ (erg $s^{-1}$)', fontsize=fontsize, labelpad=labelpad )
                        
    def showfit(self, syntax=None, filt='r', funit='mag',
                x_pred=None, step = 1, quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is not None:
                if 'fitcls' not in self.dset[syntax][objid].__dict__:
                    if self.verbose: print ('Warning: no fitcls for %s'%objid)
                    continue
                if filt not in self.dset[syntax][objid].fitcls:
                    if self.verbose: print ('Warning: no filter %s for %s fitcls'%(filt,objid))
                    continue
                z = self.dset[syntax][objid].z
                x, y, y1, y2, w = self.dset[syntax][objid].fitcls[filt].predict(x_pred=x_pred, 
                                        step=step, clobber=False, returnv=True, quant=quant)
                if funit == 'mag':
                    mm = self.dset[syntax][objid].flux_to_mag(y, dflux=None)
                    self.ax.plot( x/(1+z), mm-self.data[objid].dm,**kwargs )
                elif funit == 'flux':
                    if 'fpeak' not in self.dset[syntax][objid].__dict__:
                        if self.verbose: print ('Warning: no fmax for %s'%(objid))
                        continue
                    fmax = self.dset[syntax][objid].fpeak['fit'][filt][0]
                    self.ax.plot( x/(1+z), y/fmax,**kwargs )
    
    def showgp(self, syntax=None, filt='r', funit='mag',
                x_pred=None, step = 1, quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is not None:
                if 'gpcls' not in self.dset[syntax][objid].__dict__:
                    if self.verbose: print ('Warning: no fitcls for %s'%objid)
                    continue
                if filt not in self.dset[syntax][objid].gpcls.f_pred:
                    if self.verbose: print ('Warning: no filter %s for %s gpcls'%(filt,objid))
                    continue
                z = self.dset[syntax][objid].z
                t0 = self.dset[syntax][objid].t0
                __ = np.where(self.dset[syntax][objid].gpcls.f_pred==filt)
                x = self.dset[syntax][objid].gpcls.x_pred[__]
                y = self.dset[syntax][objid].gpcls.y_pred[__]
                ye = self.dset[syntax][objid].gpcls.y_prede[__]
                if funit == 'mag':
                    mm = self.dset[syntax][objid].flux_to_mag(y, dflux=None)
                    self.ax.plot( (x-t0)/(1+z), mm-self.data[objid].dm,**kwargs )
                elif funit == 'flux':
                    if 'fpeak' not in self.dset[syntax][objid].__dict__:
                        if self.verbose: print ('Warning: no fmax for %s'%(objid))
                        continue
                    fmax = self.dset[syntax][objid].fpeak['GP'][filt][0]
                    self.ax.plot( (x-t0)/(1+z), y/fmax,**kwargs )
                
    def showpl(self, syntax=None, filt='r', x1=-50, x2=0, limit=.9, **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is not None:
                if 'plcls' not in self.dset[syntax][objid].__dict__:
                    if self.verbose: print ('Warning: no plcls for %s'%objid)
                    continue
                _forder, nf = dict(), 0                
                for __filt in central_wavelengths:
                    if __filt in self.dset[syntax][objid].kwargs['pl_bands']:
                        _forder[__filt] = nf
                        nf += 1
                if filt in _forder:
                    nf = _forder[filt]                                    
                    theta = self.dset[syntax][objid].filter_samples(self.dset[syntax][objid].plcls.samples, 
                                                self.dset[syntax][objid].plcls.lnprob, limit=limit)[0]
                    adjust_samp = [theta[0],theta[1+4*nf],theta[2+4*nf],theta[3+4*nf]]
                    t_post = np.linspace(theta[0], x2, 1000)                
                    t_pre = np.linspace(x1, theta[0], 1000)
                    model_flux = powerlaw_post_baseline(t_post,
                         adjust_samp[2]*10**(-adjust_samp[3]), adjust_samp[0], adjust_samp[3])
                    self.ax.plot(t_pre, np.zeros_like(t_pre), **kwargs)
                    self.ax.plot(t_post, model_flux, **kwargs)
                    
    def showarnett(self, syntax=None, x_pred=None, step = 1, 
                   quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is not None:
                if 'arnettcls' not in self.dset[syntax][objid].__dict__:
                    if self.verbose: print ('Warning: no arnettcls for %s'%objid)
                    continue
                texp = self.dset[syntax][objid].texp
                x, y, y1, y2, w = self.dset[syntax][objid].arnettcls.predict(x_pred=x_pred, 
                                        step=step, clobber=False, returnv=True, quant=quant)
                if len(self.dset[syntax][objid].arnettcls.get_par()[0]) == 3:
                    x -= texp[1]
                self.ax.plot( x,y,**kwargs )
        
    def showtail(self, syntax=None, x_pred=None, step = 1, 
                   quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is not None:
                if 'tailcls' not in self.dset[syntax][objid].__dict__:
                    if self.verbose: print ('Warning: no tailcls for %s'%objid)
                    continue
                texp = self.dset[syntax][objid].texp
                x, y, y1, y2, w = self.dset[syntax][objid].tailcls.predict(x_pred=x_pred, 
                                        step=step, clobber=False, returnv=True, quant=quant)
                self.ax.plot( x,y,**kwargs )
                
    def showvelocity(self, syntax=None, quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is not None:
                if 'specls' not in self.dset[syntax][objid].__dict__:
                    if self.verbose: print ('Warning: no specls for %s'%objid)
                    continue
                x,y=[],[]
                for phase in self.dset[syntax][objid].specls:
                    _spec = self.dset[syntax][objid].specls[phase]
                    if _spec.velocity is None: continue
                    x.append( float(_spec.phase) )
                    y.append( float(_spec.velocity[1]) )
                self.ax.plot( x,y, **kwargs )
                
class snobject(object):

    # Static version info
    version = 1.0
    
    def __init__(self, objid, z=None, ra=None, dec=None,
                 mkwebv=None, hostebv=None, sntype=None, dm=None,
                 jdpeak=None, fig=None, ax=None, ax1=None, ax2=None,
                 ax3=None, ax4=None, **kwargs):
        # check dir
        check_dir()
        
        # ----- read default parameters ----- #        
        defkwargs = read_kwargs(snobject_pars, 'snobject')
        
        # ----- read meta ----- #
        self.kwargs = read_kwargs(kwargs, 'snobject')
        for _key in defkwargs: self.kwargs.setdefault(_key, defkwargs[_key])        
        
        # set meta infos
        self.objid  = objid
        self.ra     = ra
        self.dec    = dec
        if z is not None:
            self.z  = z # rest frame phase 
        else:
            self.z  = 0       # phase
        if mkwebv is not None:
            self.mkwebv = mkwebv
        else:
            self.mkwebv = 0  # ignore milky ebv
        if hostebv is not None:
            self.hostebv = hostebv
        else:           
            self.hostebv = 0        # ignore host ebv
        self.sntype = sntype
        self.dm     = dm
        if len(self.kwargs['set_tpeak_method']) == 0: assert jdpeak is not None
        if jdpeak is not None:
            assert jdpeak > 2400000.5, '!!! Error: use jd (instead of mjd)'
            self.t0 = jdpeak  # zero point thoughout the analysis                                  
        else:
            self.t0 = 0       # if None, try to get it via GP since it was needed for fits/pl/arnett...
        self.texp   = None    # explosion epoch        
        self.fpeak  = dict()  # peak flux in uJy for different bands       
        
        # plots
        self.fig = fig
        self.ax  = ax  # flux
        self.ax1 = ax1 # spectra
        self.ax2 = ax2 # mag        
        self.ax3 = ax3 # color
        self.ax4 = ax4 # luminosity
            
    def test_run(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        #self.init_fig()
        #self.config_ztfquery()
        #self.parse_coo(verbose=True, deg=True)
        #self.mjd_now(jd=False)
        if len(kwargs['lctype']) == 0 or 'ztffp' in kwargs['lctype']:
            self.get_fp_ztf()
            if not 'lc' in self.__dict__:
                self.query_fp_ztf()
                sys.exit('manual input ztf forced photometry file...')
            elif not 'ztffp' in list(self.lc['source']):
                print (  self.lc['source']  )
                input('...')
                self.query_fp_ztf()
                sys.exit('manual input ztf forced photometry file...')                 
        if len(kwargs['lctype']) == 0 or 'atlasfp' in kwargs['lctype']:
            self.get_fp_atlas(binDays=kwargs['tdbin'])
            if 'lc' not in self.__dict__:
                self.query_fp_atlas()
                self.get_fp_atlas(binDays=kwargs['tdbin'])
            elif not 'atlasfp' in list(self.lc['source']):
                self.query_fp_atlas()
                self.get_fp_atlas(binDays=kwargs['tdbin'])

        #if len(kwargs['spectype']) == 0 or 'marshal' in kwargs['spectype']:
        #    self.query_spectra(source='marshal')
        #self.query_spectra(source='fritz')        
        self.get_local_spectra()
        
        # after all lightcurve injection
        #if len(kwargs['set_texp_method'])==0: self.set_texp_midway()
        
        #self.clip_lc()
        #self.run_gp()
        #self.run_fit('multiband_early')
        #self.run_fit('multiband_main')
        self.run_fit('specline')
        
        #self.calc_colors()
        #self.lyman_bol()
        #self.run_fit('bol_main')        
        #self.set_texp_bol_main()        
        #print (self._color_at('g', 'r', 10,interpolation='fit'))
        #print (self._absmag_at_list('g', [2,5], interpolation='gp'))
        #print (self._absmag_at('r', 10, 'fit',))
        #self.est_hostebv_with_c10()
        #self.calc_colors()
        #self.lyman_bol()
        #self.bb_colors()
        #self.bb_bol()
        #print (self.mbolbb)
        #self.on_sntype()
        #self.check_lc()        
        #self.est_host_c10()
        #sys.exit (self.texp)
        #self.plot()
        #self.fig, self.ax1=plt.subplots(1,1)
        fig, ax=plt.subplots(1,1)
        #self._ax1()
        #self.savefig()
        #self.showfig()
        self.show_corner(ax, index=1, gp=False,
                    engine='specline', model=None, source=None)
        sys.exit()
        
    def run(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])         
        t_start = time.time()
        
        ''' parse local data via objid '''
        if len(kwargs['lctype']) == 0 or 'ztffp' in kwargs['lctype']:
            self.get_fp_ztf()
            if not 'lc' in self.__dict__:
                self.query_fp_ztf()
                sys.exit('manual input ztf forced photometry file...')
            elif not 'ztffp' in list(self.lc['source']):
                print (  self.lc['source']  )
                input('...')
                self.query_fp_ztf()
                sys.exit('manual input ztf forced photometry file...')                 
        if len(kwargs['lctype']) == 0 or 'atlasfp' in kwargs['lctype']:
            self.get_fp_atlas(binDays=kwargs['tdbin'])
            if 'lc' not in self.__dict__:
                self.query_fp_atlas()
                self.get_fp_atlas(binDays=kwargs['tdbin'])
            elif not 'atlasfp' in list(self.lc['source']):
                self.query_fp_atlas()
                self.get_fp_atlas(binDays=kwargs['tdbin'])
        
        #self.get_fp_ztf()  # obtain forced phot lc       
        self.clip_lc()  # remove lc outliers        
        
        ''' other infos '''        
        self.check_lc() # check lc quality
        self.on_sntype() # check sn type        
        if kwargs['verbose']:
            print("Prepare data in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        ''' intepolation'''       
        self.run_gp()     # Gaussian process
        if kwargs['verbose']:
            print("Run GP in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.set_texp_midway() # try to set texp
        
        self.run_fit()   # sn-like analytic functions
        if kwargs['verbose']:
            print("Run Fitting in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()

        self.check_lc()
        self.explosion_pl()  # power law for first light
        if kwargs['verbose']:
            print("Run early power law fit in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
                    
        ''' g-r '''
        self.calc_colors()  # g-r colour epochs
        if kwargs['verbose']:
            print("calc colour in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
        
        self.est_host_c10()   # use g-r compared to tpl for host ebv
        if kwargs['verbose']:
            print("est host ebv in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        if 1 in kwargs['bolopt']:
            self.lyman_bol()  # calculate luminosity from g,r with Lyman bol correction
            if kwargs['verbose']:
                print("calc lyman bol in = {:.2f} s".format(time.time() - t_start))
                t_start = time.time()
                
        if 2 in kwargs['bolopt']:
            self.bb_colors()  # calculate luminosity from BB
            if kwargs['verbose']:
                print("calc colours for BB in = {:.2f} s".format(time.time() - t_start))
                t_start = time.time()
                
            self.bb_bol()
            if kwargs['verbose']:
                print("calc colours for BB in = {:.2f} s".format(time.time() - t_start))
                t_start = time.time()
            
        self.shock_fit()
        if kwargs['verbose']:
            print("fit shock model in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()

        self.arnett_fit() # fit lums around peak to arnett model
        if kwargs['verbose']:
            print("fit arnett model in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.tail_fit()   # fit lums at tail to gamma leakage model
        if kwargs['verbose']:
            print("fit tail model in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()

        self.joint_fit()
        if kwargs['verbose']:
            print("fit joint model in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        ''' go to spectral part '''
        self.get_local_spectra()  # obtain local spectra
        if kwargs['verbose']:
            print("get local spectra in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.meas_specline()  # measure specific lines, fit photos v to break arnett degenaracy
        if kwargs['verbose']:
            print("fit spec modelline in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.plot()   # plot everything
        if kwargs['verbose']:
            print("plot in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
    def config_ztfquery(self):        
        if len(keypairs) > 0:
            assert 'fritz' in keypairs.keys()
            ztfquery.io.set_account('fritz', token=keypairs['fritz']['token'], force=True)
            assert 'marshal' in keypairs.keys()
            ztfquery.io.set_account('marshal', username=keypairs['marshal']['usr'],
                                    password=keypairs['marshal']['pwd'])
            ztfquery.io.test_irsa_account()
        else:
            print ('check keypairs')
            
    def parse_coo(self, verbose=True, deg=True):
        assert self.ra is not None and self.dec is not None
        degunits = False
        try:
            float(self.ra)
            degunits = True
            if verbose: print ('%s use unit deg'%self.ra)
        except:
            if verbose: print ('%s use unit hourangle'%self.ra)
        if degunits:
            c = coordinates.SkyCoord('%s %s'%(self.ra, self.dec), unit=(u.deg, u.deg))
        else:
            c = coordinates.SkyCoord('%s %s'%(self.ra, self.dec), unit=(u.hourangle, u.deg))            
        if deg: return c.ra.deg, c.dec.deg
        else:   return c.to_string('hmsdms').split()
        
    def mjd_now(self, jd=False):
        if jd: return Time.now().jd
        else:  return Time.now().mjd

    def add_lc(self, df, source=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        ck1 = 'mag' in df.keys() and 'emag' in df.keys() and 'jdobs' in df.keys() and 'filter' in df.keys()
        ck2 = 'flux' in df.keys() and 'eflux' in df.keys() and 'jdobs' in df.keys() and 'filter' in df.keys()
        assert ck1 or ck2, 'rename columns to mag/emag/jdobs/filter or flux/eflux/jdobs/filter'
        df['source'] = source
        if not 'lc' in self.__dict__:
            self.lc = df
        else:
            self.lc = self.lc.append(df, ignore_index=True).drop_duplicates()

    def _add_flux(self, df, zp, sigma):        
        if 'limmag' in df.keys() and not math.isnan(df['limmag'].to_list()[0]):            
            df = df.query('limmag<99')            
            df['flux'], df['eflux'] = self.mag_to_flux(df['mag'],
                                limmag=df['limmag'], sigma=sigma, zp=zp)
        elif 'emag' in df.keys() and not math.isnan(df['emag'].to_list()[0]):           
            df['flux'], df['eflux'] = self.mag_to_flux(df['mag'],
                                magerr=df['emag'], sigma=sigma, zp=zp)            
        else:
            df['flux'] = self.mag_to_flux(df['mag'], sigma=sigma, zp=zp)
        return df
    
    def add_flux(self, zp=23.9, source=None, **kwargs):
        '''
        zp = 23.9 for micro Jansky to AB mag        
        zp = 48.6 for ergs/s/cm2/Hz to AB mag
        e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        '''                       
        assert 'lc' in self.__dict__        
        assert 'mag' in self.lc.keys(), 'input mag'        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])

        if source is None:
            self.lc = self._add_flux(self.lc, zp, kwargs['snrt'])               
        else:            
            df1 = self.lc.query('source==@source')
            df2 = self.lc.query('source!=@source')
            __arr = []
            if len(df1) > 0:
                df1 = self._add_flux(df1, zp, kwargs['snrt'])                
                for i in df1.index: __arr.append( df1.loc[i].values )
                columns = df1.columns
            if len(df2) > 0:
                df2 = self._add_flux(df2, zp, kwargs['snrt'])                
                for i in df2.index: __arr.append( df2.loc[i].values )
                columns = df2.columns
            if len(__arr) > 0:
                self.lc = pd.DataFrame(__arr, columns=columns)
            else:
                print ('Error: no data found for source=%s'%source)

    def _add_mag(self, df, zp, sigma):
        if 'eflux' in df.keys():
            df = df.query('eflux>0')
            df['mag'], df['emag'], df['limmag'] = self.flux_to_mag(df['flux'],
                    dflux=df['eflux'], sigma=sigma, zp=zp)        
        else:
            df['mag'] = self.flux_to_mag(df['flux'], sigma=sigma, zp=zp) 
        return df
    
    def add_mag(self, zp=23.9, source=None, **kwargs):
        assert 'lc' in self.__dict__       
        assert 'flux' in self.lc.keys(), 'input flux'
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        if source is None:
            self.lc = self._add_mag(self.lc, zp, kwargs['snrt'])               
        else:            
            df1 = self.lc.query('source==@source')
            df2 = self.lc.query('source!=@source')            
            df1 = self._add_mag(df1, zp, kwargs['snrt']) 
            df2 = self._add_mag(df2, zp, kwargs['snrt'])
            __arr = []            
            for i in df1.index:            
                __arr.append( df1.loc[i].values )
            for i in df2.index:            
                __arr.append( df2.loc[i].values )                
            self.lc = pd.DataFrame(__arr, columns=self.lc.columns)  

    def get_external_phot(self, filename, source, **kwargs):
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])             
        if not os.path.exists(filename):
            print ('Error: %s not found'%f)
            return None
        _ = []
        for nn,ll in enumerate(open(filename).readlines()):
            if ll[0]=='#' or len(ll)==0:_.append(nn)
        df = pd.read_csv(filename,skiprows=_,delim_whitespace=True)
        
        self.add_lc(df, source=source, **kwargs)
        if 'mag' in df.keys() and not 'flux' in df.keys():
            self.add_mag(zp=23.9, **kwargs)
        if 'flux' in df.keys() and not 'mag' in df.keys():
            self.add_flux(zp=23.9, **kwargs)
    
    def bin_fp_atlas(self, binDays=3, resultsPath=None, outPath=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        cmd = 'python %s/plot_atlas_fp.py stack %.2f %s '%(srcpath, binDays, resultsPath)
        #if kwargs['jdmin'] is not None and kwargs['jdmax'] is not None:            
        #    cmd += '%f %f ' % (kwargs['jdmin']-2400000.5,kwargs['jdmax']-2400000.5)
        cmd += '--o %s'%(outPath)
        if kwargs['verbose']: print (cmd)
        pid = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,\
                               stderr=subprocess.PIPE)
        output,error = pid.communicate()        
        if kwargs['verbose']: print (output)
        if error: print (error)
        
    def get_fp_atlas(self, binDays=None, clobber=False, **kwargs):
        ''' get local ATLAS forced phot LCs
        or binned ATLAS forced phot LCs
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        targetdir = '%s/ForcePhot_atlas/'%LOCALSOURCE
        wdir = '%s/%s/'%(targetdir, self.objid)
        f_unbin = '%s/forcedphotometry_%s_lc.csv'%(wdir, self.objid)        
        if binDays is None:            
            if not os.path.exists(f_unbin):
                print ('Error: %s not found'%f_unbin)
                return None
            self._get_fp_atlas(f_unbin, **kwargs)
        else:            
            binDays = float(binDays)            
            f = '%s/forcedphotometry_%s_lc_atlas_fp_stacked_%.2f_days.txt'%\
                (wdir, self.objid, binDays)
            if not os.path.exists(f) or clobber:
                if not os.path.exists(f_unbin):
                    print ('Error: %s not found'%f_unbin)
                    return None
                self.bin_fp_atlas(binDays=binDays, resultsPath=f_unbin, outPath=wdir, **kwargs)
            self._get_fp_atlas(f, binned=True, **kwargs)
            
    def _get_fp_atlas(self, f, binned=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        sigma = kwargs['snrt']        
        if binned: skiprows=[0,1,2,3,4,5]
        else: skiprows=None
        df = pd.read_csv(f,skiprows=skiprows,sep = ',',).drop_duplicates().reset_index(drop=True)
        df.rename(columns={'uJy':'flux','duJy':'eflux','F':'filter',}, inplace=True)           
        df['jdobs'] = df['MJD'] + 2400000.5
        df['source'] = 'atlasfp'
        self.add_lc(df, source='atlasfp', **kwargs)
        self.add_mag(zp=23.9, **kwargs)
        
    def get_alert_ztf(self, source='marshal', **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert source in ['marshal', 'fritz']
        if source == 'marshal':
            df = marshal.get_local_lightcurves(self.objid)
        else:
            #df = fritz.download_lightcurve(self.objid, get_object=True)
            filename = fritz.FritzPhotometry._build_filename_(self.objid)
            df = fritz.FritzPhotometry.read_csv(filename).data
            df.rename(columns={'magerr':'emag',
                               'limiting_mag':'limmag'}, inplace=True)            
            for _k in df.keys():
                if _k not in ['filter','limmag','mag','emag','mjd','instrument_name']:
                    df.drop(_k, inplace=True, axis=1)            
            df['jdobs'] = df['mjd'] + 2400000.5 
            filters = np.array([''] * len(df))
            fmap = {'sdssg': 'g',
                    'sdssi': 'i',
                    'sdssr': 'r',
                    'sdssu': 'u',
                    'sdssz': 'z',
                    'uvot::b': 'B',
                    'uvot::u': 'U',
                    'uvot::uvm2': 'D',
                    'uvot::uvw1': 'A',
                    'uvot::uvw2': 'S',
                    'uvot::v': 'V',
                    'ztfg': 'g',
                    'ztfi': 'i',
                    'ztfr': 'r',}
            for _ in fmap: filters[df['filter']==_] = fmap[_]            
            df['filter'] = filters            
            mags,emags = df['mag'],df['emag']
            mags[np.isnan(mags)] = 99
            emags[np.isnan(emags)] = 99
            df['mag'] = mags
            df['emag'] = emags
        df['source'] = source        
        self.add_lc(df, source=source, **kwargs)
        self.add_flux(zp=23.9, **kwargs)
        
    def get_fp_ztf(self, seeing_cut = 7., **kwargs):
        ''' get local ZTF forced phot LCs
        ''' 
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        targetdir = '%s/ForcePhot/'%LOCALSOURCE
        sigma = kwargs['snrt']
        f = '%s/%s/forcedphotometry_%s_lc.csv'%(targetdir, self.objid, self.objid)
        if not os.path.exists(f):
            print ('Error: %s not found'%f)
            return None        
        _ = []
        for nn,ll in enumerate(open(f).readlines()):
            if ll[0]=='#' or len(ll)==0:_.append(nn)
        df = pd.read_csv(f,skiprows=_,delim_whitespace=True)
        for k in df.keys(): df = df.rename(columns={k: k[:-1]})
        
        # update df
        df = df.query('~(forcediffimflux == "NaN")').reset_index(drop=True)
        df.rename(columns={'forcediffimflux':'Fpsf',
                       'forcediffimfluxunc':'Fpsf_unc',
                       'forcediffimsnr':'Fpsf_snr',
                       'zpdiff':'zp',
                       'zpmaginpsciunc':'ezp',
                       'jd':'jdobs',
                       'forcediffimchisq':'chi2_red',
                       'sciinpseeing':'seeing'}, inplace=True)
        F0 = 10**(df['zp'].values/2.5)
        eF0 = F0 / 2.5 * np.log(10) * df['ezp'].values
        Fpsf = df['Fpsf'].values
        eFpsf = df['Fpsf_unc'].values
        # from Jy to micro Jy        
        df['flux'] = Fpsf / F0 * 3631 * 1e6
        df['eflux'] = np.sqrt( (eFpsf / F0)**2 + (Fpsf * eF0 / F0**2)**2 ) * 3631 * 1e6
        filt = df['filter']
        filterid, filters = np.zeros(len(df)), np.array([''] * len(df))
        filterid[filt=='ZTF_g']=1
        filterid[filt=='ZTF_r']=2
        filterid[filt=='ZTF_i']=3
        filters[filt=='ZTF_g']='g'
        filters[filt=='ZTF_r']='r'
        filters[filt=='ZTF_i']='i'
        df['filterid'] = filterid
        df['filter'] = filters
        df['chi2_red'] = np.array([np.float(x) for x in df['chi2_red'].values])    
        df['fcqfid'] = df['field']*10000 + df['ccdid']*100 + df['qid']*10 + df['filterid']              
        df = df.query('seeing < @seeing_cut')
        self.add_lc(df, source='ztffp', **kwargs)
        self.add_mag(zp=23.9, **kwargs)
        
    def query_fp_atlas(self, **kwargs):
        ''' ATLAS forced phot query
        https://fallingstar-data.com/forcedphot/static/apiexample.py
        '''                
        if len(keypairs) == 0:
            print ('Error: set atlas account in keypairs')
            return
        assert 'atlas' in keypairs.keys()
        BASEURL = "https://fallingstar-data.com/forcedphot"
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])

        # coordinates
        radeg, decdeg = self.parse_coo(verbose=kwargs['verbose'], deg=True)

        # query period
        try:
            mjdstart = float(kwargs['mjdstart'])
        except:
            assert self.t0 > 2400000, 'Error: specify mjdstart or t0 (in JD) first'
            mjdstart = float(self.t0) - 2400000.5 + float(kwargs['dstart'])
            if kwargs['verbose']: print ('mjdstart as %f days prior to the peak which is %.2f'%
                                         (float(kwargs['dstart']), self.t0-2400000.5))
        try:
            mjdend = float(kwargs['mjdend'])
        except:
            if len(kwargs['dend']) == 0:
                mjdend = self.mjd_now(jd=False)
                if kwargs['verbose']: print ('mjdend as current mjd %.2f'%mjdend)
            else:
                assert self.t0 > 2400000, 'Error: specify mjdend or t0 (in JD) first'
                mjdend = float(self.t0) - 2400000.5 + float(kwargs['dend'])
                if kwargs['verbose']: print ('mjdend as %f days after the peak which is %.2f'%
                                             (float(kwargs['dend']), self.t0-2400000.5))                
        targetdir = '%s/ForcePhot_atlas/'%LOCALSOURCE
        f = '%s/%s/forcedphotometry_%s_lc.csv'%(targetdir, self.objid, self.objid)
        if os.path.exists(f) and not kwargs['clobber']:
            if kwargs['verbose']: print ('file exists: %s'%f)
            return
        if os.environ.get('ATLASFORCED_SECRET_KEY'):
            token = os.environ.get('ATLASFORCED_SECRET_KEY')
            if kwargs['verbose']: print('Using stored token')
        else:
            data = {'username': keypairs['atlas']['usr'],
                    'password': keypairs['atlas']['pwd'],}
            resp = requests.post(url=f"{BASEURL}/api-token-auth/", data=data)

            if resp.status_code == 200:
                token = resp.json()['token']
                if kwargs['verbose']: 
                    print(f'Your token is {token}')
                    print('Store this by running/adding to your .zshrc file:')
                    print(f'export ATLASFORCED_SECRET_KEY="{token}"')
            else:
                if kwargs['verbose']: 
                    print(f'ERROR {resp.status_code}')
                    print(resp.text)
                return
            
        headers = {'Authorization': f'Token {token}', 'Accept': 'application/json'}
        task_url = None
        while not task_url:
            with requests.Session() as s:
                # alternative to token auth
                # s.auth = ('USERNAME', 'PASSWORD')
                resp = s.post(f"{BASEURL}/queue/", headers=headers, data={
                    'ra': radeg, 'dec': decdeg, 'mjd_min': mjdstart, 'mjd_max': mjdend, 'send_email': False})
                
                if resp.status_code == 201:  # successfully queued
                    task_url = resp.json()['url']
                    if kwargs['verbose']: print(f'The task URL is {task_url}')
                elif resp.status_code == 429:  # throttled
                    message = resp.json()["detail"]
                    if kwargs['verbose']: print(f'{resp.status_code} {message}')
                    t_sec = re.findall(r'available in (\d+) seconds', message)
                    t_min = re.findall(r'available in (\d+) minutes', message)
                    if t_sec:
                        waittime = int(t_sec[0])
                    elif t_min:
                        waittime = int(t_min[0]) * 60
                    else:
                        waittime = 10
                    if kwargs['verbose']: print(f'Waiting {waittime} seconds')
                    time.sleep(waittime)
                else:
                    if kwargs['verbose']: 
                        print(f'ERROR {resp.status_code}')
                        print(resp.text)
                    return
        
        result_url = None
        taskstarted_printed = False
        while not result_url:
            with requests.Session() as s:
                resp = s.get(task_url, headers=headers)

                if resp.status_code == 200:  # HTTP OK
                    if resp.json()['finishtimestamp']:
                        result_url = resp.json()['result_url']
                        if kwargs['verbose']:
                            print(f"Task is complete with results available at {result_url}")
                    elif resp.json()['starttimestamp']:
                        if not taskstarted_printed:
                            if kwargs['verbose']:
                                print(f"Task is running (started at {resp.json()['starttimestamp']})")
                            taskstarted_printed = True
                        time.sleep(2)
                    else:
                        if kwargs['verbose']:
                            print(f"Waiting for job to start (queued at {resp.json()['timestamp']})")
                        time.sleep(4)
                else:
                    if kwargs['verbose']: 
                        print(f'ERROR {resp.status_code}')
                        print(resp.text)
                    return

        with requests.Session() as s:
            textdata = s.get(result_url, headers=headers).text

        # if we'll be making a lot of requests, keep the web queue from being
        # cluttered (and reduce server storage usage) by sending a delete operation
        # s.delete(task_url, headers=headers).json()

        dfresult = pd.read_csv(StringIO(textdata.replace("###", "")), delim_whitespace=True)
        filepath = Path(f)  
        filepath.parent.mkdir(parents=True, exist_ok=True)
        dfresult.to_csv(f)

    def query_fp_ztf(self, get_email=True, **kwargs):
        ''' ZTF forced phot query
        via email still
        '''        
        if len(keypairs) == 0: return
        assert 'ztf' in keypairs.keys()
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        # coordinates
        radeg, decdeg = self.parse_coo(verbose=kwargs['verbose'], deg=True)

        # query period                    
        try:
            jdstart = float(kwargs['mjdstart']) + 2400000.5
        except:
            assert self.t0 > 2400000, 'Error: specify mjdstart or t0 (in JD) first'
            jdstart = float(self.t0) + float(kwargs['dstart'])
            if kwargs['verbose']: print ('jdstart as %f days prior to the peak which is %.2f'%
                                         (float(kwargs['dstart']), self.t0))                
        try:
            jdend = float(kwargs['mjdend']) + 2400000.5
        except:
            if len(kwargs['dend']) == 0:
                jdend = self.mjd_now(jd=True)
                if kwargs['verbose']: print ('jdend as current jd %.2f'%jdend)
            else:
                assert self.t0 > 2400000, 'Error: specify mjdend or t0 (in JD) first'
                jdend = float(self.t0) + float(kwargs['dend'])
                if kwargs['verbose']: print ('jdend as %f days after the peak which is %.2f'%
                                             (float(kwargs['dend']), self.t0))                
        targetdir = '%s/ForcePhot/'%LOCALSOURCE
        f = '%s/%s/forcedphotometry_%s_lc.csv'%(targetdir, self.objid, self.objid)        
        if os.path.exists(f) and not kwargs['clobber']:
            if kwargs['verbose']: print ('file exists: %s'%f)
            return
        
        if get_email:
            line = 'wget --http-user=ztffps --http-passwd=dontgocrazy!'+\
                ' -q -O log.txt '+ \
                '"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?'+\
                'ra=%.7f&dec=%.7f&jdstart=%.5f&jdend=%.5f&'%(radeg, decdeg, jdstart, jdend)+\
                'email=%s&userpass=%s"'%(keypairs['ztf']['email'], keypairs['ztf']['pwd'])
            subprocess.Popen(line, shell=True)
            print ('querying fp for %s (%.2f %.2f) from %.2f till %.2f\n->after reveive email, put content into %s' %
                   (self.objid, radeg, decdeg, jdstart, jdend, f))
            time.sleep(kwargs['waittime'])
            return
        
        # URL used to place a request        
        SendReqURL = 'https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi'

        # IPAC Auth (Don't change these)
        IUN = 'ztffps'
        IPW = 'dontgocrazy!'

        # Request Parameters to send (USER DETERMINED)
        SendRA = radeg
        SendDec = decdeg
        SendJDStart = jdstart
        SendJDEnd = jdend
        Email = keypairs['ztf']['email']
        UserPass = keypairs['ztf']['pwd']
        
        # Send request and return unformatted HTML output of get() method
        Request = requests.get( SendReqURL , auth=( IUN, IPW ),
                    params={'ra': SendRA, 'dec': SendDec, 'jdstart': SendJDStart, 'jdend': SendJDEnd,
                    'email': Email, 'userpass': UserPass})

        # Formatted (Parseable) HTML of request
        ReqSoup = BeautifulSoup(Request.text, 'html.parser')

        # For one pending job (ie: one request), check forced photometry job tables

        # Grab the table section of the HTML
        ReqTable = ReqSoup.find('table')
        # Limit to rows in said table
        for row in ReqTable.find_all('tr'):
            # Find the items in the row (omit header cells)
            cols=row.find_all('td')
            if len(cols) > 0:
                # Find RA, Dec, and JD values of request as recorded by the service
                # Use a nested for loop here if submitting multiple requests, append values to list/arr/etc
                ReqRA=float(cols[0].text.strip())
                ReqDec=float(cols[1].text.strip())
                ReqJDS=float(cols[2].text.strip())
                ReqJDE=float(cols[3].text.strip())

        # Check if request is fulfilled using a request status check and parse with beautifulsoup4
        # Note: Requests older than 30 days will not show up here
        # Iterable, as beautifulsoup can parse html columns and output lists

        TrueIfPending = True
        slp=True

        # URL for checking the status of jobs sent by user
        # Note: This webpage only updates once an hour, on the hour
        StatusURL = 'https://ztfweb.ipac.caltech.edu/cgi-bin/getForcedPhotometryRequests.cgi'

        # Open loop to periodically check the Forced Photometry job status page
        while(TrueIfPending):
   
            # Query service for jobs sent in past 30 days
            OutputStatus = requests.get(StatusURL, auth=(IUN, IPW),
                    params={'email': Email, 'userpass': UserPass, 'option': 'All recent jobs', 'action': 'Query Database'})

            # Check if job has ended and lightcurve file was created
            # Note: If an exotic error occurs and the ended field is not populated, this will go on forever
            # Table has 11 cols, reqid=0, ended=7, lc=10

            # Format HTML
            OutputSoup = BeautifulSoup(OutputStatus.text, 'html.parser')
            # Get Job information in HTML table
            OutputTable = OutputSoup.find('table')
            OutputEnded=''
            # Verify Table exists (If no requests have been sent in past 30 days or if service is updating, this will catch it)
            if OutputTable != '' and OutputTable is not None:
                # Parse Table rows
                for row in OutputTable.find_all('tr'):
                    # Parse Table entries
                    cols=row.find_all('td')
                    if len(cols) > 0:
                        # Check if values contained in a given row coorespond to the current request
                        # Use a nested for loop here to check all requests submitted if there are more than one pending
                        OutputRA=float(cols[1].text.strip())
                        OutputDec=float(cols[2].text.strip())
                        OutputJDS=float(cols[3].text.strip())
                        OutputJDE=float(cols[4].text.strip())
                        OutputEnded=cols[7].text.strip()
                        # Check if job is finished (OutputEnded)
                        # Check for equality between recorded request params and what is known to be
                        #    in the Job Status table, accounting for rounding
                        if(OutputEnded != '' and abs(OutputRA-ReqRA)<=0.00001 and abs(OutputDec-ReqDec)<=0.00001
                           and abs(OutputJDS-ReqJDS)<=0.00001 and abs(OutputJDE-ReqJDE)<=0.00001):
                            # Get end time of job and lightcurve path
                            OutputLC=cols[10].text.strip()
                            OutputEnded=cols[7].text.strip()
                            # Set open loop hooks
                            TrueIfPending = False
                            slp=False
            # Pace the script so it doesn't spam the service
            if( slp ): time.sleep(30)
            
        # If Job Status table values are not null, set the path for file to be downloaded to the path recorded in Job Status Table
        if OutputEnded != '':
            if OutputLC != '':
                ReqPath = OutputLC
            else:
                # Catch errors in processing and alert user
                print("Job is done but no lightcurve was produced, quitting")
                quit()
                
        # Set link to download lightcurve from
        ForcedPhotometryReqURL = f'https://ztfweb.ipac.caltech.edu{ReqPath}'
        
        # Set local path and filename for lightcurve (currently same as what is assigned by the service)
        LocalFilename = ForcedPhotometryReqURL.split('/')[-1]

        # Download the file with a get request
        with requests.get(ForcedPhotometryReqURL, stream=True, auth=(IUN, IPW)) as r:
            # Write to the local file in chunks in case file is large
            with open (LocalFilename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        
    def query_alert_ztf(self, source=None, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert source in ['marshal', 'fritz', None], 'source should be gm/fritz or None for both'
        if source in ['marshal',None]:
            marshal.download_lightcurve(self.objid,
                    auth=(keypairs['marshal']['usr'], keypairs['marshal']['pwd']),
                    verbose=kwargs['verbose'])
        if source in ['fritz',None]:
            fritz.download_lightcurve(self.objid, store=True, verbose=kwargs['verbose'])
        
    def calibrate_baseline(self, ax=None, key='fcqfid', source='ztffp',
                    xmin=-100, xmax=-20, ax_xlim=None, ax_ylim=None):
        assert xmin < xmax
        assert self.t0 > 2400000, 'set t0 first'
        xp = self.t0
        tb = self.lc
        if source is not None:
            tb = tb.query("source==@source")
        fcqfs = np.unique(tb[key].values)
        colors = ['limegreen', 'cyan', 'skyblue', 'purple',
                  'r', 'm', 'pink', 'gold', 'orange', 'y', 'grey']    
        baseline = dict()
        for n, fcqfid in enumerate(fcqfs):            
            ix = tb[key].values==fcqfid            
            thistime = (tb['jdobs'][ix] - 2458000)
            if ax is not None:
                ax.errorbar(thistime, tb['flux'][ix], tb['eflux'][ix], 
                         fmt='.', color=colors[n],
                         label = '%s = %s, Nobs = %d'%(key, fcqfid, np.sum(ix)))
            # baseline            
            __ = np.logical_and(tb['jdobs'][ix]>=xp+xmin, tb['jdobs'][ix]<=xp+xmax)
            if len(tb['jdobs'][ix][__]) > 0:
                xx = tb['jdobs'][ix][__] - 2458000
                yy = tb['flux'][ix][__]
                yye = tb['eflux'][ix][__]
                y0 = np.average(yy, weights=1/yye/yye)
                baseline[fcqfid] = y0
                if ax is not None:
                    ax.errorbar(xx, yy, yye, marker='o', markersize=12,
                                ls='', fillstyle='none', color=colors[n])                
                    ax.axhline(y0, color=colors[n])               
        if ax is not None:
            ax.axvline(xp-2458000, color='k', ls='-')            
            ax.axvline(xp+xmin-2458000, color='r', ls='--')
            ax.axvline(xp+xmax-2458000, color='r', ls='--')        
            ax.grid(ls=":")
            ax.set_xlabel('JD - 2458000 (days)')
            ax.set_ylabel('f')
            ax.legend(loc = 'best', fontsize=10, frameon=False)
            if ax_ylim is not None: ax.set_ylim(ax_ylim)
            if ax_xlim is not None: ax.set_xlim(ax_xlim)                        
        return baseline

    def get_external_spectra(self, filename, epoch, tel='', **kwargs):
        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if not os.path.exists(filename):
            print ('Error: %s not found'%f)
            return None
        _ = []
        for nn,ll in enumerate(open(filename).readlines()):
            if ll[0]=='#' or len(ll)==0:_.append(nn)
        df = pd.read_csv(filename,skiprows=_,delim_whitespace=True)
        assert 'wave' in df.keys() and 'flux' in df.keys()        
        if 'spec' not in self.__dict__:
            self.spec = handle_spectra(self.z, self.mkwebv, t0=self.t0, **kwargs)
        self.spec.add_spectrum(get_numpy(df['wave']),
                               get_numpy(df['flux']), tel, epoch, **kwargs)        
            
    def query_spectra(self, source=None, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert source in ['marshal', 'fritz', None]
        if source in ['marshal',None]:
            marshal.download_spectra(self.objid)            
        if source in ['fritz',None]:
            fritz.download_spectra(self.objid, get_object=True, store=True, verbose=False)   
    
    def get_local_spectra(self, source=None, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert source in ['marshal', 'fritz', None]
        if 'spec' not in self.__dict__:
            self.spec = handle_spectra(self.z, self.mkwebv, t0=self.t0, **kwargs)                      
        
        if source in ['marshal', None]:
            # spec from marshal
            dir_ = marshal.target_spectra_directory(self.objid)
            if os.path.exists( dir_ ):
                for d in os.listdir( dir_ ): 
                    objid,epoch,tel = d.split('_')[0],d.split('_')[1],d.split('_')[2]
                    try: 
                        data, header = fritz.parse_ascii(
                            open( os.path.join(dir_,d) ).read().splitlines() )                        
                    except:
                        print ( os.path.join(dir_,d) )
                        continue
                    self.spec.add_spectrum(get_numpy(data['lbda']),
                            get_numpy(data['flux']), tel, epoch, **kwargs)
                    
        if source in ['fritz', None]:
            # spec from fritz
            _spec = fritz.FritzSpectrum.from_name( self.objid,force_dl=False )
            if _spec is not None:
                if type(_spec) is list:
                    for _ in _spec: 
                        epoch, tel = _.observed_at.split('T')[0].replace('-',''), _.instrument
                        self.spec.add_spectrum(_.lbda, _.flux, tel, epoch, **kwargs)
                else:
                    _ = _spec
                    epoch, tel = _.observed_at.split('T')[0].replace('-',''), _.instrument
                    self.spec.add_spectrum(_.lbda, _.flux, tel, epoch, **kwargs)
        
    def rapid(self, light_curve, run_config):
        ''' run astrorapid
        '''
        print ( 'need to finish astrorapid' )
        return
        photflag, fd = [], True
        for _ in self.lc.index:
            if self.lc['mag'][_] == 99:
                photflag.append(0)
            elif fd:
                photflag.append(6144)
                fd = False
            else:
                photflag.append(4096)
        sys.exit( photflag )
        light_curve_info1 = (self.lc['mjd'], self.lc['flux'], self.lc['eflux'],
                             self.lc['filter'], photflag, ra, dec, self.tranid,
                             self.z, mwebv)
        
        # Classification
        light_curve_list = [light_curve_info1,]
        
        classification = Classify(known_redshift=known_redshift, class_names=self.class_names)
        predictions = classification.get_predictions(light_curve_list, return_objids=False)                
        _pred, _time = predictions                
        if _pred is None:
            self.featurelist['rapid'] = {'success':False,
                                         'cause':'No data survive selection criteria'}
            return
        else:
            self.featurelist['rapid']['prediction'] = []                
            for i,j in zip(_pred[0], _time[0]):
                ii=[]
                for _i in i: ii.append(float(_i))
                self.featurelist['rapid']['prediction'].append({'time': j, 'predict':ii})

    def run_gp(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        
        lc = self.lc
        gpfs = kwargs['gp_bands']
        lc = lc.query('filter in @gpfs')        
        if self.t0 > 0: # cut lc first
            pmin,pmax = min(kwargs['gp_fitr']), max(kwargs['gp_fitr'])
            lc = lc.query('jdobs<=@self.t0+@pmax and jdobs>=@self.t0+@pmin')        
        self.gpcls = fit_gp(np.array(lc['jdobs']), np.array(lc['flux']),
                    np.array(lc['eflux']), filters=np.array(lc['filter']))        
        self.gpcls.train(kernel=kwargs['kernel'], fix_scale=kwargs['fix_scale'],
              gp_mean=kwargs['gp_mean'], opt_routine=kwargs['gp_routine'],
              nwalkers=kwargs['nwalkers'], nsteps=kwargs['nsteps'],
              nsteps_burnin=kwargs['nsteps_burnin'], clobber=kwargs['gp_redo'],
              mcmc_h5_file='gp_%s_%s_%s'%(self.objid,kwargs['gp_routine'],kwargs['gp_mean']),
              verbose=kwargs['verbose'], datadir='%s/cache/'%LOCALSOURCE)
        self.gpcls.predict()
        self.gpcls.set_peak()
        
        # update t0 and fpeak
        if len(kwargs['set_tpeak_method'])==0 and (self.t0 ==0 or len(self.fpeak)==0):
            self.set_peak_gp(kwargs['set_tpeak_filter'])
        elif kwargs['set_tpeak_method'] == 'gp':
            self.set_peak_gp(kwargs['set_tpeak_filter'])
            
    def run_fit(self, enginetype, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        for which in kwargs['fit_methods']:            
            _which = get_pars(which)
            engine = _which['engine']
            enginename = _which['enginename']            
            if enginename != enginetype: continue
            print ('run %s %s' % (enginename, which))
            engine(self, which, enginename, **kwargs)
            
    def set_peak_gp(self, filt):
        assert filt in self.gpcls.tpeak, '%s not available by GP'%filt
        self.t0 = self.gpcls.tpeak[filt]
        self.fpeak = self.gpcls.fpeak
        
    def set_peak_bol_main(self, model_name=None, source_name=None):
        assert 'fitcls' in self.__dict__
        assert 'bol_main' in self.fitcls
        for source in self.fitcls['bol_main']:
            # mbol or mbolbb
            if source_name is not None and source != source_name: continue 
            for model in self.fitcls['bol_main'][source]:
                if model_name is not None and model != model_name: continue
                self.t0 = self.fitcls['bol_main'][source][model].tpeak
                self.fpeak[source] = self.fitcls['bol_main'][source][model].fpeak
                
    def set_peak_multiband_main(self, filt, model_name=None):        
        assert 'fitcls' in self.__dict__
        assert 'multiband_main' in self.fitcls
        for model in self.fitcls['multiband_main']:
            if model_name is not None and model != model_name: continue            
            assert filt in self.fitcls['multiband_main'][model].tpeak,\
                '%s not available by multiband_main.%s'%(filt, model)            
            self.t0 = self.fitcls['multiband_main'][model].tpeak[filt]
            self.fpeak = self.fitcls['multiband_main'][model].fpeak
            
    def set_texp_pl(self, filt, model_name=None):
        assert 'fitcls' in self.__dict__
        assert 'multiband_early' in self.fitcls
        for model in self.fitcls['multiband_early']:
            if model_name is not None and model != model_name: continue            
            assert filt in self.fitcls['multiband_early'][model].filters,\
                '%s not available by multiband_early.%s'%(filt, model)
            t1,t0,t2 = self.fitcls['multiband_early'][model].get_par(filt=filt, parname='texp')
            self.texp = [t1,t0,t2]            

    def set_texp_bol_main(self, model_name=None, source_name=None):
        assert 'fitcls' in self.__dict__
        assert 'bol_main' in self.fitcls
        for source in self.fitcls['bol_main']:
            # mbol or mbolbb
            if source_name is not None and source != source_name: continue 
            for model in self.fitcls['bol_main'][source]:
                if model_name is not None and model != model_name: continue                
                t1,t0,t2 = self.fitcls['bol_main'][source][model].get_par(parname='texp')
                self.texp = [t1,t0,t2]
    
    def set_texp_midway(self):
        assert self.t0 > 2400000, 'set t0 first'
        t0 = self.t0     
        __lc = self.lc.query('mag<99 and jdobs<@t0')
        if len(__lc) > 0:
            expt0 = min(__lc['jdobs'])            
            __lc = self.lc.query('jdobs<@expt0')
            if len(__lc)>0:
                expt1 = max(__lc['jdobs'])            
                expt = (expt1+expt0)/2.

                # convert texp from jd to phase relative to peak
                self.texp = [expt1-t0, expt-t0, expt0-t0]
            else:
                self.texp = [expt0-t0, expt0-t0, expt0-t0]
    
    def _flux_at(self, filt, phase, interpolation=None, fitmodel=0, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert interpolation in [None, 'bin', 'gp', 'fit']
        assert self.t0 > 2400000, 'set t0 first'                
        if interpolation is None or interpolation=='bin':
            # get flux without interpolation
            lc = self.lc.query('filter==@filt')
            if len(lc) > 0:
                __ = np.argmin(abs(lc['jdobs']-self.t0-phase))
                if abs(lc['jdobs'][__]-self.t0-phase)<kwargs['tdbin']:
                    return lc['flux'][__], lc['eflux'][__]
        if (interpolation is None and 'gpcls' in self.__dict__) or interpolation == 'gp':
            assert 'gpcls' in self.__dict__ 
            p_fit = [self.t0 + phase * (1+self.z)]
            xx,yy,yye,ff = self.gpcls.predict(x_pred=p_fit, returnv=True,)
            if len(yy[ np.where(ff == filt) ]) >0:
                return yy[ np.where(ff==filt) ][0], yye[ np.where(ff==filt) ][0]
        if (interpolation is None and 'fitcls' in self.__dict__) or interpolation == 'fit':
            assert 'fitcls' in self.__dict__
            if 'multiband_main' in self.fitcls:
                p_fit = [phase]
                for n, model in enumerate(self.fitcls['multiband_main']):
                    if n != fitmodel: continue
                    xx,yy,yy1,yy2,ff = self.fitcls['multiband_main'][model].predict(x_pred=p_fit, returnv=True)
                    if len(yy[ np.where(ff == filt) ])>0:
                        f = yy[ np.where(ff == filt) ][0]
                        f1 = yy1[ np.where(ff == filt) ][0]
                        f2 = yy2[ np.where(ff == filt) ][0]
                        return f, max(abs(f-f1), abs(f-f2))
        return None, None
    
    def _flux_at_list(self, filt, phaselist, interpolation=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        fluxlist, efluxlist = [], []
        for phase in phaselist:
            flux, eflux = self._flux_at(filt, phase, interpolation=interpolation, **kwargs)
            fluxlist.append(flux)
            efluxlist.append(eflux)
        return fluxlist, efluxlist
    
    def _mag_at(self, filt, phase, interpolation=None,
                corr_mkw=False, corr_host=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        flux, eflux = self._flux_at(filt, phase, interpolation=interpolation, **kwargs)
        if flux is None: return None, None
        ebv = 0
        if corr_mkw: ebv += self.mkwebv
        if corr_host: ebv += self.hostebv
        m, me, limmag = self.flux_to_mag([flux], dflux=[eflux], sigma=kwargs['snrt'], zp=23.9)
        if m < 99:
            return m[0] - ebv * Rf[filt], me[0]
        else:
            m = self.flux_to_mag([flux], dflux=None, sigma=kwargs['snrt'], zp=23.9)
            return m[0] - ebv * Rf[filt], None

    def _mag_at_list(self, filt, phaselist, interpolation=None,
                    corr_mkw=False, corr_host=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        maglist, emaglist = [], []
        for phase in phaselist:
            mag, emag = self._mag_at(filt, phase, interpolation=interpolation,
                            corr_mkw=corr_mkw, corr_host=corr_host, **kwargs)
            maglist.append(mag)
            emaglist.append(emag)
        return maglist, emaglist
        
    def _absmag_at(self, filt, phase, interpolation=None,
                corr_mkw=False, corr_host=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        mag, emag = self._mag_at(filt, phase, interpolation=interpolation,
                                 corr_mkw=corr_mkw, corr_host=corr_host, **kwargs)
        if mag is None: return None, None
        if emag is None: return mag-self.dm, None
        else: return mag-self.dm, np.sqrt(emag**2 + self.dm_error(filt)**2)

    def _absmag_at_list(self, filt, phaselist, interpolation=None,
                corr_mkw=False, corr_host=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        maglist, emaglist = [], []
        for phase in phaselist:
            mag, emag = self._absmag_at(filt, phase, interpolation=interpolation,
                            corr_mkw=corr_mkw, corr_host=corr_host, **kwargs)
            maglist.append(mag)
            emaglist.append(emag)
        return maglist, emaglist
    
    def _color_at(self, filt1, filt2, phase,
                  interpolation=None, corr_mkw=False, corr_host=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        mag1, emag1 = self._mag_at(filt1, phase, interpolation=interpolation,
                                   corr_mkw=corr_mkw, corr_host=corr_host, **kwargs)
        mag2, emag2 = self._mag_at(filt2, phase, interpolation=interpolation,
                                   corr_mkw=corr_mkw, corr_host=corr_host, **kwargs)
        if mag1 is None or mag2 is None: return None, None
        if emag1 is None or emag2 is None: return mag1-mag2, None
        else: return mag1-mag2, np.sqrt(emag1**2+emag2**2)
    
    def _rate_at(self, filt, phase1, phase2, interpolation='gp', **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        mag1, emag1 = self._mag_at(filt, phase1, interpolation=interpolation, **kwargs)
        mag2, emag2 = self._mag_at(filt, phase2, interpolation=interpolation, **kwargs)
        if mag1 is None or mag2 is None: return None, None
        if emag1 is None or emag2 is None: return mag1-mag2, None
        else: return mag1-mag2, np.sqrt(emag1**2+emag2**2)
    
    def est_hostebv_with_c10(self, interpolation='gp', **kwargs):
        '''
        estimate host ebv with color comparison approach
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        filt1, filt2 = kwargs['hostebv_bands']
        c10_temp = self.read_c10()
        sntype = self.sntype.replace('SN ','').replace('AT ','').strip()
        if not (filt1, filt2) in c10_temp:
            print ('Error: %s-%s not defined in c10 template file, check!!'%
                   (filt1, filt2))
            return
        if not sntype in c10_temp[(filt1, filt2)]:
            print ('Error: %s not defined for %s-%s in c10 template file, check!!'%
                   (self.sntype, filt1, filt2))
            return
        if '%s-%s'%(filt1,filt2) not in toebv:
            print ('Error: E(%s-%s) to E(B-V) factor is not defined, check!!'%
                   (filt1, filt2))
            return
        c10, c10e = c10_temp[(filt1, filt2)][sntype].split(',')
        c10, c10e = float(c10), float(c10e)        
        _c10, _c10e = self._color_at(filt1, filt2, 10,
                                     interpolation=interpolation, corr_mkw=True)
        #if _c10 < c10-c10e:  ehost = 0
        #elif _c10 > c10-c10e and _c10 < c10+c10e:  ehost = 0
        #else:
        #    randomc = c10-c10e+np.random.random(1)[0]*c10e*2
        #    ehost = _c10 - randomc
        if _c10 < c10:  ehost = 0
        else: ehost = _c10 - c10        
        # from E(X-Y) to E(B-V)
        self.hostebv = ehost/toebv['%s-%s'%(filt1,filt2)]

    def calc_colors(self, **kwargs):
        ''' define colors
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        assert self.t0 > 2400000, 'set t0 first'
          
        filt1, filt2 = kwargs['color_bands']        
        self.colors = dict()
        if 1 in kwargs['color_interp']:
            # bin and match g and r
            self.lc = self.bin_df(self.lc, deltah = kwargs['tdbin'] )            
            self.match_colors(**kwargs)            
            _lc = self.lc_match.query('mag_%s<99 and mag_%s<99' % (filt1, filt2))
            if len(_lc) > 0:
                self.colors[1] = [_lc['jdobs_%s'%filt1], _lc['mag_%s'%filt1], _lc['mag_%s'%filt2],
                                  _lc['emag_%s'%filt1], _lc['emag_%s'%filt2]]
        
        # or with GP/fit intepolation
        jd2,m2,mm2,me2,mme2 = [],[],[],[],[]
        jd3,m3,mm3,me3,mme3 = [],[],[],[],[]
        for _filt, _filt1 in zip([filt1, filt2], [filt2, filt1]):            
            _lc = self.lc.query('filter==@_filt and mag<99')            
            if 'fitcls' in self.__dict__ and 2 in kwargs['color_interp']: # analytic sn lc fit
                if 'multiband_main' in self.fitcls:                    
                    for model in self.fitcls['multiband_main']:
                        __lc = _lc # donot disturb GP
                        if 'xpredict' in self.fitcls['multiband_main'][model].__dict__:
                            p1,p2 = self.fitcls['multiband_main'][model].xpredict
                            __lc = _lc.query('jdobs>@self.t0+@p1 and jdobs<@self.t0+@p2')
                        p_fit = (__lc['jdobs']-self.t0) / (1+self.z)
                        xx,yy,yy1,yy2,ff = self.fitcls['multiband_main'][model].predict(x_pred=p_fit,
                                                                returnv=True, quant=kwargs['quantile'])
                        __ = np.where(ff==_filt1)
                        xx = xx[__]
                        yy = yy[__]
                        yy1 = yy1[__]
                        yy2 = yy2[__]
                        mm, mme, lm = self.flux_to_mag(yy, dflux=abs(yy1-yy2)/2, sigma=kwargs['snrt'], zp=23.9)
                        __ = np.where(mm<99)
                        if _filt == filt1:  
                            jd2 = np.append(jd2, get_numpy(__lc['jdobs'])[__])
                            m2 = np.append(m2, get_numpy(__lc['mag'])[__])
                            mm2 = np.append(mm2, mm[__])
                            me2 = np.append(me2, get_numpy(__lc['emag'])[__])
                            mme2 = np.append(mme2, mme[__])
                        else: 
                            jd2 = np.append(jd2, get_numpy(__lc['jdobs'])[__])
                            mm2 = np.append(mm2, get_numpy(__lc['mag'])[__])
                            m2 = np.append(m2, mm[__])
                            mme2 = np.append(mme2, get_numpy(__lc['emag'])[__])
                            me2 = np.append(me2, mme[__])
                        break                    
            if 'gpcls' in self.__dict__ and 3 in kwargs['color_interp']: # GP intepolation
                __lc = _lc
                p_fit = __lc['jdobs']
                xx,yy,yye,ff = self.gpcls.predict(x_pred=p_fit, returnv=True,)
                __ = np.where(ff==_filt1)
                xx = xx[__]
                yy = yy[__]
                yye = yye[__]
                mm, mme, lm = self.flux_to_mag(yy, dflux=yye, sigma=kwargs['snrt'], zp=23.9)
                __ = np.where(mm<99)
                if _filt == filt1:
                    jd3 = np.append(jd3, get_numpy(__lc['jdobs'])[__])
                    m3 = np.append(m3, get_numpy(__lc['mag'])[__])
                    mm3 = np.append(mm3, mm[__])
                    me3 = np.append(me3, get_numpy(__lc['emag'])[__])
                    mme3 = np.append(mme3, mme[__])
                else:
                    jd3 = np.append(jd3, get_numpy(__lc['jdobs'])[__])
                    mm3 = np.append(mm3, get_numpy(__lc['mag'])[__])
                    m3 = np.append(m3, mm[__])
                    mme3 = np.append(mme3, get_numpy(__lc['emag'])[__])
                    me3 = np.append(me3, mme[__])
        if len(jd2)>0: self.colors[2] = (jd2,m2,mm2,me2,mme2)
        if len(jd3)>0: self.colors[3] = (jd3,m3,mm3,me3,mme3)
            
    def lyman_bol(self, **kwargs):
        ''' 
        calculate bolometric LC with Lyman bolometric correction
        ''' 
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if not 'colors' in self.__dict__:
            if kwargs['verbose']: print ('calc colors first')
            return
        if self.sntype is None:
            if kwargs['verbose']: print ('define sn type first')
            return        
        filt1, filt2 = kwargs['color_bands']
        sntype = self.sntype.replace('SN ','').replace('AT ','').strip()
        self.mbol = dict()
        ebv = self.mkwebv + self.hostebv
        for _ in self.colors:
            jd, m1, m2, m1e, m2e = self.colors[_]
            bc, __ = BC_Lyman(m1-m2, colortype='%s-%s'%(filt1,filt2), phase='normal', sntype=sntype)
            if bc is None:
                print ('Error: check lyman BC for sntype %s, color %s-%s for normal phase'%(sntype, filt1, filt2))
                return
            mbol, embol = bc + m1 -self.dm, np.sqrt(m1e**2 + self.dm_error(filt1)**2)            
            Lbol = Mbol_to_Lbol(mbol)                        
            eLbol = abs(Mbol_to_Lbol(mbol+embol)-Lbol)
            _ck=np.isfinite(Lbol)            
            self.mbol[_] = (jd[_ck], Lbol[_ck], eLbol[_ck])

    def bb_colors(self, **kwargs):
        ''' define epochs for BB
        !!! haven't included any ebv value
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])       
        
        filters = kwargs['bb_bands']
        assert len(filters) >= 3, 'at least 3 bands needed for blackbody'
        
        self.cbb = dict()
        self.cbb[1] = dict()                
        if 1 in kwargs['bb_interp']:
            # bin and match g and r
            self.lc = self.bin_df(self.lc, deltah = kwargs['tdbin'] )
            self.match_colors(**kwargs)        
            
            _cond = ''
            for n, filt in enumerate(filters):
                assert len(self.lc_match.query('mag_%s<99' % filt)) > 0, \
                    '%s band not available, set filters correctly!' % filt
                _cond += 'mag_%s<99' %filt
                if n < len(filters)-1: _cond += ' and '                        
            _lc = self.lc_match.query(_cond)
            if len(_lc) > 0:
                self.cbb[1]['t'] = _lc['jdobs_%s'%filters[0]]
                for n, filt in enumerate(filters):
                    self.cbb[1][filt] = ( _lc['mag_%s'%filt], _lc['emag_%s'%filt] )
                    
        # or with GP/fit intepolation
        _jds = self.lc.query('mag<99')['jdobs']
        self.cbb[2] = dict()        
        self.cbb[3] = dict()
        for _filt in filters:
            if 'fitcls' in self.__dict__ and 2 in kwargs['bb_interp']: # analytic sn lc fit
                if 'multiband_main' in self.fitcls:
                    p_fit = (_jds-self.t0) / (1+self.z)
                    for model in self.fitcls['multiband_main']:                    
                        xx,yy,yy1,yy2,ff = self.fitcls['multiband_main'][model].predict(x_pred=p_fit,
                                                                returnv=True, quant=kwargs['quantile'])
                        __ = np.where(ff==_filt)
                        xx = xx[__]
                        yy = yy[__]
                        yy1 = yy1[__]
                        yy2 = yy2[__]
                        mm = self.flux_to_mag(yy, dflux=None, sigma=kwargs['snrt'], zp=23.9)
                        _mm1 = self.flux_to_mag(yy1, dflux=None, sigma=kwargs['snrt'], zp=23.9)
                        _mm2 = self.flux_to_mag(yy2, dflux=None, sigma=kwargs['snrt'], zp=23.9)
                        mme = (abs(mm-_mm1) + abs(mm-_mm2))/2                
                        self.cbb[2][_filt] = (mm,mme)
                        break
            if 'gpcls' in self.__dict__ and 3 in kwargs['bb_interp']: # GP intepolation
                p_fit = _jds
                xx,yy,yye,ff = self.gpcls.predict(x_pred=p_fit, returnv=True,)
                __ = np.where(ff==_filt)
                xx = xx[__]
                yy = yy[__]
                yye = yye[__]
                mm = self.flux_to_mag(yy, dflux=None, sigma=kwargs['snrt'], zp=23.9)
                mme = abs(self.flux_to_mag(yy, dflux=None, sigma=kwargs['snrt'], zp=23.9)
                          -self.flux_to_mag(yy+yye, dflux=None, sigma=kwargs['snrt'], zp=23.9))
                self.cbb[3][_filt] = (mm,mme)
        if len(self.cbb[2])>0: self.cbb[2]['t'] = _jds
        if len(self.cbb[3])>0: self.cbb[3]['t'] = _jds
        
    def bb_bol(self, do_Kcorr=True, ab2vega=True, **kwargs):
        ''' 
        calculate bolometric LC with diluted blackbody
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if not 'cbb' in self.__dict__:
            if kwargs['verbose']: print ('calc BB colors first')
            return        
        filters = kwargs['bb_bands']
        self.mbolbb = dict()
        ebv = self.mkwebv + self.hostebv
        SN_distance = 10**((self.dm-25)/5)*3.086e24 # in cm
        for _ in self.cbb:
            jdl, Ll, eLl, Tl, eTl, Rl, eRl, L1l, eL1l = [], [], [], [], [], [], [], [], []
            if len(self.cbb[_])==0: continue
            for __ in range(len(self.cbb[_][filters[0]][0])):
                _jd = get_numpy(self.cbb[_]['t'])[__] 
                # black body
                fluxes, efluxes, ws, bs = [], [], [], []
                for _filt in filters:
                    m, em = self.cbb[_][_filt]
                    if len(m) == 0:continue
                    _mag = get_numpy(m)[__] - ebv * Rf[_filt]
                    _mage = get_numpy(em)[__]                    
                    # If UVOT bands are in AB, need to convert to Vega
                    if _filt in ab_to_vega.keys() and ab2vega:
                        _mag -= ab_to_vega[_filt]
                    wlref = central_wavelengths[_filt]
                    fref = zp[_filt] * 1e-11
                    bandwidths = filter_width[_filt]
                    if do_Kcorr:
                        wlref /= (1+self.z)
                        fref *= (1+self.z)
                        bandwidths /= (1+self.z)
                    _flux = 4*np.pi*SN_distance**2*fref*10**(-0.4*_mag)
                    _dflux = 2.5/np.log(10) * _flux * _mage
                    fluxes.append(_flux)
                    efluxes.append(_dflux)
                    ws.append(wlref)
                    bs.append(bandwidths)                                
                fluxes, efluxes, ws, bs = np.array(fluxes), np.array(efluxes), np.array(ws), np.array(bs)
                
                # Set flux to zero at red and blue extrema matching wlref1
                #flux1 = np.insert(fluxes,0,0)
                #flux1 = np.append(flux1,0)
                
                # if any band cannot provide fluxes, just skip them
                __ = np.isnan(fluxes)
                if len(fluxes[__]) > 0: continue
                
                BBparams, covar = curve_fit( bbody, ws, fluxes/1e40, p0=(10,2),
                        sigma=efluxes/1e40, absolute_sigma=True, maxfev=kwargs['maxfev'])
                _T = BBparams[0]
                _Te = min(np.sqrt(np.diag(covar))[0],_T)
                _R = np.abs(BBparams[1])
                _Re = min(np.sqrt(np.diag(covar))[1],_R)
                Tl.append(_T)
                eTl.append(_Te)
                Rl.append(_R)
                eRl.append(_Re)                                                               
                jdl.append(_jd)

                # Get pseudobolometric luminosity by trapezoidal integration,
                #    with flux set to zero outside of observed bands
                _L = integrate.trapz(fluxes[np.argsort(ws)],ws[np.argsort(ws)])
                Ll.append(_L)
                
                # Use flux errors and bandwidths to get luminosity error
                _Le = np.sqrt(np.sum((bs*efluxes)**2))
                eLl.append(_Le)
                
                # Calculate luminosity using alternative method of Stefan-Boltzmann, and T and R from fit
                _Lbb = 4*np.pi*(_R*1e15)**2*5.67e-5*(_T*1e3)**4
                _Lbbe = _Lbb*np.sqrt((2*_Re/_R)**2+(4*_Te/_T)**2)
                L1l.append(_Lbb)
                eL1l.append(_Lbbe)                
            self.mbolbb[_] = ( jdl, Ll, eLl, Tl, eTl, Rl, eRl, L1l, eL1l )
                
    def match_colors(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])            
        self.combine_multi_obs(**kwargs)
        _lc = self.lc_match
        filters = np.unique(_lc['filter'])
        _keys = []
        for f in sorted(filters):  _keys = np.append(_keys, ['%s_%s'%(_,f) for _ in _lc.columns])
        __arr = []
        for _jd in np.unique(_lc['jdbin']):
            _lc0 = _lc.query('jdbin==@_jd')
            __arr0 = []
            for _filt in sorted(filters): 
                _lc0f = _lc0.query('filter==@_filt')
                if len(_lc0f) == 0:
                    __arr0 = np.append( __arr0, np.zeros(len(_lc.columns))+99 )
                else:
                    __arr0 = np.append( __arr0, _lc0f.values[0] )
            __arr.append (__arr0.tolist())
        _df = pd.DataFrame(__arr, columns=_keys) 
        for _ in _df.keys(): 
            try: _df[_] = pd.to_numeric(_df[_])
            except: pass
        self.lc_match = _df
        
    def combine_multi_obs(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        if not 'jdbin' in self.lc:
            _lc = self.bin_df(self.lc, deltah = kwargs['tdbin'])
        else:
            _lc = self.lc
        if True: _lc = _lc.query('mag<99')
        __arr = []
        for _jd in np.unique(_lc['jdbin']):
            _lc0 = _lc.query('jdbin==@_jd')
            unique_f, counts = np.unique(_lc0['filter'],return_counts=True)
            single_f, double_f = list(unique_f[ (counts==1)]), \
                list(unique_f[ (counts>1)]) 
            for f in single_f:
                _lc0f = _lc0.query('filter==@f')
                __arr.append( _lc0f.values[0].tolist() )
            
            if not len(double_f)==0:   
                for f in double_f:
                    _lc0f = _lc0.query('filter==@f')
                    __arr.append( self.merge_df_cols(_lc0f) )
        _df = pd.DataFrame(__arr, columns=self.lc.columns) 
        for _ in _df.keys(): 
            try: _df[_] = pd.to_numeric(_df[_])
            except: pass
        self.lc_match = _df
        
    def ext_lc(self,filt='g',num=3,**kwargs):
        ''' add LC points based on intepolations
        which can be used for instance the power law fits when no enough early data exists.
        '''
        assert self.t0 > 2400000, '!!!either input jdpeak or do GP first'
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        
        assert 'gpcls' in self.__dict__            
        rel_flux_cutoff = float(kwargs['rel_flux_cutoff'])
        rel_flux_cutfrom = float(kwargs['rel_flux_cutfrom'])
        tm = self.fpeak['GP'][filt][0]
        y1 = tm * rel_flux_cutoff
        y2 = tm * rel_flux_cutfrom
            
        __ = np.logical_and(self.gpcls.f_pred==filt, self.gpcls.x_pred<self.t0)
        xx = self.gpcls.x_pred[__]
        yy = self.gpcls.y_pred[__]
        yye = self.gpcls.y_prede[__]
        
        __ = np.logical_and(yy>=y2, yy<=y1)
        _xx, _yy, _yye = xx[__], yy[__], yye[__]
            
        num = min(num, len(_xx)-1)
        _id = np.arange(0,len(_xx),1)
        _rid = np.random.choice(_id, num, replace=False)
        
        return _xx[_rid], _yy[_rid], _yye[_rid]               
        
    def clip_lc(self, **kwargs):
        ''' Removes outlier data points using sigma-clipping '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        # remove negative flux epoch
        self.lc = self.lc.query('flux>0')
        
        # remove outliers
        if kwargs['clipsigma'] is not None:
            clipsigma = float(kwargs['clipsigma'])
            __arr = None            
            for f in np.unique(self.lc['filter']):                
                _lc = self.lc.query('filter==@f')
                outlier_mask = sigma_clip(
                    data=_lc['flux'],
                    sigma=clipsigma,
                    sigma_lower=None,
                    sigma_upper=None,
                    cenfunc='median',
                    stdfunc='std',
                ).mask
                _lc = self.clip_df(_lc, outlier_mask)
                if __arr is None: __arr = _lc
                else: __arr = __arr.append( _lc, ignore_index=True)                
            self.lc = __arr                    

    def _nepochs(self, **kwargs):
        # how many epochs of photometry in either band.
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        dets = self.lc.query('mag<99')
        ndets = dict()
        if kwargs['plot_bands'] is not None: plot_bands = kwargs['plot_bands']
        else: plot_bands = np.unique(self.lc['filter'])
        for filt in plot_bands:
            ndets[filt] = len(dets.query('filter==@filt'))
        return ndets
    
    def _ncolors(self, **kwargs):
        # how many color epochs (sampled within $\pm3$ days)
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        ncolors = dict()
        if kwargs['plot_bands'] is not None: plot_bands = kwargs['plot_bands']
        else: plot_bands = np.unique(self.lc['filter'])
        for filt in plot_bands:
            for filt1 in plot_bands:
                if filt == filt1: continue
                if len(dets.query('filter==@filt'))==0 or len(dets.query('filter==@filt1'))==0:continue
                nc = 0
                for _jd in dets.query('filter==@filt')['jdobs']:
                    if min(abs(dets.query('filter==@filt1')['jdobs'] - _jd)) <= kwargs['tdbin']:
                        nc += 1
                ncolors['%s %s'%(filt, filt1)] = nc
        return ncolors

    def _peak_accuracy(self, within=3, **kwargs):
        # Photometry available both before and after peak so that a peak can be determined?
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        peakphot = dict()        
        t0 = self.t0
        if kwargs['plot_bands'] is not None: plot_bands = kwargs['plot_bands']
        else: plot_bands = np.unique(self.lc['filter'])
        for filt in plot_bands:
            filt_lc = self.lc.query('mag<99 and filter==@filt and jdobs>@t0-@within and jdobs<@t0+@within')
            peakphot[filt] = len(filt_lc)
        return peakphot      
            
    def summarize_plot(self, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])         
        self.fig = plt.figure(num=1, clear=True, dpi=400,
                            tight_layout=False, constrained_layout=False,)                
        gs1 = self.fig.add_gridspec(nrows=3, ncols=2, left=.05,
                                        right=.95, wspace=0.05, hspace=0.3)
        self.ax = self.fig.add_subplot(gs1[0, 0]) # flux
        self.ax1 = self.fig.add_subplot(gs1[:2, 1]) # spectra
        self.ax2 = self.fig.add_subplot(gs1[1, 0]) # mag        
        self.ax3 = self.fig.add_subplot(gs1[2, 1]) # color
        self.ax4 = self.fig.add_subplot(gs1[2, 0]) # luminosity
        # plot them
        self._ax()
        self._ax1()
        self._ax2()
        self._ax3()
        self._ax4()
        
    def _ax(self, show_title=True, show_legend=True, ylabel_2right=False,
            x0=2458000, source=None, cond=None, **kwargs):  #flux LC plot                
        if self.ax is None:
            print ('Error: no ax defined, skipped flux plots')
            return
        if not 'lc' in self.__dict__:
            print ('Error: no lc defined, parse lc data first')
            return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])       
        fscale=float(kwargs['flux_scale'])        
        if kwargs['plot_bands'] is not None: plot_bands = kwargs['plot_bands']
        else: plot_bands = np.unique(self.lc['filter'])
        for filt in plot_bands:
            ''' for each filter '''        
            
            # normalize flux item in ax
            if source is None:
                __lc = self.lc.query('filter==@filt')
            else:
                __lc = self.lc.query('filter==@filt and source==@source')
            # peak flux
            fm, efm = self._flux_at(filt, 0, interpolation=None, **kwargs)            
            if fm is None: fm = __lc['flux'].max()
            self.ax.errorbar(__lc['jdobs']-x0, __lc['flux']/fm*fscale+ys[filt],
                             yerr=__lc['eflux']/fm*fscale, alpha=kwargs['alphabest'], **PROP1f[filt])
            
            # analytic sn lc fit  
            if 'fitcls' in self.__dict__:
                if 'multiband_early' in self.fitcls:
                    for model in self.fitcls['multiband_early']:
                        _model = self.fitcls['multiband_early'][model]
                        
                        # fit range
                        pmin,pmax = kwargs['multiband_early_xrangep']
                        p_fit = np.arange(pmin, pmax, .1)
                        
                        # samples
                        x,y,f=_model.predict_random(limit=kwargs['plot_mcmct'],
                                                    plotnsamples=kwargs['plot_nsamples'])
                        for xx, yy in zip(x[ np.where(f == filt) ], y[ np.where(f == filt) ]):
                            self.ax.plot(xx-x0, yy+ys[filt],
                                         alpha=kwargs['alphasample'], **PROP4[filt])
                            
                if 'multiband_main' in self.fitcls:
                    for model in self.fitcls['multiband_main']:
                        _model = self.fitcls['multiband_main'][model]
                        
                        # fit range
                        pmin,pmax = kwargs['multiband_main_xrangep']
                        p_fit = np.arange(pmin, pmax, .1)
                            
                        # samples
                        x,y,f=_model.predict_random(limit=kwargs['plot_mcmct'],
                                                    plotnsamples=kwargs['plot_nsamples'])
                        for xx, yy in zip(x[ np.where(f == filt) ], y[ np.where(f == filt) ]):
                            self.ax.plot(xx-x0, yy/fm*fscale+ys[filt],
                                         alpha=kwargs['alphasample'], **PROP2[filt])

            # Gaussian process
            if 'gpcls' in self.__dict__:
                if filt in self.gpcls.f_pred:                                        
                    __ = np.where(self.gpcls.f_pred==filt)
                    xx = self.gpcls.x_pred[__]
                    yy = self.gpcls.y_pred[__]
                    yye = self.gpcls.y_prede[__]
                    
                    if kwargs['gp_plotr'] is not None:
                        pmin,pmax = min(kwargs['gp_plotr']), max(kwargs['gp_plotr'])                    
                        __ = np.logical_and(xx>self.t0+pmin, xx<self.t0+pmax)                    
                        xx = xx[__]
                        yy = yy[__]
                        yye = yye[__] 
                   
                    self.ax.plot(xx-x0, yy/fm*fscale+ys[filt],
                                 alpha=kwargs['alphabest'], **PROP3[filt])                    
                    self.ax.fill_between(xx-x0, (yy-yye) /fm*fscale+ys[filt],
                            (yy+yye) /fm*fscale+ys[filt], alpha=kwargs['alphasample'], **PROP3[filt])
                    
        #  subplot settings            
        if kwargs['ax_xlim'] is not None:
            x1,x2 = min(kwargs['ax_xlim']), max(kwargs['ax_xlim'])
            self.ax.set_xlim([x1*(1+self.z)+self.t0-x0, x2*(1+self.z)+self.t0-x0])
        if kwargs['ax_ylim'] is not None:
            self.ax.set_ylim(kwargs['ax_ylim'])
        if show_title:
            self.ax.set_title('%s\n%s z=%.2f'%
                                  (self.objid,self.sntype,self.z),fontsize=12)                        
        if self.texp is not None:
            self.ax.axvline(self.texp[1]*(1+self.z)+self.t0-x0, ls='--', label='texp', color='k')
            self.ax.axvline(self.texp[0]*(1+self.z)+self.t0-x0, ls=':', color='k')
            self.ax.axvline(self.texp[2]*(1+self.z)+self.t0-x0, ls=':', color='k')               
        self.ax.set_xlabel('JD - %s (d)'%x0,fontsize=12)
        if ylabel_2right:
            ax=self.ax.twinx()
            ax.set_ylabel('Flux + $\mathrm{offset}$',fontsize=12)
            ax.set_yticks([])
        else:
            self.ax.set_ylabel('Flux + $\mathrm{offset}$',fontsize=12)
        if show_legend: self.ax.legend(fontsize=10, frameon=False)
        
    def _ax1(self, show_title=False, show_legend=False,
             source=None, **kwargs):  #spectra plot
        if self.ax1 is None:
            print ('Error: no ax1 defined, skipped spectra plots')
            return
        if not 'spec' in self.__dict__:
            print ('Error: no spec defined, parse spectral data first')
            return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])

        cw, region = None, None
        if 'cw' in self.__dict__: cw, region = self.cw, self.region
        for _ in self.spec.data:
            spec   = self.spec.data[_]['data']
            phase  = float(self.spec.data[_]['phase'])
            #if phase < kwargs['specfit_phase'][0] or phase > kwargs['specfit_phase'][1]:
            #    continue 
            source = '_'.join(_.split())
            ys     = self.spec.data[_]['ys']
            
            wave, flux = spec._norm_spectrum(region=region, stype='flat', **kwargs)             
            self.ax1.step(wave, flux + ys, **PROP10['data'])
            if region is not None:
                self.ax1.text(max(region), ys, '%s' % source, color='k')
            else:
                self.ax1.text(max(wave), ys, '%s' % source, color='k')
                
            # peaks
            if 'specpeaks' in self.__dict__:
                if source in self.specpeaks:
                    findp = self.specpeaks[source]                    
                    if 1 in findp:
                        for _w in findp[1]:
                            __ = np.argmin(abs(wave-_w))
                            self.ax1.plot(_w, flux[__] + ys, color='orange', marker='|', markersize=20)
                    if -1 in findp:
                        for _w in findp[-1]:
                            __ = np.argmin(abs(wave-_w))
                            self.ax1.plot(_w, flux[__] + ys, color='cyan', marker='|', markersize=20)
                            
            # fits            
            if 'fitcls' in self.__dict__:
                if 'specline' in self.fitcls:
                    if source in self.fitcls['specline']:                        
                        for model in self.fitcls['specline'][source]:
                            _model = self.fitcls['specline'][source][model]
                            if _model is None: continue
                            x,y,f=_model.predict_random(limit=kwargs['plot_mcmct'],
                                                        plotnsamples=kwargs['plot_nsamples'])
                            for xx, yy in zip(x,y):                                
                                self.ax1.plot(xx, yy+ys, **PROP10['bindata'])
        if cw is not None: self.ax1.axvline(cw, color='k', ls='--')
        
        #  subplot settings        
        if region is not None:                           
            _ticks = np.arange(min(self.region),max(self.region)+100,100)
        else:
            _ticks = np.arange(3000, 9000, 1000)        
        self.ax1.set_xticks( _ticks )
        self.ax1.set_xticklabels(['%.1f'%_ for _ in _ticks], rotation=30)
        if cw is not None:
            ax1 = self.ax1.twiny()
            ax1.set_xlim( self.region )
            _tickl = ['%.1f'%_ for _ in (_ticks/cw-1)*299792.458/1000.]
            ax1.set_xticks(_ticks)
            ax1.set_xticklabels(_tickl, rotation=30)
            ax1.set_xlabel('velocity ($10^{3}~km~s^{-1}$)', fontsize=12)
        self.ax1.set_xlabel('Rest Wavelength ($\AA$)', fontsize=12)        
        self.ax1.set_yticks([])
        self.ax1.tick_params("both", direction='in',
                             right=True, labelright=True, left=False, labelleft=False)
        if region is not None:
            self.ax1.set_xlim([min(self.region),max(self.region)])            
        if show_title:
            self.ax1.set_title('%s\n%s z=%.2f'%
                               (self.objid,self.sntype,self.z),fontsize=12)                    
        if show_legend: self.ax1.legend(fontsize=10, frameon=False)
        
    def _ax2(self, show_title=False, show_legend=False, ylabel_2right=False,
             corr_mkw=False, corr_host=False, source=None, **kwargs):  #mag LC plot
        if self.ax2 is None:
            print ('Error: no ax2 defined, skipped mag plots')
            return
        if not 'lc' in self.__dict__:
            print ('Error: no lc defined, parse lc data first')
            return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])       
        ebv = 0
        if corr_mkw: ebv += self.mkwebv
        if corr_host: ebv += self.hostebv
        if kwargs['plot_bands'] is not None: plot_bands = kwargs['plot_bands']
        else: plot_bands = np.unique(self.lc['filter'])
        for filt in plot_bands:
            ''' for each filter '''        
            
            # mag item in ax
            # set t to rest frame phase relative to peak 
            if source is None:
                __lc = self.lc.query('filter==@filt and mag<99')
            else:
                __lc = self.lc.query('filter==@filt and mag<99 and source==@source')
            self.ax2.errorbar((__lc['jdobs']-self.t0)/(1+self.z), __lc['mag']-self.dm-Rf[filt]*ebv, 
                            yerr=np.sqrt(__lc['emag']**2+self.dm_error(filt)**2),
                            alpha=kwargs['alphabest'], **PROP1[filt])
            # upper limits
            if source is None:
                __lc = self.lc.query('filter==@filt and mag==99 and jdobs<@self.t0')
            else:
                __lc = self.lc.query('filter==@filt and mag==99 and jdobs<@self.t0 and source==@source')
            if len(__lc) > 0:
                self.ax2.plot((__lc['jdobs']-self.t0)/(1+self.z), __lc['limmag']-self.dm-Rf[filt]*ebv,
                              alpha=kwargs['alphasample'], **PROP1l[filt])
                
            # analytic sn lc fit            
            if 'fitcls' in self.__dict__:
                if 'multiband_early' in self.fitcls:
                    for model in self.fitcls['multiband_early']:
                        _model = self.fitcls['multiband_early'][model]
                        
                        # fit range
                        pmin,pmax = kwargs['multiband_early_xrangep']
                        if self.texp is not None: pmin = max(min(self.texp), pmin)
                        p_fit = np.arange(pmin, pmax, .1)
                        
                        # samples
                        x,y,f=_model.predict_random(limit=kwargs['plot_mcmct'],
                                                    plotnsamples=kwargs['plot_nsamples'])
                        for xx, yy in zip(x[ np.where(f == filt) ], y[ np.where(f == filt) ]):
                            m = self.flux_to_mag(yy, dflux=None, sigma=kwargs['snrt'], zp=23.9)
                            self.ax2.plot((xx-self.t0)/(1+self.z), m-self.dm-Rf[filt]*ebv,
                                          alpha=kwargs['alphasample'], **PROP2[filt])
                
                if 'multiband_main' in self.fitcls:
                    for model in self.fitcls['multiband_main']:
                        _model = self.fitcls['multiband_main'][model]
                        
                        # fit range
                        pmin,pmax = kwargs['multiband_main_xrangep']
                        if self.texp is not None: pmin = max(min(self.texp), pmin)
                        p_fit = np.arange(pmin, pmax, .1)
                            
                        # samples
                        x,y,f=_model.predict_random(limit=kwargs['plot_mcmct'],
                                                    plotnsamples=kwargs['plot_nsamples'])
                        for xx, yy in zip(x[ np.where(f == filt) ], y[ np.where(f == filt) ]):
                            m = self.flux_to_mag(yy, dflux=None, sigma=kwargs['snrt'], zp=23.9)
                            self.ax2.plot((xx-self.t0)/(1+self.z), m-self.dm-Rf[filt]*ebv,
                                         alpha=kwargs['alphasample'], **PROP2[filt])
                                                        
            # Gaussian process
            if 'gpcls' in self.__dict__:
                if filt in self.gpcls.f_pred:                                        
                    __ = np.where(self.gpcls.f_pred==filt)
                    xx = self.gpcls.x_pred[__]
                    yy = self.gpcls.y_pred[__]
                    yye = self.gpcls.y_prede[__]
                    
                    if kwargs['gp_plotr'] is not None:
                        pmin,pmax = min(kwargs['gp_plotr']), max(kwargs['gp_plotr'])
                        if self.texp is not None: pmin = max(min(self.texp), pmin)
                        __ = np.logical_and(xx>self.t0+pmin, xx<self.t0+pmax)                    
                        xx = xx[__]
                        yy = yy[__]
                        yye = yye[__]                     
                    m,me,limmag = self.flux_to_mag(yy, dflux=yye, sigma=kwargs['snrt'], zp=23.9)
                    __=np.where(m<99)
                    self.ax2.plot((xx[__]-self.t0)/(1+self.z), m[__]-self.dm-Rf[filt]*ebv,
                                 alpha=kwargs['alphabest'], **PROP3[filt])                    
                    #self.ax2.fill_between((xx-self.t0)/(1+self.z), m-self.dm-Rf[filt]*ebv-me,
                    #            m-self.dm-Rf[filt]*ebv+me, alpha=kwargs['alphasample'], **PROP3[filt])
                    #self.ax2.errorbar((xx[__]-self.t0)/(1+self.z), m[__]-self.dm-Rf[filt]*ebv,
                    #            yerr=me[__], alpha=kwargs['alphabest'], **PROP3[filt])
                    
        #  subplot settings        
        if kwargs['ax_xlim'] is not None:
            x1,x2 = min(kwargs['ax_xlim']), max(kwargs['ax_xlim'])
            self.ax2.set_xlim([x1/(1+self.z), x2/(1+self.z)])                     
        if kwargs['ax2_ylim'] is not None:
            self.ax2.set_ylim(kwargs['ax2_ylim'])
        else:
            self.ax2.invert_yaxis()
        if show_title:
            self.ax2.set_title('%s\n%s z=%.2f'%
                               (self.objid,self.sntype,self.z),fontsize=12)                        
        if self.texp is not None:
            self.ax2.axvline(self.texp[1]/(1+self.z), ls='--', label='texp', color='k')
            self.ax2.axvline(self.texp[0]/(1+self.z), ls=':', color='k')
            self.ax2.axvline(self.texp[2]/(1+self.z), ls=':', color='k')
        self.ax2.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=12)
        if ylabel_2right:
            ax22=self.ax2.twinx()
            ax22.set_ylabel('M$_{abs}$ (mag)',fontsize=12)
            ax22.set_yticks([])
        else:
            self.ax2.set_ylabel('M$_{abs}$ (mag)',fontsize=12)
        if show_legend: self.ax2.legend(fontsize=10, frameon=False)
                
    def _ax3(self, show_title=False, show_legend=False, ylabel_2right=True,
             corr_mkw=False, corr_host=False, source=None, **kwargs):  #mag color plot
        if self.ax3 is None:
            print ('Error: no ax3 defined, skipped color plots')
            return
        if not 'colors' in self.__dict__:
            print ('Error: no colors defined, do self.calc_colors() first')
            return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])       
        ebv = 0
        if corr_mkw: ebv += self.mkwebv
        if corr_host: ebv += self.hostebv                                    
        for _ in self.colors:
            xx = (self.colors[_][0]-self.t0)/(1+self.z)
            yy = self.colors[_][1]-self.colors[_][2]
            yye = np.sqrt(self.colors[_][3]**2 + self.colors[_][4]**2)            
            yy += ebv * (Rf['r'] - Rf['g'])
            self.ax3.errorbar(xx, yy, yerr=yye, **PROP5['gr%s'%_])
            if _==1: # set limits
                self.ax3.set_ylim([min(self.colors[1][1]-self.colors[1][2])-.1,
                                   max(self.colors[1][1]-self.colors[1][2])+.1])
                
        #  subplot settings        
        if kwargs['ax_xlim'] is not None:
            x1,x2 = min(kwargs['ax_xlim']), max(kwargs['ax_xlim'])
            self.ax3.set_xlim([x1/(1+self.z), x2/(1+self.z)])                     
        if kwargs['ax3_ylim'] is not None:
            self.ax3.set_ylim(kwargs['ax3_ylim'])        
        if show_title:
            self.ax3.set_title('%s\n%s z=%.2f'%
                               (self.objid,self.sntype,self.z),fontsize=12)                        
        if self.texp is not None :
            self.ax3.axvline(self.texp[1]/(1+self.z), ls='--', label='texp', color='k')
            self.ax3.axvline(self.texp[0]/(1+self.z), ls=':', color='k')
            self.ax3.axvline(self.texp[2]/(1+self.z), ls=':', color='k')
        self.ax3.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=12)
        f1, f2 = kwargs['color_bands']
        if ylabel_2right:
            ax33=self.ax3.twinx()
            ax33.set_ylabel('%s-%s (mag)'%(f1,f2),fontsize=12)
            ax33.set_yticks([])
        else:
            self.ax3.set_ylabel('%s-%s (mag)'%(f1,f2),fontsize=12)
        if show_legend: self.ax3.legend(fontsize=10, frameon=False)

    def _ax4(self, show_title=False, show_legend=False, ylabel_2right=False,
             logscale=True, source=None, **kwargs):  #luminosity plot
        if self.ax4 is None:
            print ('Error: no ax4 defined, skipped luminosity plots')
            return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if 'mbol' in self.__dict__:            
            for _ in self.mbol:                    
                xx = (get_numpy(self.mbol[_][0])-self.t0)/(1+self.z)
                yy = get_numpy(self.mbol[_][1])
                yye = get_numpy(self.mbol[_][2]) 
                self.ax4.errorbar(xx,yy,yerr=yye,**PROP6['lyman_data_%s'%_])
                    
        if 'mbolbb' in self.__dict__:                
            for _ in self.mbolbb:                    
                xx = (get_numpy(self.mbolbb[_][0])-self.t0)/(1+self.z)
                yy = get_numpy(self.mbolbb[_][1])
                yye = get_numpy(self.mbolbb[_][2]) 
                self.ax4.errorbar(xx,yy,yerr=yye,**PROP6['bb_data_%s'%_])
                
                #yy = get_numpy(self.mbolbb[_][7])
                #yye = get_numpy(self.mbolbb[_][8]) 
                #self.ax4.errorbar(xx,yy,yerr=yye,**PROP6['bb_data_%sf'%_])
                
        # fits
        if 'fitcls' in self.__dict__:
            if 'bol_main' in self.fitcls:
                for source in self.fitcls['bol_main']:                    
                    for model in self.fitcls['bol_main'][source]:
                        _model = self.fitcls['bol_main'][source][model]
                        
                        # fit range
                        pmin,pmax = kwargs['bol_main_xrangep']
                        p_fit = np.arange(pmin, pmax, .1)
                        
                        # samples                       
                        x,y,f=_model.predict_random(limit=kwargs['plot_mcmct'],
                                                    plotnsamples=kwargs['plot_nsamples'])
                        for xx, yy in zip(x, y):                            
                            self.ax4.plot((xx-self.t0)/(1+self.z), yy,
                                          alpha=kwargs['alphasample'], **PROP6['fit'])
                            
            if 'bol_tail' in self.fitcls:
                for source in self.fitcls['bol_tail']:                    
                    for model in self.fitcls['bol_tail'][source]:
                        _model = self.fitcls['bol_tail'][source][model]
                        
                        # fit range
                        pmin,pmax = kwargs['bol_tail_xrangep']
                        p_fit = np.arange(pmin, pmax, .1)
                        
                        # samples                       
                        x,y,f=_model.predict_random(limit=kwargs['plot_mcmct'],
                                                    plotnsamples=kwargs['plot_nsamples'])
                        for xx, yy in zip(x, y):                            
                            self.ax4.plot((xx-self.t0)/(1+self.z), yy,
                                          alpha=kwargs['alphasample'], **PROP7['fit'])
                            
        #  subplot settings
        if logscale: self.ax4.set_yscale('log')
        if kwargs['ax_xlim'] is not None:
            x1,x2 = min(kwargs['ax_xlim']), max(kwargs['ax_xlim'])
            self.ax4.set_xlim([x1/(1+self.z), x2/(1+self.z)])                     
        if kwargs['ax4_ylim'] is not None:
            self.ax4.set_ylim(kwargs['ax4_ylim']) 
        if show_title:
            self.ax4.set_title('%s\n%s z=%.2f'%
                               (self.objid,self.sntype,self.z),fontsize=12)                        
        if self.texp is not None :
            self.ax4.axvline(self.texp[1]/(1+self.z), ls='--', label='texp', color='k')
            self.ax4.axvline(self.texp[0]/(1+self.z), ls=':', color='k')
            self.ax4.axvline(self.texp[2]/(1+self.z), ls=':', color='k')
        self.ax4.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=12)
        f1, f2 = kwargs['color_bands']
        if ylabel_2right:
            ax44=self.ax.twinx()
            ax44.set_ylabel('L$_{bol}$ (erg $s^{-1}$)',fontsize=12)
            ax44.set_yticks([])
        else:
            self.ax4.set_ylabel('L$_{bol}$ (erg $s^{-1}$)',fontsize=12) 
        if show_legend: self.ax4.legend(fontsize=10, frameon=False)        
                
    def show_corner(self, ax, index=0, gp=False,
                    engine=None, model=None, source=None, **kwargs):
        if gp:
            assert 'gpcls' in self.__dict__, 'Error: no GP had been done...'
        else:
            assert 'fitcls' in self.__dict__, 'Error: no fittings had been done...'            
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        if gp: # gp
            figpath = '{}_gpcorner'.format(self.objid)
            fignames = self.gpcls.save_corner(figpath)            
        else:  # fit
            models, modelnames, enginenames, sourcenames = self.get_model(
                engine=engine,model=model, source=source
            )
            if len(models) == 0:
                print ('No models found')
                return            
            _model, _modelname, _enginename, _sourcename = models[index],\
                modelnames[index], enginenames[index], sourcenames[index]
            if _sourcename is not None:
                figpath = '{}_fitcorner_{}_{}'.format(self.objid, _modelname, _sourcename)
            else:
                figpath = '{}_fitcorner_{}'.format(self.objid, _modelname)                
            fignames = _model.save_corner(figpath)       
        if len(fignames) == 0:
            print ('No Contours found')
            return
        plt.figure(constrained_layout=True, figsize=kwargs['figsize'], dpi=400)
        print ('show corner plot %s'%fignames[0])
        img = mpimg.imread(fignames[0])
        ax.imshow(img)        
        ax.axis("off")        

    def get_model(self, engine=None, model=None, source=None):
        assert 'fitcls' in self.__dict__
        models, modelnames, enginenames, sourcenames = [], [], [], []
        for engine_name in self.fitcls:
            if engine is not None and engine != engine_name: continue
            if 'bol_' in engine_name or engine_name == 'specline':                
                for source_name in self.fitcls[engine_name]:
                    if source is not None and source != source_name: continue
                    for model_name in self.fitcls[engine_name][source_name]:
                        if model is not None and model != model_name: continue
                        models.append( self.fitcls[engine_name][source_name][model_name] )
                        modelnames.append( model_name )
                        enginenames.append( engine_name )
                        sourcenames.append( source_name )
            else:
                for model_name in self.fitcls[engine_name]:
                    if model is not None and model != model_name: continue
                    models.append( self.fitcls[engine_name][model_name] )
                    modelnames.append( model_name )
                    enginenames.append( engine_name )
                    sourcenames.append( None )
        return models, modelnames, enginenames, sourcenames
    
    def savefig(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.fig is None:
            print ('Error: no fig to be used to savefig')
            return
        self.fig.set_size_inches(kwargs['figsize'][0], kwargs['figsize'][1])                
        if '%s' in kwargs['figpath']:
            figpath = '{}/plots/{}'.format(LOCALSOURCE, kwargs['figpath']%self.objid)
        else:
            figpath = '{}/plots/{}'.format(LOCALSOURCE, kwargs['figpath'])
        self.fig.savefig(figpath, dpi=400, bbox_inches='tight')
        if kwargs['verbose']: print ('saved fig to %s'%figpath)
        
    def showfig(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        if '%s' in kwargs['figpath']:                    
            figpath = '{}/plots/{}'.format(LOCALSOURCE, kwargs['figpath']%self.objid)
        else:
            figpath = '{}/plots/{}'.format(LOCALSOURCE, kwargs['figpath'])
        if not os.path.exists(figpath):
            print ('savefig first')
            return
        plt.figure(constrained_layout=True, figsize=kwargs['figsize'], dpi=400)
        img = mpimg.imread(figpath)
        imgplot = plt.imshow(img)
        ax = plt.gca()
        ax.axis("off")
    
    def dm_error(self, filt):
        #assert filt in Rf, 'Error: define R for other filters first'
        if filt not in Rf: R = 0
        else: R=Rf[filt]
        
        # mkw error
        mkav = R * self.mkwebv
        mkav_err = mkav*0.15

        # z error
        if self.z > .1:
            c=3e5
            vpel=150.
            dz=vpel/(c*self.z)*self.z                         
            dmag=dz*5/np.log(10)*((1+self.z)/(self.z*(1+0.5*self.z)))
        else:
            dmag=0
            
        # H0 error
        deltaH0=3.
        H0=70
        distmoduncertainty=deltaH0/0.461/H0
    
        # total for error
        return np.sqrt(mkav_err**2+dmag**2+distmoduncertainty**2)        
        
    @staticmethod
    def clip_df(df, mask):
        __arr = dict()
        for kk in df.keys():  __arr[kk] = df[kk][~mask]   
        return pd.DataFrame(data=__arr)

    @staticmethod
    def read_c10():
        filename = '%s/data/c10_template.txt' % srcpath
        assert os.path.exists(filename)
        c10 = dict()
        for ll in open(filename).readlines():
            if ll[0] == '#': continue
            if ll.split()[0] == 'color':
                sntypes = ll.split()[1:]
            else:
                assert '-' in ll.split()[0], ll.split()[0]
                assert len(ll.split()[0].split('-')) == 2
                c1, c2 = ll.split()[0].split('-')
                if (c1,c2) not in c10: c10[(c1,c2)] = dict()
                for n, i in enumerate(ll.split()[1:]):                
                    c10[(c1,c2)][sntypes[n]] = i
        return c10
    
    @staticmethod
    def mag_to_flux(mag, magerr=None, limmag=None, sigma=5., units='zp', zp=23.9, wavelength=None):
        ''' with flux/mag zerop '''
        mag = get_numpy(mag)
        zp = get_numpy(zp)                          
        if units not in ['zp', 'phys']:
            raise ValueError("units must be 'zp' or 'phys'")
        elif units == 'zp':
            if zp is None: raise ValueError("zp must be float or array if units == 'zp'")
            flux = 10**(-(mag-zp)/2.5)
        else:
            if wavelength is None: raise ValueError("wavelength must be float or array if units == 'phys'")
            flux = 10**(-(mag+2.406)/2.5) / wavelength**2                
        __ = np.where( mag==99 )
        if limmag is None and magerr is None:
            flux[__] = min(flux)
            return flux
        elif limmag is not None:
            limmag = get_numpy(limmag)
            dflux = 10**((limmag - zp)/(-2.5))/sigma           
        elif magerr is not None:        
            magerr = get_numpy(magerr)
            dflux = np.abs(flux*(-magerr/2.5*np.log(10)))        
        flux[__] = dflux[__]
        return flux, dflux
        
    @staticmethod    
    def flux_to_mag(flux, dflux=None, sigma=5., units='zp', zp=23.9, wavelength=None):
        """ Converts fluxes (erg/s/cm2/A) into AB magnitudes with flux/mag zerop """
        
        flux = get_numpy(flux)
        zp = get_numpy(zp)                   
        if units not in ['zp', 'phys']:
            raise ValueError("units must be 'zp' or 'phys'")
        elif units == 'zp':
            if zp is None: raise ValueError("zp must be float or array if units == 'zp'")
            wavelength = 1.
        else:
            if wavelength is None: raise ValueError("wavelength must be float or array if units == 'phys'")
            zp = -2.406
        mag = -2.5*np.log10(flux*wavelength**2) + zp
        if dflux is None: return mag
        
        dflux = get_numpy(dflux)
        snr = flux/dflux
        limmag = -2.5*np.log10(sigma*dflux) + zp
        dmag = np.abs( -2.5/np.log(10) * dflux / flux )
        
        __ = np.where( snr < sigma )
        mag[__] = 99
        dmag[__] = 99
        return mag, dmag, limmag                        
    
    @staticmethod
    def bin_df(df, deltah = 1., xkey='jdobs', fkey='filter'):        
        
        jds = np.array(df[xkey])
        if len(jds)==0: return
        _jdbin = np.arange(min(jds), max(jds), float(deltah))
        if len(_jdbin)==0: return
    
        jdl = []
        for jd in jds:
            _id = np.argmin(abs(jd - _jdbin))
            jdl.append(_jdbin[_id])
        jdl = np.array(jdl)
    
        unique_jd, counts = np.unique(jdl,return_counts=True)
        single_jd, double_jd = list(unique_jd[ (counts==1)]), \
            list(unique_jd[ (counts>1)]) 
        __arr = []
        if not len(double_jd)==0:   
            for jd in double_jd:
                bDoubleJD = ( jdl==jd )
                __id = np.arange(len(jdl))[bDoubleJD]
                for __df0 in df.values[__id]: __arr.append(np.append(__df0, jd))
            
        for jd in single_jd:
            bSingleJD = ( jdl==jd )
            __id = np.arange(len(jdl))[bSingleJD]
            __arr.append( np.append(df.values[__id][0], jd) )
        
        return pd.DataFrame(__arr, columns=np.append(df.columns, 'jdbin'))

    @staticmethod
    def merge_df_cols(__df):
        '''
        will bug when combining columns from differents sources, i.e. ZTF FP r and marshal r
        '''
        rows = __df.shape[0]
        __df1, __n = __df.values[0], 1
        for __nk, _kk in enumerate(range(rows)):
            if __nk==0:continue
            __df1 += np.array(__df.values[_kk])
            __n += 1
        _arr = [__df2[0:int(len(__df2)/__n)] if type(__df2) is str else float(__df2)/__n for __df2 in __df1]
        return _arr

    def query_lasiar(self):
        urllink = 'https://lasair.roe.ac.uk/object/%s/json/'%self.objid
        with urllib.request.urlopen(urllink) as url:
            data = json.loads(url.read().decode())
            gid = data['crossmatches'][0]['catalogue_object_id']
            gc = data['crossmatches'][0]['catalogue_table_name']
            sep = float(data['crossmatches'][0]['separationArcsec'])
            try: gz = float(data['crossmatches'][0]['photoZ'])
            except: gz = -99            
            try: gg = float(data['crossmatches'][0]['g'])
            except: gg = -99
            if gc == 'PS1':
                gid = 'PS1 %s'%(gid)
            else: 
                try:gid = data['objectData']['annotation'].split('</a>')[0].split('<a')[1].split('>')[1].strip()
                except: gid = '' 
        return '%s & %.2f & %.2f & %.3f \\\\ \n'%(gid, sep, gg, gz)
    

def main():    
    description="run analysis on a list of SNe"    
    parser = argparse.ArgumentParser(description=description,\
         formatter_class=argparse.ArgumentDefaultsHelpFormatter)   
    parser.add_argument('-s', '--syntax', help='meta syntax',
                        default='type in ["SN Ib", "SN Ic"]')     
    parser.add_argument('-z', '--ztfid', help='specify ztfid of sn')
    parser.add_argument('-i', '--iauid', help='specify iauid of sn')
    parser.add_argument('-t', '--sntype', help='specify type of sn')
    parser.add_argument("-v", "--verbose",dest="verbose",default=False,\
                        action="store_true",help='Enable progress report')
    parser.add_argument("-c", "--clobber",dest="clobber",default=False,\
                        action="store_true",help='Redo analysis')
    args = parser.parse_args()  
    pars = vars(args)    
    
    # run
    if pars['ztfid'] is not None:  pars.update(syntax = 'ZTFID =="%s"'%pars['ztfid'])
    if pars['iauid'] is not None:  pars.update(syntax = 'IAUID =="%s"'%pars['iauid'])
    if pars['sntype'] is not None: pars.update(syntax = 'type in ["%s"]' %pars['sntype'])    
    cls = snelist( **pars )
    cls.run(axes=None)
