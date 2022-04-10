#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/ztfanalysis.py
# Author            : syang <sheng.yang@astro.su.se>
# Date              : 12.11.2019
# Last Modified Date: 11.02.2022
# Last Modified By  : syang <sheng.yang@astro.su.se>

from __future__ import print_function
description="run analysis on SNe based on ZTF data"
import os, sys, re, math, glob, warnings, logging, requests,\
       emcee, corner, random, shlex, subprocess, argparse, time, json
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab as pl
import matplotlib.image as mpimg
from scipy.optimize import curve_fit
from scipy import integrate
from scipy.integrate import simps
from astropy.cosmology import Planck13 as cosmo
from astropy.time import Time
from astroquery.ned import Ned
import astropy.units as u
from astropy import coordinates
from astropy.stats.sigma_clipping import sigma_clip
from joblib import dump, load
from multiprocessing import Pool, cpu_count
from io import StringIO
from pathlib import Path
from collections import OrderedDict

from .gaussian_process import fit_gp
from .model_fitters import fit_model
from .filters import *
from .models import *
from .corner_hack import corner_hack
from .functions import Mbol_to_Lbol, get_numpy, bbody
from .specline_fits import handle_spectrum
from .pbar import get_progress_bar
from . import __path__
srcpath = __path__[0]

__all__ = ('ztfsingle', 'ztfmultiple', 'plotter')

if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser(description=description,\
         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dir', dest='datadir', help='directory to put data', default='./')
    parser.add_argument('--meta', dest='metafile', help='meta file', default='rcf_query.txt')
    parser.add_argument('--prior', dest='parfile', help='prior info file', default='priors.txt')   
    parser.add_argument('--data', dest='datafile', help='store products to local file', default='data.npz')
    parser.add_argument('-s', '--syntax', help='meta syntax', default='type in ["SN Ib", "SN Ic"]')     
    parser.add_argument('-z', '--ztfid', help='specify ztfid of sn')
    parser.add_argument('-i', '--iauid', help='specify iauid of sn')
    parser.add_argument('-t', '--sntype', help='specify type of sn')
    parser.add_argument('-n', '--n', dest='nsn', help='specify number index of sn', type=int)   
    parser.add_argument("-v", "--verbose",dest="verbose",default=False,\
                        action="store_true",help='Enable progress report')
    parser.add_argument("-c", "--clobber",dest="clobber",default=False,\
                        action="store_true",help='Redo analysis')
    args = parser.parse_args()        

class ztfmultiple(object):

    # Static version info
    version = 1.0
    
    def __init__(self, logger=None, **kwargs):                
        """ initialize """
        defkwargs = {
            'datadir'   : '/Users/yash0613/Library/CloudStorage/Box-Box/ztf_data/',
                                  # where to put and read products
            'metafile': 'rcf_query.txt',            
            'parfile' : 'priors.txt',
            'datafile' : 'tmp/%s_data.clf',
            'syntax'  : 'type in ["SN Ib", "SN Ic"]',
            'ztfid'   : None,
            'iauid'   : None,
            'sntype'  : None,
            'nsn'     : None,                                  
            'verbose' : True,
            'clobber' : False,                   
            'sep'     : ',',
            'skiprows': [0,1],
        }
        
        # ----- define logger ----- #
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        # ----- read meta ----- #
        self.kwargs = kwargs        
        for _key in defkwargs:
            self.kwargs.setdefault(_key, defkwargs[_key])
        
        # ----- define products ----- #
        if not 'data' in self.__dict__: self.data = {}
        if not 'meta' in self.__dict__: self.meta = {}
        if not 'params' in self.__dict__: self.params = {}
        
    def parse_meta(self, force=False, **kwargs):
        if len(self.meta) > 0 and not force: return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        metafile = '%s/%s' % (kwargs['datadir'], kwargs['metafile'])
        if not os.path.exists(metafile): return
        meta = pd.read_csv(metafile,sep=kwargs['sep'], skiprows=kwargs['skiprows']).drop_duplicates()
        if kwargs['syntax'] is not None: meta = meta.query(kwargs['syntax'])        
        self.meta = meta
        if kwargs['verbose']:  print ( 'meta %i objs'%( len(meta)) )

    def format_meta(self, Rv=3.1, idkey='ZTFID',
                    sortkey='peakabs', numkey=['redshift'],
                    zkey='redshift', avkey='A_V', ebvkey='ebv', **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        assert 'meta' in self.__dict__
        if idkey in self.meta.keys():
            self.meta = self.meta.set_index(idkey)        
        if sortkey in self.meta.keys():
            self.meta = self.meta.sort_values(sortkey, ascending=False)
        #####
        try:  self.meta['redshift']["ZTF20acpzaoa"]=0.03
        except: pass
        #####
        for k in numkey:  self.meta[k] = pd.to_numeric(self.meta[k])
        if zkey is not None:
            self.meta['dist'] = cosmo.luminosity_distance( self.meta[zkey] ).value
            self.meta['dm'] = 5*np.log10(self.meta['dist'])  + 25
        assert ebvkey in self.meta.keys() or avkey in self.meta.keys()
        if ebvkey in self.meta.keys():
            self.meta['ebv'] = self.meta[ebvkey]
        else:
            self.meta['ebv'] = self.meta[avkey] / Rv

    def parse_meta_info(self, ztfid, key=None, test=False):
        _meta = self.meta.query('ZTFID==@ztfid')
        if len(_meta) != 0:
            if test:
                iauid,ra,dec,z,dist,dm,mkwebv,sntype,jdpeak = _meta['IAUID'][0],\
                    _meta['RA'][0],_meta['Dec'][0],_meta['redshift'][0],_meta['dist'][0],\
                    _meta['dm'][0],_meta['ebv'][0], _meta['type'][0], _meta['peakt'][0] + 2458000            
                return iauid,ra,dec,z,dist,dm,mkwebv,sntype,jdpeak
            if key in _meta:  return _meta[key][0]
            else: print ('select key from %s'%_meta.keys())
        
    def parse_params(self, force=False, **kwargs):
        if len(self.params) > 0 and not force: return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        parfile = '%s/%s' % (kwargs['datadir'], kwargs['parfile'])
        if not os.path.exists(parfile): return
        for line in open(parfile).readlines():
            if line[0]=='#' or len(line.strip())==0: continue
            obj = line.split()[0]
            assert 'ZTF' in obj            
            if obj not in self.params: self.params[obj] = dict()
            for _ in line.split()[1:]: 
                self.params[obj][_.split('=')[0].strip()]=_.split('=')[1].strip()
    
    def load_data(self, ztfid, force=False, **kwargs):
        if ztfid in self.data and not force: return
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        datafile = '%s/%s' % (kwargs['datadir'], kwargs['datafile']%ztfid )        
        if os.path.exists(datafile):
            self.data[ztfid] = load(datafile)
            
    def save_data(self, ztfid, force=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        datafile = '%s/%s' % (kwargs['datadir'], kwargs['datafile']%ztfid)
        if kwargs['verbose']: print ('saved cache')
        if not os.path.exists(datafile) or force:
            dump(self.data[ztfid], datafile)
            
    def run(self, axes=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        self.parse_meta()
        self.format_meta()  
        self.parse_params()       

        with get_progress_bar(kwargs['verbose'], len(self.meta.index)) as pbar:            
            for i, ztfid in enumerate(self.meta.index):                            
                iauid,ra,dec,z,dist,dm,mkwebv,sntype,jdpeak = self.parse_meta_info(ztfid,test=True)
                if kwargs['nsn'] is not None and i != kwargs['nsn']:continue
                if kwargs['ztfid'] is not None and ztfid != kwargs['ztfid']:continue
                if kwargs['iauid'] is not None and iauid != kwargs['iauid']:continue
                if kwargs['sntype'] is not None and sntype != kwargs['sntype']:continue
                
                self.load_data(ztfid)
                if not ztfid in self.data or kwargs['clobber']:
                    par = dict()
                    if ztfid in self.params: par = self.params[ztfid]
                    
                    # for each object                         
                    self.data[ztfid] = ztfsingle(ztfid, iauid=iauid, z=z,
                            mkwebv=mkwebv, sntype=sntype, dm=dm, jdpeak=jdpeak,
                            axes=axes, **par)                                
                    try:
                        self.data[ztfid].run()
                    except:
                        self.data[ztfid].run0()
                    self.save_data(ztfid,force=True)
                pbar.update(1)
                
class ztfsingle(object):

    # Static version info
    version = 1.0
    
    def __init__(self, ztfid, iauid=None, z=None, ra=None, dec=None,
                 mkwebv=None, sntype=None, dm=None, jdpeak=None,
                 logger=None, axes=None, **kwargs):
        """ initialize """   
        defkwargs = {
            ###   general parameters ###
            'targetdir' : "/Users/yash0613/Library/CloudStorage/Box-Box/ztf_data/",
                                  # set ztfquery ZTFDATA dir
                                  #    from where to take phot/spec file         
            'datadir'   : '/Users/yash0613/Library/CloudStorage/Box-Box/ztf_data/tmp/',
                                  # where to put and read products,
                                  #   e.g. mc samples, plots, etc
            'sigma'     : 5,      # signal noise ratio for detections
            'fsigma'    : 1,      # signal noise ratio for fittings
            'color_thre': 1.,     # cross match g and r within how many days for color epoch
            'c10temp'   : 1.,     # which template to use for g-r comparison                        
            'clipsigma' : None,   # for LC clipping
            'verbose'   : False,   # show detailed processing informations            
            'jdthre'    : 1,      # threshold to define colour epoch
            'peaksep'   : 3,      # threshold to define peak accuracy
            'copt'      : [1,2,3],   # how to make color
                                  #  1: epochs with both filters (g and r) within *color_thre* hours
                                  #  2: epochs with one filter, and the other from fits
                                  #  3. ...                                   from GP
            'cfilt1'    : 'g',    # filt1 for colors
            'cfilt2'    : 'r',    # filt2 for colors
            'cfilters'  : ['g','r','i'],  # filters to make blackbody
            'bolopt'    : [1,],  # how to make bolometric lc
                                  #  1. use Lyman bolometric corrections, when g/r or g/i available
                                  #  2. diluted blackbody fits, when at least 3 bands available
            'lyman_copt': [1,2,3],  # how to make lyman
                                  #  1: epochs with all filters (g/r/i) within *color_thre* hours
                                  #  2: epochs with one filter, and the others from fits
                                  #  3. ...                                    from GP
            'bb_copt'   : [1,2,3],   # how to make blackbody
                                  #  1: epochs with all filters (g/r/i) within *color_thre* hours
                                  #  2: epochs with one filter, and the others from fits
                                  #  3. ...                                    from GP
            'plot_bands'  : ['r', 'g', 'i'],  # which bands to show
            'jdmin'      : None,
            'jdmax'      : None,
            ###   early power law parameters ###
            'pl_type'   : 0,  # 0: donot do power law fits
                              # 1: when opt with mc and cache file exists, read samples
                              # 2: when opt with mc and cache file exists, redo mc for new samples
                              # 3: do but not show
            'pl_style'   : 1,  # style=1 fit for multiple filters
                               # style=2 fit for ztf forced phot only
                               # style=3 fit for ztf forced phot only, fixed alpha=2
            'pl_cache'  : 'pl_%s_%s',  # will cache all samplings to datadir,
                                          #  with ZTFname, routine as %s_%s, donot necessarily need suffix
            'pl_cache1' : 'pl1_%s_%s', # will cache best and a  random number of samplings for plot purpose
            'pl_bands'  : ['r', 'g'],  # power law fits in which bands
            'pl_routine': 'mcmc', #  Which technic to be used to realize optimization.
                                     #   Possible choices are:'mcmc', 'lm', 'trf', 'dogbox'.
                                     #   for power law, I suggest to use mcmc!!!
            'tm_pl'       : 50,      # min phase for power law relative to peak
            'flux_scale'  : 100,     # normalize flux peak           
            'rel_flux_cutoff' :.4,   # fit power from baseline up to which level of max
            'rel_flux_cutfrom':.1,   # check early data down to which level
            'pl_plotmax' :  5,     # if plot power law, from baseline up to how many days relative to peak
                                    #   if None, will not show power law samplings
            ###   Gaussian process parameters ###  
            'gp_type'   : 1,          # similar to power law
            'gp_cache'  : 'gp_%s_%s_%s',
            'gp_bands'  : ['r', 'g'],           
            'gp_routine': 'minimize', #  Possible choices are: 'minimize', 'mcmc', 'leastsq'.
            'gp_mean'   : 'bazin',  # Mean y_data function.
                                    #  Possible choices are: 'mean', 'gaussian', 'bazin', 'villar'.
            'kernel'    : 'matern32', # Kernel to be used with the gaussian process. 
                                      #  Possible choices are: 'matern52', 'matern32', 'squaredexp'.
            'fix_scale' : True,   # If fix default gaussian process param            
            'gp_fitr'   : [-60, 120], # gp fit range, relative to peak
            'gp_plotr'  : [-60, 120], # gp plot range, relative to peak
                                      #   if None, will not plot gp samplings
            ### analytic fit on peak parameters ###
            'fit_type'   : 1,
            'fit_cache'  : 'fit_%s_%s_%s',
            'fit_cache1'  : 'fit1_%s_%s_%s', 
            'fit_bands' : ['r', 'g'],           
            'fit_routine' : 'trf', # opt: 'mcmc', 'trf', 'dogbox'
            'fit_method'  : 'bazin',            
            'fit_fitr'    : [-30, 40],
            'fit_plotr'   : [-30, 40],            
            ### Arnett fitting parameters ###
            'Arnett_type'    : 1,
            'Arnett_style'   : 3,  # style=1 fit Mni and taum, v is needed to break degenracy
                                   # style=2 fit v as well, together with Mni, Ek and Mej
                                   # style=3 fit Mni, taum, and texp as well
                                   # style=4 fit Arnett (Mni, taum, texp) + radioactive tail (Mni, t0)
                                   # style=5 fit Arnett (Mni, taum, texp) + Piro SBO (Me, Re, Ee, texp)
                                   # style=6 fit Arnett + SBO + radioactive tail
            'Arnett_cache'  : 'Arnett_%s_%s_%i',
            'Arnett_cache1'  : 'Arnett1_%s_%s_%i',
            'Arnett_routine' : 'trf',            
            'Arnett_fitr'  : [-30, 40],
            'Arnett_plotr' : [-30, 40],
            'Arnett_copt'   : [1,2],  # which *copt* should be used to fit Arnett models
            'Arnett_bolopt' : 1,  # which *bolopt* should be used to fit Arnett models
            ### Tail fitting parameters ###
            'Tail_type'    : 1,
            'Tail_style'   : 1,  # style=1 fit Mni and t0 in the range of *Tail_fitr*
                                 # style=2 fit Mni, t0 and the start of tail as well,
            'Tail_cache'   : 'tail_%s_%s_%s',
            'Tail_cache1'   : 'tail1_%s_%s_%s',
            'Tail_routine' : 'trf',           
            'Tail_fitr'    : [60, 120], 
            'Tail_plotr'   : [60, 120],
            'Tail_copt'    : [1,3],  # which *copt* should be used to fit tail models
            'Tail_bolopt'   : 1,  # which *bolopt* should be used to fit Tail models
            ### for mcmc ###
            'ncores'    : None,     # how many cores to run multi-processing
            'nwalkers'  : 30,     # number of walkers
            'nsteps'    : 20000,    # number of MC steps
            'nsteps_burnin' : 5000,    # number of MC burn in steps        
            'thin_by'   : 1,  # If you only want to store and yield every thin_by samples in the chain,
                              # set thin_by to an integer greater than 1. When this is set,
                              # iterations * thin_by proposals will be made.            
            'maxfev'    : 20000,  # for scipy
                                  # The maximum number of calls to the function.          
            'emcee_burnin' : True,
            'use_emcee_backend' : True,
            'plotsamples'  : 8,  # how many random MC samples to be plotted
            'quantile'   : [0.16, 0.50, 0.84],  # use 50 percentile as mean
                                                # 1 sigma: 68% -> 16% - 84%
            'scipysamples' : 100,
            'sguess'    : True,
            ### for spectra fits ###
            'bin_method' : 'savgol',
            'bin_size' : 9,            
            'pfactor'   : 20, # a parameter needed to decide peaks
            'spec_type'    : 1,
            'spec_cache'   : 'spec_%s_%s_%s',
            'spec_cache1'   : 'spec_%s_%s_%s',
            'spec_routine' : 'trf',
            'spec_fitr'    : [6000, 20000], # unit in km/s
            'spec_plotr'   : [6000, 20000],
            'spec_guess_red'  : True,
            ### for plotter ###            
            'figsize' : (8, 10), # figure size
            'figpath' : '%s.png', # figure path
            'alphabest'   : .6, # alpha p for best sample
            'alphasample' : .2, # alpha p for other sample
            'ax_ylim'   : None,  # ylim
            'ax1_ylim'   : None,
            'ax2_ylim'   : None,
            'ax3_ylim'   : None,
            'ax4_ylim'   : None,
            'ax_xlim'    : [-50,120],  # x lim, x is phase relative to peak,
                                       #  only working when t0 avaibale
        }
                 
        # ----- define logger ----- #
        #if logger is None:
        #    logging.basicConfig(level = logging.INFO)
        #    self.logger = logging.getLogger(__name__)
        #else:
        #    self.logger = logger  

        # ----- read meta ----- #
        self.kwargs = dict()  # read from parfile
        for _key in kwargs: self.kwargs[_key] = eval( kwargs[_key] )
        
        for _key in defkwargs: # read default
            self.kwargs.setdefault(_key, defkwargs[_key])
        
        # set to ZTF data dir
        os.environ["ZTFDATA"] = self.kwargs['targetdir']
        global bts, marshal, fritz
        from ztfquery import bts, marshal, fritz
        
        self.ztfid  = ztfid
        self.iauid  = iauid
        self.ra  = ra
        self.dec  = dec
        if z is not None: self.z = z # rest frame phase 
        else:  self.z = 0       # phase
        if mkwebv is not None: self.mkwebv = mkwebv
        else:  self.mkwebv = 0  # ignore milky ebv
        self.hostebv = 0        # initially with no host
        self.sntype = sntype
        self.dm     = dm
        if jdpeak is not None:
            self.t0 = jdpeak  # zero point thoughout the analysis                                  
        else:
            self.t0 = 0  # if None, try to get it via GP since it was needed for fits/pl/arnett...
        self.tpeak  = dict()  # peak epoch relative to self.t0
        self.fpeak  = dict()  # peak Fmcmc for GP or fit
        self.axes   = axes    # if axes is None, will init axes
                              # if axes is False, will not make plots
                              # otherwise, make axes = [ax,ax1,ax2, ...], where ax is defined by user              

    def config_ztfquery(self, fritz_token='733be155-a0c2-41c6-b994-d877fb4a0088',
                        marshal_username='saberyoung', marshal_password='Ys_19900615', ):
        ztfquery.io.set_account('fritz', token=fritz_token, force=True)
        ztfquery.io.set_account('marshal', username=marshal_username, password=marshal_password)
        ztfquery.io.test_irsa_account()

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

    def add_lc(self, df, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert 'mag' in df.keys() and 'emag' in df.keys() and 'jdobs' in df.keys()
        if not 'lc' in self.__dict__:
            self.lc = df
        else:
            self.lc = self.lc.append(df, ignore_index=True).drop_duplicates()
        jdmin, jdmax = kwargs['jdmin'], kwargs['jdmax']
        if jdmin is not None: self.lc = self.lc.query('jdobs>=@jdmin')
        if jdmax is not None: self.lc = self.lc.query('jdobs<=@jdmax')
        
    def bin_fp_ztf(self, binDays=3, resultsPath=None, outPath=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        cmd = 'python %s/plot_atlas_fp.py stack %.2f %s '%(srcpath, binDays, resultsPath)
        if kwargs['jdmin'] is not None and kwargs['jdmax'] is not None:            
            cmd += '%f %f ' % (kwargs['jdmin'],kwargs['jdmax'])
        cmd += '--o %s'%(outPath)                
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
        targetdir = '%s/ForcePhot_atlas/'%kwargs['targetdir']
        wdir = '%s/%s/'%(targetdir, self.ztfid)
        f_unbin = '%s/forcedphotometry_%s_lc.csv'%(wdir, self.ztfid)        
        if binDays is None:            
            if not os.path.exists(f_unbin):
                print ('Error: %s not found'%f_unbin)
                return None
            self._get_fp_atlas(f_unbin, **kwargs)
        else:            
            binDays = float(binDays)            
            f = '%s/forcedphotometry_%s_lc_atlas_fp_stacked_%.2f_days.txt'%\
                (wdir, self.ztfid, binDays)
            if not os.path.exists(f) or clobber:
                if not os.path.exists(f_unbin):
                    print ('Error: %s not found'%f_unbin)
                    return None
                self.bin_fp_ztf(binDays=binDays, resultsPath=f_unbin,
                                outPath=wdir, **kwargs)
            self._get_fp_atlas(f, binned=True, **kwargs)
            
    def _get_fp_atlas(self, f, binned=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        sigma = kwargs['sigma']
        if binned:
            df = pd.read_csv(f,sep = ',',skiprows=[0,1,2,3,4,5]).drop_duplicates().reset_index(drop=True)
        else:
            df = pd.read_csv(StringIO('\n'.join(open(f).readlines()).replace("###", "")), delim_whitespace=True)
        df.rename(columns={'uJy':'flux','duJy':'eflux','F':'filter','MJD':'mjd',}, inplace=True)
        # calc mag
        mags, sigmamags, limmags = [], [], []
        f, fe = [], []
        for flux, fluxerr in zip(df['flux'], df['eflux']):
            snr = flux/fluxerr                         
            if flux > 0.0 and snr>sigma:
                magpsf = -2.5*np.log10(flux) + 23.9
                sigmamagpsf = abs(-2.5/np.log(10) * fluxerr / flux)
            else:
                magpsf = 99
                sigmamagpsf = 99            
            limmag = -2.5*np.log10(sigma*fluxerr) + 23.9
            mags.append(magpsf)
            sigmamags.append(sigmamagpsf)
            limmags.append(limmag)
            F0 = 10**(23.9/2.5)           
            f.append( flux / F0 )
            fe.append( fluxerr / F0 )            
        df['mag'] = mags
        df['emag'] = sigmamags
        df['limmag'] = limmags        
        df['jdobs'] = df['mjd'] + 2400000.5
        df['Fratio'] = f
        df['Fratio_unc'] = fe
        self.add_lc(df, **kwargs)
        
    def get_alert_ztf(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        df = marshal.get_local_lightcurves(self.ztfid)
        self.add_lc(df, **kwargs)
    
    def get_fp_ztf(self, **kwargs):
        ''' get local ZTF forced phot LCs
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        targetdir = '%s/ForcePhot/'%kwargs['targetdir']
        sigma = kwargs['sigma']
        f = '%s/%s/forcedphotometry_%s_lc.csv'%(targetdir, self.ztfid, self.ztfid)
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
        Fratio = Fpsf / F0
        eFratio2 = (eFpsf / F0)**2 + (Fpsf * eF0 / F0**2)**2
        eFratio = np.sqrt(eFratio2)
        df['Fratio'] = Fratio
        df['Fratio_unc'] = eFratio
        
        filt = df['filter']
        filterid = np.zeros(len(df))
        filterid[filt=='ZTF_g']=1
        filterid[filt=='ZTF_r']=2
        filterid[filt=='ZTF_i']=3
        df['filterid'] = filterid    
        df['chi2_red'] = np.array([np.float(x) for x in df['chi2_red'].values])    
        df['fcqfid'] = df['field']*10000 + df['ccdid']*100 + df['qid']*10 + df['filterid']
    
        # calc mag
        mags, sigmamags, limmags, filters = [], [], [], []
        for filtro, zpdiff, flux, fluxerr, snr in zip(df['filter'], df['zp'], 
                    df['Fpsf'], df['Fpsf_unc'], df['Fpsf_snr']):
            if flux > 0.0 and snr>sigma:
                magpsf = -2.5*np.log10(flux) + zpdiff
                sigmamagpsf = abs(-2.5/np.log(10) * fluxerr / flux)
            else:
                magpsf = 99
                sigmamagpsf = 99
            limmag = -2.5*np.log10(sigma*fluxerr) + zpdiff
            mags.append(magpsf)
            sigmamags.append(sigmamagpsf)
            limmags.append(limmag)
            filters.append(filtro.replace('ZTF_',''))    
        df['mag'] = mags
        df['emag'] = sigmamags
        df['limmag'] = limmags
        df['filter'] = filters
        self.add_lc(df, **kwargs)
        
    def query_fp_atlas(self, clobber=False, verbose=True,
            mjdstart=None, mjdend=None, atlas_username='saberyoung', 
            atlas_passord='Ys19900615', **kwargs):
        ''' ATLAS forced phot query
        https://fallingstar-data.com/forcedphot/static/apiexample.py
        '''
        BASEURL = "https://fallingstar-data.com/forcedphot"
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        radeg, decdeg = self.parse_coo(verbose=verbose, deg=True)        
        try: mjdstart = float(mjdstart)
        except:
            if verbose: print ('specify mjdstart')
            return
        try: mjdend = float(mjdend)
        except:
            mjdend = self.mjd_now(jd=False)
            if verbose: print ('use current mjd')            
        targetdir = '%s/ForcePhot_atlas/'%kwargs['targetdir']        
        f = '%s/%s/forcedphotometry_%s_lc.csv'%(targetdir, self.ztfid, self.ztfid)
        if os.path.exists(f) and not clobber:
            if verbose: print ('file exists: %s'%f)
            return
        if os.environ.get('ATLASFORCED_SECRET_KEY'):
            token = os.environ.get('ATLASFORCED_SECRET_KEY')
            if verbose: print('Using stored token')
        else:
            data = {'username': atlas_username,
                    'password': atlas_passord}
            resp = requests.post(url=f"{BASEURL}/api-token-auth/", data=data)

            if resp.status_code == 200:
                token = resp.json()['token']
                if verbose: 
                    print(f'Your token is {token}')
                    print('Store this by running/adding to your .zshrc file:')
                    print(f'export ATLASFORCED_SECRET_KEY="{token}"')
            else:
                if verbose: 
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
                    if verbose: print(f'The task URL is {task_url}')
                elif resp.status_code == 429:  # throttled
                    message = resp.json()["detail"]
                    if verbose: print(f'{resp.status_code} {message}')
                    t_sec = re.findall(r'available in (\d+) seconds', message)
                    t_min = re.findall(r'available in (\d+) minutes', message)
                    if t_sec:
                        waittime = int(t_sec[0])
                    elif t_min:
                        waittime = int(t_min[0]) * 60
                    else:
                        waittime = 10
                    if verbose: print(f'Waiting {waittime} seconds')
                    time.sleep(waittime)
                else:
                    if verbose: 
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
                        if verbose: print(f"Task is complete with results available at {result_url}")
                    elif resp.json()['starttimestamp']:
                        if not taskstarted_printed:
                            if verbose: print(f"Task is running (started at {resp.json()['starttimestamp']})")
                            taskstarted_printed = True
                        time.sleep(2)
                    else:
                        if verbose: print(f"Waiting for job to start (queued at {resp.json()['timestamp']})")
                        time.sleep(4)
                else:
                    if verbose: 
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

    def query_fp_ztf(self, query=False, ztf_email='sheng.yang@astro.su.se',
            ztf_userpass='augj975', clobber=False, verbose=True, jdstart=None,
            jdend=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        radeg, decdeg = self.parse_coo(verbose=verbose, deg=True)        
        try: jdstart = float(jdstart)
        except:
            if verbose: print ('specify jdstart')
            return
        try: jdend = float(jdend)
        except:
            jdend = self.mjd_now(jd=True)
            if verbose: print ('use current jd')
        targetdir = '%s/ForcePhot/'%kwargs['targetdir']        
        f = '%s/%s/forcedphotometry_%s_lc.csv'%(targetdir, self.ztfid, self.ztfid)        
        if os.path.exists(f) and not clobber:
            if verbose: print ('file exists: %s'%f)
            return
        c = coordinates.SkyCoord('%s %s'%(self.ra, self.dec), unit=(u.hourangle, u.deg))        
        line = 'wget --http-user=ztffps --http-passwd=dontgocrazy!'+\
            ' -q -O log.txt '+ \
            '"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?'+\
            'ra=%.7f&dec=%.7f&jdstart=%.5f&jdend=%.5f&'%(radeg, decdeg, jdstart, jdend)+\
            'email=%s&userpass=%s"'%(ztf_email, ztf_userpass)
        if query: subprocess.Popen(line, shell=True)
        else:  print (self.ztfid, c.ra.deg, c.dec.deg, jdstart, jdend, line)
            
    def query_alert_ztf(self, source=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert source in ['marshal', 'fritz',None]
        if source in ['marshal',None]:
            marshal.download_lightcurve(self.ztfid)         
        if source in ['fritz',None]:
            fritz.download_lightcurve(self.ztfid)
            
    def add_flux(self, **kwargs):
        assert 'lc' in self.__dict__        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        lc = self.lc
        if 'Fratio' in self.lc:
            # to flux Jy        
            lc['flux'], lc['flux_unc'] = self.mag_to_flux(lc['mag'], magerr=lc['emag'])
            lcdet = lc.query('mag<99')
            maxflux = max(lcdet['flux'])            
            lc['Fmcmc'] = maxflux*lc['Fratio']/max(lcdet['Fratio'])           
            lc['Fmcmc_unc'] = lc['Fmcmc']*lc['Fratio_unc']/lc['Fratio']
        else:            
            # to flux Jy
            ''' !!! make sure no upper limit mag numbers here !!!'''
            lc['flux'], lc['flux_unc'] = self.mag_to_flux(lc['mag'], magerr=lc['emag'])
            lc['Fmcmc'] = lc['flux']
            lc['Fmcmc_unc'] = lc['flux_unc']
        self.lc = lc.query('Fmcmc_unc>0')

    def query_spectra(self, source=None, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert source in ['marshal', 'fritz',None]
        if source in ['marshal',None]:
            marshal.download_spectra(self.ztfid)            
        if source in ['fritz',None]:
            fritz.download_spectra(self.ztfid, get_object=True, store=True, verbose=False)

    def query_tns(self, tns_botid=131335, tns_botname = "kinder_bot",
                  tns_api='603cc592cbb2b9ac8d63fcf5dfbe18e0bd981982',
                  **kwargs):
        TNS         = "sandbox.wis-tns.org"
        url_tns_api = "https://" + TNS + "/api/get"
        search_url  = url_tns_api + "/search"
        tns_marker = 'tns_marker{"tns_id": "' + str(tns_botid) + '", "type": "bot", "name": "' + tns_botname + '"}'        
        headers = {'User-Agent': tns_marker}
        json_file = OrderedDict(search_obj)
        search_data = {'api_key': TNS_API_KEY, 'data': json.dumps(json_file)}
        response = requests.post(search_url, headers = headers, data = search_data)
        
    def get_local_spectra(self,**kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        self.spec = dict()        
        
        # spec from marshal
        dir_ = marshal.target_spectra_directory(self.ztfid)
        if os.path.exists( dir_ ):
            for d in os.listdir( dir_ ): 
                ztfid,epoch,tel = d.split('_')[0],d.split('_')[1],d.split('_')[2]
                try: 
                    data, header = fritz.parse_ascii( open( os.path.join(dir_,d) ).read().splitlines() )
                    self.spec['%s %s'%(epoch,tel)] = {'source':'GM', 'epoch':epoch,
                    'instru': tel, 'w': get_numpy(data['lbda']),'f': get_numpy(data['flux'])}
                except:  print ( os.path.join(dir_,d) )
        # spec from fritz
        _spec = fritz.FritzSpectrum.from_name( self.ztfid )
        if _spec is not None:
            if type(_spec) is list:
                for _ in _spec: 
                    epoch, tel = _.observed_at.split('T')[0].replace('-',''), _.instrument
                    self.spec['%s %s'%(epoch,tel)] = {'source':'fritz', 'epoch':epoch,
                                'instru': tel, 'w':_.lbda, 'f':_.flux} 
            else:
                _ = _spec
                epoch, tel = _.observed_at.split('T')[0].replace('-',''), _.instrument
                self.spec['%s %s'%(epoch,tel)] = {'source':'fritz',
                        'epoch':epoch, 'instru': tel, 'w':_.lbda, 'f':_.flux} 
        # spec from TNS
        dir_ = os.path.join(kwargs['targetdir'],"TNS/spectra",self.ztfid)
        if os.path.exists( dir_ ):
            for d in os.listdir( dir_ ): 
                ztfid,epoch,tel = d.split('_')[0],d.split('_')[1],d.split('_')[2]
                data, header = fritz.parse_ascii( open( os.path.join(dir_,d) ).read().splitlines() )
                self.spec['%s %s'%(epoch,tel)] = {'source':'TNS',
                    'epoch':epoch, 'instru': tel, 
                    'w': get_numpy(data['lbda']), 'f': get_numpy(data['flux'])}
        
    def run(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        self.init_fig()
        t_start = time.time()
        
        ''' parse local data via ztfid '''    
        self.get_fp_ztf()  # obtain forced phot lc              
        self.add_flux() # add Fmcmc item
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

        self.set_t0_withgp()
        #self.set_t0_withfit()
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
            
        self.arnett_fit() # fit lums around peak to arnett model
        if kwargs['verbose']:
            print("fit arnett model in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.tail_fit()   # fit lums at tail to gamma leakage model
        if kwargs['verbose']:
            print("fit tail model in = {:.2f} s".format(time.time() - t_start))
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

    def run0(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        self.init_fig()
        t_start = time.time()
        
        ''' parse local data via ztfid '''        
        self.get_local_forced_lightcurves()  # obtain forced phot lc              
        self.add_flux() # add Fmcmc item        
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

        self.set_t0_withgp()
        #self.set_t0_withfit()
        self.set_texp_midway() # try to set texp
        
        self.run_fit()   # sn-like analytic functions
        if kwargs['verbose']:
            print("Run Fitting in = {:.2f} s".format(time.time() - t_start))
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
                  
    def taum_to_EM(self, vpeak=10, vpeake=1, vopt=1, **kwargs):
        ''' taum and v to Mej and Ek
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        assert 'sntype' in self.__dict__, 'on_sntype() first ...'
        assert 'arnettcls' in self.__dict__, 'run Arnett fitting first ...'
        assert vpeak < 100 and vpeak > 1, 'unit of vpeak is 1000 km/s'

        if vopt == 1:
            # following Dessart 16, convert photospheric velocity (vej) to characteristic velicity (vm)
            #  vej is O 7772 line velocity for SNe Ic, and He 5886 for Ib, at the peak epoch.
            #  https://academic.oup.com/mnras/article/458/2/1618/2589109
            if self.sntype in ['SN Ib', 'SN IIb']:
                vm, dvm = (vpeak-2.64) / 0.765, vpeake / 0.765
            elif self.sntype in ['SN Ic', 'SN Ic-BL']:
                vm, dvm = (vpeak-2.99) / 0.443, vpeake / 0.443
            else:
                if kwargs['verbose']: print ('check if Dessart et al suitable for your sn type')
                return
            func = Arnett_mej_ek
        elif vopt == 2:
            # E/M = 3/10*V**2
            vm, dvm, func = vpeak, vpeake, Arnett_mej_ek_1            
            
        __ = self.arnettcls.get_par(filt=None, quant=kwargs['quantile'])               
        mni = __[1][0]
        mnierr = self.penc_to_errors(mni, __[0][0], __[2][0])
        taum = __[1][1]        
        taumerr = self.penc_to_errors(taum, __[0][1], __[2][1])                
        return func(taum, vm, taumerr=taumerr, vejerr=dvm)
        
    def meas_specline(self, **kwargs):
        ''' fits on bolometric tail
        '''        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        assert 'spec' in self.__dict__, 'parse local spectra first'
        assert self.t0 > 2400000, 'set t0 first'
        
        if kwargs['spec_type'] == 0: return
        elif kwargs['spec_type'] == 1: clobber = False
        elif kwargs['spec_type'] == 2: clobber = True
        elif kwargs['spec_type'] == 3: clobber = False
        else: return
        
        self.specls = dict()
        for n, k in enumerate(sorted(self.spec)[::-1]):
            epoch = self.spec[k]['epoch']
            t = '%s-%s-%s'%(epoch[:4],epoch[4:6],epoch[6:])
            at = Time(t, format='isot', scale='utc')
            phase = (at.jd-self.t0)/(1+self.z)
            phase = '%.2f' % phase
            self.specls[phase] = handle_spectrum(self.spec[k]['w'],
                                   self.spec[k]['f'],
                                   z = self.z,
                                   ebv = self.mkwebv,
                                   ax = self.ax1,
                                   cw = self.cw,
                                   fnorm = 1,
                                   ys = n,
                                   region = self.region,
                                   instru = self.spec[k]['instru'],
                                   epoch = epoch,
                                   source = self.spec[k]['source'],                                   
                                   phase = phase,
                                   clobber = clobber,
                                   mcmc_h5_file=kwargs['spec_cache']%\
                                      (self.ztfid, phase, kwargs['spec_routine']),**kwargs) 
            self.specls[phase].run()        
            
        #if self.vpeak is None: # exponential fits
        if False:
            p0=(10,-15,-.4, 7)
            bounds=((0,-20,-1,0),(40,0,0,15))
            BBparams, covar = curve_fit(exp,pl,vl,sigma=vel,absolute_sigma=True,
                                        maxfev=20000,p0=p0,bounds=bounds)
            perr = np.sqrt(np.diag(covar))
            self.vpeak = exp(0, *BBparams)
            
            sigma=kwargs['fsigma']
            nsteps=kwargs['nsteps_burnin']
            samples = np.ndarray(shape=(nsteps,4), dtype=float)
            self.vpeake = 0.            
            for _n, (_1,_2) in enumerate(zip(BBparams + perr*sigma, BBparams - perr*sigma)):                
                __ = np.random.uniform(low=_1, high=_2, size=(nsteps))
                for _nn in range(nsteps):  samples[_nn, _n] = __[_nn]
            for _theta in samples:                
                model = exp(0, *_theta)
                self.vpeake = max(self.vpeake, abs(self.vpeak-model))
            
    def tail_fit(self,**kwargs):
        ''' fits on bolometric tail
        '''        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if kwargs['Tail_type'] == 0: return
        elif kwargs['Tail_type'] == 1: clobber = False
        elif kwargs['Tail_type'] == 2: clobber = True
        elif kwargs['Tail_type'] == 3: clobber = False
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
        
        if kwargs['Tail_style'] == 1:
            fit_mean='tail'        
        elif kwargs['Tail_style'] == 2:
            fit_mean='tail_fitt'       
        else: return
        
        self.tailcls = fit_model(xx[__],yy[__],yye[__],filters=None)        
        self.tailcls.train(t_fl=18, opt_routine=kwargs['Tail_routine'],
                    fit_mean=fit_mean, nwalkers=kwargs['nwalkers'],
                    nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'],ncores=ncores,
                    thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
                    mcmc_h5_file=kwargs['Tail_cache']%
                           (self.ztfid, kwargs['Tail_routine'], kwargs['Tail_style']),
                    emcee_burnin=kwargs['emcee_burnin'], p0=None, datadir=kwargs['datadir'],
                    use_emcee_backend=kwargs['use_emcee_backend'], scipysamples=kwargs['scipysamples'],
                    clobber=clobber, verbose=kwargs['verbose'], sguess=kwargs['sguess'])
        #self.tailcls.predict(quant=kwargs['quantile'])
        self.tail_theta, self.tail_rtheta = self.get_random_samples(self.tailcls.samples,
                    self.tailcls.lnprob, cachefile=kwargs['Tail_cache1']%
                            (self.ztfid, kwargs['Tail_routine'], kwargs['Tail_style']),                    
                    clobber=clobber, verbose=kwargs['verbose'], datadir=kwargs['datadir'],)
         
    def arnett_fit(self,**kwargs):
        ''' power law on early lc
        currently explosion epoch is not a fitted parameter
        so do power law fits first
        '''        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
                
        if kwargs['Arnett_type'] == 0: return
        elif kwargs['Arnett_type'] == 1: clobber = False
        elif kwargs['Arnett_type'] == 2: clobber = True
        elif kwargs['Arnett_type'] == 3: clobber = False
        else: return

        if kwargs['Arnett_bolopt'] == 1:
            assert 'mbol' in self.__dict__, print ('constrcut bolometric lc first')
            mbol = self.mbol
        elif kwargs['Arnett_bolopt'] == 2:
            assert 'mbolbb' in self.__dict__, print ('constrcut bolometric lc first')
            mbol = self.mbolbb
        else:
            if kwargs['verbose']: print ('skip tail fit')
            return 
        assert self.t0 > 2400000, 'set t0 first'
        
        ncores = kwargs['ncores']
        if ncores is None: ncores = cpu_count() - 1
        
        xx, yy, yye = [], [], []
        for _ in kwargs['Arnett_copt']:
            if _ in mbol:
                xx = np.append(xx, get_numpy(mbol[_][0]))
                yy = np.append(yy, get_numpy(mbol[_][1]))
                yye = np.append(yye, get_numpy(mbol[_][2]))

        # Arnett fits        
        xx = (xx-self.t0)/(1+self.z)
        if kwargs['Arnett_style'] == 1:
            fit_mean='arnett'
            assert 'texp' in self.__dict__
            xx -= self.texp[1]            
        elif kwargs['Arnett_style'] == 2:
            fit_mean='arnett_fitv'
            assert 'texp' in self.__dict__
            xx -= self.texp[1]            
        elif kwargs['Arnett_style'] == 3:
            fit_mean='arnett_fitt'
        elif kwargs['Arnett_style'] == 4:
            fit_mean='arnett_tail'
        else: return
        p1, p2 = min(kwargs['Arnett_fitr']), max(kwargs['Arnett_fitr'])
        __ = np.logical_and(xx>=p1, xx<=p2)
        
        self.arnettcls = fit_model(xx[__],yy[__],yye[__],filters=None)         
        self.arnettcls.train(t_fl=18, opt_routine=kwargs['Arnett_routine'],
                    fit_mean=fit_mean, nwalkers=kwargs['nwalkers'],
                    nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'], ncores=ncores,
                    thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
                    mcmc_h5_file=kwargs['Arnett_cache']%
                             (self.ztfid, kwargs['Arnett_routine'], kwargs['Arnett_style']),
                    emcee_burnin=kwargs['emcee_burnin'], datadir=kwargs['datadir'],
                    use_emcee_backend=kwargs['use_emcee_backend'], clobber=clobber,
                    verbose=kwargs['verbose'], p0=None, bounds=None, scipysamples=kwargs['scipysamples'],
                    sguess=kwargs['sguess'])
        self.arnett_theta, self.arnett_rtheta = self.get_random_samples(self.arnettcls.samples,
                    self.arnettcls.lnprob, cachefile=kwargs['Arnett_cache1']%
                            (self.ztfid, kwargs['Arnett_routine'], kwargs['Arnett_style']),                    
                    clobber=clobber, verbose=kwargs['verbose'],datadir=kwargs['datadir'],)
        if kwargs['Arnett_style'] == 3:
            # define texp here
            _, _1, _2 = self.arnettcls.get_par(filt=None, quant=kwargs['quantile'])               
            self.texp = [_[2], _1[2], _2[2]]
            
    def set_texp_midway(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
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
                
    def est_host_c10(self, **kwargs):
        ''' 
        should include milky way ebv, and due the additional color to host
        '''    
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        if kwargs['cfilt1'] != 'g' and kwargs['cfilt2'] != 'r':
            if kwargs['verbose']:
                print ('define bc function here in est_host_c10() for your filters first')
            return
        filt1, filt2 = kwargs['cfilt1'], kwargs['cfilt2']
        if 'c10_temp' not in self.__dict__:
            if kwargs['verbose']:
                print ('define template g-r at 10 rest frame days post peak, otherwise assume no host ebv')
            return
        assert self.t0 > 2400000, 'set t0 first'
        
        # get g-r at 10 d
        # from analytic sn lc fit
        if 'fitcls' in self.__dict__ and filt1 in self.fitcls and filt2 in self.fitcls:
            p_fit = [10]
            xx,yy,yy1,yy2,_ = self.fitcls[filt1].predict(x_pred=p_fit,
                                returnv=True, quant=kwargs['quantile'])
            m1 = self.flux_to_mag(yy, dflux=None)[0] - self.mkwebv * Rf[filt1]
            xx,yy,yy1,yy2,_ = self.fitcls[filt2].predict(x_pred=p_fit,
                                returnv=True, quant=kwargs['quantile'])
            m2 = self.flux_to_mag(yy, dflux=None)[0] - self.mkwebv * Rf[filt2]
            gr10 = m1 - m2
        # or from GP intepolation        
        elif 'gpcls' in self.__dict__: 
            p_fit = [self.t0 + 10 * (1+self.z)]
            xx,yy,yye,ff = self.gpcls.predict(x_pred=p_fit, returnv=True,)            
            _yy = yy[ np.where(ff==filt1) ]            
            m1 = self.flux_to_mag(_yy, dflux=None)[0] - self.mkwebv * Rf[filt1]            
            _yy = yy[ np.where(ff==filt2) ]
            m2 = self.flux_to_mag(_yy, dflux=None)[0] - self.mkwebv * Rf[filt2]
            gr10 = m1 - m2
        else:
            gr10 = None
        
        if gr10 is not None:
            if gr10 < self.c10_temp[0]-self.c10_temp[1]:
                egrhost = 0
            elif gr10 > self.c10_temp[0]-self.c10_temp[1] and gr10 < self.c10_temp[0]+self.c10_temp[1]:
                egrhost = 0
            else:
                egrhost = gr10 - self.c10_temp[0]
            self.hostebv = egrhost/0.901
            
    def lyman_bol(self, **kwargs):
        ''' 
        should correct with both host and milky ebv      
        will only use BCg
        '''        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if not 'colors' in self.__dict__:
            if kwargs['verbose']: print ('calc colors first')
            return
        if self.sntype is None or self.sntype not in ['SN Ib', 'SN Ibn', 'SN IIb',
                        'SN Ic', 'SN Ic-BL', 'SN II', 'SN Ic-BL?']:
            if kwargs['verbose']: print ('define sn type first')
            return
        if kwargs['cfilt1'] != 'g' and kwargs['cfilt2'] != 'r':
            if kwargs['verbose']:
                print ('define bc function here in lyman_bol() for your filters first')
            return
        filt1, filt2 = kwargs['cfilt1'], kwargs['cfilt2']
        self.mbol = dict()
        ebv = self.mkwebv + self.hostebv
        for _ in self.colors:
            jd, g, r, ge, re = self.colors[_]
            '''
            !!! for SE SNe, if for type II, use sntype=2
            and for cooling phase, use phase=2
            '''
            if self.sntype in ['SN Ib', 'SN Ic', 'SN Ibn', 'SN Ic-BL', 'SN Ic-BL?', 'SN IIb']: sntype = 1
            elif self.sntype in ['SN II',]: sntype = 2            
            BCg, __ = BC_Lyman(g - r + ebv * (Rf['r'] - Rf['g']) , phase=1, sntype=sntype)
            mbol, embol = BCg + g -self.dm, np.sqrt(ge**2 + self.dm_error('g')**2)
            Lbol = Mbol_to_Lbol(mbol)
            eLbol = abs(Mbol_to_Lbol(mbol+embol)-Lbol)
            _ck=np.isfinite(Lbol)
            self.mbol[_] = (jd[_ck], Lbol[_ck], eLbol[_ck])

    def calc_colors(self, **kwargs):
        ''' define g-r epochs / filt1 - filt2
        !!! haven't included any ebv value
        '''        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])       
        assert self.t0 > 2400000

        filt1, filt2 = kwargs['cfilt1'], kwargs['cfilt2']
        self.colors = dict()
        
        if 1 in kwargs['copt']:
            # bin and match g and r
            self.lc = self.bin_df(self.lc, deltah = kwargs['color_thre'] )            
            self.match_colors(**kwargs)        

            _lc = self.lc_match.query('mag_%s<99 and mag_%s<99' % (filt1, filt2))
            if len(_lc) > 0:
                self.colors[1] = [_lc['jdobs_%s'%filt1], _lc['mag_%s'%filt1], _lc['mag_%s'%filt2],
                                  _lc['emag_%s'%filt1], _lc['emag_%s'%filt2]]
        
        # or with GP/fit intepolation
        jd2,m2,mm2,me2,mme2 = [],[],[],[],[]
        jd3,m3,mm3,me3,mme3 = [],[],[],[],[]
        for _filt, _filt1 in zip([filt1, filt2], [filt2, filt1]):
            __lc = self.lc.query('filter==@_filt and mag<99')            
            if 'fitcls' in self.__dict__ and 2 in kwargs['copt']: # analytic sn lc fit
                p_fit = (__lc['jdobs']-self.t0) / (1+self.z)
                __ = np.logical_and(p_fit > min(kwargs['fit_fitr']),
                                    p_fit < max(kwargs['fit_fitr']),)
                if _filt1 in self.fitcls: 
                    xx,yy,yy1,yy2,_ = self.fitcls[_filt1].predict(x_pred=p_fit[__],
                                    returnv=True, quant=kwargs['quantile'])
                    mm = self.flux_to_mag(yy, dflux=None)
                    _mm1 = self.flux_to_mag(yy1, dflux=None)
                    _mm2 = self.flux_to_mag(yy1, dflux=None)
                    mme = (abs(mm-_mm1) + abs(mm-_mm2))/2
                    if _filt == filt1:  
                        jd2 = np.append(jd2, __lc['jdobs'][__])
                        m2 = np.append(m2, __lc['mag'][__])
                        mm2 = np.append(mm2, mm)
                        me2 = np.append(me2, __lc['emag'][__])
                        mme2 = np.append(mme2, mme)
                    else: 
                        jd2 = np.append(jd2, __lc['jdobs'][__])
                        mm2 = np.append(mm2, __lc['mag'][__])
                        m2 = np.append(m2, mm)
                        mme2 = np.append(mme2, __lc['emag'][__])
                        me2 = np.append(me2, mme)
                    
            if 'gpcls' in self.__dict__ and 3 in kwargs['copt']: # GP intepolation
                p_fit = __lc['jdobs']
                xx,yy,yye,ff = self.gpcls.predict(x_pred=p_fit, returnv=True,)
                __ = np.where(ff==_filt1)
                xx = xx[__]
                yy = yy[__]
                yye = yye[__]
                mm = self.flux_to_mag(yy, dflux=None)
                mme = abs(self.flux_to_mag(yy, dflux=None)-
                          self.flux_to_mag(yy+yye, dflux=None))                          
                if _filt == filt1:
                    jd3 = np.append(jd3, __lc['jdobs'])
                    m3 = np.append(m3, __lc['mag'])
                    mm3 = np.append(mm3, mm)
                    me3 = np.append(me3, __lc['emag'])
                    mme3 = np.append(mme3, mme)
                else:
                    jd3 = np.append(jd3, __lc['jdobs'])
                    mm3 = np.append(mm3, __lc['mag'])
                    m3 = np.append(m3, mm)
                    mme3 = np.append(mme3, __lc['emag'])
                    me3 = np.append(me3, mme)
        if len(jd2)>0: self.colors[2] = (jd2,m2,mm2,me2,mme2)
        if len(jd3)>0: self.colors[3] = (jd3,m3,mm3,me3,mme3)

    def bb_bol(self, do_Kcorr=True, **kwargs):
        ''' 
        should correct with both host and milky ebv      
        will only use BCg
        '''        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if not 'cbb' in self.__dict__:
            if kwargs['verbose']: print ('calc BB colors first')
            return        
        filters = kwargs['cfilters']        
        self.mbolbb = dict()
        ebv = self.mkwebv + self.hostebv
        SN_distance = 10**((self.dm-25)/5)*3.086e24 # in cm
        for _ in self.cbb:
            jdl, Ll, eLl, Tl, eTl, Rl, eRl, L1l, eL1l = [], [], [], [], [], [], [], [], []
            for __ in range(len(self.cbb[_][filters[0]][0])):
                _jd = get_numpy(self.cbb[_]['t'])[__]                                
                jdl.append(_jd)
                # black body
                fluxes, efluxes, ws, bs = [], [], [], []
                for _filt in filters:
                    m, em = self.cbb[_][_filt]
                    _mag = get_numpy(m)[__] - ebv * Rf[_filt]
                    _mage = get_numpy(em)[__]
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

    def bb_colors(self, **kwargs):
        ''' define epochs for BB
        !!! haven't included any ebv value
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])       
        assert self.t0 > 2400000
        
        filters = kwargs['cfilters']
        assert len(filters) >= 3, 'at least 3 bands needed for blackbody'

        self.cbb = dict()
        self.cbb[1] = dict()                
        if 1 in kwargs['bb_copt']:
            # bin and match g and r
            self.lc = self.bin_df(self.lc, deltah = kwargs['color_thre'] )
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
            if 'fitcls' in self.__dict__ and 2 in kwargs['bb_copt']: # analytic sn lc fit
                p_fit = (_jds-self.t0) / (1+self.z)                
                __ = np.logical_and(p_fit > min(kwargs['fit_fitr']),
                                    p_fit < max(kwargs['fit_fitr']),) 
                xx,yy,yy1,yy2,_ = self.fitcls[_filt].predict(x_pred=p_fit[__],
                                    returnv=True, quant=kwargs['quantile'])
                mm = self.flux_to_mag(yy, dflux=None)
                _mm1 = self.flux_to_mag(yy1, dflux=None)
                _mm2 = self.flux_to_mag(yy1, dflux=None)
                mme = (abs(mm-_mm1) + abs(mm-_mm2))/2                
                self.cbb[2][_filt] = (mm,mme)
                    
            if 'gpcls' in self.__dict__ and 3 in kwargs['bb_copt']: # GP intepolation
                p_fit = _jds
                xx,yy,yye,ff = self.gpcls.predict(x_pred=p_fit, returnv=True,)
                __ = np.where(ff==_filt)
                xx = xx[__]
                yy = yy[__]
                yye = yye[__]
                mm = self.flux_to_mag(yy, dflux=None)
                mme = abs(self.flux_to_mag(yy, dflux=None)-
                          self.flux_to_mag(yy+yye, dflux=None))                          
                self.cbb[3][_filt] = (mm,mme)
        if len(self.cbb[2])>0: self.cbb[2]['t'] = _jds
        if len(self.cbb[3])>0: self.cbb[3]['t'] = _jds
                
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
            _lc = self.bin_df(self.lc, deltah = kwargs['color_thre'])
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
    
    def explosion_pl(self,**kwargs):
        ''' power law on early lc
        '''        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if kwargs['pl_type'] == 0: return
        elif kwargs['pl_type'] == 1: clobber = False
        elif kwargs['pl_type'] == 2: clobber = True
        elif kwargs['pl_type'] == 3: clobber = False
        else: return
                
        ncores = kwargs['ncores']
        if ncores is None: ncores = cpu_count() - 1

        assert self.t0 > 2400000, '!!!either input jdpeak or do GP first'
        # I need peak flux here
        assert len(self.fpeak) > 0, '!!!Error: do GP/fitting first...'
        
        tm_pl = kwargs['tm_pl']        
        rel_flux_cutoff = float(kwargs['rel_flux_cutoff'])
        fscale=float(kwargs['flux_scale'])
        
        df = self.lc.query('jdobs>@self.t0-@tm_pl and jdobs<@self.t0 and filter in %s'%
                           kwargs['pl_bands'])        
        if len(df)==0:  return
        
        # t: rest frame phase relative to r band peak
        t = (get_numpy(df['jdobs'])-self.t0)/(1+self.z)        
        f = get_numpy(df['Fmcmc'])
        f_unc = get_numpy(df['Fmcmc_unc'])
        f_zp = np.zeros_like(f)
        f_zp_unc = np.zeros_like(f_unc)
        bands = np.array(['']*len(t))
        for filt in kwargs['pl_bands']:            
            this_chip = np.where(df['filter'] == filt)
            if 'fit' in self.fpeak:
                tm = self.fpeak['fit'][filt][0]
            else:
                tm = self.fpeak['GP'][filt][0]
            f_zp[this_chip] = f[this_chip]/tm*fscale
            f_zp_unc[this_chip] = f_unc[this_chip]/tm*fscale
            bands[this_chip] = filt
        __ = np.where(f_zp < rel_flux_cutoff * fscale)
        t = t[__]
        f_zp = f_zp[__]
        f_zp_unc = f_zp_unc[__]
        bands = bands[__]
        
        self.plcls = fit_model(t,f_zp,f_zp_unc,filters=bands)        
        self.plcls.train(t_fl=18, opt_routine=kwargs['pl_routine'],
                    fit_mean='powerlaw', nwalkers=kwargs['nwalkers'],
                    nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'], ncores=ncores,
                    thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
                    mcmc_h5_file=kwargs['pl_cache']%(self.ztfid, kwargs['pl_routine']),
                    emcee_burnin=kwargs['emcee_burnin'], p0=None, datadir=kwargs['datadir'],
                    use_emcee_backend=kwargs['use_emcee_backend'], scipysamples=kwargs['scipysamples'],
                    clobber=clobber, verbose=kwargs['verbose'], sguess=kwargs['sguess']) 
        self.pl_theta, self.pl_rtheta = self.get_random_samples(self.plcls.samples,
                    self.plcls.lnprob, cachefile=kwargs['pl_cache1']%
                                (self.ztfid, kwargs['pl_routine']),
                    clobber=clobber, verbose=kwargs['verbose'],datadir=kwargs['datadir'],)
        
        # explosion epoch
        _, _1, _2 = self.plcls.get_par(filt=None, quant=kwargs['quantile'])               
        self.texp = [_[0], _1[0], _2[0]]

    def corner(self, figpath=None, show=True, which='fit', limit=0, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        assert which in ['fit', 'arnett', 'arnett_fitv', 'arnett_fitt', 'arnett_tail',
                         'gp', 'pl', 'tail', 'tail_fitt']
        quantiles = [kwargs['quantile'][0], kwargs['quantile'][-1]]
        
        if which == 'fit' and 'fitcls' in self.__dict__:            
            for filt in self.fitcls:
                samples = self.filter_samples(self.fitcls[filt].samples, self.fitcls[filt].lnprob, limit=limit)
                
                if kwargs['fit_method'] == 'bazin':
                    theta = np.zeros((np.shape(samples)[0],3))
                    for n, alpha  in enumerate([1,2,3]):                                                  
                        theta[:,n] = samples[:,alpha]
                    paramsNames = [r'$t_\mathrm{0,%s}$'%filt,
                                   r'$\tau_\mathrm{fall,%s}$'%filt, r'$\tau_{rise,%s}$'%filt]

                cfig = corner_hack(theta, labels=paramsNames,
                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                           show_titles=True, quantiles=quantiles,
                           title_fmt=".2f", title_kwargs={'fontsize':16},
                           plot_datapoints=True, plot_contours=True)
                if show: plt.show()
                
        if which == 'pl':           
            samples = self.filter_samples(self.plcls.samples, self.plcls.lnprob, limit=limit)
            
            _forder, nf = dict(), 0                
            for __filt in central_wavelengths:
                if __filt in kwargs['pl_bands']:
                    _forder[__filt] = nf
                    nf += 1
            nf = int( (samples.shape[1] - 1) / 4. )
            theta = np.zeros((np.shape(samples)[0],1+nf))
            theta[:,0] = samples[:,0]
            paramsNames = [r'$t_\mathrm{fl}$']
            for filt in kwargs['pl_bands']:
                nf = _forder[filt]                
                theta[:,1+nf] = samples[:,3+4*nf]
                paramsNames += [r"$\alpha_%s$"%filt]                
            cfig = corner_hack(theta, labels=paramsNames,
                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                           show_titles=True, quantiles=quantiles,
                           title_fmt=".2f", title_kwargs={'fontsize':16},
                           plot_datapoints=True, plot_contours=True)
            if show: plt.show()
            
        if which == 'arnett':
            samples = self.filter_samples(self.arnettcls.samples, self.arnettcls.lnprob, limit=limit)
            
            theta = np.zeros((np.shape(samples)[0],3))
            for n, alpha  in enumerate([0,1]):
                theta[:,n] = samples[:,alpha]
            theta[:,2] = Arnett_taum(samples[:,1])            
            paramsNames = [r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'M$_{ej}^{3/4}$ E$_{kin}^{-1/4}$']
            
            cfig = corner_hack(theta, labels=paramsNames,
                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                           show_titles=True, quantiles=quantiles,
                           title_fmt=".2f", title_kwargs={'fontsize':16},
                           plot_datapoints=True, plot_contours=True)
            
            if show: plt.show()

        if which == 'arnett_fitv':
            samples = self.filter_samples(self.arnettcls.samples, self.arnettcls.lnprob, limit=limit)
            
            theta = np.zeros((np.shape(samples)[0],4))
            for n, alpha  in enumerate([0,1,2]):
                theta[:,n] = samples[:,alpha]
            theta[:,3] = np.sqrt(2 * samples[:,1] / samples[:,2])                     
            paramsNames = [r'$M_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'v$_{m}$',]
            
            cfig = corner_hack(theta, labels=paramsNames,
                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                           show_titles=True, quantiles=quantiles,
                           title_fmt=".2f", title_kwargs={'fontsize':16},
                           plot_datapoints=True, plot_contours=True)
            
            if show: plt.show()

        if which == 'arnett_fitt':
            samples = self.filter_samples(self.arnettcls.samples, self.arnettcls.lnprob, limit=limit)
            
            theta = np.zeros((np.shape(samples)[0],3))
            for n, alpha  in enumerate([0,1,2]):
                theta[:,n] = samples[:,alpha]                    
            paramsNames = [r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'$t_\mathrm{fl}$',]
            
            cfig = corner_hack(theta, labels=paramsNames,
                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                           show_titles=True, quantiles=quantiles,
                           title_fmt=".2f", title_kwargs={'fontsize':16},
                           plot_datapoints=True, plot_contours=True)
            
            if show: plt.show()

        if which == 'arnett_tail':
            samples = self.filter_samples(self.arnettcls.samples, self.arnettcls.lnprob, limit=limit)
            
            theta = np.zeros((np.shape(samples)[0],5))
            for n, alpha  in enumerate([0,1,2,3,4]):
                theta[:,n] = samples[:,alpha]                    
            paramsNames = [r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'$t_{0}$', r'$t_\mathrm{fl}$', r'$t_{turn}$']
            
            cfig = corner_hack(theta, labels=paramsNames,
                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                           show_titles=True, quantiles=quantiles,
                           title_fmt=".2f", title_kwargs={'fontsize':16},
                           plot_datapoints=True, plot_contours=True)
            
            if show: plt.show()

        if which == 'tail':
            #samples = self.tailcls.samples
            samples = self.filter_samples(self.tailcls.samples, self.tailcls.lnprob, limit=limit)
            
            theta = np.zeros((np.shape(samples)[0],2))
            for n, alpha  in enumerate([0,1]):
                theta[:,n] = samples[:,alpha]
            paramsNames = [r'$M_\mathrm{Ni}$', r'$t_{0}$']
        
            cfig = corner_hack(theta, labels=paramsNames,
                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                           show_titles=True, quantiles=quantiles,
                           title_fmt=".2f", title_kwargs={'fontsize':16},
                           plot_datapoints=True, plot_contours=True)
            
            if show: plt.show()

        if which == 'tail_fitt':
            #samples = self.tailcls.samples
            samples = self.filter_samples(self.tailcls.samples, self.tailcls.lnprob, limit=limit)
            
            theta = np.zeros((np.shape(samples)[0],3))
            for n, alpha  in enumerate([0,1,2]):
                theta[:,n] = samples[:,alpha]                
            paramsNames = [r'$M_\mathrm{Ni}$', r'$t_{0}$', r'$t_{s}$',]
        
            cfig = corner_hack(theta, labels=paramsNames,
                           label_kwargs={'fontsize':16}, ticklabelsize=13,
                           show_titles=True, quantiles=quantiles,
                           title_fmt=".2f", title_kwargs={'fontsize':16},
                           plot_datapoints=True, plot_contours=True)
            
            if show: plt.show()
            
        if figpath is not None:
            cfig.savefig('{}/{}.png'.format(kwargs['datadir'], figpath), bbox_inches='tight')
                
    def run_gp(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
             
        if kwargs['gp_type'] == 0: return
        elif kwargs['gp_type'] == 1: clobber = False
        elif kwargs['gp_type'] == 2: clobber = True
        elif kwargs['gp_type'] == 3: clobber = False
        else: return
        
        if self.t0 > 0: # cut lc first
            pmin,pmax = min(kwargs['gp_fitr']), max(kwargs['gp_fitr'])
            '''
            !!! will cut lc permanently depending on GP fit range !!!
            '''
            self.lc = self.lc.query('jdobs<=@self.t0+@pmax and jdobs>=@self.t0+@pmin')
        
        lc = self.lc
        self.gpcls = fit_gp(np.array(lc['jdobs']), np.array(lc['Fmcmc']),
                    np.array(lc['Fmcmc_unc']), filters=np.array(lc['filter']))        
        self.gpcls.train(kernel=kwargs['kernel'], fix_scale=kwargs['fix_scale'],
              gp_mean=kwargs['gp_mean'], opt_routine=kwargs['gp_routine'],
              nwalkers=kwargs['nwalkers'], nsteps=kwargs['nsteps'],
              nsteps_burnin=kwargs['nsteps_burnin'], clobber=clobber,
              mcmc_h5_file=kwargs['gp_cache']%(self.ztfid,kwargs['gp_routine'],kwargs['gp_mean']),
              verbose=kwargs['verbose'], datadir=kwargs['datadir'])
        self.gpcls.predict()        
        
        for filt in np.unique(lc['filter']):
            __ = np.where(self.gpcls.f_pred==filt)
            xx = self.gpcls.x_pred[__]
            yy = self.gpcls.y_pred[__]
            yye = self.gpcls.y_prede[__]
            if 'GP' not in self.tpeak: self.tpeak['GP'] = dict()
            self.tpeak['GP'][filt] = xx[np.argmax(yy)]            
            if 'GP' not in self.fpeak: self.fpeak['GP'] = dict()
            self.fpeak['GP'][filt] = ( max(yy), yye[np.argmax(yy)] )
        
    def set_t0_withgp(self, filt='r', **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert 'GP' in self.tpeak, '!!! do GP first'                 
        self.t0 = self.tpeak['GP'][filt]

    def set_t0_withfit(self, filt='r', **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert 'fit' in self.tpeak, '!!! do fitting first'
        self.t0 = self.tpeak['fit'][filt]
        
    def run_fit(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        if kwargs['fit_type'] == 0: return
        elif kwargs['fit_type'] == 1: clobber = False
        elif kwargs['fit_type'] == 2: clobber = True
        elif kwargs['fit_type'] == 3: clobber = False
        else: return
        
        assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
        ncores = kwargs['ncores']
        if ncores is None: ncores = cpu_count() - 1
        
        pmin,pmax = min(kwargs['fit_fitr']), max(kwargs['fit_fitr'])        
        lc_df = self.lc.query('jdobs<=@self.t0+@pmax and jdobs>=@self.t0+@pmin')

        self.fitcls = dict()
        self.fit_theta = dict()
        self.fit_rtheta = dict()
        for filt in kwargs['fit_bands']:
            lc = lc_df.query('filter==@filt')
            if len(lc) <= 5: continue
            # fit on rest frame phase relative to r band peak            
            self.fitcls[filt] = fit_model((np.array(lc['jdobs'])-self.t0)/(1+self.z),
                    np.array(lc['Fmcmc']), np.array(lc['Fmcmc_unc']), filters=None)            
            self.fitcls[filt].train(t_fl=18, opt_routine=kwargs['fit_routine'],
                    fit_mean=kwargs['fit_method'], nwalkers=kwargs['nwalkers'],
                    nsteps=kwargs['nsteps'], nsteps_burnin=kwargs['nsteps_burnin'], ncores=ncores,
                    thin_by=kwargs['thin_by'], maxfev=kwargs['maxfev'],
                    mcmc_h5_file=kwargs['fit_cache']%(self.ztfid, filt, kwargs['fit_routine']),
                    emcee_burnin=kwargs['emcee_burnin'], p0=None, datadir=kwargs['datadir'],
                    use_emcee_backend=kwargs['use_emcee_backend'], scipysamples=kwargs['scipysamples'],
                    clobber=clobber, verbose=kwargs['verbose'],sigma=kwargs['fsigma'], sguess=kwargs['sguess'])
            self.fitcls[filt].predict(quant=kwargs['quantile'])
            self.fit_theta[filt], self.fit_rtheta[filt] = self.get_random_samples(self.fitcls[filt].samples,
                    self.fitcls[filt].lnprob, cachefile=kwargs['fit_cache1']%
                            (self.ztfid, filt, kwargs['fit_routine']),                    
                    clobber=clobber, verbose=kwargs['verbose'],datadir=kwargs['datadir'],)
            # for fits
            xx = self.fitcls[filt].x_pred
            yy = self.fitcls[filt].y_pred            
            __ = np.argmax(yy)
            yye = self.penc_to_errors(yy[__], self.fitcls[filt].y_pred1[__], self.fitcls[filt].y_pred2[__])                        
            
            if 'fit' not in self.tpeak: self.tpeak['fit'] = dict()
            self.tpeak['fit'][filt] = (xx[np.argmax(yy)]*(1+self.z))+self.t0
            if 'fit' not in self.fpeak: self.fpeak['fit'] = dict()
            self.fpeak['fit'][filt] = ( max(yy), yye )     

    @staticmethod
    def penc_to_errors(mean, perc_low, perc_up):        
        return abs(max(mean-perc_low, perc_up-mean))
    
    def clip_lc(self, **kwargs):
        ''' Removes outlier data points using sigma-clipping '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if kwargs['clipsigma'] is not None:
            clipsigma = float(kwargs['clipsigma'])
            __arr = None            
            for f in np.unique(self.lc['filter']):                
                _lc = self.lc.query('filter==@f')
                outlier_mask = sigma_clip(
                    data=_lc['Fmcmc'],
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
            
    def check_lc(self, **kwargs):
        ''' check Data quality '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        # how many epochs of photometry in either band.        
        dets = self.lc.query('mag<99')
        self.ndets = dict()
        for filt in kwargs['plot_bands']:
            self.ndets[filt] = len(dets.query('filter==@filt'))
            
        # how many color epochs (sampled within $\pm3$ days)
        self.ncolors = dict()
        for filt in kwargs['plot_bands']:
            for filt1 in kwargs['plot_bands']:
                if filt == filt1: continue
                if len(dets.query('filter==@filt'))==0 or len(dets.query('filter==@filt1'))==0:continue
                nc = 0
                for _jd in dets.query('filter==@filt')['jdobs']:
                    if min(abs(dets.query('filter==@filt1')['jdobs'] - _jd)) <= kwargs['jdthre']:
                        nc += 1
                self.ncolors['%s %s'%(filt, filt1)] = nc

        # Photometry available both before and after peak so that a peak can be determined 
        # to within an accuracy of $\pm 3$ days.
        self.peakphot = dict()
        peaksep = kwargs['peaksep']
        t0 = self.t0
        for filt in kwargs['plot_bands']:   
            filt_lc = self.lc.query('mag<99 and filter==@filt and jdobs>@t0-@peaksep and jdobs<@t0+@peaksep')
            self.peakphot[filt] = len(filt_lc)

        # how many early points
        self.earlypoints = dict()
        if 't0' in self.__dict__:
            level1 = float(kwargs['rel_flux_cutfrom'])
            level2 = float(kwargs['rel_flux_cutoff'])
            t0 = self.t0
            for filt in kwargs['plot_bands']:
                filt_lc = self.lc.query('mag<99 and filter==@filt and jdobs<@t0')
                if len(filt_lc)==0: continue
                if 'fit' in self.fpeak and filt in self.fpeak['fit']:
                    fmax = self.fpeak['fit'][filt][0]
                elif 'GP' in self.fpeak and filt in self.fpeak['GP']:
                    fmax = self.fpeak['GP'][filt][0]
                else:
                    fmax = max(filt_lc['Fmcmc'])
                fl1, fl2 = fmax*level1, fmax*level2            
                self.earlypoints[filt] = len(filt_lc.query('Fmcmc >= @fl1 and Fmcmc <= @fl2'))
        
    def on_sntype(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        ''' decide which region to measure v '''
        if self.sntype in ['SN Ib', 'SN IIb']:                           
            self.cw, self.region = 5876, [5300,6000]
            if kwargs['c10temp'] == 1:
                # c10 temp from Stritzinger et al 2018, table 2
                self.c10_temp = (0.790, 0.048)
            else:
                # Taddia et al 2015
                self.c10_temp = (0.64, 0.13)         
        elif self.sntype in ['SN Ic',]: 
            self.cw, self.region = 7774, [7200,7900]
            if kwargs['c10temp'] == 1:
                # c10 temp from Stritzinger et al 2018, table 2
                self.c10_temp = (0.851, 0.018)
            else:
                # Taddia et al 2015
                self.c10_temp = (0.64, 0.13)
        elif self.sntype in ['SN Ic-BL','SN Ic-BL?']: 
            self.cw, self.region = 7774, [6900,7900]            

    def init_fig(self, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        if self.axes is None:                     
            self.fig = plt.figure(num=1, clear=True, dpi=400,
                            tight_layout=False, constrained_layout=False,)                
            gs1 = self.fig.add_gridspec(nrows=3, ncols=2, left=.05,
                                        right=.95, wspace=0.05, hspace=0.3)
            self.ax = self.fig.add_subplot(gs1[0, 0]) # flux
            self.ax1 = self.fig.add_subplot(gs1[:2, 1]) # spectra
            self.ax2 = self.fig.add_subplot(gs1[1, 0]) # mag        
            self.ax3 = self.fig.add_subplot(gs1[2, 1]) # color
            self.ax4 = self.fig.add_subplot(gs1[2, 0]) # luminosity            
        elif not self.axes:
            self.fig, self.ax, self.ax1, self.ax2, self.ax3, self.ax4 = None, None, None, None, None, None
        else:
            assert len( self.axes ) == 6
            self.fig, self.ax, self.ax1, self.ax2, self.ax3, self.ax4 = self.axes
        
    def plot(self, show_title=True, **kwargs):        
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        ''' 
        - no extinction in ax1, flux plot
        - correct for milky way extinction in ax2, abs mag plot
        - correct for milky way extinction in ax3, color curve plot
        - correct for milky way and host extinction in ax4, bolometric luminosity plot
        '''        
        level = kwargs['rel_flux_cutfrom']
        if 'dm_rerr' in self.__dict__:  dm_rerr, dm_gerr = self.dm_rerr, self.dm_gerr
        else:                           dm_rerr, dm_gerr = 0, 0            
        fscale=float(kwargs['flux_scale'])
        allbands = kwargs['plot_bands']
        
        showbands, dmes, fms = [], [], []
        for filt in allbands:
            if len(self.lc.query('filter==@filt'))==0:continue
            try:
                dmes.append( self.dm_error(filt) )  # generate dm error
            except:
                print ('Warning: no dme for %s'%filt)
                dmes.append( 0 )
            if 'fit' in self.fpeak and filt in self.fpeak['fit']:
                fms.append( self.fpeak['fit'][filt][0] )            
            elif 'GP' in self.fpeak and filt in self.fpeak['GP']:
                fms.append( self.fpeak['GP'][filt][0] )
            else:
                fms.append( max(self.lc.query('filter==@filt')['Fmcmc']) )                 
            showbands.append(filt)
            
        t0 = self.t0
        axt0 = 2450000 # jd zero point in flux subplot        
        for filt, _dme, fm in zip(showbands, dmes, fms):
            ''' for each filter '''            
            if 'ax' not in self.__dict__ and 'ax2' not in self.__dict__ and self.ax is None and self.ax2 is None:
                continue
            if filt in Rf:  R = Rf[filt]
            else:   R = 3.1            
            if 'fitcls' in self.__dict__ and kwargs['fit_type'] != 3 and filt in self.fitcls:
                # analytic sn lc fit                

                # fit peaks
                pmin,pmax = min(kwargs['fit_plotr']), max(kwargs['fit_plotr'])                
                p_fit = np.arange(pmin, pmax, .1)

                # best sample                         
                theta, random_theta = self.fit_theta[filt], self.fit_rtheta[filt]
                
                # on ax
                model_flux = bazin(p_fit, *theta[:5])
                if 'ax' in self.__dict__ and self.ax is not None:
                    self.ax.plot(p_fit*(1+self.z) + t0-axt0,
                      model_flux/fm*fscale+ys[filt], alpha=kwargs['alphabest'], **PROP2[filt])
                    
                # on ax2
                __ = np.where(model_flux > fm*level)
                model_mag = self.flux_to_mag(model_flux[__], dflux=None)
                if  'ax2' in self.__dict__ and self.ax2 is not None:
                    self.ax2.plot(p_fit[__], model_mag-self.dm-R*self.mkwebv,
                                  alpha=kwargs['alphabest'], **PROP2[filt])

                # if show peak?
                #if 'ax' in self.__dict__ and self.ax is not None:
                #    __ = np.argmax( model_flux )
                #    self.ax.axvline(p_fit[__] *(1+self.z) + t0 - axt0, **PROP2p[filt])                    

                # other samples
                for theta_samp in random_theta:
                    # ax
                    model_flux = bazin(p_fit, *theta_samp[:5])
                    if 'ax' in self.__dict__ and self.ax is not None:
                        self.ax.plot(p_fit*(1+self.z) + t0 - axt0,
                            model_flux/fm*fscale+ys[filt], alpha=kwargs['alphasample'], **PROP2[filt])
                    # ax2
                    __ = np.where(model_flux > fm*level)
                    model_mag = self.flux_to_mag(model_flux[__], dflux=None)
                    if 'ax2' in self.__dict__ and self.ax2 is not None:
                        self.ax2.plot(p_fit[__], model_mag-self.dm-R*self.mkwebv,
                                      alpha=kwargs['alphabest'], **PROP2[filt])
                            
            if 'gpcls' in self.__dict__ and kwargs['gp_type'] != 3: # Gaussian process
                if filt in self.gpcls.f_pred:                                        
                    __ = np.where(self.gpcls.f_pred==filt)
                    xx = self.gpcls.x_pred[__]
                    yy = self.gpcls.y_pred[__]
                    yye = self.gpcls.y_prede[__] * .5

                    if kwargs['gp_plotr'] is not None:
                        pmin,pmax = min(kwargs['gp_plotr']), max(kwargs['gp_plotr'])                    
                        __ = np.logical_and(xx>t0+pmin, xx<t0+pmax)                    
                        xx = xx[__]
                        yy = yy[__]
                        yye = yye[__] 
                    # ax
                    if 'ax' in self.__dict__ and self.ax is not None:
                        self.ax.plot(xx - axt0, yy/fm*fscale+ys[filt],
                                     alpha=kwargs['alphabest'], **PROP3[filt])
                        self.ax.fill_between(xx - axt0,
                            (yy-yye) /fm*fscale+ys[filt], (yy+yye) /fm*fscale+ys[filt],
                            alpha=kwargs['alphasample'], **PROP3[filt])
                    # ax2
                    __ = np.where(yy > fm*level)
                    model_mag, model_mage = self.flux_to_mag(yy[__], dflux=yye[__])
                    if 'ax2' in self.__dict__ and self.ax2 is not None:
                        self.ax2.plot((xx[__]-t0)/(1+self.z), model_mag-self.dm-R*self.mkwebv,
                                      alpha=kwargs['alphabest'], **PROP3[filt])
                        self.ax2.fill_between((xx[__]-t0)/(1+self.z),
                                model_mag-model_mage-self.dm-R*self.mkwebv,
                                model_mag+model_mage-self.dm-R*self.mkwebv,
                                alpha=kwargs['alphasample'], **PROP3[filt])
            
            if 'plcls' in self.__dict__ and kwargs['pl_type'] != 3: # power law on explosion            
                _forder, nf = dict(), 0                
                for __filt in central_wavelengths:
                    if __filt in kwargs['pl_bands']:
                        _forder[__filt] = nf
                        nf += 1
                if filt in _forder:
                    nf = _forder[filt]                                    
                    theta, random_theta = self.pl_theta, self.pl_rtheta                
                    adjust_samp = [theta[0],theta[1+4*nf],theta[2+4*nf],theta[3+4*nf]]
                    t_post = np.linspace(theta[0], kwargs['pl_plotmax'], 1000)                
                    t_pre = np.linspace((min(self.lc['jdobs']) - t0) / (1+self.z), theta[0], 1000)
                    model_flux = powerlaw_post_baseline(t_post,
                    adjust_samp[2]*10**(-adjust_samp[3]), adjust_samp[0], adjust_samp[3])
                    if 'ax' in self.__dict__ and self.ax is not None:
                        self.ax.plot(t_post*(1+self.z) + t0 - axt0,
                        model_flux + ys[filt], alpha=kwargs['alphabest'], **PROP4[filt])
                        self.ax.plot(t_pre*(1+self.z) + t0 - axt0,
                        np.zeros_like(t_pre) + ys[filt], alpha=kwargs['alphabest'], **PROP4[filt])
                                            
                    for theta_samp in random_theta:
                        nf = _forder[filt]                    
                        adjust_samp = [theta_samp[0],theta_samp[1+4*nf],
                                   theta_samp[2+4*nf],theta_samp[3+4*nf]]                    
                        # ax
                        t_post = np.linspace(theta_samp[0], kwargs['pl_plotmax'], 1000)
                        t_pre = np.linspace((min(self.lc['jdobs']) - t0) / (1+self.z), theta_samp[0], 1000)
                        model_flux = powerlaw_post_baseline(t_post,
                            adjust_samp[2]*10**(-adjust_samp[3]), adjust_samp[0], adjust_samp[3])
                        if 'ax' in self.__dict__ and self.ax is not None:
                            self.ax.plot(t_post*(1+self.z) + t0 - axt0,
                                model_flux + ys[filt], alpha=kwargs['alphasample'], **PROP4[filt])
                            self.ax.plot(t_pre*(1+self.z) + t0 - axt0,
                                np.zeros_like(t_pre) + ys[filt], alpha=kwargs['alphasample'], **PROP4[filt])
            
            # data points:                  
            ''' ax flux plot '''
            # normalize flux item in ax
            if 'ax' in self.__dict__ and self.ax is not None:
                __lc = self.lc.query('filter==@filt')
                self.ax.errorbar(__lc['jdobs']-axt0, __lc['Fmcmc']/fm*fscale+ys[filt],
                             yerr=__lc['Fmcmc_unc']/fm*fscale, alpha=kwargs['alphabest'], **PROP1f[filt])
                
            ''' ax2 mag plot '''
            if 'ax2' in self.__dict__ and self.ax2 is not None:
            # detections            
                __lc = self.lc.query('filter==@filt and mag<99 and emag<99')           
                # set t to rest frame phase relative to peak                
                self.ax2.errorbar((__lc['jdobs']-t0)/(1+self.z), __lc['mag']-self.dm-R*self.mkwebv, 
                      yerr=np.sqrt(__lc['emag']**2+_dme**2), alpha=kwargs['alphabest'], **PROP1[filt])
                # upper limits
                __lc = self.lc.query('filter==@filt and mag==99 and jdobs<@t0')
                if len(__lc) > 0:
                    self.ax2.plot((__lc['jdobs']-t0)/(1+self.z), __lc['limmag']-self.dm-R*self.mkwebv,
                                  alpha=kwargs['alphasample'], **PROP1l[filt])

        # colors
        if 'colors' in self.__dict__ and 'ax3' in self.__dict__ and self.ax3 is not None:            
            for _ in kwargs['copt']:
                if _ in self.colors:
                    xx = (self.colors[_][0]-t0)/(1+self.z)
                    yy = self.colors[_][1]-self.colors[_][2]
                    yye = np.sqrt(self.colors[_][3]**2 + self.colors[_][4]**2)
                    ''' correct for milky way extinction for g-r '''
                    yy += self.mkwebv * (Rf['r'] - Rf['g'])
                    self.ax3.errorbar(xx, yy, yerr=yye, **PROP5['gr%s'%_])
                    if _==1:
                        self.ax3.set_ylim([min(self.colors[1][1]-self.colors[1][2])-.1,
                                           max(self.colors[1][1]-self.colors[1][2])+.1])
            # template color
            if 'c10_temp' in self.__dict__:
                self.ax3.errorbar(10, self.c10_temp[0], yerr=self.c10_temp[1], marker='+', color='g')
                
        # bolometric LC
        if 'mbol' in self.__dict__ and 'ax4' in self.__dict__ and self.ax4 is not None and 1 in kwargs['bolopt']:            
            for _ in kwargs['copt']:
                if _ in self.mbol:                    
                    xx = (get_numpy(self.mbol[_][0])-t0)/(1+self.z)
                    yy = get_numpy(self.mbol[_][1])
                    yye = get_numpy(self.mbol[_][2]) 
                    self.ax4.errorbar(xx,yy,yerr=yye,**PROP6['lyman_data_%s'%_])
                    
        if 'mbolbb' in self.__dict__ and 'ax4' in self.__dict__ and self.ax4 is not None and 2 in kwargs['bolopt']:
            for _ in kwargs['bb_copt']:
                if _ in self.mbolbb:                    
                    xx = (get_numpy(self.mbolbb[_][0])-t0)/(1+self.z)
                    yy = get_numpy(self.mbolbb[_][1])
                    yye = get_numpy(self.mbolbb[_][2]) 
                    self.ax4.errorbar(xx,yy,yerr=yye,**PROP6['bb_data_%s'%_])
                    
                    yy = get_numpy(self.mbolbb[_][7])
                    yye = get_numpy(self.mbolbb[_][8]) 
                    self.ax4.errorbar(xx,yy,yerr=yye,**PROP6['bb_data_%sf'%_])
                
        if 'arnettcls' in self.__dict__ and 'ax4' in self.__dict__ and self.ax4 is not None:                 
            texp = self.texp[1] 
            theta, random_theta = self.arnett_theta, self.arnett_rtheta
            
            t_ = np.linspace(min(kwargs['Arnett_plotr']), max(kwargs['Arnett_plotr']), 1000)
            t_ -= texp
            __ = np.where(t_ > 1)
            t_ = t_[__]
            model_flux = Arnett_fit(t_, *theta[:2])            
            self.ax4.plot(t_+texp, model_flux, alpha=kwargs['alphabest'], **PROP6['fit'], label='Arnett fits')
            for theta_samp in random_theta:            
                model_flux = Arnett_fit(t_, *theta_samp[:2])
                self.ax4.plot(t_+texp, model_flux, alpha=kwargs['alphasample'], **PROP6['fit'])

        if 'tailcls' in self.__dict__ and 'ax4' in self.__dict__ and self.ax4 is not None:                         
            texp = self.texp[1]           
            theta, random_theta = self.tail_theta, self.tail_rtheta
            
            t_ = np.linspace(min(kwargs['Tail_plotr']), max(kwargs['Tail_plotr']), 1000)
            t_ -= texp
            __ = np.where(t_ > 1)
            t_ = t_[__]
            model_flux = tail_nickel(t_, *theta[:2])            
            self.ax4.plot(t_+texp, model_flux, alpha=kwargs['alphabest'], **PROP7['fit'], label='Tail fits')
            for theta_samp in random_theta:              
                model_flux = tail_nickel(t_, *theta_samp[:2])                
                self.ax4.plot(t_+texp, model_flux, alpha=kwargs['alphasample'], **PROP7['fit'])

        # spectra
        if 'specls' in self.__dict__ and 'ax1' in self.__dict__ and self.ax1 is not None:
            for phase in self.specls:
                _spec = self.specls[phase]
                if _spec.feature_fit is None: continue
                _spec.ax.step( _spec.rest_wave, _spec.norm_flux+_spec.ys, **PROP8['data'])
                _spec.ax.step( _spec.rest_wave, _spec.bin_flux+_spec.ys, **PROP8['bindata'])
                _spec.ax.text(max(_spec.region), _spec.ys, '%s %s d' % (_spec.instru, _spec.phase), color='k')
                _spec.ax.axvline(_spec.cw, color='k', ls='--')
                
                if 'peaks' in _spec.__dict__:
                    xln,yln = _spec.peaks
                    if xln is not None:
                        _spec.ax.plot(xln, yln+_spec.ys,'r|',markersize=8, markeredgewidth=2)

                if _spec.feature_fit is not None:                    
                    w, fit, fit1, fit2 = _spec.feature_fit.predict(w_pred=_spec.rest_wave, quant=kwargs['quantile'])
                    if w is not None: 
                        _spec.ax.plot( w, fit+_spec.ys, **PROP8['fit'])
                        _spec.ax.plot( w, fit1+_spec.ys, **PROP8['fit'])
                        _spec.ax.plot( w, fit2+_spec.ys, **PROP8['fit'])
            
                    q, q1, q2 = _spec.feature_fit.get_par(quant=kwargs['quantile'])
                    if q is not None:
                        _spec.ax.vlines(x=q[1], ymin=_spec.ys+_spec.fnorm*.2,ymax=_spec.ys+_spec.fnorm*.8,color='k',ls='--',)
                        _spec.ax.vlines(x=q1[1]-q1[2]/2., ymin=_spec.ys+_spec.fnorm*.2,ymax=_spec.ys+_spec.fnorm*.8,color='g',ls='--',)
                        _spec.ax.vlines(x=q2[1]+q2[2]/2., ymin=_spec.ys+_spec.fnorm*.2,ymax=_spec.ys+_spec.fnorm*.8,color='g',ls='--',)
            
        #  other settings
        x1, x2 = None, None
        if kwargs['ax_xlim'] is not None and 't0' in self.__dict__:
            x1,x2 = min(kwargs['ax_xlim']), max(kwargs['ax_xlim'])
        if 'ax' in self.__dict__ and self.ax is not None and show_title:
            try:  self.ax.set_title('%s %s\n%s z=%.2f'%
                (self.iauid,self.ztfid,self.sntype,self.z),fontsize=12)
            except:
                self.ax.set_title('%s'%(self.ztfid), fontsize=12)
        if 'ax' in self.__dict__ and self.ax is not None:
            if kwargs['ax_ylim'] is not None:
                self.ax.set_ylim(kwargs['ax_ylim'])
            if x1 is not None and self.t0>2450000:
                self.ax.set_xlim([self.t0+x1-axt0,self.t0+x2-axt0])
            if 'texp' in self.__dict__ :                            
                self.ax.axvline(t0-axt0+self.texp[1], ls='--', label='texp', color='k')
                self.ax.axvline(t0-axt0+self.texp[1], ls=':', color='k')
                self.ax.axvline(t0-axt0+self.texp[1], ls=':', color='k')
            if 't0' in self.__dict__ and self.t0>2450000:
                self.ax.axvline(t0-axt0, ls='--', label='t0', color='r')                
            self.ax.set_xlabel('JD - %s (d)'%axt0,fontsize=12)
            self.ax.set_ylabel('Flux + $\mathrm{offset}$',fontsize=12)
            self.ax.legend(fontsize=10, frameon=False)

        if 'ax1' in self.__dict__ and self.ax1 is not None and 'region' in self.__dict__:
            ax1 = self.ax1.twiny()            
            _ticks = np.arange(min(self.region),max(self.region)+100,100)
            self.ax1.set_xticks( _ticks )
            self.ax1.set_xticklabels(_ticks, rotation=30)
            ax1.set_xlim( self.region )
            _tickl = ['%.1f'%_ for _ in (_ticks/self.cw-1)*299792.458/1000.]
            ax1.set_xticks(_ticks)
            ax1.set_xticklabels(_tickl, rotation=30)
            self.ax1.set_xlabel('Rest Wavelength ($\AA$)', fontsize=12)
            ax1.set_xlabel('velocity ($10^{3}~km~s^{-1}$)', fontsize=12)
            self.ax1.set_yticks([])
            self.ax1.tick_params("both", direction='in',
                right=True, labelright=True, left=False, labelleft=False)
            self.ax1.set_xlim([min(self.region),max(self.region)])
         
        if 'ax2' in self.__dict__ and self.ax2 is not None:
            self.ax2.invert_yaxis()
            #if 'texp' in self.__dict__ :                             
            #    self.ax2.axvline(self.texp[1]/(1+self.z), ls=':', label=None, color='k')
            if kwargs['ax2_ylim'] is not None:
                self.ax2.set_ylim(kwargs['ax2_ylim'])
            if x1 is not None and self.t0>2450000:
                self.ax2.set_xlim([x1/(1+self.z),x2/(1+self.z)])
            self.ax2.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=12)
            self.ax2.set_ylabel('M$_{abs}$ (mag)',fontsize=12)
            self.ax2.legend(fontsize=10, frameon=False)
                       
        if 'ax3' in self.__dict__ and self.ax4 is not None:
            self.ax3.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=12)
            ax33=self.ax3.twinx()
            ax33.set_ylabel('g-r (mag)',fontsize=12,labelpad=24)
            ax33.set_yticks([])
            if kwargs['ax3_ylim'] is not None:
                self.ax3.set_ylim(kwargs['ax3_ylim'])
            if x1 is not None and self.t0>2450000:
                self.ax3.set_xlim([x1/(1+self.z),x2/(1+self.z)])
            self.ax3.tick_params("both", direction='in',
                    right=True, labelright=True, left=False, labelleft=False)
            self.ax3.legend(fontsize=10, frameon=False)
            
        if 'ax4' in self.__dict__ and self.ax4 is not None:
            self.ax4.set_yscale('log')
            self.ax4.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=12)
            self.ax4.set_ylabel('L$_{bol}$ (erg $s^{-1}$)',fontsize=12)    
            
            if kwargs['ax4_ylim'] is not None:
                self.ax4.set_ylim(kwargs['ax4_ylim'])
            if x1 is not None and self.t0>2450000:
                self.ax4.set_xlim([x1/(1+self.z),x2/(1+self.z)])                        
            self.ax4.legend(fontsize=10, frameon=False)
            
        if 'fig' in self.__dict__ and self.fig is not None:
            self.fig.set_size_inches(kwargs['figsize'][0], kwargs['figsize'][1])
            #self.fig.tight_layout()
            #self.fig.show()
            if kwargs['figpath'] is not None:
                if '%s' in kwargs['figpath']:
                    self.fig.savefig('{}/{}'.format(kwargs['datadir'], kwargs['figpath']%self.ztfid),
                                     dpi=400, bbox_inches='tight')
                else:
                    self.fig.savefig('{}/{}'.format(kwargs['datadir'], kwargs['figpath']),
                                     dpi=400, bbox_inches='tight')

    def show(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        figpath = '{}/{}'.format(kwargs['datadir'], kwargs['figpath']%self.ztfid)

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
        c=3e5
        vpel=150.
        dz=vpel/(c*self.z)*self.z                         
        dmag=dz*5/np.log(10)*((1+self.z)/(self.z*(1+0.5*self.z)))
    
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
    def mag_to_flux(mag, magerr=None, units='zp', zp=25., wavelength=None):
        mag = get_numpy(mag)            
        magerr = get_numpy(magerr)           
        ''' with flux/mag zerop '''
        if units not in ['zp', 'phys']:
            raise ValueError("units must be 'zp' or 'phys'")
        elif units == 'zp':
            if zp is None: raise ValueError("zp must be float or array if units == 'zp'")
            flux = 10**(-(mag-zp)/2.5)
        else:
            if wavelength is None: raise ValueError("wavelength must be float or array if units == 'phys'")
            flux = 10**(-(mag+2.406)/2.5) / wavelength**2                
        __ = np.where( mag==99 )
        flux[__] = min(flux)        
        if magerr is None: return flux        
        dflux = np.abs(flux*(-magerr/2.5*np.log(10)))
        dflux[__] = 0.1
        return flux, dflux
        
    @staticmethod
    def flux_to_mag(flux, dflux=None, units='zp', zp=25., wavelength=None):
        flux = get_numpy(flux)        
        dflux = get_numpy(dflux)          
        """ Converts fluxes (erg/s/cm2/A) into AB magnitudes with flux/mag zerop """
        if units not in ['zp', 'phys']: raise ValueError("units must be 'zp' or 'phys'")
        elif units == 'zp':
            if zp is None: raise ValueError("zp must be float or array if units == 'zp'")
            wavelength = 1.
        else:
            if wavelength is None: raise ValueError("wavelength must be float or array if units == 'phys'")
            zp = -2.406                    
        mag = -2.5*np.log10(flux*wavelength**2) + zp
        if dflux is None: return mag
        __ = np.where( flux/dflux < 1 )
        mag[__] = np.nan        
        dmag = np.abs( -2.5/np.log(10) * dflux / flux )
        dmag[__] = np.nan   
        return mag, dmag

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
        rows = __df.shape[0]
        __df1, __n = __df.values[0], 1
        for __nk, _kk in enumerate(range(rows)):
            if __nk==0:continue
            __df1 += np.array(__df.values[_kk])
            __n += 1
        _arr = [__df2[0:int(len(__df2)/__n)] if type(__df2) is str else float(__df2)/__n for __df2 in __df1]
        return _arr

    def get_random_samples(self, samples, lnprob, limit=0.5, datadir=None,
                           cachefile=None, clobber=False, verbose=False, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])        
        
        if cachefile is not None and datadir is not None:
            cache = '%s/%s.npz'%(datadir, cachefile)
        else:
            cache = None
            
        if cache is not None:
            if os.path.exists(cache) and not clobber:
                if verbose: print ('reload from %s'%cache)
                return np.load(cache)['theta'], np.load(cache)['theta_samp']
            
        # best sample
        theta = samples[np.argmax(lnprob)]
        
        # other samples
        _samples = self.filter_samples(samples, lnprob, limit=limit)
        nsamples = min(len(_samples), kwargs['plotsamples'])
        
        theta_samp = []        
        for samp_num in np.random.choice(range(len(_samples)),nsamples, replace=False):
            theta_samp.append( _samples[samp_num] )
            
        if cache is not None:
            np.savez(cache, theta=theta, theta_samp=theta_samp)
            if verbose: print ('generate %s'%cache)
            
        return theta, theta_samp

    def filter_samples(self, samples, lnprob, limit=0., cache=None, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        # other samples
        thre = min(lnprob) + (max(lnprob) - min(lnprob)) * limit        
        theta_samp = []
        for nn, samp_num in enumerate(samples):           
            if lnprob[nn] >= thre: theta_samp.append( samp_num )        
        return np.array(theta_samp)        
    
class plotter:
    
    # Static version info
    version = 1.0
    
    def __init__(self, ztfmultiple, verbose=False):                
        self.data = ztfmultiple.data
        self.meta = ztfmultiple.meta
        self.verbose = verbose
        if not 'dset' in self.__dict__: self.dset = {}
        
    def init_fig(self, dpi=100, figsize=[6,6], 
                 tight_layout=False, constrained_layout=False):
        self.fig, self.ax = plt.subplots(1, 1, num=1, clear=True, 
                dpi=dpi, tight_layout=tight_layout, 
                constrained_layout=constrained_layout,
                figsize=figsize)
        
    def init_hist_axes(self, pad=0.1, labelbottom=False, labelleft=False):     
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        assert self.ax is not None
        # create new axes on the right and on the top of the current axes
        divider = make_axes_locatable(self.ax)
        # below height and pad are in inches
        self.ax_histx = divider.append_axes("top", 1.2, pad=pad, sharex=self.ax)
        self.ax_histy = divider.append_axes("right", 1.2, pad=pad, sharey=self.ax)
        # make some labels invisible
        self.ax_histx.xaxis.set_tick_params(labelbottom=labelbottom)
        self.ax_histy.yaxis.set_tick_params(labelleft=labelleft)

    def add_subset(self, syntax=None):
        _meta = self.meta
        if syntax is not None: _meta = _meta.query(syntax)
        self.dset[syntax] = {ztfid:self.data[ztfid] if ztfid in list(_meta.index) else None for ztfid in self.data}
        
    def show2d(self, k1, k2, syntax=None, fontsize=12, labelpad=12, **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        x, y = [], []
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                x.append(self.dset[syntax][ztfid].__dict__[k1])
                y.append(self.dset[syntax][ztfid].__dict__[k2])
        self.ax.plot( x, y, **kwargs)
        self.ax.set_xlabel(k1, fontsize=fontsize, labelpad=labelpad )
        self.ax.set_ylabel(k2, fontsize=fontsize, labelpad=labelpad )
        self.ax.set_xlim([min(x)/1.1, max(x)*1.1])
        self.ax.set_ylim([min(y)/1.1, max(y)*1.1])
        
    def add_hist(self, k1, k2, syntax=None, nbinx = 10, nbiny = 10, 
                 xticks=None, yticks=None, **kwargs):        
        assert 'ax_histx' in self.__dict__ and 'ax_histy' in self.__dict__, 'init_hist_axes() first'
        self.add_subset(syntax=syntax)
        x, y = [], []
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                x.append(self.dset[syntax][ztfid].__dict__[k1])
                y.append(self.dset[syntax][ztfid].__dict__[k2])        
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
    
    def showlc(self, syntax=None, cond='mag<99 and filter=="r"', 
               fontsize=12, labelpad=12, **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                lc = self.dset[syntax][ztfid].lc
                t0 = self.dset[syntax][ztfid].t0
                z = self.dset[syntax][ztfid].z
                dm = self.dset[syntax][ztfid].dm
                if cond is not None: lc = lc.query(cond)
                self.ax.plot( (lc['jdobs']-t0)/(1+z), lc['mag']-dm, **kwargs )
        self.ax.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$', 
                           fontsize=fontsize, labelpad=labelpad )
        self.ax.set_ylabel('M$_{abs}$ (mag)', fontsize=fontsize, labelpad=labelpad )
        
    def showcolor(self, syntax=None, copt=[1], fontsize=12, labelpad=12, **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                if 'colors' not in self.dset[syntax][ztfid].__dict__:
                    if self.verbose: print ('Warning: no colors for %s'%ztfid)
                    continue
                colors = self.dset[syntax][ztfid].colors
                t0 = self.dset[syntax][ztfid].t0
                z = self.dset[syntax][ztfid].z
                for _copt in copt:
                    if _copt in colors:
                        jd,g,r,ge,re = colors[_copt]
                        self.ax.plot( (jd-t0)/(1+z), g-r, **kwargs )
        self.ax.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$', 
                           fontsize=fontsize, labelpad=labelpad )
        self.ax.set_ylabel('g-r (mag)', fontsize=fontsize, labelpad=labelpad )
        
    def showlcbol(self, syntax=None, copt=[1], fontsize=12, labelpad=12, **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                if 'mbol' not in self.dset[syntax][ztfid].__dict__:
                    if self.verbose: print ('Warning: no mbol for %s'%ztfid)
                    continue
                mbol = self.dset[syntax][ztfid].mbol
                t0 = self.dset[syntax][ztfid].t0
                z = self.dset[syntax][ztfid].z
                for _copt in copt:
                    if _copt in mbol:
                        self.ax.plot( (mbol[_copt][0]-t0)/(1+z), mbol[_copt][1], **kwargs )
        self.ax.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',
                           fontsize=fontsize, labelpad=labelpad )
        self.ax.set_ylabel('L$_{bol}$ (erg $s^{-1}$)', fontsize=fontsize, labelpad=labelpad )     
        
    def showfit(self, syntax=None, filt='r', funit='mag',
                x_pred=None, step = 1, quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                if 'fitcls' not in self.dset[syntax][ztfid].__dict__:
                    if self.verbose: print ('Warning: no fitcls for %s'%ztfid)
                    continue
                if filt not in self.dset[syntax][ztfid].fitcls:
                    if self.verbose: print ('Warning: no filter %s for %s fitcls'%(filt,ztfid))
                    continue
                z = self.dset[syntax][ztfid].z
                x, y, y1, y2, w = self.dset[syntax][ztfid].fitcls[filt].predict(x_pred=x_pred, 
                                        step=step, clobber=False, returnv=True, quant=quant)
                if funit == 'mag':
                    mm = self.dset[syntax][ztfid].flux_to_mag(y, dflux=None)
                    self.ax.plot( x/(1+z), mm-self.data[ztfid].dm,**kwargs )
                elif funit == 'flux':
                    if 'fpeak' not in self.dset[syntax][ztfid].__dict__:
                        if self.verbose: print ('Warning: no fmax for %s'%(ztfid))
                        continue
                    fmax = self.dset[syntax][ztfid].fpeak['fit'][filt][0]
                    self.ax.plot( x/(1+z), y/fmax,**kwargs )
    
    def showgp(self, syntax=None, filt='r', funit='mag',
                x_pred=None, step = 1, quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                if 'gpcls' not in self.dset[syntax][ztfid].__dict__:
                    if self.verbose: print ('Warning: no fitcls for %s'%ztfid)
                    continue
                if filt not in self.dset[syntax][ztfid].gpcls.f_pred:
                    if self.verbose: print ('Warning: no filter %s for %s gpcls'%(filt,ztfid))
                    continue
                z = self.dset[syntax][ztfid].z
                t0 = self.dset[syntax][ztfid].t0
                __ = np.where(self.dset[syntax][ztfid].gpcls.f_pred==filt)
                x = self.dset[syntax][ztfid].gpcls.x_pred[__]
                y = self.dset[syntax][ztfid].gpcls.y_pred[__]
                ye = self.dset[syntax][ztfid].gpcls.y_prede[__]
                if funit == 'mag':
                    mm = self.dset[syntax][ztfid].flux_to_mag(y, dflux=None)
                    self.ax.plot( (x-t0)/(1+z), mm-self.data[ztfid].dm,**kwargs )
                elif funit == 'flux':
                    if 'fpeak' not in self.dset[syntax][ztfid].__dict__:
                        if self.verbose: print ('Warning: no fmax for %s'%(ztfid))
                        continue
                    fmax = self.dset[syntax][ztfid].fpeak['GP'][filt][0]
                    self.ax.plot( (x-t0)/(1+z), y/fmax,**kwargs )
                
    def showpl(self, syntax=None, filt='r', x1=-50, x2=0, limit=.9, **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                if 'plcls' not in self.dset[syntax][ztfid].__dict__:
                    if self.verbose: print ('Warning: no plcls for %s'%ztfid)
                    continue
                _forder, nf = dict(), 0                
                for __filt in central_wavelengths:
                    if __filt in self.dset[syntax][ztfid].kwargs['pl_bands']:
                        _forder[__filt] = nf
                        nf += 1
                if filt in _forder:
                    nf = _forder[filt]                                    
                    theta = self.dset[syntax][ztfid].filter_samples(self.dset[syntax][ztfid].plcls.samples, 
                                                self.dset[syntax][ztfid].plcls.lnprob, limit=limit)[0]              
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
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                if 'arnettcls' not in self.dset[syntax][ztfid].__dict__:
                    if self.verbose: print ('Warning: no arnettcls for %s'%ztfid)
                    continue
                texp = self.dset[syntax][ztfid].texp
                x, y, y1, y2, w = self.dset[syntax][ztfid].arnettcls.predict(x_pred=x_pred, 
                                        step=step, clobber=False, returnv=True, quant=quant)
                if len(self.dset[syntax][ztfid].arnettcls.get_par()[0]) == 3:
                    x -= texp[1]
                self.ax.plot( x,y,**kwargs )
        
    def showtail(self, syntax=None, x_pred=None, step = 1, 
                   quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                if 'tailcls' not in self.dset[syntax][ztfid].__dict__:
                    if self.verbose: print ('Warning: no tailcls for %s'%ztfid)
                    continue
                texp = self.dset[syntax][ztfid].texp
                x, y, y1, y2, w = self.dset[syntax][ztfid].tailcls.predict(x_pred=x_pred, 
                                        step=step, clobber=False, returnv=True, quant=quant)
                self.ax.plot( x,y,**kwargs )
                
    def showvelocity(self, syntax=None, quant=[0.16, 0.50, 0.84], **kwargs):
        assert 'ax' in self.__dict__, 'init_fig() first'
        self.add_subset(syntax=syntax)
        for ztfid in self.dset[syntax]:
            if self.dset[syntax][ztfid] is not None:
                if 'specls' not in self.dset[syntax][ztfid].__dict__:
                    if self.verbose: print ('Warning: no specls for %s'%ztfid)
                    continue
                x,y=[],[]
                for phase in self.dset[syntax][ztfid].specls:
                    _spec = self.dset[syntax][ztfid].specls[phase]
                    if _spec.velocity is None: continue
                    x.append( float(_spec.phase) )
                    y.append( float(_spec.velocity[1]) )
                self.ax.plot( x,y, **kwargs )
        
if __name__ == "__main__":
    
    cls = ztfmultiple(**vars(args))    
    cls.run(axes=None)
