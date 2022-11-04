#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/sne.py
# Author            : syang <sheng.yang@astro.su.se>
# Date              : 12.11.2019
# Last Modified Date: 11.02.2022
# Last Modified By  : syang <sheng.yang@astro.su.se>

from __future__ import print_function
import os, sys, io, requests, emcee, corner, random, urllib, urllib.request,\
    shlex, glob, shutil, subprocess, time, json, math, re, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit
from scipy import integrate
from scipy.integrate import simps
from scipy.interpolate import interp1d
from astropy.cosmology import Planck13 as cosmo
from astropy.time import Time
import astropy.units as u
from astropy import coordinates
from joblib import dump, load
from io import StringIO
from pathlib import Path
from collections import OrderedDict
from bs4 import BeautifulSoup
import extinction

# self import
from sdapy.gaussian_process import fit_gp
from sdapy.model_fitters import get_pars
from sdapy.filters import *
from sdapy.functions import *
from sdapy.constants import telluric_lines, line_location, line_forsn
from sdapy.specline_handler import handle_spectra, handle_spectrum
from sdapy.pbar import get_progress_bar
from sdapy import __path__, corner_hack, read_default, models, constants

# all
__all__ = ('snelist', 'snobject')

# sdapy source path
srcpath = __path__[0]

# ztfquery source path
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def print_logo(returns=False):
    ''' Print HAFFET Logo.

    Parameters
    ----------
    returns:         `bool`
        if ``returns`` is True, return string, otherwise print in terminal.
    '''
    if returns: _s = ''
    for l in open('%s/data/logo.txt' % srcpath, encoding='utf-8').readlines():
        if not returns: print (l.strip())    
        else: _s += l
    if returns: return _s

def check_dir(check=False, askifdirmissed=False, update=False):
    ''' Check HAFFET data directory, if any folders or files missing, will copy into the data directory.
    
    Parameters
    ----------
    check:         `bool`
        if True, will only check if file existed or not, will return the checking result as a bool
    askifdirmissed :         `bool`
        if any directory missed, ask the permission to make them or just make them directly.
    update:         `bool`
        if True, will replace all data files in the data directory, with those defualts ones, in
        https://github.com/saberyoung/HAFFET/tree/master/sdapy/data
    '''    
    if askifdirmissed: answ = None
    else:  answ = 'Y'
    if not os.path.isdir( LOCALSOURCE ):
        if check: return False
        if answ is None:
            answ = input('as set, %s was the folder used to deal with data however not exists, create directory?(Y/N)'%LOCALSOURCE)
        if answ.strip() in ['y', 'Y']:
            try: os.mkdir( LOCALSOURCE )
            except: raise OSError("Can't create destination directory (%s)!" % (LOCALSOURCE))
        else:
            sys.exit('exist')
    for filestype in ['fritz', 'marshal', 'ForcePhot', 'ForcePhot_atlas', 'oac', 'cache', 'plots', 'images']:
        destdir = '%s/%s'%(LOCALSOURCE, filestype)
        if not os.path.isdir( destdir ):
            if check: return False
            if answ is None:
                answ = input('as set, %s was the folder used to deal with data however some sub folder not exists, create them?(Y/N)'%(LOCALSOURCE))
            if answ.strip() in ['y', 'Y']:
                try:  os.mkdir( destdir )
                except: raise OSError("Can't create destination directory (%s)!" % (destdir))
            else: sys.exit()
    for filename in glob.glob('%s/data/*' % srcpath):
        if not filename.endswith('.txt'): continue
        f = os.path.basename(filename)
        filename1 = '%s/%s' % (LOCALSOURCE, f)
        if not os.path.exists( filename1 ) and check: return False
        if update:
            try:    shutil.copy(filename, filename1)
            except: raise OSError("Can't copy file %s!" % (f))
        elif not os.path.exists( filename1 ):
            try:    shutil.copy(filename, filename1)
            except: raise OSError("Can't copy file %s!" % (f))
    if check: return  True
    
def read_kwargs(kwargs, clsname, check=False, verbose=False):
    ''' Read default parameters.

    Parameters
    ----------
    kwargs :         `dict`
        parameters
    clsname :         `str`
        class name: 

        - snelist

        - snobject
    
    check :         `bool`
        check if given paranames are what we needed
    verbose :         `bool`
        show detailed processing informations
    '''
    assert clsname in ['snelist', 'snobject']
    partype = read_default.get_parameters(return_type=True)[clsname]
    _kwargs = dict()
    for _key in kwargs:
        if not _key in partype:
            if check: raise ValueError(f'Unknown keyword %s'%_key)
            elif verbose: print (f'Unknown keyword %s'%_key)
            continue
        if partype[_key] == 'eval':
            _kwargs[_key] = eval( str(kwargs[_key]) )
        elif partype[_key] == 'str':
            _kwargs[_key] = str(kwargs[_key])
    return _kwargs

class snelist(object):
    """snelist: define a list of *snobject*, handle their data and fittings, aimed for a population study
    
    See Also
    --------
    snobject
    
    Notes
    ----------
    Take careful of meta table of *snelist*, especially their types.
    """
    
    version = 1.0
    ''' Static version info '''
        
    parlist = ['objid', 'aliasid', 'z', 'ra', 'dec', 'mkwebv', 'hostebv', 'dm',
               'tpeak', 'fpeak', 'texp', 'type', 'mag', 'Mag', 'deltam', 'color',
               'mni', 'alpha', 'taum', 'mej', 'ekin', 'trise', 'tfall', 'vexp']
    ''' parameter list '''
    
    def __init__(self, fig=None, ax=None, errors='ignore', updatepar=False, **kwargs):
        """ initialize *snelist* 
        
        Parameters
        ----------
        fig :         `matplotlib.subplot`
            used for histogram/scatter or other population plots
        ax :         `matplotlib.axes`
            used for histogram/scatter or other population plots
        errors :      `str`
              Control raising of exceptions on invalid data for provided dtype.
        
              - ``raise`` : allow exceptions to be raised
              - ``ignore``  : suppress exceptions. On error return original object
              - ``warning`` : suppress exceptions. On error print it as warning        
        updatepar :    `bool`
            if update parameter, when reading them
        kwargs :     `Keyword Arguments`
            see https://github.com/saberyoung/HAFFET/blob/master/sdapy/data/default_par.txt,
            **snelist** part
        
        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()
        >>> a
        <sdapy.snerun.snelist object at 0x7fd6fb805f60>
        """
        if not check_dir(check=True): sys.exit ('check data dir first')            
        
        # ----- read default parameters ----- #
        snelist_pars = read_default.get_parameters(keylist='snelist')['snelist']
        defkwargs = read_kwargs(snelist_pars, 'snelist')
        
        # ----- read meta ----- #
        self.kwargs = read_kwargs(kwargs, 'snelist')
        for _key in defkwargs: self.kwargs.setdefault(_key, defkwargs[_key])        
        
        # ----- define products ----- #
        if not 'data' in self.__dict__: self.data = {}
        if not 'meta' in self.__dict__: self.meta = None
        if not 'params' in self.__dict__: self.params = {}
        
        self.fig = fig        
        self.ax = ax
        self.errors = errors
        self.dset = dict()
        self.plotargs = dict()
        self.pars = []
        self.updatepar = updatepar
        
    def read_kwargs(self, **kwargs):
        """ Define a proper way to read and update optional parameters
        
        Parameters
        ----------                
        kwargs :     `Keyword Arguments`
              optional parameters       
        
        Returns
        ---------- 
        kwargs   :   `Keyword Arguments`
              optional parameters   
        
        See Also
        --------
        snelist.__init__

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()
        >>> a.read_kwargs()
        {'source': 'BTS', 'syntax': "type in ['Ib', 'Ic']", 'verbose': True, 'clobber': False, 'idkey': 'objid', 'idkey1': 'alias', 'use_alias': True, 'sortkey': 'peakm', 'rakey': 'ra', 'deckey': 'dec', 'zkey': 'z', 'distkey': 'dist', 'dmkey': 'dm', 'mkwebvkey': 'ebv', 'hostebvkey': 'hostebv', 'mkwavkey': 'A_V', 'hostavkey': 'host_AV', 'typekey': 'type', 'peaktkey': 'peakt', 'rv': 3.1}
        """        
        if not self.updatepar: _kwargs = self.kwargs.copy()
        for _key in kwargs: 
            if _key in self.kwargs:
                if self.updatepar: self.kwargs[_key] = kwargs[_key]
                else: _kwargs[_key] = kwargs[_key]
            else:
                if self.errors == 'warning':
                    print(f'!!! Error: Unknown keyword %s for snobject, skipped'%_key)
                elif self.errors == 'ignore':
                    pass
                else:
                    raise ValueError(f'Unknown keyword %s'%_key)                
        if self.updatepar: return self.kwargs
        else: return _kwargs
    
    def parse_meta(self, withnew='skip', metafile=None, **kwargs):
        """ Read a meta table from local
        
        Parameters
        ----------                
        withnew :        `bool`
              if meta already exists, and a new meta table comimg:
                use [new] meta instead of the original one
                or [skip] the new one 
                or [merge] them together
        source :      `str`
              which source for metadata: OAC or BTS
        metafile :      `str`
              if you want to use your own meta, set source to None
              and input a existing table metafile
        syntax  :        `str`
              syntax used to make subset of meta table, e.g.
              type in ["SN Ib", "SN Ic"], which will only parse all SNe Ibcs
        verbose :        `bool`
              show detailed running informations   
        idkey :        `str`
              object ID column name, e.g. ZTF name or IAU name
        idkey1    :    `str`
              meta key1, could be IAU or other externel survey name
        sortkey :    `str`
              the column used to sort self.meta table                                
        use_alias :  `bool`
              use aliasid instead of objid? For ZTF objects, it's much easier to use ZTFID instead of IAUID, however by default, the objid is IAUID. Set *use_alias* as True can use ZTFID as objid throughout the following process        
        rakey     :  `str`
              meta key for ra
        deckey    :  `str`
              meta key for dec
        zkey      :  `str`
              meta key for redshift
        distkey   :  `str`
              meta key for distance
        dmkey     :  `str`
              meta key for distance module
        mkwebvkey :  `str`
              meta key for milky way ebv
        hostebvkey:  `str`
              meta key for host ebv
        mkwavkey  :   `str`
              meta key for milky way Av
        hostavkey :    `str`
              meta key for host Av
        typekey   :    `str`
              meta key for SN type
        peaktkey  :   `str`
             meta key for jd at peak        
        
        See Also
        --------
        snelist.parse_meta_one, snelist.parse_meta_all

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()
        >>> a.parse_meta(withnew='skip', source='BTS', metafile=None, syntax='type in ["Ib","Ic"]')
        meta 154 objs
        >>> print (a.meta)
                      alias           ra          dec       peakt    peakfilt  ...   vissave vispeak30    visnow          b    A_V
        objid                                                                  ...
        ZTF21abmlldj  SN2021uvy  00:29:30.87  +12:06:21.0  1449.86    r        ...   5.177543  8.253165  5.257744 -50.406787  0.185
        '''
        [154 rows x 20 columns]
        """ 
        kwargs = self.read_kwargs(**kwargs)
        source = kwargs['source']        
        assert source is not None or metafile is not None, 'input a meta table'
        assert withnew in ['new', 'skip', 'merge']
        if self.meta is not None and withnew == 'skip': return
        
        if metafile is None and source is not None:        
            metafilelist = {
                'BTS' : 'bts_meta.txt',
                'OAC' : 'oac_meta.txt',
            }
            assert source in metafilelist.keys()
            metafile = metafilelist[source]                
            if os.path.isfile( metafile ):
                print (">> WARNING: using local meta file (%s)"%metafile)            
            else:
                metafile = '%s/%s' % (LOCALSOURCE, metafile)        
        if not os.path.exists(metafile):
            print ('!!! Warning: no meta parsed, check if %s exsist' % metafile)
            return 
        # read meta as pandas.dataframe        
        tab, nheaders = '', None
        for n,ll in enumerate(open(metafile, encoding='utf-8').readlines()):        
            if ll[0] == '#' or len(ll.strip()) == 0: continue            
            ll = remove_quote(ll)
            if nheaders is None: nheaders = len(ll.split(','))
            nn = len(ll.split(','))
            if nn != nheaders:
                if kwargs['verbose']: print ('Wrong format: %s, -> Skipped' % ll)
                continue
            tab += ll   
        meta = pd.read_csv(io.StringIO(tab), sep=',').drop_duplicates()
        if len(meta) == 0:
            print ('!!! Warning: no meta data parsed')
            return    
        # format meta keys
        if not kwargs['use_alias']:  id1, id2 = kwargs['idkey'], kwargs['idkey1']
        else:  id1, id2 = kwargs['idkey1'], kwargs['idkey']
        
        if source == 'BTS':                 
            meta.rename(
                columns={
                    'IAUID'    : id1,
                    'ZTFID'    : id2,
                    'peakabs'  : kwargs['sortkey'],
                    'RA'       : kwargs['rakey'],
                    'Dec'      : kwargs['deckey'],
                    'redshift' : kwargs['zkey'],
                    'dist'     : kwargs['distkey'],
                    'dm'       : kwargs['dmkey'],
                    'ebv'      : kwargs['mkwebvkey'],
                    'hostebv'  : kwargs['hostebvkey'],
                    'A_V'      : kwargs['mkwavkey'],
                    'host_AV'  : kwargs['hostavkey'],
                    'type'     : kwargs['typekey'],
                    'peakt'    : kwargs['peaktkey'],
                },
                inplace=True,
            )
            meta['peakt'] += 2458000
        elif source == 'OAC':
            # OAC alias is complicated, so will always use event column as index
            meta.rename(
                columns={
                    'event'    : kwargs['idkey'],
                    'alias'    : kwargs['idkey1'],
                    'maxabsmag': kwargs['sortkey'],
                    'ra'       : kwargs['rakey'],
                    'dec'      : kwargs['deckey'],
                    'redshift' : kwargs['zkey'],
                    'lumdist'  : kwargs['distkey'],
                    'dm'       : kwargs['dmkey'],
                    'ebvebv'   : kwargs['mkwebvkey'],
                    'hostebv'  : kwargs['hostebvkey'],
                    'A_V'      : kwargs['mkwavkey'],
                    'host_AV'  : kwargs['hostavkey'],
                    'claimedtype' : kwargs['typekey'],
                    'maxdate'  : kwargs['peaktkey'],
                },
                inplace=True,
            )        
        
        # format ids
        assert kwargs['idkey'] in meta.keys()
        
        if kwargs['idkey1'] in meta:
            meta = meta.astype({kwargs['idkey']: str, kwargs['idkey1']: str})
        else:
            meta = meta.astype({kwargs['idkey']: str})
        meta = meta.set_index(kwargs['idkey'])

        # format type
        _type = []
        if kwargs['typekey'] in meta.keys():        
            for index in meta.index:            
                sntype = type_clean(meta[kwargs['typekey']][index])            
                _type.append(sntype)
        else:
            _type = '-' * len(meta.index)
        meta[kwargs['typekey']] = _type
        
        # cut meta if syntax available
        if len(kwargs['syntax']) != 0: meta = meta.query(kwargs['syntax'])        

        # sort meta              
        if len(kwargs['sortkey']) != 0 and kwargs['sortkey'] in meta.keys():
            meta = meta.sort_values(kwargs['sortkey'], ascending=False)

        # define meta
        if self.meta is None or withnew == 'new':
            self.meta = meta
        else:           
            self.meta = self.meta.append(meta, ignore_index=False).drop_duplicates()            
        if kwargs['verbose']:
            print ( 'meta %i objs'%( len(meta)) )        
                
    def parse_meta_one(self, idkey, objid, key):
        """ Obtain value with object ID and a meta key
        
        Parameters
        ----------                
        idkey :        `str`
              object ID column name, e.g. ZTF name or IAU name
        objid :       `str`
              object ID string
        key   :       `str`
              a column key of self.meta       

        Returns
        ---------- 
        meta   :   `pandas.dataframe`
        
        See Also
        --------
        snelist.parse_meta, snelist.parse_meta_all

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()
        >>> a.parse_meta(withnew='skip', source='BTS', metafile=None, syntax='type in ["Ib","Ic"]')
        meta 154 objs
        >>> a.parse_meta_one('objid', 'ZTF21abmlldj', 'z')
        '0.0944'
        """        
        _meta = self.meta.query('%s==@objid'%idkey)        
        if len(_meta) == 0: return
        if key in _meta:  return _meta[key][0]
        return

    def parse_meta_all(self, kwargs, objid):
        """ properly read a list of meta infomations from self.meta, i.e. coordinates self.ra, self.dec,
        redshift self.z, distance self.dist, distance module self.dm, mkily way extinction self.mkwebv,
        host galaxy extinction self.hostebv, type self.sntype and peak time self.jdpeak.
        If ra dec missed, user should manully input lightcurve later instead of build-in sources.
        If redshift missed, will make all analysis in obervational frame instead of rest frame.
        If distance missed, will calculate it from redshift with a standard cosmology (astropy.cosmology.Planck13).
        If milky way E(B-V) missed, will check if A_V available, if not, make sure you had dustmaps.sfd installed,
        and SFD dust map is downloaded propoerly with dustmaps, Otherwise will ignore milky way extinction.
        If host galaxy E(B-V) missed, will check if A_V available, otherwise temporarily assign 0 to host E(B-V),
        which can be updated later from colour comparison or Na Id fittings.
        Type is used by colour comparison (which template colour should be compared) and line measurements 
        (which line should be fitted), if missed, will make trouble in these 2 parts.
        jdpeak is can be decided by *snobject* in many ways, but a prior input is important to guess the JD range
        to query photometry.
        (https://github.com/saberyoung/HAFFET/blob/master/sdapy/data/default_par.txt)
        
        Parameters
        ----------
        kwargs :     `Keyword Arguments`
              optional parameters, will use  `idkey`,  `rakey`,  `deckey`,  `zkey`, 
                    `distkey`,  `dmkey`,  `mkwebvkey`,  `mkwavkey`,  `hostebvkey`, 
                    `hostavkey`,  `typekey`,  `peaktkey`,  `idkey`,  and `rv`
        objid :       `str`
              object ID string         

        See Also
        --------
        snelist.parse_meta, snelist.parse_meta_one

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()
        >>> a.parse_meta(withnew='skip', source='BTS', metafile=None, syntax='type in ["Ib","Ic"]')
        meta 154 objs
        >>> a.parse_meta_all(a.read_kwargs(), 'ZTF21abmlldj') 
        >>> a.z
        0.0944
        >>> a.ra
        '00:29:30.87'
        >>> a.dm
        38.25085889208277
        """        
        self.aliasid = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['idkey1']))
        self.ra = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['rakey']))
        self.dec = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['deckey']))
        self.z = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['zkey']))
        self.dist = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['distkey']))
        self.dm = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['dmkey']))
        self.mkwebv = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['mkwebvkey']))
        mkwav = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['mkwavkey']))
        self.hostebv = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['hostebvkey']))
        hostav = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['hostavkey']))
        self.sntype = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['typekey']))
        self.jdpeak = str_clean(self.parse_meta_one(kwargs['idkey'], objid, kwargs['peaktkey']))
        
        if self.ra is None or not is_coordinate(self.ra):  self.ra = None
        if self.dec is None or not is_coordinate(self.dec):  self.dec = None
        if not is_number(self.z): self.z = 0
        self.z = float(self.z)
        if not is_number(self.dist):
            #if kwargs['verbose']: print ('Warning: %s distance set with redshift and standard cosmology' % (objid))
            if self.z > 0:
                self.dist = cosmo.luminosity_distance( self.z ).value
            else:
                self.dist = 0
        self.dist = float(self.dist)
        if not is_number(self.dm):
            if self.dist > 0:
                self.dm = 5*np.log10(self.dist) + 25
            else:
                self.dm = 0
        self.dm = float(self.dm)
        if is_number(self.mkwebv): 
            self.mkwebv = float(self.mkwebv)
        elif is_number(mkwav): 
            self.mkwebv = float(mkwav) / kwargs['rv']            
        elif self.ra is not None and self.dec is not None:
            try:
                from dustmaps.sfd import SFDQuery
                sfd = SFDQuery()                
                coords = coordinates.SkyCoord('%s %s'%(self.ra, self.dec), unit=(u.hourangle, u.deg))                    
                self.mkwebv = sfd(coords)
                #if kwargs['verbose']: print ('Warning: for %s, query mkw ebv via dustmaps.sfd' % (objid))
            except:
                #if kwargs['verbose']: print ('Warning: %s no mkw ebv and av was found, set mkw ebv=0' % (objid))
                self.mkwebv = 0
        else:
            #if kwargs['verbose']: print ('Warning: %s no mkw ebv and av was found, set mkw ebv=0' % (objid))
            self.mkwebv = 0
        if is_number(self.hostebv): 
            self.hostebv = float(self.hostebv)
        elif is_number(hostav):            
            self.hostebv = float(hostav) / kwargs['rv']
        else:
            #if kwargs['verbose']: print ('Warning: %s no host ebv and av was found, set host ebv=0' % (objid))
            self.hostebv = 0
        if is_number(self.jdpeak): 
            self.jdpeak = float(self.jdpeak)
        else:
            self.jdpeak = 0
        
    def parse_params(self, clobber=False, verbose=True, parfile='individual_par.txt'):       
        """ Besides the general parameter settings, for SNe with peculiar properties, a specific parameter is
        sometimes needed, and *parse_params* can read a text file (individual_par.txt) that includes all special
        settings for particular SNe.
        
        Parameters
        ----------                
        clobber :        `bool`
              if meta already exists, reload or skip
        verbose :        `bool`
              Enable progress report
        parfile :      `str`
              parameter filename, file folder is os.getenv('ZTFDATA',"./Data/")        
        
        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()
        >>> a.params
        {}
        >>> a.parse_params(clobber=False, parfile='individual_par.txt')
        >>> a.params
        {'ZTF19abqxibq': {'fit_redo': True, 'fit_methods': ['bazin1', 'arnett_fit_taum', 'gauss'], 'bol_main_xrange': [-30, 20], 'bol_tail_xrange': [40, 120], 'verbose': True}, 'ZTF20aajcdad': {'fit_redo': False, 'plot_bands': ['g', 'r', 'i', 'o', 'c'], 'verbose': True, 'bol_main_xrange': [-20, 40], 'bol_tail_xrange': [40, 120]}}
        """ 
        if len(self.params) > 0 and not clobber: return
        if os.path.isfile( parfile ):
            if verbose: print (">> WARNING: using local parameter file (%s)"%parfile)            
        else:
            parfile = '%s/%s' % (LOCALSOURCE, parfile)
        if not os.path.exists(parfile):
            if verbose: print ('!!! Warning: no prior information found, check if %s exsist' % parfile)
            return
        partype = read_default.get_parameters(return_type=True)['snobject']
        for line in open(parfile, encoding='utf-8').readlines():
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
                    
    def load_data(self, objid, datafile='%s_data.clf', reloadclass=False, **kwargs):
        """ for each object, load their *snobject* classes if they're cached before.
        
        Parameters
        ----------
        objid :       `str`
              object ID string
        clobber :        `bool`
              if meta already exists, reload or skip
        datafile :        `str`
              cached file name
        verbose :        `bool`
              show detailed running informations        
        reloadclass  :     `bool`
              if ``reloadclass`` True, will only include previous data, and reload functions.
              Otherwise, will reload everything that were cached.
        
        Returns
        ---------- 
        flag   :  `bool`
              if cachefile exists return True, otehriwse False
        
        See Also
        --------
        snelist.save_data

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()        
        >>> a.load_data('ZTF20aajcdad')
        True
        >>> a.data
        {'ZTF20aajcdad': <sdapy.snerun.snobject object at 0x7fb9eb5d9f28>}
        """
        kwargs = self.read_kwargs(**kwargs)
        if objid in self.data and not kwargs['clobber']: return            
        datafile = '%s/cache/%s' % (LOCALSOURCE, datafile%objid ) 
        if os.path.exists(datafile):            
            if reloadclass:
                cls = load(datafile)
                self.data[objid] = snobject(objid)                
                for k in ['objid', 'aliasid', 'ra', 'dec', 'z', 'mkwebv', 'hostebv',
                          'sntype', 'dm', 't0', 'texp', 'vexp', 'fpeak', 'fig', 'ax',
                          'ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'lc', 'gpcls', 'fitcls',
                          'colors', 'mbol', 'mbolbb', 'mbolspec', 'spec']:
                    if k not in cls.__dict__: continue
                    self.data[objid].__dict__[k] = cls.__dict__[k]
            else:
                try:    
                    self.data[objid] = load(datafile)
                except EOFError:                    
                    print ('Error: %s (%f) empty file? remove and redo it?' % (datafile, os.path.getsize(datafile)))
                    #os.remove(datafile)
                    return False
            return True
        return False
    
    def save_data(self, objid, datafile='%s_data.clf', **kwargs):
        """ for each object, save their *snobject* classes to local cached files.
        
        Parameters
        ----------
        objid :       `str`
              object ID string
        clobber :        `bool`
              if meta already exists, reload or skip
        datafile :        `str`
              cached file name
        verbose :        `bool`
              show detailed running informations        
        
        Returns
        ---------- 
        flag   :  `bool`
              if saved cache return True, otehrwise False

        See Also
        --------
        snelist.load_data    

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()        
        >>> a.save_data('ZTF20aajcdad')
        False
        >>> a.save_data('ZTF20aajcdad', clobber=True)
        saved cache
        save: /Users/yash0613/Desktop/haffet_data//cache/ZTF20aajcdad_data.clf
        True
        """        
        kwargs = self.read_kwargs(**kwargs)
        datafile = '%s/cache/%s' % (LOCALSOURCE, datafile%objid)
        if kwargs['verbose']: print ('saved cache')
        if not os.path.exists(datafile) or kwargs['clobber']:
            if kwargs['verbose']: print (  'save:', datafile )
            dump(self.data[objid], datafile)
            return True
        return False    
    
    def run(self, ax=None, ax1=None, ax2=None, ax3=None, ax4=None, debug=False, **kwargs):
        """ get a list of SNe, for each SN, define a dedicated *snobject*, and run snobject.run() for
        all of them.
        
        Parameters
        ----------
        ax/ax1/ax2/ax3/ax4 :       `str`
            matplotlib.axes, used for histogram/scatter or other population plots        
        verbose :        `bool`
              show detailed running informations        
        clobber :        `bool`
              if cached file (from *save_data*) exists, redo snobject.run() or just read the cach
        
        See Also
        --------
        snobject.run

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()   
        >>> a.run(syntax='type in ["Ib","Ic"]')
        meta 154 objs
        3%|██▋                                                                              | 5/154 [00:00<00:18,  7.93it/s]
        >>> a.data
        {'ZTF20aajcdad': <sdapy.snerun.snobject object at 0x7fb9eb5d9f28>, 'ZTF21abmlldj': <sdapy.snerun.snobject object at 0x7fb9ec009780>, 'ZTF19abqxibq': <sdapy.snerun.snobject object at 0x7fb9ed674cc0>, ... ...}
        """
        kwargs_snelist = self.read_kwargs(**kwargs)
        self.parse_meta(withnew='skip', metafile=None, **kwargs_snelist)
        self.parse_params(clobber=kwargs_snelist['clobber'], verbose=kwargs_snelist['verbose'])
        
        failed_objs = []
        with get_progress_bar(True, len(self.meta.index)) as pbar:
            for i, objid in enumerate(self.meta.index):               
                # parse meta infos
                self.parse_meta_all(kwargs_snelist, objid)
                
                # load data
                loaded = self.load_data(objid, clobber=kwargs_snelist['clobber'])
                
                if not loaded or kwargs_snelist['clobber']:
                    # read specific kwargs from file
                    par = dict()
                    if objid in self.params: par = self.params[objid]
                    
                    # read specific kwargs from input
                    par = {**par, **kwargs}                    
                    
                    # for each object                         
                    self.data[objid] = snobject(objid, aliasid=self.aliasid, z=self.z,
                        ra=self.ra, dec=self.dec, mkwebv=self.mkwebv, hostebv=self.hostebv,
                        sntype=self.sntype, dm=self.dm, jdpeak=self.jdpeak,
                        ax=ax, ax1=ax1, ax2=ax2, ax3=ax3, ax4=ax4, **par)
                    
                    # run
                    if debug: self.data[objid].run()
                    else:
                        try:  self.data[objid].run()
                        except:
                            failed_objs.append(objid)
                            continue
                    
                    # save data
                    saved = self.save_data(objid, clobber=True)
                    
                pbar.update(1)
        
        print ('\n'*10 + '#'*5 + ' Report:')
        if len(failed_objs) == 0:
            print (' %i objects, all done' % len(self.meta.index))
        else:
            print (' %i objects, %i done, %i failed: %s' %
                (len(self.meta.index), len(self.meta.index)-len(failed_objs), len(failed_objs), failed_objs))
            
    def add_subset(self, syntax=None, astype={"z": float}, **kwargs):
        """ add a data subset, and corresponding plotting kwargs
        
        Parameters
        ----------
        syntax  :        `str`
              syntax used to make subset of meta table, e.g.
              type in ["SN Ib", "SN Ic"], which will only parse all SNe Ibcs        
        astype  :    `dict`
              force the syntax column to be a number, or other types for query
        nbinx  :        `int`, `list`
              histogram bins in x axis
        fontsize  :        `str`
              figure label font size
        labelpad  :        `str`
              figure label pad        
        label   :       `str`
              label for the subset
        color  :        `str`
              color for plotting the subset
        ls     :         `str`
              linestyle for plotting the subset
        marker  :         `str`
              marker for plotting the subset
        markersize  :         `str`
              markersize for plotting the subset
        fillstyle   :         `str`
              fillstyle for plotting the subset

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()   
        >>> a.run(syntax='type in ["Ib"]')
        meta 72 objs
        3%|██▋                                                                              | 5/154 [00:00<00:18,  7.93it/s]
        >>> print (a.dset)
        {}
        >>> a.add_subset('z>.01')
        >>> print (a.dset)
        {'z>.01': {'ZTF21abmlldj': <sdapy.snerun.snobject object at 0x7fd11c9efe80>, 'ZTF19abqxibq': <sdapy.snerun.snobject object at 0x7fd11cbf10b8>, 'ZTF19abcegvm': <sdapy.snerun.snobject object at 0x7fd11dcbdfd0>, 'ZTF19abztknu': <sdapy.snerun.snobject object at 0x7fd11dcc19e8>, 'ZTF20acvebcu': <sdapy.snerun.snobject object at 0x7fd11b4c7be0>, 'ZTF19acgjosf': <sdapy.snerun.snobject object at 0x7fd11cb619b0>, ... ...}}
        """
        _meta = self.meta        
        if syntax is not None:
            if astype is not None:
                for k in astype: assert k in _meta.keys(), 'astype key should be from %s'%_meta.keys()
                _meta = _meta.astype(astype, copy=True, errors='ignore')            
            _meta = _meta.query(syntax)
        else: syntax = 'all'
        self.dset[syntax] = {objid:self.data[objid] if objid in list(_meta.index) else None for objid in self.data}
        self.plotargs[syntax] = kwargs
        
    def add_parameter(self, **kwargs):
        ''' add a parameter
        
        Parameters
        ----------
        parname  :        `str`
              a parameter name from self.parlist
        filt1   :         `str`
              used when e.g. parname is mag, or color
        filt2   :         `str`
              used only when e.g. parname is color, together with filt1
        quant   :         `list`
              used to get fitting samples, e.g. 5th to 95th percentiles as the error
        peakwith  :         `str`
              how to decide the peak
        expwith   :         `str`
              how to decide the explosion epoch
        phase1   :         `str`
              used when parname is detlam, together with phase2
        phase2   :         `str`
              used when parname is detlam, together with phase1
        interpolation  : `str`
              what interpolation to be used when e.g. calculate magnitudes, colours, etc
        corr_mkw   :   `str`
              when calculating absulte mag, if correct milky way extinction if there're any
        corr_host  :   `str`
              when calculating absulte mag, if correct host galaxy extinction if there're any
        index   :   `int`
              if there're multple models can provide, e.g. the nickel mass mni, which one to use
        engine  :   `str`
              which engine models to be used
        model   :   `str`
              which models to be used
        source  :   `str`
              which source to be used

        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snelist()   
        >>> print (a.pars)
        []
        # add r band magnitude at peak
        >>> a.add_parameter(parname='mag', filt1='r', phase1=0)
        # add g-r at 10 day past peak
        >>> a.add_parameter(parname='color', filt1='g', filt2='r', phase1=10)
        # add nickel mass from Arnett model fits
        >>> a.add_parameter(parname='mni', model='Arnett', quant=[.05,.5,.95])
        >>> print (a.pars)
        [{'parname': 'mag', 'filt1': 'r', 'phase1': 0}, 
        {'parname': 'color', 'filt1': 'g', 'filt2': 'r', 'phase1': 10}, 
        {'parname': 'mni', 'model': 'Arnett', 'quant': [0.05, 0.5, 0.95]}]
        '''
        assert 'parname' in kwargs        
        assert kwargs['parname'] in self.parlist
        if not kwargs in self.pars: self.pars.append(kwargs)
        
    def table(self, syntax='all', tablepars=None, tablename=None, style='latex'):
        """ create a latex table for a subset of SNe
        
        Parameters
        ----------
        syntax  :       `str`
              syntax used to make subset of meta table, e.g.
              type in ["SN Ib", "SN Ic"], which will only parse all SNe Ibcs                
        tablename  :        `str`
              output file name 
        tablepars  :        `list`
              which parameter to be included
        style  :        `str`
              table style
        
        Notes
        ---------
        Need to run snelist.add_subset to add a dataset, and snelist.add_parameter
        to add a list of parameters first.
        
        Returns
        ---------- 
        table   :  `str`
              if tablename is None, return a table, otherwise, store the table to tablename
        
        See Also
        ----------
        snelist.add_subset, snelist.add_parameter
        """
        if len(self.dset) == 0:
            print ('Warning: no subset defined !!!')
            return
        if len(self.pars) == 0:
            print ('Warning: no parameter defined !!!')
            return
        if not syntax in self.dset:
            print ('Warning: syntax %s not found in dset !!!'%syntax)
            return
        assert style in ['latex', 'space', 'comma']
        if tablepars is None: ns = len(self.pars)-1
        else: ns = len(tablepars) -1
        outtable = ''
        for n, k in enumerate(self.pars):
            tn = self.get_par(list(self.dset[syntax].keys())[0], returnname=True, **k)
            if tablepars is not None and tn not in tablepars: continue
            if n < ns:
                if style == 'latex': outtable += '%s & ' % tn
                elif style == 'space': outtable += '%s ' % tn
                else: outtable += '%s, ' % tn
            else:
                if style == 'latex': outtable += '%s \\\\ \n' % tn
                else: outtable += '%s \n' % tn
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue            
            for n, k in enumerate(self.pars):
                vs = self.get_par(objid, **k)                
                if vs is not None: 
                    v, vlow, vup = vs
                    _str = self.format_par(v, vlow, vup, digits=3, style=style)
                else:
                    _str = '-'
                if n < ns:
                    if style == 'latex': outtable += '%s & ' % _str
                    elif style == 'space': outtable += '%s ' % _str
                    else: outtable += '%s, ' % _str
                else:
                    if style == 'latex': outtable += '%s \\\\ \n' % _str
                    else: outtable += '%s \n' % _str
        if tablename is None: return outtable
        t = open(tablename, 'w')
        t.write( outtable )
        t.close()
        
    def get_par(self, objid, returnname=False, **kwargs):
        """ get parameter of one SN
        
        Parameters
        ----------
        objid :       `str`
              object ID string        
        returnname :     `bool`
              get number or name of parameter
        parname  :        `str`
              a parameter name from self.parlist
        filt1   :         `str`
              used when e.g. parname is mag, or color
        filt2   :         `str`
              used only when e.g. parname is color, together with filt1
        quant   :         `list`
              used to get fitting samples, e.g. 5th to 95th percentiles as the error
        peakwith  :         `str`
              how to decide the peak
        expwith   :         `str`
              how to decide the explosion epoch
        phase1   :         `str`
              used when parname is detlam, together with phase2
        phase2   :         `str`
              used when parname is detlam, together with phase1
        interpolation  : `str`
              what interpolation to be used when e.g. calculate magnitudes, colours, etc
        corr_mkw   :   `str`
              when calculating absulte mag, if correct milky way extinction if there're any
        corr_host  :   `str`
              when calculating absulte mag, if correct host galaxy extinction if there're any
        index   :   `int`
              if there're multple models can provide, e.g. the nickel mass mni, which one to use
        engine  :   `str`
              which engine models to be used
        model   :   `str`
              which models to be used
        source  :   `str`
              which source to be used
        """
        parname = kwargs['parname']
        filt1 = kwargs.get('filt1', 'g')
        filt2 = kwargs.get('filt2', 'r')
        quant = kwargs.get('quant', [.05,.5,.95])
        peakwith = kwargs.get('peakwith', 'gp')
        expwith = kwargs.get('expwith', 'pl')
        phase1 = kwargs.get('phase1', 0)
        phase2 = kwargs.get('phase2', 15)
        interpolation = kwargs.get('interpolation', None) 
        corr_mkw = kwargs.get('corr_mkw', False)
        corr_host = kwargs.get('corr_host', False)
        index = kwargs.get('index', 0)
        engine = kwargs.get('engine', None)
        model = kwargs.get('model', None)
        source = kwargs.get('source', None)
        
        if parname == 'z':
            # bestfit, lower limit, upper limit
            if returnname: return 'z'
            return self.data[objid].z, None, None
        if parname == 'objid':
            if returnname: return 'objid'
            return self.data[objid].objid, None, None
        if parname == 'aliasid':
            if returnname: return 'aliasid'
            return self.data[objid].aliasid, None, None
        if parname == 'ra':
            if returnname: return 'Ra'
            return self.data[objid].parse_coo(deg=True)[0], None, None
        if parname == 'dec':
            if returnname: return 'Dec'
            return self.data[objid].parse_coo(deg=True)[1], None, None
        if parname == 'mkwebv':
            if returnname: return 'E(B-V)$_{MKW}$'
            return self.data[objid].mkwebv, None, None
        if parname == 'hostebv':
            if returnname: return 'E(B-V)$_{host}$'
            return self.data[objid].hostebv, None, None
        if parname == 'dm':
            if returnname: return 'DM_%s' % filt1            
            return self.data[objid].dm, self.data[objid].dm_error(filt1), None
        if parname in ['tpeak', 'fpeak']:
            if peakwith == 'gp':
                if returnname: return parname + '_GP_%s'%filt1
                if not 'gpcls' in self.data[objid].__dict__: return                
                t0, fpeak = self.data[objid].set_peak_gp(
                    set_tpeak_filter=filt1, returnv=True
                )
            elif peakwith == 'bol':
                if returnname: return parname + '_bol'
                if not 'fitcls' in self.data[objid].__dict__: return
                if not 'bol_main' in self.data[objid].fitcls: return                
                t0, fpeak = self.data[objid].set_peak_bol_main(
                    model_name=model, source_name=source, returnv=True
                )
            elif peakwith == 'multiband':
                if returnname: return parname + '_multimain'
                if not 'fitcls' in self.data[objid].__dict__: return
                if not 'multiband_main' in self.data[objid].fitcls: return
                
                t0, fpeak = self.data[objid].set_peak_multiband_main(
                    set_tpeak_filter=filt1, model_name=model, returnv=True
                )
            else:
                print ('Error: wrong peakwith value')
                return
            if t0 is None: return            
            if parname == 'tpeak':  return t0, None, None
            else: return fpeak[filt1][0], fpeak[filt1][1], None
        if parname == 'texp':            
            if expwith == 'pl':
                if returnname: return 'texp_pl_%s'%filt1                
                if not 'fitcls' in self.data[objid].__dict__: return
                if not 'multiband_early' in self.data[objid].fitcls: return
                
                return self.data[objid].set_texp_pl(
                    set_texp_filter = filt1, model_name=model, returnv=True
                )
            elif expwith == 'bol':
                if returnname: return 'texp_bol'
                if not 'fitcls' in self.data[objid].__dict__: return
                if not 'bol_main' in self.data[objid].fitcls: return                
                return self.data[objid].set_texp_bol_main(
                    model_name=model, source_name=source, returnv=True
                )
            elif expwith == 'mid':
                if returnname: return 'texp_mid'
                return self.data[objid].set_texp_midway(returnv=True)
            else:
                print ('Error: wrong expwith value')        
        if parname == 'type':
            if returnname: return 'Type'
            return self.data[objid].sntype, None, None
        if parname == 'mag':
            if returnname: return 'mag_%s_%s'%(filt1, phase1) 
            if not 'lc' in self.data[objid].__dict__: return            
            m, dm = self.data[objid]._mag_at(
                filt1, phase1, interpolation=interpolation, index=index,
                corr_mkw=corr_mkw, corr_host=corr_host, quantile=quant,
                snrt=3, tdbin=1,
            )
            return m, dm, None
        if parname == 'Mag':
            if returnname: return 'Mag_%s_%s'%(filt1, phase1)
            if not 'lc' in self.data[objid].__dict__: return            
            m, dm = self.data[objid]._absmag_at(
                filt1, phase1, interpolation=interpolation, index=index,
                corr_mkw=corr_mkw, corr_host=corr_host, quantile=quant,
                snrt=3, tdbin=1,
            ) 
            return m, dm, None
        if parname == 'deltam':
            if returnname: return '$\deltam$_%s_%s-%s$'%(filt1,phase1,phase2)
            if not 'lc' in self.data[objid].__dict__: return            
            m, dm = self.data[objid]._rate_at(
                filt1, phase1, phase2, interpolation=interpolation, index=index,
                corr_mkw=corr_mkw, corr_host=corr_host, quantile=quant,
                snrt=3, tdbin=1,
            )
            return m, dm, None
        if parname == 'color':
            if returnname: return '%s-%s_%s'%(filt1,filt2,phase1)
            if not 'lc' in self.data[objid].__dict__: return            
            m, dm = self.data[objid]._color_at(
                filt1, filt2, phase1, interpolation=interpolation, index=index,
                corr_mkw=corr_mkw, corr_host=corr_host, quantile=quant,
                snrt=3, tdbin=1,
            )
            return m, dm, None                              
        if parname in ['mni', 'alpha', 'taum', 'mej', 'ekin', 'trise', 'tfall']:
            if returnname: return parname
            _dict = self.data[objid].get_par(                            
                parname, index=index, engine=engine, model=model,
                source=source, filt=filt1, k_opt=None, quantile=quant,
            )
            if _dict is None: return
            elif len(_dict) == 0: return
            v = list(_dict.values())[0]
            return v[1], v[0], v[2]
        if parname == 'vexp':
            if returnname: return parname
            if not 'fitcls' in self.data[objid].__dict__: return
            if not 'specv_evolution' in self.data[objid].fitcls: return
            v = self.data[objid].set_vexp(
                model_name=model, returnv=True, quantile=quant,
            )
            if v is None: return
            return v[1], v[0], v[2]
        return
    
    @staticmethod
    def format_par(v, vlow, vup, digits=3, style='latex'):
        """ make parameter into latex format
        
        Parameters
        ----------
        v    :       `float`
              best fit value of one parameter   
        vup  :       `float`
              upper limit of one parameter   
        vlow  :       `float`
              lower limit of one parameter   
        digits  :        `int`
              number digits
        style  :        `str`
              table style        
        """ 
        if vlow is None and vup is None: # just value
            if is_number(v):
                return '%s' % round(v, digits)
            else:
                return v
        elif vup is None: # errorbars
            if is_number(v):
                return '%s (%s)' % ( round(v, digits), round(abs(vlow), digits) )
            else:
                return '%s %s' % ( v, vlow)
        else: # upper limit and lower limits
            if is_number(v):
                if style == 'latex':
                    return '%s$_{-%s}^{%s}$' % ( round(v, digits), round(abs(v-vlow), digits), round(abs(v-vup), digits) )
                else:
                    return '%s (-%s, %s)' % ( round(v, digits), round(abs(v-vlow), digits), round(abs(v-vup), digits) )                
            else:
                return '%s %s %s' % ( v, vlow, vup)
            
    def show1d(self, index=0, style='pdf'):
        """ 1D histograms plot for one parameter
        
        Parameters
        ----------
        index    :       `int`
              which par to be used       
        style    :       `str`
              plot style: [pdf] or [cdf]
        """
        assert self.ax is not None
        if len(self.dset) == 0:
            print ('Warning: no subset defined !!!')
            return
        if len(self.pars) == 0:
            print ('Warning: no parameter defined !!!')
            return
        for syntax in self.dset:            
            x = []            
            for objid in self.dset[syntax]:
                if self.dset[syntax][objid] is None: continue                        
                _x = self.get_par(objid, **self.pars[index])
                if _x is None or not is_number(_x[0]):
                    #print ('can not get %s parameter for %s'%(self.pars[index], objid))
                    continue
                if is_number(_x[0]): x.append( float(_x[0]) )
                else: x.append(_x[0])
            if len(x) == 0:
                print ('no data parsed')
                return
            if is_number(x[0]):
                nbin = self.plotargs[syntax].get('nbin', 10)
                binwidthx = (np.max(np.abs(x))-np.min(np.abs(x)))/nbin                
                binx = np.arange(np.min(np.abs(x)), np.max(np.abs(x))+binwidthx, binwidthx)
                if style == 'pdf':
                    cumulative=False
                elif style == 'cdf':
                    cumulative=True
                else:
                    return                
                self.ax.hist(x, bins=binx, histtype='step', cumulative=cumulative,
                             label=self.plotargs[syntax].get('label', None),
                             color=self.plotargs[syntax].get('color', None),
                             ls='-' #self.plotargs[syntax].get('ls', '-'),
                )
            else:            
                labels, counts = np.unique(x,return_counts=True)                
                self.ax.bar(labels, counts, 
                            label=self.plotargs[syntax].get('label', None),
                            color=self.plotargs[syntax].get('color', None),
                            ls='-' #self.plotargs[syntax].get('ls', '-'),
                )
        self.ax.set_xlabel(self.pars[index]['parname'], fontsize=self.plotargs[syntax].get('fontsize', 12),
                           labelpad=self.plotargs[syntax].get('labelpad', 12))
        self.ax.set_ylabel('counts', fontsize=self.plotargs[syntax].get('fontsize', 12),
                           labelpad=self.plotargs[syntax].get('labelpad', 12))
        self.ax.legend()
        
    def shownd(self, syntax='all', **kwargs):
        """ contour plots for the bestfit value of all parameters
        
        Parameters
        ----------        
        syntax    :       `int`
              which subset for scatter
        quant   :         `list`
              used to get fitting samples, e.g. 5th to 95th percentiles as the error        
        """
        if not syntax in self.dset:
            print ('Warning: syntax %s not found in dset !!!'%syntax)
            return
        quant = kwargs.get('quant', [.05,.5,.95])        
        if len(self.dset) == 0:
            print ('Warning: no subset detected !!!')
        if len(self.pars) == 0:
            print ('Warning: no parameter defined !!!')
            return        
        xl, parname = [], []
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue
            if len(parname) == 0:                    
                for _kwargs in self.pars:
                    if _kwargs['parname'] == 'type': continue
                    parname.append(self.get_par(objid, returnname=True, **_kwargs))
            _xl, flag = [], True
            for _kwargs in self.pars:
                if _kwargs['parname'] == 'type': continue
                _x = self.get_par(objid, **_kwargs)
                if _x is None: continue
                elif not is_number(_x[0]): continue
                _xl.append( float(_x[0]) )
            if len(_xl) == len(parname): xl.append( _xl )
        if len(xl) == 0 or len(parname) == 0:
            print ('no data parsed')
            return
        fig = corner_hack.corner_hack(xl, labels=parname,
                        label_kwargs={'fontsize':16}, ticklabelsize=13,
                        show_titles=True, quantiles=quant,
                        title_fmt=".2f", title_kwargs={'fontsize':16},
                        plot_datapoints=True, plot_contours=True)     
        return fig
    
    def show2d(self, index1, index2):
        """ 2D scatter plot for two parameter
        
        Parameters
        ----------
        index1    :       `int`
              index for column1 name for x axis
        index2    :       `int`
              index for column2 name for y axis          
        """
        assert self.ax is not None
        if len(self.dset) == 0:
            print ('Warning: no subset detected !!!')
        if len(self.pars) == 0:
            print ('Warning: no parameter defined !!!')
            return
        xlim,ylim=None,None
        # histograms        
        self.init_hist_axes(pad=0.1, labelbottom=False, labelleft=False)
        
        for syntax in self.dset:            
            x, y, xe1, xe2, ye1, ye2 = [], [], [], [], [], []
            for objid in self.dset[syntax]:
                if self.dset[syntax][objid] is None: continue
                _x = self.get_par(objid, **self.pars[index1])
                _y = self.get_par(objid, **self.pars[index2])
                if _x is None or _y is None:  continue 
                if not is_number(_x[0]) or not is_number(_y[0]): continue                
                # append
                x.append( float(_x[0]) )
                y.append( float(_y[0]) )
                if _x[1] is not None: xe1.append( float(_x[1]) )
                else: xe1.append( 0 )
                if _x[2] is not None: xe2.append( float(_x[2]) )
                else: xe2.append( 0 )
                if _y[1] is not None: ye1.append( float(_y[1]) )
                else: ye1.append( 0 )
                if _y[2] is not None: ye2.append( float(_y[2]) )
                else: ye2.append( 0 )
            if len(x) == 0:
                print ('no data parsed')
                return
            xerr = (xe1,xe2)
            yerr = (ye1,ye2)            
            self.ax.errorbar( x, y, xerr=xerr, yerr=yerr,
                        marker=self.plotargs[syntax].get('marker', 'o'),
                        markersize=self.plotargs[syntax].get('markersize', 4),
                        fillstyle=self.plotargs[syntax].get('fillstyle', 'none'),
                        label=self.plotargs[syntax].get('label', None),
                        color=self.plotargs[syntax].get('color', None),
                        ls=self.plotargs[syntax].get('ls', '-') )
            
            # histograms
            self.add_hist(x, y, 
                    nbinx = self.plotargs[syntax].get('nbin', 10),
                    nbiny = self.plotargs[syntax].get('nbin', 10),
                    color = self.plotargs[syntax].get('color', None),
                    xticks=None, yticks=None)
            
            if xlim is None:  xlim = [min(x)/1.1, max(x)*1.1]
            else:  xlim = [min(min(xlim), min(x)/1.1), max(max(xlim), max(x)/1.1)]
            
            if ylim is None:  ylim = [min(y)/1.1, max(y)*1.1]
            else:  ylim = [min(min(ylim), min(y)/1.1), max(max(ylim), max(y)/1.1)]
            
        # limit
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        
        self.ax.set_xlabel(self.pars[index1]['parname'],
                           fontsize=self.plotargs[syntax].get('fontsize', 12),
                           labelpad=self.plotargs[syntax].get('labelpad', 12))
        self.ax.set_ylabel(self.pars[index2]['parname'],
                           fontsize=self.plotargs[syntax].get('fontsize', 12),
                           labelpad=self.plotargs[syntax].get('labelpad', 12))
        self.ax.legend()
        
    def init_hist_axes(self, pad=0.1, labelbottom=False, labelleft=False):
        """ create 2 subplots as histograms for the scatter plots 
        
        Parameters
        ----------
        pad    :       `float`
              distance between scatter plot and the histograms
        labelbottom    :       `str`
              bottom histogram label name
        labelletf  :        `str`
              left histogram label name       
        """
        assert self.ax is not None
        # create new axes on the right and on the top of the current axes
        divider = make_axes_locatable(self.ax)
        # below height and pad are in inches
        self.ax_histx = divider.append_axes("top", 1.2, pad=pad, sharex=self.ax)
        self.ax_histy = divider.append_axes("right", 1.2, pad=pad, sharey=self.ax)
        # make some labels invisible
        self.ax_histx.xaxis.set_tick_params(labelbottom=labelbottom)
        self.ax_histy.yaxis.set_tick_params(labelleft=labelleft)
        
    def add_hist(self, x, y, nbinx = 10, nbiny = 10, 
                 xticks=None, yticks=None, **kwargs):
        """ make histograms
        
        Parameters
        ----------
        x    :       `float`
              table column name for x axis
        y    :       `str`
              table column name for y axis        
        nbinx  :        `int`, `list`
              histogram bins for x axis
        nbiny  :        `int`, `list`
              histogram bins for y axis
        xticks  :        `int`, `list`
              ticks for x axis
        yticks  :        `int`, `list`
              ticks for y axis
        kwargs :     `Keyword Arguments`
              **matplotlib.ax.hist** kwargs
        """
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

    def showax(self, syntax=None, show_data=True, show_fits=True, show_gp=True,
            show_fit_error=True, interpolation=None, index=0, plot_bands=None): 
        """ make flux plot for a large set of SNe
        
        Parameters
        ----------   
        syntax  :        `str`
              syntax used to make subset of meta table, e.g.
              type in ["SN Ib", "SN Ic"], which will only parse all SNe Ibcs   

        See Also
        ----------       
        snobject._ax
        """
        assert self.ax is not None
        if not syntax in self.dset:
            print ('Warning: syntax %s not found in dset !!!'%syntax)
            return        
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue
            self.dset[syntax][objid].ax = self.ax
            self.dset[syntax][objid]._ax(
                show_title=False, show_legend=False, ylabel_2right=False,
                show_data=show_data, show_fits=show_fits, show_gp=show_gp,
                show_fit_error=show_fit_error, show_texp=False,
                interpolation=interpolation, index=index, snr_thre=3,
                plot_bands=plot_bands, ax_xstyle='rp', ax_ystyle='norm',
                color=self.plotargs[syntax].get('color', 'k'),
                marker=self.plotargs[syntax].get('marker', 'o'),
                markersize=self.plotargs[syntax].get('markersize', 12),
                label=self.plotargs[syntax].get('label', None),
                ls=self.plotargs[syntax].get('ls', 'none'),
                fillstyle=self.plotargs[syntax].get('fillstyle', None),
                fontsize=self.plotargs[syntax].get('fontsize', 12),
                verbose=False,
            )
            
    def showax2(self, syntax=None, show_data=True, show_fits=True,
                show_gp=True, show_fit_error=True, plot_bands=None):
        """ make mag plot for a large set of SNe
        
        Parameters
        ----------   
        syntax  :        `str`
              syntax used to make subset of meta table, e.g.
              type in ["SN Ib", "SN Ic"], which will only parse all SNe Ibcs   
        
        See Also
        ----------       
        snobject._ax2
        """        
        assert self.ax is not None
        if not syntax in self.dset:
            print ('Warning: syntax %s not found in dset !!!'%syntax)
            return
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue
            self.dset[syntax][objid].ax2 = self.ax
            self.dset[syntax][objid]._ax2(
                show_title=False, show_legend=False, ylabel_2right=False,
                show_data=show_data, show_fits=show_fits, show_gp=show_gp,
                show_fit_error=show_fit_error, show_texp=False,
                interpolation=None, index=0, show_limits=False,
                plot_bands=plot_bands, ax2_xstyle='rp', ax2_ystyle='abs',
                color=self.plotargs[syntax].get('color', 'k'),
                marker=self.plotargs[syntax].get('marker', 'o'),
                markersize=self.plotargs[syntax].get('markersize', 12),
                label=self.plotargs[syntax].get('label', None),
                ls=self.plotargs[syntax].get('ls', 'none'),
                fillstyle=self.plotargs[syntax].get('fillstyle', None),
                fontsize=self.plotargs[syntax].get('fontsize', 12),
                corr_mkw=True, corr_host=True, verbose=False,
            )
        self.ax.invert_yaxis()
        
    def showax3(self, syntax=None, show_data=True, show_fits=False,
        show_gp=False, show_fit_error=True):
        """ make colour plot for a large set of SNe
        
        Parameters
        ----------   
        syntax  :        `str`
              syntax used to make subset of meta table, e.g.
              type in ["SN Ib", "SN Ic"], which will only parse all SNe Ibcs   
        
        See Also
        ----------       
        snobject._ax3
        """      
        assert self.ax is not None
        if not syntax in self.dset:
            print ('Warning: syntax %s not found in dset !!!'%syntax)
            return
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue
            self.dset[syntax][objid].ax3 = self.ax
            self.dset[syntax][objid]._ax3(                
                show_title=False, show_legend=False, ylabel_2right=False,
                show_data=show_data, show_fits=show_fits, show_gp=show_gp,
                show_fit_error=show_fit_error, show_texp=False, ax3_xstyle='rp',
                color=self.plotargs[syntax].get('color', 'k'),
                marker=self.plotargs[syntax].get('marker', 'o'),
                markersize=self.plotargs[syntax].get('markersize', 12),
                label=self.plotargs[syntax].get('label', None),
                ls=self.plotargs[syntax].get('ls', 'none'),
                fillstyle=self.plotargs[syntax].get('fillstyle', None),
                fontsize=self.plotargs[syntax].get('fontsize', 12),
                verbose=False, color_interp=['bin'],
                corr_mkw=True, corr_host=True,
            )
        
    def showax4(self, syntax=None, show_data=True, show_fits=True, show_fit_error=True):
        """ make bolometric lc plot for a large set of SNe
        
        Parameters
        ----------   
        syntax  :        `str`
              syntax used to make subset of meta table, e.g.
              type in ["SN Ib", "SN Ic"], which will only parse all SNe Ibcs   
        
        See Also
        ----------       
        snobject._ax4
        """      
        assert self.ax is not None
        if not syntax in self.dset:
            print ('Warning: syntax %s not found in dset !!!'%syntax)
            return
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue
            self.dset[syntax][objid].ax4 = self.ax
            self.dset[syntax][objid]._ax4(
                show_title=False, make_bol=['lyman'], show_legend=False,
                ylabel_2right=False, show_data=show_data, show_fits=show_fits,
                show_fit_error=show_fit_error, show_texp=False, ax4_xstyle='rp',
                color=self.plotargs[syntax].get('color', 'k'),
                marker=self.plotargs[syntax].get('marker', 'o'),
                markersize=self.plotargs[syntax].get('markersize', 12),
                label=self.plotargs[syntax].get('label', None),
                ls=self.plotargs[syntax].get('ls', 'none'),
                fillstyle=self.plotargs[syntax].get('fillstyle', None),
                fontsize=self.plotargs[syntax].get('fontsize', 12),
                verbose=False,
            )

    def showax6(self, syntax=None, show_data=True, show_fits=True, show_fit_error=True):
        """ velocity evolution plot for a large set of SNe
        
        Parameters
        ----------   
        syntax  :        `str`
              syntax used to make subset of meta table, e.g.
              type in ["SN Ib", "SN Ic"], which will only parse all SNe Ibcs   
        
        See Also
        ----------       
        snobject._ax4
        """      
        assert self.ax is not None
        if not syntax in self.dset:
            print ('Warning: syntax %s not found in dset !!!'%syntax)
            return
        for objid in self.dset[syntax]:
            if self.dset[syntax][objid] is None: continue
            self.dset[syntax][objid].ax6 = self.ax
            self.dset[syntax][objid]._ax6(
                show_title=False, show_legend=False,
                show_data=show_data, show_fits=show_fits,
                show_fit_error=show_fit_error,               
                color=self.plotargs[syntax].get('color', 'k'),
                marker=self.plotargs[syntax].get('marker', 'o'),
                markersize=self.plotargs[syntax].get('markersize', 12),
                label=self.plotargs[syntax].get('label', None),
                ls=self.plotargs[syntax].get('ls', 'none'),
                fillstyle=self.plotargs[syntax].get('fillstyle', None),
                fontsize=self.plotargs[syntax].get('fontsize', 12),
                verbose=False,
            )
            
class snobject(object):
    """snobject: define *snobject* for one SN, handle its data and fittings.
    
    See Also
    --------
    snelist    
    """
    
    version = 1.0
    ''' Static version info '''
    
    urllist = {
        'AL'  : 'https://api.astrocats.space/',
        'SN'  : 'https://api.sne.space/',
        'TDE' : 'https://api.tde.space/',
        'KN'  : 'https://api.kilonova.space/',
        'FT'  : 'https://api.faststars.space/'
    }
    ''' different url origins for the Open Astronomical Catalog (*OAC*) catalog '''
    
    keys2query_lc = {
        'time'              : True,
        'e_time'            : False,
        'e_lower_time'      : False,
        'e_upper_time'      : False,
        'u_time'            : False,
        'instrument'        : False,
        'host'              : False,
        'includeshost'      : False,
        'magnitude'         : True,
        'e_magnitude'       : True,
        'e_lower_magnitude' : False,
        'e_upper_magnitude' : False,
        'band'              : True,
        'system'            : False,
        'zeropoint'         : False,
        'upperlimit'        : False,
        'upperlimitsigma'   : False
    }
    '''
    keys to query LCs from *OAC* catalog,    
    if value True, key mush be included by transient, otherwise, could be empty.
    Description: https://github.com/astrocatalogs/schema
    '''
    
    keys2query_spec = [
        'time', 'instrument', 'data' 
    ]
    '''
    keys to query spectra from *OAC* catalog
    '''    
    def __init__(self, objid, errors='ignore', aliasid=None, z=None, ra=None, dec=None,
                 mkwebv=None, hostebv=None, sntype=None, dm=None, jdpeak=None, fig=None,
                 ax=None, ax1=None, ax2=None, ax3=None, ax4=None, ax5=None,
                 updatepar=False, **kwargs):
        """ initialize *snelist* 
        
        Parameters
        ----------
        objid :       `str`
              object ID string   
        errors :      `str`
              Control raising of exceptions on invalid data for provided dtype.

              - ``raise`` : allow exceptions to be raised
              - ``ignore``  : suppress exceptions. On error return original object
              - ``warning`` : suppress exceptions. On error print it as warning
        aliasid :       `str`
              other name of object
        z :       `float`
              redshift
        ra :       `str`
              R.A. in ``hh:mm:ss``
        dec :      `str`
              Dec. in ``dd:mm:dd``
        mkwebv :    `float`
              milky way extinction, E(B-V)
        hostebv :   `float`
              host galaxy extinction, E(B-V)
        sntype :    `str`
              transient/supernova type
        dm :       `float`
              distance module
        jdpeak :    `float`
              julian date of peak epoch
        fig :         `matplotlib.subplot`
              for figures
        ax  :         `matplotlib.axes`
              for flux lightcurve plot
        ax1 :         `matplotlib.axes`
              for spectra plot
        ax2 :         `matplotlib.axes`
              for magnitude lightcurve plot
        ax3 :         `matplotlib.axes`
              for colour plot
        ax4 :         `matplotlib.axes`
              for luminosoty plot       
        ax4 :         `matplotlib.axes`
              for SED plot
        updatepar :    `bool`
              if update parameter, when reading them
        kwargs :     `Keyword Arguments`
              see https://github.com/saberyoung/HAFFET/blob/master/sdapy/data/default_par.txt,
              **snobject** part
        
        Notes
        -----        
        Using ``sdapy`` instead of ``haffet``
        
        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snobject('sss')
        >>> a
        <sdapy.snerun.snobject object at 0x7fa4c8284780>
        """
        if not check_dir(check=True): sys.exit ('check data dir first')            
        
        # auth
        global keypairs, keypairs_type
        keypairs = read_default.get_keypairs()
        keypairs_type = read_default.get_keypairs(return_type=True)
        
        # ----- read default parameters ----- #
        snobject_pars = read_default.get_parameters(keylist='snobject')['snobject']        
        defkwargs = read_kwargs(snobject_pars, 'snobject')
        
        # ----- read meta ----- #
        self.kwargs = read_kwargs(kwargs, 'snobject')
        for _key in defkwargs: self.kwargs.setdefault(_key, defkwargs[_key])        
        
        # set meta infos
        objid = name_clean(objid)
        if objid is None and aliasid is not None: objid = name_clean(aliasid)
        if objid is None:
            print ('Error: can not set the name properly, check %s (%s)' % (objid, aliasid))
            return            
        self.objid  = objid
        self.errors = errors
        self.aliasid= aliasid
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
        if jdpeak is not None:
            self.t0 = jdpeak  # zero point thoughout the analysis                                  
        else:
            self.t0 = 0       # if None, try to get it via GP since it was needed for fits/pl/arnett...
        self.texp   = None    # explosion epoch
        self.vexp   = None    # photospheric velocity, unit in 1000 km/s
        self.fpeak  = dict()  # peak flux in uJy for different bands       
        
        # plots
        self.fig = fig
        self.ax  = ax  # flux
        self.ax1 = ax1 # spectra
        self.ax2 = ax2 # mag        
        self.ax3 = ax3 # color
        self.ax4 = ax4 # luminosity
        self.ax5 = ax5 # sed
        self.updatepar = updatepar
        
    def read_kwargs(self, **kwargs):
        """ Define a proper way to read and update optional parameters
        
        Parameters
        ----------                
        kwargs :     `Keyword Arguments`
              optional parameters       
        
        Notes
        -----        
        The defined ``errors`` parameter decide how to react when getting undefined parameters
        
        See Also
        --------
        snobject.__init__, snelist.read_kwargs
        
        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snobject(sss)
        >>> a.read_kwargs()
        >>> print (a.kwargs)
        {'jdkey': 'jdobs', 'magkey': 'mag', 'emagkey': 'emag', 'limmagkey': 'limmag', ... ...}
        """
        if not self.updatepar: _kwargs = self.kwargs.copy()
        for _key in kwargs: 
            if _key in self.kwargs:
                if self.updatepar: self.kwargs[_key] = kwargs[_key]
                else: _kwargs[_key] = kwargs[_key]
            else:
                if self.errors == 'warning':
                    print(f'!!! Error: Unknown keyword %s for snobject, skipped'%_key)
                elif self.errors == 'ignore':
                    pass
                else:
                    raise ValueError(f'Unknown keyword %s'%_key)                
        if self.updatepar: return self.kwargs
        else: return _kwargs
        
    def run(self, **kwargs):
        """ Running a lot of actions for one SN, according to default settings.
        
        Parameters
        ----------
        kwargs :     `Keyword Arguments`
              optional parameters      
              see https://github.com/saberyoung/HAFFET/blob/master/sdapy/data/default_par.txt,
        
        Examples
        --------
        >>> from sdapy import snerun
        >>> a = snerun.snobject(objid='ZTF20aajcdad', aliasid='SN2020bcq',
        ...      z=0.0186, dm=34.6, mkwebv=0.01387, hostebv=0, sntype='Ib', 
        ...     ra='13:26:29.65', dec='+36:00:31.1', jdpeak=2458888.02)
        >>> a.run()
        
        See Also
        --------
        snelist.run
        """         
        kwargs = self.read_kwargs(**kwargs)         
        t_start = time.time()
        
        ''' parse local data via objid '''
        # photometry
        if len(kwargs['lctype']) == 0 or 'ztffp' in kwargs['lctype']:
            self.get_fp_ztf(**kwargs)
        if len(kwargs['lctype']) == 0 or 'atlasfp' in kwargs['lctype']:
            self.get_fp_atlas(**kwargs)
        if len(kwargs['lctype']) == 0 or 'marshal' in kwargs['lctype']:
            self.get_alert_ztf(source='marshal', **kwargs)
        if len(kwargs['lctype']) == 0 or 'fritz' in kwargs['lctype']:
            self.get_alert_ztf(source='fritz', **kwargs)
        if len(kwargs['lctype']) == 0 or 'oac' in kwargs['lctype']:
            self.get_oac(which='photometry', **kwargs)
        if 'lc' not in self.__dict__:
            print ('no lc found in local, try query them first...')            
            return
        
        # spectra
        if len(kwargs['spectype']) == 0 or 'marshal' in kwargs['spectype']:            
            self.get_local_spectra(source='marshal', **kwargs)
        if len(kwargs['spectype']) == 0 or 'fritz' in kwargs['spectype']:            
            self.get_local_spectra(source='fritz', **kwargs)
        if len(kwargs['spectype']) == 0 or 'oac' in kwargs['spectype']:            
            self.get_oac(which='spectra', **kwargs)
        if kwargs['verbose']:
            print("Get all Data in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        ''' intepolation'''       
        self.run_gp(**kwargs)     # Gaussian process
        if kwargs['verbose']:
            print("Run GP in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
                
        if len(kwargs['set_texp_method'])==0: # set texp with detections
            self.set_texp_midway()
        
        self.clip_lc(**kwargs)  # remove lc outliers                  
        if kwargs['verbose']:
            print("Clip data in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.bin_lc(**kwargs)  # remove lc outliers                  
        if kwargs['verbose']:
            print("Bin phorometric data in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.run_fit('multiband_main',**kwargs)   # sn-like analytic functions
        if kwargs['verbose']:
            print("Run multiband Main peak fitting in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()

        self.run_fit('multiband_early',**kwargs)  # power law for first light
        if kwargs['verbose']:
            print("Run multiband early fitting in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()        

        ''' build sed for bolometric '''                
        self.calc_colors(**kwargs)  # calculate colours
        if kwargs['verbose']:
            print("Calculate colours in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
        
        self.est_hostebv_with_colours(**kwargs)   # compare colours to template for host ebv
        if kwargs['verbose']:
            print("estimate host ebv with colours in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.lyman_bol(**kwargs)  # calculate luminosity from Lyman BC correction
        if kwargs['verbose']:
            print("Calculate bolometric LCs with lyman BC in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()        
        
        self.bb_colors(**kwargs)  # match epochs for BB
        if kwargs['verbose']:
            print("Match epochs for SED with BB in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
        
        self.run_fit('sed', **kwargs)  # fit BB on multiband photometry or spectra
        if kwargs['verbose']:
            print("Run SED fittings in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
        
        self.bb_bol(fastsedfitting=False, **kwargs)   # calculate luminosity from BB or spectra
        if kwargs['verbose']:
            print("Make lumnisoty LCs from SED fittings in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()                    
        
        ''' bolometric fit'''                       
        self.run_fit('bol_early', **kwargs)  # fit lums in early phases to shock cooling tail model
        if kwargs['verbose']:
            print("Run bolometric early model fittings in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
        
        self.run_fit('bol_main', **kwargs)  # fit lums around peak to arnett model
        if kwargs['verbose']:
            print("Run bolometric main peak model fittings in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.run_fit('bol_tail',**kwargs)    # fit lums at tail to gamma leakage model
        if kwargs['verbose']:
            print("Run bolometric tail model fittings in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.run_fit('bol_full',**kwargs)    # fit full range lums to some joint models
        if kwargs['verbose']:
            print("Run bolometric full model fittings in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
        
        ''' go to spectral fitting part '''        
        if 'spec' not in self.__dict__:
            if kwargs['verbose']: print("Skip spectral fittings since no spectra found")
            return
        
        self.run_fit('specline',**kwargs)  # measure specific lines, fit photos v to break arnett degenaracy
        if kwargs['verbose']:
            print("fit spec modelline in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()
            
        self.run_fit('specv_evolution',**kwargs)  # line velocity evolution
        if kwargs['verbose']:
            print("fit spec evolution in = {:.2f} s".format(time.time() - t_start))
            t_start = time.time()        
            
    def parse_coo(self, verbose=False, deg=True, hpx=False, nside=128):
        """ handle SN coordinate
        
        Parameters
        ----------                
        verbose :        `bool`
              show detailed running informations      
        deg :        `bool`
              return degrees, or hourse: ``hms``, ``dms``
        hpx :        `bool`
              if deg is False, return healpix index or not
        nside :        `bool`
              the resolution of healpix index

        Returns
        ----------

        - if ``deg`` = *True*:
        
           * ra:   `float`
        
           * dec:  `float`
        
        - if ``hpx`` = *True*:
        
           index:   `int`

        - otherwise:

           coordinate : `str`
        """
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
        if deg:
            return c.ra.deg, c.dec.deg
        elif hpx:
            import healpy as hp
            return hp.pixelfunc.ang2pix(nside,np.radians(-c.dec.deg+90.),np.radians(360.-c.ra.deg))
        else:
            return c.to_string('hmsdms').split()
        
    def mjd_now(self, jd=False):
        """ Get current juliand date via astropy.time
        
        Parameters
        ----------                
        jd :        `bool`
              get Julian date or midified Julian date
        
        Returns
        ----------

        - if ``jd`` = *True*:
        
           jd    :   `float`
        
        - if ``jd`` = *False*:
        
           mjd   :   `float`
        """
        if jd: return Time.now().jd
        else:  return Time.now().mjd
        
    def add_lc(self, df, source=None):
        """ Add a lightcurve into self.lc, with a given source name
        
        Parameters
        ----------                
        df :        `panda.dataframe`
              lightcurve data
        source :    `str`
              source name, e.g. ztffp for ZTF forced photometry, or myfavour for you own data

        Notes
        ---------- 
        ``df`` should include: ``jdobs``, ``filter``, ``mag``/``emag`` or ``flux``/``eflux``
        """ 
        # format lc keys:
        ck1 = 'mag' in df.keys() and 'emag' in df.keys() and 'jdobs' in df.keys() and 'filter' in df.keys()
        ck2 = 'flux' in df.keys() and 'eflux' in df.keys() and 'jdobs' in df.keys() and 'filter' in df.keys()
        assert ck1 or ck2, 'rename lc columns to make sure either mag/emag/jdobs/filter or flux/eflux/jdobs/filter are available'
        # format filters
        fs = []            
        for index in df.index:            
            fs.append( filter_clean(df['filter'][index]) )           
        df['filter'] = fs
        df = df.query('filter != "-"')        
        # save
        df['source'] = source
        if not 'lc' in self.__dict__: self.lc = df
        else: self.lc = self.lc.append(df, ignore_index=True).drop_duplicates()        
        
    def _add_flux(self, df, zp, sigma):
        """ (static) add **flux** and **eflux** column, based on **mag**, **emag**, and/or **limmag** column.
        
        Parameters
        ----------                
        df :        `panda.dataframe`
              lightcurve data
        zp :     `float`
              zeropoint to convert magnitude to flux
        sigma  :    `float`
              SNR threshold to distinguish detection/limit
        """ 
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
    
    def add_flux(self, source=None, **kwargs):
        ''' add **flux** and **eflux** column, based on **mag**, **emag**, and/or **limmag** column.

        Parameters
        ----------                
        source :    `str`
              source name, if None will add flux column for all data
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit        
        '''                       
        assert 'lc' in self.__dict__        
        assert 'mag' in self.lc.keys(), 'input mag'        
        kwargs = self.read_kwargs(**kwargs)
        
        if source is None:
            self.lc = self._add_flux(self.lc, kwargs['zp'], kwargs['snrt'])               
        else:
            df1 = self.lc.query('source==@source')
            df2 = self.lc.query('source!=@source')                        
            __arr = []
            if len(df1) > 0:                
                df1 = self._add_flux(df1, kwargs['zp'], kwargs['snrt'])                
                for i in df1.index: __arr.append( df1.loc[i].values )
                columns = df1.columns
            if len(df2) > 0:
                #df2 = self._add_flux(df2, kwargs['zp'], kwargs['snrt'])                
                for i in df2.index: __arr.append( df2.loc[i].values )
                columns = df2.columns            
            if len(__arr) > 0:
                self.lc = pd.DataFrame(__arr, columns=columns)
            else:
                print ('Error: no data found for source=%s'%source)

    def _add_mag(self, df, zp, sigma):
        """ (static) add **mag** and/or *limmag*/**emag** column, based on **flux**, **eflux** column.
        
        Parameters
        ----------                
        df :        `panda.dataframe`
              lightcurve data
        zp :     `float`
              zeropoint to convert flux to magnitude
        sigma  :    `float`
              SNR threshold to distinguish detection/limit
        """
        if 'eflux' in df.keys():
            df = df.query('eflux>0')
            df['mag'], df['emag'], df['limmag'] = self.flux_to_mag(df['flux'],
                    dflux=df['eflux'], sigma=sigma, zp=zp)        
        else:
            df['mag'] = self.flux_to_mag(df['flux'], sigma=sigma, zp=zp) 
        return df
    
    def add_mag(self, source=None, **kwargs):
        ''' add **mag** and/or *limmag*/**emag** column, based on **flux**, **eflux** column.

        Parameters
        ----------                
        source :    `str`
              source name, if None will add flux column for all data
        zp :     `float`
              zeropoint to convert flux to magnitude.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit        
        '''                       
        assert 'lc' in self.__dict__       
        assert 'flux' in self.lc.keys(), 'input flux'
        kwargs = self.read_kwargs(**kwargs)        
    
        if source is None:
            self.lc = self._add_mag(self.lc, kwargs['zp'], kwargs['snrt'])               
        else:
            df1 = self.lc.query('source==@source')
            df2 = self.lc.query('source!=@source')            
            df1 = self._add_mag(df1, kwargs['zp'], kwargs['snrt']) 
            df2 = self._add_mag(df2, kwargs['zp'], kwargs['snrt'])
            assert len(df1.columns) == len(df2.columns)
            __arr = []            
            for i in df1.index:            
                __arr.append( df1.loc[i].values )
            for i in df2.index:            
                __arr.append( df2.loc[i].values )                      
            self.lc = pd.DataFrame(__arr, columns=df1.columns)  
        
    def get_external_phot(self, filename, source='myfavour', **kwargs):
        ''' Parse and add user defined photometric data to ``snobject.lc``
        
        Parameters
        ----------                
        filename :   `str`
              path of the photometric file
        source :    `str`
              source name, if None will add flux column for all data        
        magkey    :  `str`
              meta key for AB magnitude
        emagkey   :  `str`
              meta key for magnitude error
        limmagkey :   `str`
              meta key for limiting magnitude
        fluxkey   :  `str`
              meta key for flux, unit in 1e-6 Jy
        efluxkey  :  `str`
              meta key for flux error
        filterkey :   `str`
              meta key for filters
        jdkey     :   `str`
              meta key for julian date
        zp :     `float`
              zeropoint to convert flux to magnitude.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit  
        
        Notes
        ----------
        After reading filename into ``pandas.dataframe``, ``snobject.add_lc`` will be called to add the dataframe into 
        ``snobject.lc``, and if ``mag`` or ``flux`` item missed, ``snobject.add_flux`` or ``snobject.add_mag`` will be called correspondingly.
        
        See Also
        ----------
        snobject.add_flux, snobject.add_mag, snobject.add_lc
        ''' 
        kwargs = self.read_kwargs(**kwargs)             
        if not os.path.exists(filename):
            print ('Error: %s not found'%filename)
            return None
        _ = []
        for nn,ll in enumerate(open(filename, encoding='utf-8').readlines()):
            if ll[0]=='#' or len(ll)==0:_.append(nn)
        df = pd.read_csv(filename,skiprows=_,delim_whitespace=True)

        # format lc keys
        df.rename(
            columns={
                kwargs['magkey']   : 'mag',
                kwargs['emagkey']  : 'emag',
                kwargs['limmagkey']: 'limmag',
                kwargs['fluxkey']  : 'flux',
                kwargs['efluxkey'] : 'eflux',
                kwargs['jdkey']    : 'jdobs',
                kwargs['filterkey']: 'filter',              
            },
            inplace=True,
        )
        self.add_lc(df, source=source)
        if 'mag' in df.keys() and not 'flux' in df.keys():            
            self.add_flux(source=source, **kwargs)
        if 'flux' in df.keys() and not 'mag' in df.keys():            
            self.add_mag(source=source, **kwargs)
        
    def bin_fp_atlas(self, binDays=3, resultsPath=None, outPath=None, verbose=False):
        ''' Bin ATLAS forced photometry with Dave's code,
        https://gist.github.com/thespacedoctor/86777fa5a9567b7939e8d84fd8cf6a76
        
        Parameters
        ----------                
        binDays :   `int`
              days bin
        resultsPath :   `str`
              path to store the binning plots
        outPath :   `str`
              path to store the binned data
        verbose:   `bool`
              show detailed running informations
        '''        
        cmd = 'python %s/plot_atlas_fp.py stack %.2f %s '%(srcpath, binDays, resultsPath)        
        cmd += '--o %s'%(outPath)
        if verbose: print (cmd)
        pid = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        output,error = pid.communicate()
        if verbose: print (output)
        if error: print (error)
        
    def get_fp_atlas(self, **kwargs):        
        ''' Parse and add local ATLAS forced (binned or not binned) photometric data to ``snobject.lc`` (with source name as **atlasfp**).
        
        Parameters
        ----------                
        tdbin :   `int`
              days bin
        clobber :   `bool`
              if cached file exists, redo or read the cache        
        zp :     `float`
              zeropoint to convert flux to magnitude.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit     
        
        Notes
        ----------

        - Download ATLAS forced photometry file in advance, with snobject.query_fp_atlas.

        - The ATLAS forced photometry service is open public, however in order to use it, one need to register first. And the user information should be put correctly in the auth.txt file.
        
        - The ATLAS photometry file will be read as ``pandas.dataframe``, and added to ``snobject.lc`` with ``snobject.add_lc``. Then ``snobject.add_flux`` or ``snobject.add_mag`` will be called if ``mag`` or ``flux`` items are missing.
        
        See Also
        ---------
        snobject.query_fp_atlas, snobject._get_fp_atlas, snobject.add_lc, snobject.add_mag, snobject.add_flux
        '''
        kwargs = self.read_kwargs(**kwargs)
        targetdir = '%s/ForcePhot_atlas/'%LOCALSOURCE
        wdir = '%s/%s/'%(targetdir, self.objid)
        f_unbin = '%s/forcedphotometry_%s_lc.csv'%(wdir, self.objid)        
        if kwargs['tdbin'] is None:            
            if not os.path.exists(f_unbin):
                print ('Error: %s not found'%f_unbin)
                return None            
            self._get_fp_atlas(f_unbin, **kwargs)
        else: 
            f = '%s/forcedphotometry_%s_lc_atlas_fp_stacked_%.2f_days.txt'%\
                (wdir, self.objid, kwargs['tdbin'])
            if not os.path.exists(f) or kwargs['clobber']:
                if not os.path.exists(f_unbin):
                    print ('Error: %s not found'%f_unbin)
                    return None
                self.bin_fp_atlas(binDays=kwargs['tdbin'],\
                    resultsPath=f_unbin, outPath=wdir, verbose=kwargs['verbose'])                        
            self._get_fp_atlas(f, binned=True, **kwargs)
            
    def _get_fp_atlas(self, f, binned=False, **kwargs):
        ''' (static) parse local ATLAS forced (binned or not binned) photometric data to ``snobject.lc``
        
        Parameters
        ----------                
        f :   `str`
              ATLAS file name
        binned :   `bool`
              if file is binned data or not     
        zp :     `float`
              zeropoint to convert flux to magnitude.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   
        
        See Also
        ---------
        snobject.add_lc, snobject.add_mag
        '''        
        kwargs = self.read_kwargs(**kwargs)
        sigma = kwargs['snrt']        
        if binned:
            df = pd.read_csv(f,skiprows=[0,1,2,3,4,5],sep = ',',).drop_duplicates().reset_index(drop=True)
        else:
            if '###' in open(f, encoding='utf-8').readlines()[0]:
                df = pd.read_csv(f,delim_whitespace=True).drop_duplicates().reset_index(drop=True)
                df.rename(columns={'###MJD':'MJD',}, inplace=True)
            else:
                df = pd.read_csv(f,sep = ',',).drop_duplicates().reset_index(drop=True)
        df.rename(columns={'uJy':'flux','duJy':'eflux','F':'filter',}, inplace=True)        
        df['jdobs'] = df['MJD'] + 2400000.5
        self.add_lc(df, source='atlasfp')
        self.add_mag(source='atlasfp', **kwargs)
        
    def get_oac(self, which='photometry', **kwargs):
        ''' Parse and add local Open Astronomical Catalog data to ``snobject.lc`` (with source name as ``oac``)
        
        Parameters
        ----------                
        which :   `str`
              which to parse: photometry or spectra        
        verbose   :   `bool`
              Enable progress report
        clobber   :   `bool`
              Redo analysis
        
        if which = photometry :
        
           zp :     `float`
              zeropoint to convert flux to magnitude.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
           snrt  :    `float`
              SNR threshold to distinguish detection/limit  

        if which = spectra :
        
           rv     :   `float`
              E(B-V) to AV, default is 3.1               
           bin_method  :  `str`
              method used to bin spectrum
           bin_size   :  `float`
              size used to bin spectrum, in AA
           savgol_order  : `float`
              polynomial order for savgol filter
           continuum_method  :  `str`
              The function type for continumm fitting, valid functions are "scalar", "linear", "quadratic", "cubic", "poly", and "exponential"
           continuum_degree  :  `float`
              degree of polynomial when method="poly", for continuum fitting
        
        Notes
        ----------

        - Download OAC photometry or spectra file in advance, with snobject.query_oac.
        
        - The OAC photometry file will be read as ``pandas.dataframe``, and added to ``snobject.lc`` with ``snobject.add_lc``. Then ``snobject.add_flux`` or ``snobject.add_mag`` will be called if ``mag`` or ``flux`` items are missing.        

        - The OAC spectra file will be read and pre-processed by specline_handler.handle_spectra
        
        - The OAC data are open public, and no registrations are needed in advance.
        
        See Also
        ----------   
        snobject.query_oac, snobject.add_lc, snobject.add_flux, snobject.add_mag, specline_handler.handle_spectra
        ''' 
        kwargs = self.read_kwargs(**kwargs)
        targetdir = '%s/oac/'%LOCALSOURCE        
        f = '%s/%s/%s_%s.csv'%(targetdir, self.objid, self.objid, which)
        if not os.path.exists(f):
            print ('Error: %s not found'%f)
            return None
        with open(f, 'r', encoding='utf-8') as ff:
            data = json.load(ff)[self.objid][which]
            
            if which == 'photometry':
                df = pd.DataFrame(data, columns = self.keys2query_lc.keys())
                df.rename(columns={'magnitude':'mag', 'e_magnitude':'emag', 'band': 'filter'}, inplace=True)
                df = df.astype({"time": float, "mag": float, 'emag':float})
                df['jdobs'] = df['time'] + 2400000.5 
                self.add_lc(df, source='oac')
                self.add_flux(source='oac', **kwargs)
            else:                
                for _data in data:                        
                    mjd, instru = _data[0], _data[1]                    
                    if len(_data[2][0]) == 2:
                        df = pd.DataFrame(_data[2], columns = ['wave', 'flux'])
                    elif len(_data[1][0]) == 3:
                        df = pd.DataFrame(_data[2], columns = ['wave', 'flux', 'eflux'])
                    else:
                        print ('skipped %s %s %s, check' % (self.objid, mjd, instru))
                        continue
                    df = df.astype({"wave": float, "flux": float})
                    if 'spec' not in self.__dict__:
                        self.spec = handle_spectra(self.z, self.mkwebv, t0=self.t0, **kwargs)
                        self.spec.add_spectrum(get_numpy(df['wave']), get_numpy(df['flux']), instru, mjd, **kwargs)         
                
    def query_oac(self, db='AL', which='photometry', **kwargs):
        ''' Query OAC data, via OACAPI (https://github.com/astrocatalogs/OACAPI)
        
        Parameters
        ----------                
        db :   `str`
        
              which catalog to query: 
              
              - **AL**  : https://api.astrocats.space/

              - **SN**  : https://api.sne.space/

              - **TDE** : https://api.tde.space/

              - **KN**  : https://api.kilonova.space/

              - **FT**  : https://api.faststars.space/

        which :   `str`

              what data to query: 

              - **photometry**

              - **spectra**
        
        clobber :    `bool`
              if file already exists, re-query it or not
        verbose:   `bool`
              show detailed running informations

        See Also
        ----------
        snobject.oac_phot_url, snobject.keys2query_spec, snobject.urllist, snobject.keys2query_lc
        '''
        kwargs = self.read_kwargs(**kwargs)       
        targetdir = '%s/oac/'%LOCALSOURCE        
        f = '%s/%s/%s_%s.csv'%(targetdir, self.objid, self.objid, which)
        if os.path.exists(f) and not kwargs['clobber']:
            if kwargs['verbose']: print ('file exists: %s'%f)
            return
        if not os.path.isdir( '%s/%s'%(targetdir, self.objid) ):
            os.mkdir( '%s/%s'%(targetdir, self.objid) )        
        assert db in self.urllist.keys()
        assert which in ['photometry', 'spectra']
        
        # base url                                
        url = '%s%s/%s/' % (self.urllist[db], self.objid, which)
        url = url.replace(' ','')
        if which == 'photometry':
            url = self.oac_phot_url(url)
        else:
            for _nkey, _key in enumerate(self.keys2query_spec):
                if _nkey == len(self.keys2query_spec)-1: url += '%s'%_key
                else: url += '%s+'%_key
                
        # complete url
        if kwargs['verbose']: print('OACAPI query URL: %s' % url)
        
        # download       
        urllib.request.urlretrieve(url, f)
        
    def oac_phot_url(self, url):
        ''' (static) make url to query OAC
        
        Parameters
        ----------                
        url :   `str`
              basic url

        See Also
        -----------
        snobject.keys2query_lc
        '''
        # query keys
        for _nkey, _key in enumerate(self.keys2query_lc):
            if _nkey == len(self.keys2query_lc)-1: url += '%s'%_key
            else: url += '%s+'%_key
        
        # constrains
        _init = True
        for _key in self.keys2query_lc:
            if self.keys2query_lc[_key]:
                if _init:
                    url += '?'
                    _init = False
                else:
                    url += '&'
                url += '%s'%_key
        return url
    
    def get_alert_ztf(self, source='marshal', **kwargs):
        ''' Parse and add local ZTF alert photometic data to ``snobject.lc``.
        
        Parameters
        ----------                
        source :   `str`
              which one to parse: from Growth [marshal], or [fritz]                      
        zp :     `float`
            zeropoint to convert flux to magnitude.
            zp = 23.9 for micro Jansky to AB mag
            zp = 48.6 for ergs/s/cm2/Hz to AB mag
            e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
            SNR threshold to distinguish detection/limit  
        
        Notes
        ----------
        
        - This function is dealing with ZTF alter photometry file from private brokers, i.e. the Growth marshal and the fritz, which are dedicated for the internal ZTF users. In order to use it, one need to define their user information in the auth.txt file.
        
        - Download ZTF alter photometry in advance, with snobject.query_alert_ztf.
        
        - The ZTF alert photometry file will be read as ``pandas.dataframe``, and added to ``snobject.lc`` with ``snobject.add_lc``. Then ``snobject.add_flux`` or ``snobject.add_mag`` will be called if ``mag`` or ``flux`` items are missing.        

        See Also
        ----------   
        snobject.add_lc, snobject.add_flux, snobject.add_mag, snobject.query_alert_ztf
        '''
        kwargs = self.read_kwargs(**kwargs)
        assert source in ['marshal', 'fritz']
        if source == 'marshal':
            from ztfquery import marshal
            df = marshal.get_local_lightcurves(self.objid)
        else:
            from ztfquery import fritz
            #df = fritz.download_lightcurve(self.objid, get_object=True)
            filename = fritz.FritzPhotometry._build_filename_(self.objid)
            df = fritz.FritzPhotometry.read_csv(filename).data
            df.rename(columns={'magerr':'emag', 'limiting_mag':'limmag'}, inplace=True)            
            for _k in df.keys(): # drop other columns
                if _k not in ['filter','limmag','mag','emag','mjd','instrument_name']:
                    df.drop(_k, inplace=True, axis=1)            
            df['jdobs'] = df['mjd'] + 2400000.5         
        mags,emags = df['mag'],df['emag']
        mags[np.isnan(mags)] = 99
        emags[np.isnan(emags)] = 99
        df['mag'] = mags
        df['emag'] = emags     
        self.add_lc(df, source=source)
        self.add_flux(source=source, **kwargs)
        
    def get_fp_ztf(self, seeing_cut = 7., **kwargs):
        ''' Parse and add local ZTF forced photometic data to ``snobject.lc``.
        
        Parameters
        ----------                
        seeing_cut :   `float`
              there's a ``seeing`` column in ztf forced photometry file,
              which can be used to remove epochs with quite poor observing conditions.
        zp :     `float`
            zeropoint to convert flux to magnitude.
            zp = 23.9 for micro Jansky to AB mag
            zp = 48.6 for ergs/s/cm2/Hz to AB mag
            e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
            SNR threshold to distinguish detection/limit  
        
        Notes
        ----------
        
        - Download ZTF forced photometry in advance, with snobject.query_fp_ztf.
        
        - The ZTF forced photometry service is open public, however in order to use it, one need to register first. And the user information should be put correctly in the auth.txt file.

        - The ZTF forced photometry file will be read as ``pandas.dataframe``, and added to ``snobject.lc`` with ``snobject.add_lc``. Then ``snobject.add_flux`` or ``snobject.add_mag`` will be called if ``mag`` or ``flux`` items are missing.        
        
        See Also
        ----------   
        snobject.add_lc, snobject.add_flux, snobject.add_mag, snobject.query_fp_ztf
        '''
        kwargs = self.read_kwargs(**kwargs)
        targetdir = '%s/ForcePhot/'%LOCALSOURCE
        sigma = kwargs['snrt']
        f = '%s/%s/forcedphotometry_%s_lc.csv'%(targetdir, self.objid, self.objid)
        if not os.path.exists(f):
            print ('Error: %s not found'%f)
            return None        
        _ = []
        for nn,ll in enumerate(open(f, encoding='utf-8').readlines()):
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
        self.add_lc(df, source='ztffp')
        self.add_mag(source='ztffp', **kwargs)
        
    def query_fp_atlas(self, **kwargs):
        ''' Qeury ATLAS forced photometry, see https://fallingstar-data.com/forcedphot/static/apiexample.py.
        
        Parameters
        ----------                        
        verbose :        `bool`
              show detailed running informations        
        clobber :        `bool`
              if ATLAS forced phot file exists, re-downloadit or not       
        mjdstart  :    `float`
              start julian date to query
        dstart  :    `float`
              if **mjdstart** is None, how many days prior than **t0** to query
        mjdend  :    `float`
              end julian date to query
        dend    :    `float`
              if **mjdend** is None, how many days later than **t0** to query        
        
        Notes
        ----------
        The ATLAS forced photometry service is open public, however in order to use it, one need to register first. And the user information should be put correctly in the auth.txt file.
        
        See Also
        ----------   
        snobject.query_fp_ztf, snobject.get_fp_atlas, snobject.parse_coo, snobject.mjd_now
        '''
        if len(keypairs) == 0:
            print ('Error: set atlas account in keypairs')
            return
        assert 'atlas' in keypairs.keys()
        BASEURL = "https://fallingstar-data.com/forcedphot"
        kwargs = self.read_kwargs(**kwargs)

        # coordinates
        radeg, decdeg = self.parse_coo(verbose=kwargs['verbose'], deg=True)

        # query period
        if kwargs['mjdstart'] is not None:
            mjdstart = float(kwargs['mjdstart'])
        else:
            assert self.t0 > 2400000, 'Error: specify mjdstart or t0 (in JD) first'
            mjdstart = float(self.t0) - 2400000.5 + float(kwargs['dstart'])
            if kwargs['verbose']: print ('mjdstart as %f days prior to the peak which is %.2f'%
                                         (float(kwargs['dstart']), self.t0-2400000.5))
        if kwargs['mjdend'] is not None:
            mjdend = float(kwargs['mjdend'])
        else:
            if kwargs['dend'] is None:
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
        ''' Qeury ZTF forced photometry, see documentation: https://irsa.ipac.caltech.edu/data/ZTF/docs/forcedphot.pdf,
        and the web query: https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi.
        
        Parameters
        ----------                        
        get_email :        `bool`
              ZTF force phot query service is time consuming. 
              One can set ``get_email`` as *True*, to receive an email including data file
              instead of query via API that would stuck the GUI for a while. 
        verbose :        `bool`
              show detailed running informations        
        clobber :        `bool`
              if ZTF forced phot file exists, re-downloadit or not       
        mjdstart  :    `float`
              start julian date to query
        dstart  :    `float`
              if **mjdstart** is None, how many days prior than **t0** to query
        mjdend  :    `float`
              end julian date to query
        dend    :    `float`
              if **mjdend** is None, how many days later than **t0** to query
        
        Notes
        ----------
        The ZTF forced photometry service is open public, however in order to use it, one need to register first. And the user information should be put correctly in the auth.txt file.
        
        See Also
        ----------   
        snobject.query_fp_atlas, snobject.get_fp_ztf, snobject.parse_coo, snobject.mjd_now
        '''
        if len(keypairs) == 0: return
        assert 'ztf' in keypairs.keys()
        kwargs = self.read_kwargs(**kwargs)
        
        # coordinates
        radeg, decdeg = self.parse_coo(verbose=kwargs['verbose'], deg=True)

        # query period                    
        if kwargs['mjdstart'] is not None:
            jdstart = float(kwargs['mjdstart']) + 2400000.5
        else:
            assert self.t0 > 2400000, 'Error: specify mjdstart or t0 (in JD) first'
            jdstart = float(self.t0) + float(kwargs['dstart'])
            if kwargs['verbose']: print ('jdstart as %f days prior to the peak which is %.2f'%
                                         (float(kwargs['dstart']), self.t0))                
        if kwargs['mjdend'] is not None:
            jdend = float(kwargs['mjdend']) + 2400000.5
        else:
            if kwargs['dend'] is None:
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
            time.sleep(10)
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
        ''' Qeury ZTF alert photometry via ztfquery, see https://github.com/MickaelRigault/ztfquery.        
        
        Parameters
        ----------                        
        source :        `str`
              **fritz**, or **marshal**, or None for both
        verbose :        `bool`
              show detailed running informations        

        Notes
        ----------        
        This function is dealing with ZTF alter photometry file from private brokers, i.e. the Growth marshal and the fritz, which are dedicated for the internal ZTF users. In order to use it, one need to define their user information in the auth.txt file.              
        
        See Also
        ----------   
        snobject.get_alert_ztf, ztfquery
        '''
        kwargs = self.read_kwargs(**kwargs)
        assert source in ['marshal', 'fritz', None], 'source should be gm/fritz or None for both'
        if source in ['marshal',None]:
            from ztfquery import marshal
            marshal.download_lightcurve(self.objid,
                    auth=(keypairs['marshal']['usr'], keypairs['marshal']['pwd']),
                    verbose=kwargs['verbose'])
        if source in ['fritz',None]:
            from ztfquery import fritz
            fritz.download_lightcurve(self.objid, store=True, verbose=kwargs['verbose'])
        
    def correct_baseline(self, baseline, key='fcqfid', source='ztffp', **kwargs):
        ''' Correct baseline for ZTF, when targets was covered by multiple CCDs.
        
        Parameters
        ----------                        
        baseline :        `dictionary`
              data returned from **self.calibrate_baseline**
        key :        `str`
              which column to distinguish photometry between different fields/ccds
        source :     `str`
              which source of lc to be corrected
        zp :     `float`
            zeropoint to convert flux to magnitude.
            zp = 23.9 for micro Jansky to AB mag
            zp = 48.6 for ergs/s/cm2/Hz to AB mag
            e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
            SNR threshold to distinguish detection/limit

        Notes
        ----------
        This function was only tested with ZTF forced photometry.
        
        See Also
        ----------
        snobject.calibrate_baseline
        '''
        kwargs = self.read_kwargs(**kwargs)
        __arr = []
        columns = self.lc.columns
        lc = self.lc.query('source==@source')
        lc0 = self.lc.query('source!=@source')         
        if len(lc0) > 0:
            for i in lc0.index: __arr.append( lc0.loc[i].values )                         
        for fid in baseline:
            baselevel = baseline[fid]            
            _lc = lc.query('%s==@fid'%key)            
            _lc['flux'] -= baselevel            
            _m = []
            for __m, __f in zip(_lc['mag'], _lc['flux']):
                if __m < 99:
                    __m = self.flux_to_mag(
                        __f, dflux=None, sigma=kwargs['snrt'],
                        units='zp', zp=kwargs['zp'], wavelength=None
                    )
                _m.append(__m)
            _lc['mag'] = _m
            for i in _lc.index: __arr.append( _lc.loc[i].values )               
        self.lc = pd.DataFrame(__arr, columns=columns)         
        
    def calibrate_baseline(self, ax=None, key='fcqfid', source='ztffp',
                        xmin=-100, xmax=-20, ax_xlim=None, ax_ylim=None):
        ''' Calculate ZTF forced photometry baseline from different fields with different filters.
        
        Parameters
        ----------              
        ax   :          `matplotlib.axes`
              matplotlib subplot, is None will not show        
        key :        `str`
              which column to distinguish photometry between different fields/ccds
        source :     `str`
              which source of lc to be corrected
        xmin   :     `float`
              left side of LC to be used to calibrate the baseline
        xmax   :     `float`
              right side of LC to be used to calibrate the baseline
        ax_xlim  :    `range`
              x axis range for the plot
        ax_ylim  :    `range`
              y axis range for the plot
        
        Returns
        -----------
        baseline :    `dict`
              estimated baseline flux levels
        '''
        assert xmin < xmax
        assert self.t0 > 2400000, 'set t0 first'
        xp = self.t0
        tb = self.lc
        if source is not None:
            tb = tb.query("source==@source")
        assert len(tb) > 0
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
            else: ax.set_xlim([xp+xmin-20-2458000, xp+xmax+20-2458000])
        return baseline

    def get_external_spectra(self, filename, epoch, tel='', **kwargs):
        ''' Parse and add user defined photometric data to ``snobject.spec``.
        
        Parameters
        ----------              
        filename   :          `str`
              filename.
              All lines starting with **#** will be skipped.
              The first line should be the header, at least including: ``wave`` and ``flux`` columns,
              while in the following lines are the data. All should be seperated with spaces.
        epoch :        `str`
              astropy.time, that will be used to calculate the phase
        tel :     `str`
              telescope/instrument
        rv     :   `float`
              E(B-V) to AV, default is 3.1               
        bin_method  :  `str`
              method used to bin spectrum
        bin_size   :  `float`
              size used to bin spectrum, unit in ``A``
        savgol_order  : `float`
              polynomial order for ``savgol`` filter
        continuum_method  :  `str`        
              The function type for continumm fitting, valid functions are:
        
              - ``scalar``

              - ``linear`` 
             
              - ``quadratic``

              - ``cubic``
        
              - ``poly``
        
              - ``exponential``
        
        continuum_degree  :  `float`
              degree of polynomial when ``continuum_method`` = **poly**, for continuum fitting
        
        Notes
        ------------        
        The spectra file will be read and pre-processed by specline_handler.handle_spectra        

        See Also
        ------------
        specline_handler.handle_spectra, snobject.get_local_spectra
        '''
        kwargs = self.read_kwargs(**kwargs)
        if not os.path.exists(filename):
            print ('Error: %s not found'%f)
            return None
        _ = []
        for nn,ll in enumerate(open(filename, encoding='utf-8').readlines()):
            if ll[0]=='#' or len(ll)==0:_.append(nn)
        df = pd.read_csv(filename,skiprows=_,delim_whitespace=True)
        assert 'wave' in df.keys() and 'flux' in df.keys()        
        if 'spec' not in self.__dict__:
            self.spec = handle_spectra(self.z, self.mkwebv, t0=self.t0, **kwargs)
        self.spec.add_spectrum(get_numpy(df['wave']), get_numpy(df['flux']), tel, epoch, **kwargs)        
        
    def query_spectra(self, source=None, **kwargs):
        ''' Qeury ZTF spectra via ztfquery, see https://github.com/MickaelRigault/ztfquery.
        
        Parameters
        ----------                        
        source :        `str`
              **fritz**, or **marshal**, or None for both
        verbose :        `bool`
              show detailed running informations        
        
        Notes
        ----------        
        This function is dealing with ZTF spectral file from private brokers, i.e. the Growth marshal and the fritz, which are dedicated for the internal ZTF users. In order to use it, one need to define their user information in the auth.txt file.
        
        See Also
        ----------   
        snobject.get_local_spectra, ztfquery
        '''
        kwargs = self.read_kwargs(**kwargs)
        assert source in ['marshal', 'fritz', None]
        if source in ['marshal',None]:
            from ztfquery import marshal
            marshal.download_spectra(self.objid)            
        if source in ['fritz',None]:
            from ztfquery import fritz
            fritz.download_spectra(self.objid, get_object=True, store=True, verbose=kwargs['verbose'])   
        
    def get_local_spectra(self, source=None, **kwargs):
        ''' Parse amd add local ZTF fritz/marshal spectra to ``snobject.spec``.
        
        Parameters
        ----------                        
        source :        `str`
              **fritz**, or **marshal**, or ``None`` for both
        verbose :        `bool`
              show detailed running informations         
        rv     :   `float`
              E(B-V) to AV, default is 3.1
        bin_method  :  `str`
              method used to bin spectrum
        bin_size   :  `float`
              size used to bin spectrum, in ``A``
        savgol_order  : `float`
              polynomial order for savgol filter
        continuum_method  :  `str`        
              The function type for continumm fitting, valid functions are:
        
              - ``scalar``

              - ``linear`` 
             
              - ``quadratic``

              - ``cubic``
        
              - ``poly``
        
              - ``exponential``
        
        continuum_degree  :  `float`
              degree of polynomial when ``continuum_method`` = **poly**, for continuum fitting

        Notes
        -------------
        - Download ZTF spectra file in advance, with snobject.query_spectra.
        
        - The ZTF spectra file will be read and pre-processed by specline_handler.handle_spectra
        
        - This function is dealing with ZTF private spectra file from private brokers, i.e. the Growth marshal and the fritz, which are dedicated for the internal ZTF users. In order to use it, one need to define their user information in the auth.txt file.        
        
        See Also
        ------------
        specline_handler.handle_spectra
        '''      
        kwargs = self.read_kwargs(**kwargs)
        assert source in ['marshal', 'fritz', None]
        if 'spec' not in self.__dict__:
            self.spec = handle_spectra(self.z, self.mkwebv, t0=self.t0, **kwargs)                      
        
        if source in ['marshal', None]:
            from ztfquery import marshal, fritz
            # spec from marshal
            dir_ = marshal.target_spectra_directory(self.objid)
            if os.path.exists( dir_ ):
                for d in os.listdir( dir_ ): 
                    objid,epoch,tel = d.split('_')[0],d.split('_')[1],d.split('_')[2]
                    try: 
                        data, header = fritz.parse_ascii(
                            open( os.path.join(dir_,d), encoding='utf-8' ).read().splitlines() )                        
                    except:
                        if kwargs['verbose']:
                            print ( 'Error: failed to parse spectra from %s' % os.path.join(dir_,d) )
                        continue
                    self.spec.add_spectrum(get_numpy(data['lbda']), get_numpy(data['flux']), tel, epoch, **kwargs)
                    
        if source in ['fritz', None]:
            from ztfquery import fritz
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
        
    def rapid(self, ax=None, returnv=False, **kwargs):
        ''' Run astrorapid codes (a Deep learning classifier to distiinguish LCs between different SNe type) on ``snobject.lc``.
        
        Parameters
        ----------                        
        ax   :          `matplotlib.axes`
              matplotlib subplot, is None will not show          
        clobber :        `bool`
              if astrorapid plot exists, re-do or load it
        verbose :        `bool`
              show detailed running informations    
        returnv  :        `bool`
               if return the Machine learning scores (i.e. ``ml``), or store them in the class
        
        Notes
        ----------
        
        This is only available when ``astrorapid`` packagte is correcly installed.        
        
        See Also
        ----------
        snobject.rapid_plot, snobject.parse_coo
        ''' 
        kwargs = self.read_kwargs(**kwargs)        
        try:
            import astrorapid
            from astrorapid.classify import Classify
        except:
            print ( 'need to install astrorapid' )            
            return
        class_names = ('Pre-explosion', 'SNIa-norm', 'SNIbc', 'SNII',
                       'SNIa-91bg', 'SNIa-x', 'Kilonova', 'SLSN-I', 'TDE')
        class_color = astrorapid.classify.CLASS_COLOR
        if 'ml' in self.__dict__ and not kwargs['clobber']:
            if ax is not None:
                self.rapid_plot( ax, class_names, class_color, figdir='%s/plots/'%LOCALSOURCE )
            return        
        ra, dec = self.parse_coo(verbose=kwargs['verbose'], deg=True)
        lc = self.lc.query('filter in ["r","i","g"]')        
        photflag, fd = [], True
        for _ in lc.index:
            if lc['mag'][_] == 99:
                photflag.append(0)
            elif fd:
                photflag.append(6144)
                fd = False
            else:
                photflag.append(4096)
        # flux zp from 23.9 (uJy) to 26.2
        light_curve_info1 = (lc['jdobs'].to_numpy(), lc['flux'].to_numpy()*10.,
                             lc['eflux'].to_numpy()*10., lc['filter'].to_numpy(),
                             photflag, ra, dec, self.objid, self.z, self.mkwebv)        
        # Classification
        light_curve_list = [light_curve_info1,]
        classification = Classify(known_redshift=True, class_names=class_names)                
        predictions = classification.get_predictions(light_curve_list, return_objids=False)
        _pred, _time = predictions                
        if _pred is None:
            ml = None            
        else:
            ml = []                
            for i,j in zip(_pred[0], _time[0]):
                ii=[]
                for _i in i: ii.append(float(_i))
                ml.append({'time': j, 'predict':ii})        
        if returnv: return ml
        self.ml = ml
        if ax is not None:
            self.rapid_plot( ax, class_names, class_color, figdir='%s/plots/'%LOCALSOURCE )
        
    def rapid_plot(self, ax, class_names, class_color, figdir='./'):
        ''' (static) Make astrorapid plot.
        
        Parameters
        ----------                        
        ax   :          `matplotlib.axes`
              matplotlib subplot, is None will not show          
        class_names :     `list`
              different transient classification type names
        class_color  :          `matplotlib.axes`
              colors of them
        figdir   :          `matplotlib.axes`
              path to store the figure
        
        Notes
        -------------
        working when ``snobject.ml`` available, and is not **None**.
        
        See Also
        -------------
        snobject.rapid
        ''' 
        assert 'ml' in self.__dict__
        assert self.ml is not None
        plt.rcParams["text.usetex"] = False                
        for classnum, classname in enumerate(class_names):
            tt, pp = [], []
            for _epoch in self.ml:
                phase = _epoch['time']
                pred = _epoch['predict']                
                tt.append(phase)
                pp.append(pred[classnum])
            tt, pp = np.array(tt), np.array(pp)
            ax.plot(tt-self.t0, pp, '-', label=classname, color=class_color[classname], linewidth=3)
        ax.legend(frameon=True, fontsize=8)
        ax.set_xlabel("Days since trigger (rest frame)")        
        ax.set_ylabel("Class Probability")       
        ax.set_ylim(0, 1)        
        ax.grid(True)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune='upper'))
        plt.tight_layout()
        savename = 'rapid_classification_%s.png'%self.objid
        plt.savefig(os.path.join(figdir, savename), dpi=400, bbox_inches='tight')
        
    def run_gp(self, source=None, **kwargs):
        ''' Run Gaussian Process interpolation via george package on multiband lightcurves, i.e. ``snobject.lc``.
        
        Parameters
        ----------                        
        source :        `str`
              which source of ``snobject.lc`` to be used
        gp_fit    :    `bool`
              if fit GP or not
        gp_redo   :    `bool`
              if cached file exists, redo or read
        gp_bands  :    `list`
              which bands for GP, None, or a list. e.g. ['g', 'r']
        gp_routine:    `str`
              Possible choices are: 

              - ``minimize``

              - ``mcmc``

              - ``leastsq``
        
        gp_mean   :    `str`
              Mean y_data function. Possible choices are: 

              - ``mean``
         
              - ``gaussian``
        
              - ``bazin``
        
              - ``villar``
        
        gp_xrange :  `list`
              gp fit range, relative to peak
        gp_xrangep:  `list`
              gp plot range, relative to peak. if None, will not plot gp samplings
        verbose  :      `list`
              show warnings or not
        kernel    :  `str`
              Kernel to be used with the gaussian process. Possible choices are: 

              - ``matern52``

              - ``matern32``

              - ``squaredexp``
        
        fix_scale :  `bool`
              If fix default gaussian process param        
        nwalkers  :  `int`
              number of walkers
        nsteps    :  `int`
              number of MC steps
        nsteps_burnin : `int`
              number of MC burn in steps               
        set_tpeak_method :  `str`
              set t0 with peak from:
        
              - [gp] Gaussian process

              - [fit] multiband_main/Bazin fits

              - [bol] bol_main/Arnett fits

              - leave blanckt as with input jdpeak
        
        set_tpeak_filter :  `str`
              set t0 with peak of which band
        
        Notes
        ---------------
        
        - This function is using ``sdapy.gaussian_process.fit_gp``, which was afterwards stored as ``snobject.gpcls``.
        
        - This function is using a constant kernel for wavelengthes and a demanded kernel for fluxes.
        
        See Also
        ---------------
        snobject.set_peak_gp, snobject.run_fit
        ''' 
        kwargs = self.read_kwargs(**kwargs)        
        if not kwargs['gp_fit']: return
        
        lc = self.lc
        if kwargs['gp_bands'] is not None: gpfs = kwargs['gp_bands']
        else: gpfs = np.unique(lc['filter'])
        lc = lc.query('filter in @gpfs')
        if source is not None: lc = lc.query('source==@source')       
        
        if self.t0 > 0: # cut lc first
            pmin,pmax = min(kwargs['gp_xrange']), max(kwargs['gp_xrange'])
            lc = lc.query('jdobs<=@self.t0+@pmax and jdobs>=@self.t0+@pmin')        
        self.gpcls = fit_gp(np.array(lc['jdobs']), np.array(lc['flux']),
                    np.array(lc['eflux']), filters=np.array(lc['filter']))        
        self.gpcls.train(
            kernel=kwargs['kernel'], fix_scale=kwargs['fix_scale'],
            gp_mean=kwargs['gp_mean'], opt_routine=kwargs['gp_routine'],
            nwalkers=kwargs['nwalkers'], nsteps=kwargs['nsteps'],
            nsteps_burnin=kwargs['nsteps_burnin'], clobber=kwargs['gp_redo'],
            mcmc_h5_file='gp_%s_%s_%s'%(self.objid,kwargs['gp_routine'],kwargs['gp_mean']),
            verbose=kwargs['verbose'], datadir='%s/cache/'%LOCALSOURCE,
            t0=self.t0, timedilation=1, xpredict=None, #kwargs['gp_xrangep'],
        )
        self.gpcls.predict()
        self.gpcls.set_peak()
        
        # update t0 and fpeak
        if len(kwargs['set_tpeak_method'])==0 and (self.t0 ==0 or len(self.fpeak)==0):
            self.set_peak_gp(**kwargs)
        elif kwargs['set_tpeak_method'] == 'gp':
            self.set_peak_gp(**kwargs)
            
    def run_fit(self, enginetype, source=None, **kwargs):
        ''' Run model fittings on varies of data from different engines.
        
        Parameters
        ----------                        
        enginetype :        `str`
              which engine to be used, e.g. ``multiband_main``, ``bol_main``, etc.
              Check all engines at https://github.com/saberyoung/HAFFET/tree/master/sdapy/engines,
              or one can define their own engine.
        source :       `str`
              **optional**, which source to be used for fitting, e.g. for ``bol_main`` engine, if there're multiple
              bolometric LCs available, fit on which of them, or for ``sed`` engine, which epochs to be fitted, etc.        
        fit_methods :   `str`
              which fitting model to be used, e.g. bazin, gauss, etc
        fit_redo    :   `bool`
              if False by default, when cached file exists, read samples from cached sample file. for mcmc, 
              if given nsteps is larger than that of cached sample, continue fits to nsteps, otherwise return samples. 
              if True, when cached file exists, redo everything for new fitting samples        
        
        for engine ``x`` :
        
        x_type   :  `str`
               what data to fit
        x_xrange :  `list`
               x range to fit the data        
        x_xrangep :  `list`
               x range to predict the data  
        x_yrange :  `list`
               y range to fit the data
        x_bands  : `list`
               data from which bands for the fit
        x_routine: `str`
               Which technic to be used to realize optimization. Possible choices are: minimize, mcmc, leastsq

        if ``routine`` is **mcmc** :

        ncores    :  `int`
               how many cores to run multi-processing
        nwalkers  :   `int`
               number of walkers
        nsteps    :   `int`  
               number of MC steps
        nsteps_burnin : `int` 
               number of MC burn in steps        
        thin_by   :   `int` 
               If you only want to store and yield every thin_by samples in the chain, set thin_by to an integer greater than 1. When this is set, iterations * thin_by proposals will be made.              
        emcee_burnin : `bool`
               if use emcee to burnin chains
        use_emcee_backend : `bool`
               if use emcee backend
        quantile   :  `list` 
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        
        if ``routine`` is **minimize** :

        maxfev    :  `int` 
               for scipy. The maximum number of calls to the function.    
        scipysamples : `int` 
               generated sampling numbers for scipy approach
        quantile   :  `list` 
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        
        Examples
        --------
        Define a ``snobject`` object
        
        >>> a = snerun.snobject(objid='ZTF20aajcdad', aliasid='SN2020bcq',
        ...        z=0.0186, dm=34.6, mkwebv=0.01387, hostebv=0, sntype='Ib', 
        ...        ra='13:26:29.65', dec='+36:00:31.1', jdpeak=2458888.02)    

        Run ``bz`` model on **g**, **r**, and **i** band data prepared by the ``multiband_main`` engine, with ``minimize`` routine.
        
        >>> a.run_fit('multiband_main', fit_methods=['bz'],
        ...       multiband_main_bands=['g','r','i'], multiband_main_routine='minimize')
        
        Run ``pl`` model on **g**, **r**, and **i** band data prepared by the ``multiband_early`` engine, with ``mcmc`` routine.
        For the ``mcmc``, run maximum 5000 steps for burnin process, and after burnin, continue 30000 steps as maximum for samplings.
        
        >>> a.run_fit('multiband_early', fit_methods=['pl'],
        ...       multiband_early_bands=['g','r','i'], multiband_early_routine='mcmc',
        ...       nsteps_burnin=5000, nsteps=30000) 
        
        Run ``bb`` model on SED data prepared by ``snobject.bb_colors`` and the ``sed`` engine, with ``minimize`` routine.
        ``fit_redo`` will remove the samplings obtained before.
        
        >>> a.run_fit('sed', fit_methods=['bb'], sed_routine='minimize', fit_redo=True)
        
        Run ``arnett_fit_taum`` model on bolometric LCs that were prepared by the ``bol_main`` engine, with ``minimize`` routine.
        The bolometric LCs for model fittings are selected from -18 to 40 rest frame days relative to peak epoch. 
        After generating a number of samplings, the LC realization would be from -18 to 80 days.        
        
        >>> a.run_fit('bol_main', fit_methods=['arnett_fit_taum'], fit_redo=True,
        ...       bol_main_xrange=[-18, 40],  bol_main_xrangep=[-18, 80],  bol_main_routine='minimize')
        
        Run ``gauss`` model on spectral line minima prepared by the ``specline`` engine.
        
        >>> a.run_fit('specline', fit_methods=['gauss'])
        
        Run ``exp`` model on spectral line minima prepared by the ``specv_evolution`` engine.
        
        >>> a.run_fit('specv_evolution', fit_methods=['exp'])        
        
        Notes
        ------------
        kwargs :     `Keyword Arguments`
            see more parameters in https://github.com/saberyoung/HAFFET/blob/master/sdapy/data/default_par.txt,
            the **snobject** part.
        
        ``kwargs`` are different for different engines:
        
        - for multi_main and multi_early engine:
           set_tpeak_method :  `str`
              set t0 with peak from [gp] Gaussian process or [fit] multiband_main/Bazin fits or [bol] bol_main/Arnett fits or leave blanckt as with input jdpeak
           set_tpeak_filter :  `str`
              set t0 with peak of which band
           set_texp_method  :   `str`
              set texp with [fit] power law fits or [bolmain] bolometric LC fits or [bolearly] shock cooling fits or leave blanckt as with middle of last nondetection to the first detection
           set_texp_filter  :  `str`
              set texp with peak of which band
        
        - for spectral engines :
        spec_snr          :  `int` 
               if no error from spectra, add noise level of snr
        bin_method        :  `str`      
               method used to bin spectrum
        bin_size          :  `int`  
               size used to bin spectrum, in AA
        savgol_order      :  `int` 
               polynomial order for savgol filter
        continuum_method  :  `str`       
               The function type for continumm fitting, valid functions are "scalar", "linear", "quadratic", "cubic", "poly", and "exponential"
        continuum_degree  :  `int`       
               degree of polynomial when method="poly", for continuum fitting
        pfactor           :  `int`         
               threshold used for peak detection
        sn_line           :  `str`           
               which line to fit, e.g. 'H~$\alpha$', 'He~5876$\AA$', 'O~7774$\AA$'
        specfit_phase     :  `list` 
               phase range to decide which spectra should be fitted        
        v_p               :  `float` 
               guessed velocity (unit: 1e3 km/s)
        v_bounds          :  `list` 
               guessed velocity range (unit: 1e3 km/s)        

        - for sed engines:
           sed_phaserange     : `list` 

           sed_abscal_bands   : `list` 
               which photometry bands to be used to absolute calibrate the spectra
           sed_abscal_method  : `str` 

               * [cal]ibrate
        
               * [mangle] the spetra
        '''        
        kwargs = self.read_kwargs(**kwargs)        
        for which in kwargs['fit_methods']:                   
            _which = get_pars(which)
            engine = _which['engine']
            enginename = _which['enginename']            
            if enginename != enginetype: continue
            engine(self, which, enginename, sourcename=source, **kwargs)
            
    def set_peak_gp(self, returnv=False, **kwargs):
        ''' Set peak time and fluxes with GP interpolation.
        
        Parameters
        ----------                        
        set_tpeak_filter :  `str`
                set t0 with peak of which band
        returnv  :        `bool`
               if return the peak infos, or store them in the class
        verbose :        `bool`
              show detailed running informations 
        
        Notes
        -----------
        Only working when GP on ``set_tpeak_filter`` band LC had been done.
        
        See Also
        -----------
        snobject.set_peak_bol_main, snobject.set_peak_multiband_main
        '''        
        assert 'gpcls' in self.__dict__
        kwargs = self.read_kwargs(**kwargs)
        filt = kwargs['set_tpeak_filter']
        if not filt in self.gpcls.tpeak:
            if kwargs['verbose']:
                print ('Warning: by default %s is used to set t0, however not available'%filt)
            # use the one closet instead
            cw0 = central_wavelengths[filt]
            fs, ds = [], []
            for f in np.unique(self.lc['filter']):
                fs.append( f )
                ds.append( abs(central_wavelengths[f]-cw0) )
            filt = fs[np.argmin(ds)]
            if kwargs['verbose']:
                print ('       use %s (whose central wavelength is closest) peak as t0 instead'%filt)
        if returnv: return self.gpcls.tpeak[filt], self.gpcls.fpeak
        self.t0 = self.gpcls.tpeak[filt]
        self.fpeak = self.gpcls.fpeak
        
    def set_peak_bol_main(self, model_name=None, source_name=None, returnv=False):
        ''' Set peak time and fluxes with bolometric LC fittings.
        
        Parameters
        ----------                        
        model_name :        `str`
               optional, which bolometric model to be used
        source_name :        `str`
               optional, which source to be used, e.g. mbol or mbolbb
        returnv  :        `bool`
               if return the peak infos, or store them in the class

        Notes
        -----------
        Only working when model fittings with ``bol_main`` engine had been done.
        
        See Also
        -----------
        snobject.set_peak_gp, snobject.set_peak_multiband_main
        ''' 
        assert 'fitcls' in self.__dict__
        assert 'bol_main' in self.fitcls
        for source in self.fitcls['bol_main']:
            # mbol or mbolbb
            if source_name is not None and source != source_name: continue 
            for model in self.fitcls['bol_main'][source]:
                if model_name is not None and model != model_name: continue
                if returnv:
                    return self.fitcls['bol_main'][source][model].tpeak, \
                        self.fitcls['bol_main'][source][model].fpeak
                self.t0 = self.fitcls['bol_main'][source][model].tpeak
                self.fpeak[source] = self.fitcls['bol_main'][source][model].fpeak
                
    def set_peak_multiband_main(self, model_name=None, returnv=False, **kwargs):
        ''' Set peak time and fluxes with multiband LC fittings.
        
        Parameters
        ----------                        
        set_tpeak_filter :  `str`
            set t0 with peak of which band
        model_name :        `str`
            optional, which multiband model to be used
        returnv  :        `bool`
            if return the peak infos, or store them in the class

        Notes
        -----------
        Only working when model fittings with ``multiband_main`` engine had been done.
        
        See Also
        -----------
        snobject.set_peak_gp, snobject.set_peak_bol_main
        '''        
        assert 'fitcls' in self.__dict__
        assert 'multiband_main' in self.fitcls
        kwargs = self.read_kwargs(**kwargs)
        filt = kwargs['set_tpeak_filter']
        for model in self.fitcls['multiband_main']:
            if model_name is not None and model != model_name: continue            
            assert filt in self.fitcls['multiband_main'][model].tpeak,\
                '%s not available by multiband_main.%s'%(filt, model)
            if returnv:
                return self.fitcls['multiband_main'][model].tpeak[filt], \
                    self.fitcls['multiband_main'][model].fpeak
            self.t0 = self.fitcls['multiband_main'][model].tpeak[filt]
            self.fpeak = self.fitcls['multiband_main'][model].fpeak
        
    def set_vexp(self, model_name=None, returnv=False, **kwargs):
        ''' Return peak velocity from models with ``specv_evolution`` engine.
        
        Parameters
        ----------         
        model_name :        `str`
               optional, which multiband model to be used
        quantile   :  `list` 
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors  
        returnv  :        `bool`
               if return the peak infos, or store them in the class
        
        Notes
        -----------
        Only working when model fittings with ``specv_evolution`` engine had been done.        
        '''
        kwargs = self.read_kwargs(**kwargs)
        assert 'fitcls' in self.__dict__
        assert 'specv_evolution' in self.fitcls
        for model in self.fitcls['specv_evolution']:
            if model_name is not None and model != model_name: continue             
            x, y, y1, y2, _ = self.fitcls['specv_evolution'][model].predict(x_pred=0, returnv=True, quant=kwargs['quantile'])
            if returnv: return float(y1),float(y),float(y2)
            self.vexp = [float(y1),float(y),float(y2)]
        
    def set_texp_pl(self, model_name=None, returnv=False, **kwargs):
        ''' Set the explosion epoch with ``multiband_early`` LC fittings.
        
        Parameters
        ----------                        
        set_texp_filter  :  `str`
              set texp with peak of which band
        model_name :        `str`
               optional, which multiband model to be used
        returnv  :        `bool`
               if return the peak infos, or store them in the class

        Notes
        -----------
        Only working when model fittings with ``multiband_early`` engine had been done. 
        
        See Also
        -------------
        snobject.set_texp_bol_main, snobject.set_texp_midway
        ''' 
        assert 'fitcls' in self.__dict__
        assert 'multiband_early' in self.fitcls
        kwargs = self.read_kwargs(**kwargs)
        filt = kwargs['set_texp_filter']
        for model in self.fitcls['multiband_early']:
            if model_name is not None and model != model_name: continue            
            assert filt in self.fitcls['multiband_early'][model].filters,\
                '%s not available by multiband_early.%s'%(filt, model)
            t1,t0,t2 = self.fitcls['multiband_early'][model].get_par(filt=filt, parname='texp')
            if returnv: return t1,t0,t2
            self.texp = [t1,t0,t2]            
        
    def set_texp_bol_main(self, model_name=None, source_name=None, returnv=False):
        ''' Set the explosion epoch with ``bol_main`` LC fittings.
        
        Parameters
        ----------                        
        model_name :        `str`
               optional, which bolometric model to be used
        source_name :        `str`
               optional, which source to be used, e.g. mbol or mbolbb
        returnv  :        `bool`
               if return the peak infos, or store them in the class

        Notes
        -----------
        Only working when model fittings with ``bol_main`` engine had been done. 
        And there's the ``texp`` parameter as free parameter in the model.
        
        See Also
        -------------
        snobject.set_texp_pl, snobject.set_texp_midway
        '''
        assert 'fitcls' in self.__dict__
        assert 'bol_main' in self.fitcls
        for source in self.fitcls['bol_main']:
            # mbol or mbolbb
            if source_name is not None and source != source_name: continue 
            for model in self.fitcls['bol_main'][source]:
                if model_name is not None and model != model_name: continue                
                t1,t0,t2 = self.fitcls['bol_main'][source][model].get_par(parname='texp')
                if returnv: return t1,t0,t2
                self.texp = [t1,t0,t2]
    
    def set_texp_midway(self, returnv=False):
        ''' Set the explosion epoch with the first detection and the last non-detection of ``sbobject.lc``.
        
        Parameters
        ----------         
        returnv  :        `bool`
               if return the peak infos, or store them in the class
        
        Note
        -------------
        This function is relying on the ``snrt`` parameter, which decide the detections and upper limits.
        If ``snrt`` is too low, e.g. 1 sigma, there would be some precursor or fake datapoints in the very early part of LC included.
        
        See Also
        -------------
        snobject.set_texp_bol_main, snobject.set_texp_pl
        '''
        if self.t0 <= 2400000: return
        if not 'lc' in self.__dict__: return
        t0 = self.t0     
        __lc = self.lc.query('mag<99 and jdobs<@t0')
        if len(__lc) > 0:
            expt0 = min(__lc['jdobs'])            
            __lc = self.lc.query('jdobs<@expt0')
            if len(__lc)>0:
                expt1 = max(__lc['jdobs'])            
                expt = (expt1+expt0)/2.

                # convert texp from jd to phase relative to peak
                if returnv: return expt1-t0, expt-t0, expt0-t0
                self.texp = [expt1-t0, expt-t0, expt0-t0]
            else:
                if returnv: return expt0-t0, expt0-t0, expt0-t0
                self.texp = [expt0-t0, expt0-t0, expt0-t0]
        
    def _flux_at(self, filt, phase, interpolation=None, index=0, **kwargs):
        ''' Estimate flux at a given epoch.
        
        Parameters
        ----------         
        filt :        `str`
               filter
        phase :       `float`
               rest frame phase (days) relative to t0
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        
        Returns
        ---------- 
        flux   :   `float`        
        flux error  :   `float`

        See Also
        ---------- 
        sbobject._flux_at_list
        '''
        kwargs = self.read_kwargs(**kwargs)
        
        assert interpolation in [None, 'bin', 'gp', 'fit']
        if self.t0 > 2400000: pass
        else:
            print ('set t0 correctly')
            return None, None
        if interpolation is None or interpolation=='bin':
            # get flux without interpolation
            lc = self.lc.query('filter==@filt')
            if len(lc) > 0:
                # rest frame
                __ = np.argmin(abs( (lc['jdobs'].to_numpy()-self.t0)/(1+self.z) - phase) )                
                if abs(lc['jdobs'].to_numpy()[__]-self.t0-phase)<kwargs['tdbin']:
                    return lc['flux'].to_numpy()[__], lc['eflux'].to_numpy()[__]
        if (interpolation is None and 'gpcls' in self.__dict__) or interpolation == 'gp':
            if not 'gpcls' in self.__dict__: return None, None
            p_fit = [self.t0 + phase * (1+self.z)]
            xx,yy,yye,ff = self.gpcls.predict(x_pred=p_fit, clobber=kwargs['clobber'], returnv=True)
            if len(yy[ np.where(ff == filt) ]) >0:
                return yy[ np.where(ff==filt) ][0], yye[ np.where(ff==filt) ][0]
        if (interpolation is None and 'fitcls' in self.__dict__) or interpolation == 'fit':
            if not 'fitcls' in self.__dict__: return None, None
            if 'multiband_main' in self.fitcls:
                p_fit = [phase]
                for n, model in enumerate(self.fitcls['multiband_main']):
                    if n != index: continue
                    xx,yy,yy1,yy2,ff = self.fitcls['multiband_main'][model].predict(
                        x_pred=p_fit, returnv=True,  quant=kwargs['quantile']
                    )
                    if len(yy[ np.where(ff == filt) ])>0:
                        f = yy[ np.where(ff == filt) ][0]
                        f1 = yy1[ np.where(ff == filt) ][0]
                        f2 = yy2[ np.where(ff == filt) ][0]
                        return f, max(abs(f-f1), abs(f-f2))
        return None, None
    
    def _flux_at_list(self, filt, phaselist, interpolation=None, index=0, **kwargs):
        ''' Estimate fluxes for a list of epochs.
        
        Parameters
        ----------         
        filt :        `str`
               filter
        phaselist :       `list`
               rest frame phases (days) relative to t0
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        
        Returns
        ---------- 
        flux list   :   `list`        
        flux error list  :   `list`
        
        See Also
        ---------- 
        sbobject._flux_at
        '''
        kwargs = self.read_kwargs(**kwargs)
        fluxlist, efluxlist = [], []
        for phase in phaselist:
            flux, eflux = self._flux_at(filt, phase, interpolation=interpolation, index=index, **kwargs)
            fluxlist.append(flux)
            efluxlist.append(eflux)
        return fluxlist, efluxlist
    
    def _mag_at(self, filt, phase, interpolation=None, index=0, **kwargs):
        ''' Estimate apparent magnitude at a given epoch.
        
        Parameters
        ----------         
        filt :        `str`
               filter
        phase :       `float`
               rest frame phase (days) relative to t0
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   
        
        Returns
        ---------- 
        mag   :   `float`        
        mag error  :   `float`

        See Also
        ---------- 
        sbobject._flux_at, sbobject._mag_at_list, sbobject.flux_to_mag
        '''
        kwargs = self.read_kwargs(**kwargs)
        corr_mkw = kwargs['corr_mkw']
        corr_host= kwargs['corr_host']
        flux, eflux = self._flux_at(filt, phase, interpolation=interpolation, index=index, **kwargs)
        if flux is None: return None, None
        ebv = 0
        if corr_mkw: ebv += self.mkwebv
        if corr_host: ebv += self.hostebv
        m, me, limmag = self.flux_to_mag(
            [flux], dflux=[eflux], sigma=kwargs['snrt'], zp=kwargs['zp'],
        )
        if m < 99:
            return m[0] - ebv * Rf[filt], me[0]
        else:
            m = self.flux_to_mag([flux], dflux=None, sigma=kwargs['snrt'], zp=kwargs['zp'])
            return m[0] - ebv * Rf[filt], None

    def _mag_at_list(self, filt, phaselist, interpolation=None, index=0, **kwargs):
        ''' Estimate apparent magnitudes from a list of epochs.
        
        Parameters
        ----------         
        filt :        `str`
               filter
        phaselist :       `list`
               rest frame phases (days) relative to t0
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   

        Returns
        ---------- 
        mag list  :   `list`        
        mag error list :   `list`

        See Also
        ---------- 
        sbobject._flux_at_list, sbobject._mag_at
        '''
        kwargs = self.read_kwargs(**kwargs)        
        maglist, emaglist = [], []
        for phase in phaselist:
            mag, emag = self._mag_at(filt, phase, interpolation=interpolation, index=index, **kwargs)
            maglist.append(mag)
            emaglist.append(emag)
        return maglist, emaglist
        
    def _absmag_at(self, filt, phase, interpolation=None, index=0, **kwargs):
        ''' Estimate absolute magnitude at a given epoch.
        
        Parameters
        ----------         
        filt :        `str`
               filter
        phase :       `float`
               rest frame phase (days) relative to t0
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   

        Returns
        ---------- 
        mag   :   `float`        
        mag error  :   `float`

        See Also
        ---------- 
        sbobject._flux_at, sbobject._mag_at
        '''
        kwargs = self.read_kwargs(**kwargs)        
        mag, emag = self._mag_at(filt, phase, interpolation=interpolation, index=index, **kwargs)
        if mag is None: return None, None
        if emag is None: return mag-self.dm, None
        else: return mag-self.dm, np.sqrt(emag**2 + self.dm_error(filt)**2)

    def _absmag_at_list(self, filt, phaselist, interpolation=None, index=0, **kwargs):
        ''' Estimate absolute magnitudes from a list of phases.
        
        Parameters
        ----------         
        filt :        `str`
               filter
        phaselist :       `list`
               rest frame phases (days) relative to t0
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   

        Returns
        ---------- 
        mag list  :   `list`        
        mag error list :   `list`
        
        See Also
        ---------- 
        sbobject._flux_at, sbobject._absmag_at, sbobject._mag_at
        '''
        kwargs = self.read_kwargs(**kwargs)         
        maglist, emaglist = [], []
        for phase in phaselist:
            mag, emag = self._absmag_at(filt, phase, interpolation=interpolation, index=index, **kwargs)
            maglist.append(mag)
            emaglist.append(emag)
        return maglist, emaglist
    
    def _color_at(self, filt1, filt2, phase, interpolation=None, index=0, **kwargs):
        ''' Estimate colour at a given phase.
        
        Parameters
        ----------         
        filt1 :        `str`
               filter1
        filt2 :        `str`
               filter2
        phase :       `float`
               rest frame phase (days) relative to t0
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   

        Returns
        ---------- 
        color   :   `float`        
        color error  :   `float`
        
        See Also
        ---------- 
        sbobject._flux_at, sbobject._mag_at
        '''
        kwargs = self.read_kwargs(**kwargs)          
        mag1, emag1 = self._mag_at(filt1, phase, interpolation=interpolation, index=index, **kwargs)
        mag2, emag2 = self._mag_at(filt2, phase, interpolation=interpolation, index=index, **kwargs)
        if mag1 is None or mag2 is None: return None, None
        if emag1 is None or emag2 is None: return mag1-mag2, None
        else: return mag1-mag2, np.sqrt(emag1**2+emag2**2)
    
    def _rate_at(self, filt, phase1, phase2, interpolation=None, index=0, **kwargs):
        ''' Estimate the magnitude different of two phases on a given filter.
        
        Parameters
        ----------         
        filt :        `str`
               filter
        phase1 :       `float`
               rest frame phase 1 (days) relative to t0
        phase2 :       `float`
               rest frame phase 2 (days) relative to t0
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits        
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   

        Returns
        ---------- 
        rate   :   `float`        
        rate error  :   `float`
        
        See Also
        ---------- 
        sbobject._flux_at, sbobject._mag_at
        '''
        kwargs = self.read_kwargs(**kwargs)
        mag1, emag1 = self._mag_at(filt, phase1, interpolation=interpolation, index=index, **kwargs)
        mag2, emag2 = self._mag_at(filt, phase2, interpolation=interpolation, index=index, **kwargs)
        if mag1 is None or mag2 is None: return None, None
        if emag1 is None or emag2 is None: return mag1-mag2, None
        else: return mag1-mag2, np.sqrt(emag1**2+emag2**2)
        
    def est_hostebv_with_colours(self, index=0, returnv=False, **kwargs):
        '''Estimate host ebv with color comparison approach.
        
        Parameters
        ----------  
        hostebv_bands   :        `list`
               two bands of the colour                 
        hostebv_interp :     `str`
               which interpolation to calculate colour for host ebv estimation
        tdbin :     `float`
               threshold for binning
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        verbose  :    `bool`
              show process
        clobber   :   `bool`
              Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit  
        returnv  :        `bool`
               if return the ebv, or store them in the class

        Notes
        ----------
        - For IIb/Ib/Ic, we compared colours at 10d post peak to templates from Stritzinger et al 2018, table 2, 
        at https://github.com/saberyoung/HAFFET/blob/master/sdapy/data/c10_template.txt.

        - For SLSNe, we compared colours at the peak phases to Chen at al 2022: 
        **At the peak phases, the median observed color is close to zero**.       
        
        See Also
        -----------
        snobject._color_at, snobject.read_c10
        '''
        kwargs = self.read_kwargs(**kwargs)
        filt1, filt2 = kwargs['hostebv_bands']        
        interpolation= kwargs['hostebv_interp']
        
        c10_temp = self.read_c10()
        sntype = type_clean(self.sntype)
        
        if 'SLSN' in sntype:
            if (filt1 == 'g', filt2 == 'r') or (filt1 == 'r', filt2 == 'g'):
                phase, ct, cte = 0, 0, 0
            else:
                if kwargs['verbose']:
                    print ('Error: %s-%s not defined for SLSNe, check!!' % (filt1, filt2))
                return
        else:
            if not (filt1, filt2) in c10_temp:
                if kwargs['verbose']:
                    print ('Error: %s-%s not defined in c10 template file, check!!' % (filt1, filt2))
                return        
            if not sntype in c10_temp[(filt1, filt2)]:
                if kwargs['verbose']:
                    print ('Error: %s not defined for %s-%s in c10 template file, check!!'%
                           (self.sntype, filt1, filt2))
                return 
            if '%s-%s'%(filt1,filt2) not in toebv:
                if kwargs['verbose']:
                    print ('Error: E(%s-%s) to E(B-V) factor is not defined, check!!'%
                           (filt1, filt2))
                return
            ct, cte = c10_temp[(filt1, filt2)][sntype].split(',')
            ct, cte = float(ct), float(cte)
            phase = 10
            
        # correct mikly way ebv only
        _par = kwargs.copy()
        _par['corr_mkw'] = True
        _par['corr_host'] = False
        _ct, _cte = self._color_at(
            filt1, filt2, phase, interpolation=interpolation, index=index, **_par
        )
        if _ct is None:
            print ('cannot get colour at 10 days')            
            return
        if _ct < ct:  ehost = 0
        else: ehost = _ct - ct
        
        # from E(X-Y) to E(B-V)
        hostebv = ehost*toebv['%s-%s'%(filt1,filt2)]
        if returnv: return hostebv
        self.hostebv = hostebv
        
    def calc_colors(self, xpred=None, index=0, returnv=False, **kwargs):
        ''' Calculate colors for two bands.
        
        Parameters
        ----------         
        xpred     :        `list`
               colour phases, if None will use observing epochs of the color_bands
        color_bands   :        `list`
               two bands of the colour        
        color_interp :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits        
        returnv  :        `bool`
               if return the colors, or store them in the class
        bin_range       :  `list`
               for bin method, predict range
        gp_range        :  `list`
               for GP method, predict range
        fit_range       :  `list`
               for fit method, predict range
        tdbin :     `float`
               threshold for binning        
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   
        
        See Also
        ----------    
        snobject._mag_at_list
        '''
        kwargs = self.read_kwargs(**kwargs)        
        assert self.t0 > 2400000, 'set t0 first'
        if not 'lc' in self.__dict__:
            if kwargs['verbose']: print ('parse LC data for %s first' % self.objid)
            return
        filt1, filt2 = kwargs['color_bands']
        if xpred is None:
            df = self.lc.query('mag<99 and filter in [@filt1, @filt2]')
            xpred = (np.array(df['jdobs'].tolist()) - self.t0) / (1+self.z)
        else:           
            xpred = np.array(xpred)
        
        colors = dict()        
        for interpolation in kwargs['color_interp']:
            if interpolation not in colors: colors[interpolation] = dict()
            _xpred = xpred
            if kwargs['%s_range'%interpolation] is not None:
                __ = np.logical_and(
                    xpred>=min(kwargs['%s_range'%interpolation]),
                    xpred<=max(kwargs['%s_range'%interpolation]),
                )
                _xpred = xpred[__]                
            m1, m1e = self._mag_at_list(filt1, _xpred, interpolation=interpolation, index=index, **kwargs)
            m2, m2e = self._mag_at_list(filt2, _xpred, interpolation=interpolation, index=index, **kwargs)
            jds = _xpred * (1+self.z) + self.t0
            pl,m1l,m1el,m2l,m2el=[],[],[],[],[]
            for _p,_m1,_m2,_m1e,_m2e in zip(jds,m1,m2,m1e,m2e):
                if _m1 is None or _m1e is None or _m2 is None or _m2e is None: continue
                pl.append(_p)
                m1l.append(_m1)
                m1el.append(_m1e)
                m2l.append(_m2)
                m2el.append(_m2e)
            colors[interpolation]= (np.array(pl),np.array(m1l),np.array(m1el),np.array(m2l),np.array(m2el))
        if returnv: return colors
        self.colors = colors
        self.colorsbands = (filt1, filt2)
        
    def lyman_bol(self, xpred=None, index=0, returnv=False, **kwargs):
        ''' Calculate bolometric LC from colours, with Lyman bolometric correction.
        
        Parameters
        ----------
        make_bol  :       `list`
               how to make bolometric lcs: [lyman] with BC defined by Lyman approach, [bb] blackbody fits, or [spec] absolute calibrated spectra
        lyman_bands   :     `list`
               two bands of the colour        
        lyman_interp :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits       
        lyman_bc_phase  :     `str`
               which Lyman BC to be used: [normal] phase, [cool] phase
        verbose  :    `bool`
              show process
        returnv  :        `bool`
               if return the luminosities, or store them in the class
        xpred     :        `list`
               colour phases, if None will use observing epochs of the color_bands                             
        bin_range       :  `list`
               for bin method, predict range
        gp_range        :  `list`
               for GP method, predict range
        fit_range       :  `list`
               for fit method, predict range
        tdbin :     `float`
               threshold for binning        
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   
        
        Notes
        ---------
        - This function will use colours from ``snobject.calc_colors`` if it had been done, otherwise, will call ``snobject.calc_colors`` first.
        Therefore, parts of ``kwargs`` above were for ``snobject.calc_colors``.

        - ``sntype`` should be defined properly, so that can be used to find a BC function,
        at https://github.com/saberyoung/HAFFET/blob/master/sdapy/data/bc_table.txt.
        Here, we include Ib/Ic/IIb/Ic-BL for SESNe. If you assume for instance SLSN had similar BC as SESNe, you could add it in the ``BC_Lyman``,
        at https://github.com/saberyoung/HAFFET/blob/master/sdapy/functions.py.
        
        See Also
        ----------       
        snobject.calc_colors
        '''
        kwargs = self.read_kwargs(**kwargs)
        if 'lyman' not in kwargs['make_bol']: return 
        if self.sntype is None:
            if kwargs['verbose']: print ('define sn type first')
            return        
        # calculate colours
        kwargs['color_bands'] = kwargs['lyman_bands']
        kwargs['color_interp'] = kwargs['lyman_interp']
        colors = self.calc_colors(xpred=xpred, index=index, returnv=True, **kwargs)
        if colors is None: return
        
        # to bol mag
        filt1, filt2 = kwargs['lyman_bands']
        sntype = type_clean(self.sntype)
        mbols = dict()
        for _ in colors:            
            jd,m1,m1e,m2,m2e = colors[_]
            bc, __ = BC_Lyman(m1-m2, colortype='%s-%s'%(filt1,filt2), phase=kwargs['lyman_bc_phase'], sntype=sntype)
            if bc is None:
                print ('Error: check lyman BC for sntype %s, color %s-%s for %s phase'%(sntype, filt1, filt2, kwargs['lyman_bc_phase']))
                return
            mbol, embol = bc + m1 -self.dm, np.sqrt(m1e**2 + self.dm_error(filt1)**2)            
            Lbol = Mbol_to_Lbol(mbol)                        
            eLbol = abs(Mbol_to_Lbol(mbol+embol)-Lbol)            
            _ck=np.isfinite(Lbol)            
            mbols[_] = (jd[_ck], Lbol[_ck], eLbol[_ck])
        if returnv: return mbols
        self.mbol = mbols
        
    def bb_colors(self, xpred=None, index=0, returnv=False, **kwargs):
        ''' Match or interpolate multiple band LCs, for blackbody (BB) construction.
        
        Parameters
        ----------         
        xpred     :        `list`
               colour phases, if None will use observing epochs of the sed_bands
        sed_bands   :        `list`
               multiple bands of the BB
        sed_ref_bands  :   `list`
               reference filters been used to provide epochs as for bolometric LCs,
               e.g. give a list ['g','r'], will create bolometric LCs with g and r epochs, or None on all epochs that had ever been observed.
        sed_color_interp :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits               
        tdbin :     `float`
               threshold for binning
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   
        returnv  :        `bool`
               if return the SED, or store them in the class
        
        See Also
        ----------       
        snobject.calc_colors, snobject._mag_at_list
        '''
        kwargs = self.read_kwargs(**kwargs)       
        assert self.t0 > 2400000, 'set t0 first'
        
        filters = kwargs['sed_bands']
        if filters is None: filters = np.unique(self.lc['filter'])
        assert len(filters) >= 2, 'at least 2 bands needed for blackbody'
        
        if xpred is not None:
            xpred = np.array(xpred)
        elif kwargs['sed_ref_bands'] is not None:
            filters1 = kwargs['sed_ref_bands']
            df = self.lc.query('mag<99 and filter in @filters1')
            xpred = (np.array(df['jdobs'].tolist()) - self.t0) / (1+self.z)
        else:
            df = self.lc.query('mag<99 and filter in @filters')
            xpred = (np.array(df['jdobs'].tolist()) - self.t0) / (1+self.z)        
        cbb = dict()
        for interpolation in kwargs['sed_color_interp']:
            if interpolation not in cbb: cbb[interpolation] = dict()
            for f in filters:
                cbb[interpolation][f] = self._mag_at_list(f, xpred, interpolation=interpolation, index=index, **kwargs)
            cbb[interpolation]['jd'] = xpred * (1+self.z) + self.t0
        if returnv: return cbb
        self.cbb = cbb
        
    def bb_bol(self, index=0, returnv=False, fastsedfitting=True, **kwargs):
        '''  Calculate bolometric LCs with diluted blackbody fits (that were fitted with the ``snobject.run_fit`` on **sed** engine).
        
        Parameters
        ---------- 
        fastsedfitting :       `bool`
               if no SED model fittings found, fit them in minimize rountine
        index     :        `int`
               if multiple models available, which one to use
        make_bol  :       `list`
               how to make bolometric lcs: [lyman] with BC defined by Lyman approach, [bb] blackbody fits, or [spec] absolute calibrated spectra
        sed_bands   :        `list`
               multiple bands of the BB
        sed_color_interp   :  `list`
               how to make blackbody. 1: epochs with all filters within *color_thre* hours. 2: epochs with one filter, and the others from analytic fits.  3. epochs with one filter, and the other from Gaussian process
        sed_xrangep  :    `str`
               after fitting, range to reproduce the SED
        quantile :    `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        verbose  :    `bool`
               show process
        returnv  :        `bool`
               if return the luminosities, or store them in the class        
        If activate SED fitting with fastsedfitting, then more parameters are needed to run SED fittings.
        
        Notes
        ---------
        Make sure SED model fittings had been done, otherwise set ``fastsedfitting`` to fit them with *minimize* rountine.
        
        See Also
        ----------
        snobject.run_fit, snobject.bb_colors
        '''
        kwargs = self.read_kwargs(**kwargs)        
        
        # if no sed fitting performed, use minimize to run on all of them
        if 'fitcls' in self.__dict__ and 'sed' in self.fitcls:            
            pass
        elif fastsedfitting:
            _par = kwargs.copy()
            _par['sed_routine'] = 'minimize'
            _par['sed_type'] = 'bb,spec'
            if not 'bb' in _par['make_bol']: _par['make_bol'].append('bb')
            if not 'spec' in _par['make_bol']: _par['make_bol'].append('spec')
            if not 'bb' in _par['fit_methods']: _par['fit_methods'].append('bb')
            self.run_fit('sed', **_par)
            if 'fitcls' in self.__dict__ and 'sed' in self.fitcls: pass
            else: return None, None
        else:
            if kwargs['verbose']: print ('Warning: no SED fittings found')
            return None, None        
        
        if 'mbolbb' in self.__dict__ and not kwargs['clobber']:
            mbolbb = self.mbolbb                       
        else:
            mbolbb = dict()
        if 'mbolspec' in self.__dict__ and not kwargs['clobber']:
            mbolspec = self.mbolspec
        else:
            mbolspec = [[],[],[],[],[],[],[]]
        
        for source in self.fitcls['sed']:            
            _source = source.split('_')[0]            
            if _source == 'bb':
                interpolation = source.split('_')[1]
                if interpolation not in mbolbb:
                    #jdl, L, Le, T, Te, R, Re    
                    mbolbb[interpolation] = [[],[],[],[],[],[],[]]
            jd = float(source.split('_')[-1])
            if len(self.fitcls['sed'][source].keys()) == 0:continue
            model = list(self.fitcls['sed'][source].keys())[index]            
            _model = self.fitcls['sed'][source][model]
            if kwargs['sed_xrangep'] is None:
                x, y, y1, y2, w = _model.predict(
                    x_pred=None, step=1, returnv=True, quant=kwargs['quantile']
                )
            else:
                x, y, y1, y2, w = _model.predict(
                    x_pred=np.arange(min(kwargs['sed_xrangep']), max(kwargs['sed_xrangep']), 1),
                    returnv=True, quant=kwargs['quantile']
                )                                            
                
            # Get pseudobolometric luminosity by trapezoidal integration,
            #    with flux set to zero outside of observed bands                
            _L = integrate.trapz(y[np.argsort(x)],x[np.argsort(x)])
            _L1 = integrate.trapz(y1[np.argsort(x)],x[np.argsort(x)])
            _L2 = integrate.trapz(y2[np.argsort(x)],x[np.argsort(x)])
            
            # temperature and radius
            _T = _model.get_par(filt=None, quant=kwargs['quantile'], parname='t')
            _R = _model.get_par(filt=None, quant=kwargs['quantile'], parname='r')
            
            # Calculate luminosity using alternative method of Stefan-Boltzmann, and T and R from fit
            #_Lbb = 4*np.pi*(_R*1e15)**2*5.67e-5*(_T*1e3)**4
            #_Lbbe = _Lbb*np.sqrt((2*_Re/_R)**2+(4*_Te/_T)**2)
            
            if 'bb' in kwargs['make_bol'] and _source == 'bb':
                mbolbb[interpolation][0].append(jd)
                mbolbb[interpolation][1].append(_L)                
                mbolbb[interpolation][2].append((abs(_L1-_L)+abs(_L2-_L))/2.)
                mbolbb[interpolation][3].append(_T[0])
                mbolbb[interpolation][4].append((abs(_T[1]-_T[0])+abs(_T[2]-_T[0]))/2.)
                mbolbb[interpolation][5].append(_R[0])
                mbolbb[interpolation][6].append((abs(_R[1]-_R[0])+abs(_R[2]-_R[0]))/2.)
            if 'spec' in kwargs['make_bol'] and _source == 'spec':
                mbolspec[0].append(jd)
                mbolspec[1].append(_L)
                mbolspec[2].append((abs(_L1-_L)+abs(_L2-_L))/2.)
                mbolspec[3].append(_T[0])
                mbolspec[4].append((abs(_T[1]-_T[0])+abs(_T[2]-_T[0]))/2.)
                mbolspec[5].append(_R[0])
                mbolspec[6].append((abs(_R[1]-_R[0])+abs(_R[2]-_R[0]))/2.)
        if returnv: return mbolbb, mbolspec
        self.mbolbb = mbolbb
        self.mbolspec = mbolspec
        
    def bin_lc(self, returnv=False, **kwargs):
        ''' Bin multiband LCs.
        
        Parameters
        ----------             
        bin_lc:     `bool`
               bin the lcs or not
        tdbin :     `float`
               threshold to be used to bin lcs (unit in days)
        returnv  :        `bool`
               if return the LCs, or store them in the class
        
        See Also
        ----------
        snobject.bin_df, snobject.match_colors
        '''
        kwargs = self.read_kwargs(**kwargs)
        if not kwargs['bin_lc']: return
        
        # add one column, jdbin
        self.lc = self.bin_df(self.lc, deltah = kwargs['tdbin'] )

        # merge those columns with same jdbin
        _df = self.match_colors(style=1, **kwargs)
        
        if returnv: return _df
        self.lc = _df
        
    def match_colors(self, style=1, **kwargs):
        ''' Match colors.
        
        Parameters
        ----------                
        tdbin :     `float`
               threshold for binning        
        style :     `int`
               decide different output format
        
        See Also
        ----------  
        snobject.combine_multi_obs, snobject.bin_df
        '''        
        kwargs = self.read_kwargs(**kwargs)            
        _lc = self.combine_multi_obs(**kwargs)
        if style == 1: return _lc
        
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
        return _df
    
    def combine_multi_obs(self, **kwargs):
        ''' Bin and combine observations from common epoch.
        
        Parameters
        ----------                 
        tdbin :     `float`
               threshold for binning
        
        See Also
        ----------  
        snobject.match_colors, self.bin_df, self.merge_df_cols
        '''
        kwargs = self.read_kwargs(**kwargs)        
        if not 'jdbin' in self.lc:
            _lc = self.bin_df(self.lc, deltah = kwargs['tdbin'])
        else:
            _lc = self.lc        
        #if True: _lc = _lc.query('mag<99')
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
        return _df         

    def cut_lc(self, **kwargs):
        ''' Cut LC if it's too long.

        Parameters
        ----------       
        verbose :        `bool`
              show detailed running informations            
        mjdstart  :    `float`
              start julian date to query
        dstart  :    `float`
              if **mjdstart** is None, how many days prior than **t0** to query
        mjdend  :    `float`
              end julian date to query
        dend    :    `float`
              if **mjdend** is None, how many days later than **t0** to query          
        '''
        kwargs = self.read_kwargs(**kwargs)
        
        if kwargs['mjdstart'] is not None:
            mjdstart = float(kwargs['mjdstart'])
        else:
            assert self.t0 > 2400000, 'Error: specify mjdstart or t0 (in JD) first'
            mjdstart = float(self.t0) - 2400000.5 + float(kwargs['dstart'])
            if kwargs['verbose']: print ('mjdstart as %f days prior to the peak which is %.2f'%
                                         (float(kwargs['dstart']), self.t0-2400000.5))
        if kwargs['mjdend'] is not None:
            mjdend = float(kwargs['mjdend'])
        else:
            if kwargs['dend'] is None:
                mjdend = self.mjd_now(jd=False)
                if kwargs['verbose']: print ('mjdend as current mjd %.2f'%mjdend)
            else:
                assert self.t0 > 2400000, 'Error: specify mjdend or t0 (in JD) first'
                mjdend = float(self.t0) - 2400000.5 + float(kwargs['dend'])
                if kwargs['verbose']: print ('mjdend as %f days after the peak which is %.2f'%
                                             (float(kwargs['dend']), self.t0-2400000.5))
        jdmin, jdmax = mjdstart + 2400000.5, mjdend + 2400000.5
        self.lc = self.lc.query('jdobs>=@jdmin and jdobs<=@jdmax')
        
    def clip_lc(self, **kwargs):
        ''' Removes outlier data points using GP interpolation.

        Parameters
        ----------                 
        clipsigma :     `float`
               sigma for LC clipping
        
        See Also
        ----------
        snobject.run_gp
        '''
        kwargs = self.read_kwargs(**kwargs)
        if 'gpcls' not in self.__dict__: return
        
        # remove negative flux epoch
        #self.lc = self.lc.query('flux>0')        
        
        # remove outliers        
        __arr = None        
        if kwargs['clipsigma'] is not None:
            clipsigma = float(kwargs['clipsigma'])
            for f in np.unique(self.lc['filter']):
                _lc = self.lc.query('filter==@f')
                p_fit, f_fit, fe_fit = _lc['jdobs'], _lc['flux'], _lc['eflux']
                if len(p_fit) >= 3:
                    xx,yy,yye,ff = self.gpcls.predict(x_pred=p_fit, returnv=True,) 
                    __ = np.where(ff==f)                
                    xx = xx[__]
                    if len(xx) > 0: 
                        yy = yy[__]
                        yye = yye[__]
                        sigma = abs(yy-f_fit)/np.sqrt(fe_fit**2)#+yye**2)
                        where = sigma <= clipsigma                        
                        _lc.where(where, inplace = True)                      
                if __arr is None: __arr = _lc
                else: __arr = __arr.append( _lc, ignore_index=True)                
            self.lc = __arr
        
    def _nepochs(self, **kwargs):
        ''' How many epochs of photometry in each band.
        
        Parameters
        ----------                 
        plot_bands :     `list`
               photometric filters
        
        Return
        ----------                 
        ndet :     `int`
               detection number of the LC
        '''
        kwargs = self.read_kwargs(**kwargs)        
        plot_bands = kwargs['plot_bands']
        if plot_bands is None: plot_bands = np.unique(self.lc['filter'])        
        ndets = dict()
        for filt in plot_bands:
            if not filt in ndets: ndets[filt] = dict()
            dets = self.lc.query('mag<99 and filter==@filt')
            lims = self.lc.query('mag==99 and filter==@filt')
            ndets[filt]['detections'] = len(dets)
            ndets[filt]['limits'] = len(lims)
        return ndets
    
    def _ncolors(self, **kwargs):
        ''' How many color epochs.
        
        Parameters
        ----------                 
        color_bands :     `list`
               two photometric filters for colour
        tdbin :     `float`
               threshold for binning  
        
        Return
        ----------                 
        ncolors :     `int`
               number of colour epochs of the LC
        '''
        kwargs = self.read_kwargs(**kwargs)        
        color_bands = kwargs['color_bands']
        if color_bands is None:
            color_bands = []
            for f in np.unique(self.lc['filter']): color_bands.append(f)            
            if len(color_bands) <2: return
            color_bands = color_bands[:2]
        assert len(color_bands) == 2
        filt1, filt2 = color_bands
        dets = self.lc.query('mag<99')
        ncolors = 0
        if len(dets.query('filter==@filt1'))==0 or len(dets.query('filter==@filt2'))==0:
            return ncolors        
        for _jd in dets.query('filter==@filt1')['jdobs']:
            if min(abs(dets.query('filter==@filt2')['jdobs'] - _jd)) <= kwargs['tdbin']:
                ncolors += 1
        return ncolors
    
    def _peak_accuracy(self, within=3, **kwargs):
        ''' How accurate a peak can be determined, i.e. how many photometric points available within a range to the peak.
        
        Parameters
        ----------                 
        within :     `int`
               the range in days
        plot_bands :     `list`
               photometric filters
        
        Return
        ----------                 
        peakphot :     `dictionary`
               within a range, how many points in each filters
        '''                
        assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
        t0 = self.t0
        kwargs = self.read_kwargs(**kwargs)
        peakphot = dict()
        plot_bands = kwargs['plot_bands']
        if plot_bands is None: plot_bands = np.unique(self.lc['filter'])        
        for filt in plot_bands:
            filt_lc = self.lc.query('mag<99 and filter==@filt and jdobs>@t0-@within and jdobs<@t0+@within')
            peakphot[filt] = len(filt_lc)
        return peakphot      

    def _earlypoints(self, **kwargs):
        ''' How many early datapoints for power law fittings (ask texp and fpeak available).
        
        Parameters
        ----------                 
        multiband_early_bands :     `list`
               photometric filters
        multiband_early_yrange :    `list`
               phase range
        verbose   :   `bool`
               Enable progress report
        
        Return
        ----------                 
        ndet :     `int`
               detection number of the LC
        '''
        assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
        assert self.texp is not None, '!!!estimate texp first'
        assert len(self.fpeak) > 0, 'Error: set peak flux with GP or main LC fits before run early LC fits'        
        t0 = self.t0
        texp = self.texp[1]
        fpeak = self.fpeak
        
        kwargs = self.read_kwargs(**kwargs)
        fmin,fmax = min(kwargs['multiband_early_yrange']), max(kwargs['multiband_early_yrange'])
        
        earlypoints = dict()
        plot_bands = kwargs['multiband_early_bands'].copy()
        if plot_bands is None: plot_bands = np.unique(self.lc['filter'])
        
        for filt in plot_bands:
            if filt not in fpeak:
                if kwargs['verbose']: print ('peak of filter %s not available, skip'%filt)
                continue
            tm = self.fpeak[filt][0]        
            filt_lc = self.lc.query('mag<99 and filter==@filt and jdobs>@t0+@texp and jdobs<@t0 and flux<@tm*@fmax and flux>@tm*@fmin')
            earlypoints[filt] = len(filt_lc)
        return earlypoints
    
    def summarize_plot(self, ax=None, savefig=False, **kwargs):
        ''' Show a summarize plot, including flux/mag/bolometric LCs, colour curve and the spectra.
        
        Parameters
        ---------- 
        ax    :        `matplotlib.subplot`
               matplotlib plot
        savefig    :        `bool`
               save figure or not
        dpi   :        `int`
               matplotlib resolution
        figsize   :    `list`
               figure size
        figpath    :   `str`
               figure path
        
        See Also
        ---------- 
        snobject._ax, snobject._ax1, snobject._ax2, snobject._ax3, snobject._ax4, snobject.savefig, snobject.showfig
        '''        
        kwargs = self.read_kwargs(**kwargs)        
        self.fig = plt.figure(num=1, clear=True, dpi=kwargs['dpi'], tight_layout=False, constrained_layout=False,)        
        gs1 = self.fig.add_gridspec(nrows=3, ncols=2, left=.05, right=.95, wspace=0.05, hspace=0.3)
        self.ax = self.fig.add_subplot(gs1[0, 0]) # flux
        self.ax1 = self.fig.add_subplot(gs1[:2, 1]) # spectra
        self.ax2 = self.fig.add_subplot(gs1[1, 0]) # mag        
        self.ax3 = self.fig.add_subplot(gs1[2, 1]) # color
        self.ax4 = self.fig.add_subplot(gs1[2, 0]) # luminosity
        
        # plot them
        self._ax(**kwargs)                   
        self._ax1(**kwargs)
        self._ax2(**kwargs)
        self._ax3(**kwargs)
        self._ax4(**kwargs)
        if savefig:
            self.savefig(**kwargs)
            self.showfig(ax=ax, **kwargs)
        
    def _ax(self, show_title=True, show_legend=True, ylabel_2right=False,
            show_data=True, show_fits=True, show_gp=True, show_fit_error=True,
            show_texp=True, interpolation=None, index=0, color=None,
            marker=None, markersize=None, label=None, ls=None, fillstyle=None,
            yshift=None, fontsize=None, snr_thre=None, **kwargs):                
        ''' Flux LC plot.
        
        Parameters
        ----------            
        jd_x0      :   `float`
              x axis zeropoint for mag/flux LCs
        ax_xstyle     :   `str`
              x axis for flux (_ax) LCs: [rp] rest frame since peak [jd] Junlian date since jd_x0
        ax_ystyle  :   `str`
              y axis for flux (_ax) LCs: [original] flux, or [norm] normalized flux
        ax_xlim    :  `list`
              x limit
        ax_ylim    :  `list`
              y limit
        plot_sources :   `list`
              which source LCs to show
        plot_bands :   `list`
              which filters to show
        flux_scale   :   `float`
              normalize flux peak
        show_title :   `bool`
              if show title
        show_legend :    `bool`
              if show legend
        ylabel_2right :    `bool`
              if put y label to right      
        show_data :    `bool`
              if show data points
        show_fit_error :    `bool`
              if False, only show best fit, otherwise, show errors or random samplings
        show_fits :    `bool`
              if show model fittings
        show_gp  :   `bool`
              if show GP modellings
        show_texp :   `bool`
              if show explosion epochs        
        alphabest :  `float` between 0 and 1
              matplotlib alpha for best fit fitting
        alphasample  :   `float` between 0 and 1
              matplotlib alpha for random samplings or errors
        plot_mcmct  :    `float` between 0 and 1
              a threshold that select good mc samples for plotting, rangiing from 0 to 1, e.g. 0.5 means selecting samplings from the top 50 percent of all samplings relying on the likelihoods
        plot_nsamples  :   `int`
              how many random MC samples to be plotted
        multiband_early_xrangep  :   `list`
              range to reproduce the multiband_early models
        gp_xrangep  :   `list`
              range to reproduce the GP interpolations
        verbose   :   `bool`
               Enable progress report
        
        **if ax_ystyle=norm, snobject._flux_at would be used to guess the peak fluxes**:        
        tdbin :     `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis        
               
        Notes
        ---------
        Only working when ``snobject.ax`` is defined.
        
        See Also
        ----------        
        snobject.summarize_plot, snobject._flux_at
        '''
        kwargs = self.read_kwargs(**kwargs)   
        if self.ax is None:
            if kwargs['verbose']: print ('Error: no ax defined, skipped flux plots')
            return
        if not 'lc' in self.__dict__:
            if kwargs['verbose']: print ('Error %s: no lc defined, parse lc data first'%self.objid)
            return
        fscale=float(kwargs['flux_scale'])
        lc = self.lc        
        if kwargs['plot_sources'] is not None:
            source = kwargs['plot_sources']
            lc = lc.query('source==@source')
        if kwargs['plot_bands'] is not None:
            plot_bands = kwargs['plot_bands']
        else:
            plot_bands = np.unique(lc['filter'])
        if kwargs['plot_bands'] is not None:
            lc = lc.query('filter==@plot_bands')        
        labelfit, labelgp = True, True
        for nfilt, filt in enumerate(plot_bands):
            ''' for each filter '''
            if filt not in PROP1:
                if kwargs['verbose']: print ('Warning %s: skip %s since it was not registed in the filters'%(self.objid,filt))
                continue
            __lc = lc.query('filter==@filt')
            xx, yy, yye = __lc['jdobs'], __lc['flux'], __lc['eflux']
            if len(xx) == 0: continue
            
            # normalize flux item in ax                        
            if filt in self.fpeak: fm, efm = self.fpeak[filt][0], self.fpeak[filt][1]
            else: fm, efm = self._flux_at(filt, 0, interpolation=interpolation, index=index, **kwargs)
            if kwargs['ax_ystyle'] == 'norm':
                # peak flux                          
                if fm is None:
                    if kwargs['verbose']: print ('Warning %s: no peak found for %s, skipped' % (self.objid, filt))
                    continue
                if yshift is not None:
                    _yshift = yshift
                elif not filt in ys:
                    if kwargs['verbose']: print ('Warning %s: no shift defined for %s, skipped' % (self.objid, filt))
                    continue
                else:
                    _yshift = ys[filt]
                yy = yy/fm*fscale+_yshift
                yye = yye/fm*fscale
            else:
                assert kwargs['ax_ystyle'] == 'original', 'ax_ystyle should be either original or norm'
            if kwargs['ax_xstyle'] == 'rp':
                assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
                xx = (xx-self.t0)/(1+self.z)
            else:
                assert kwargs['ax_xstyle'] == 'jd', 'ax_xstyle should be either rp or jd'
                xx = xx - kwargs['jd_x0']
            
            # show data points
            if show_data:
                PROP = PROP1[filt].copy()
                if color is not None: PROP['color'] = color
                if marker is not None: PROP['marker'] = marker
                if markersize is not None: PROP['markersize'] = markersize
                if label is not None: PROP['label'] = label
                if ls is not None: PROP['ls'] = ls
                if fillstyle is not None: PROP['fillstyle'] = fillstyle
                if snr_thre is not None: # set threshold on snr
                    snr = abs(yy/yye)
                    _ = np.where(snr>=snr_thre)
                    xx,yy,yye = get_numpy(xx)[_],get_numpy(yy)[_],get_numpy(yye)[_]                
                self.ax.errorbar(
                    xx, yy, yerr=yye, alpha=kwargs['alphabest'], **PROP
                )
                
            # analytic sn lc fit  
            if 'fitcls' in self.__dict__ and show_fits:
                if 'multiband_early' in self.fitcls:
                    for model in self.fitcls['multiband_early']:
                        _model = self.fitcls['multiband_early'][model]
                        
                        # samplings                        
                        if show_fit_error:                            
                            x,y,f = _model.predict_random(
                                limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                                x_pred=None, step=1
                            )
                            for xx, yy in zip(x[ np.where(f == filt) ], y[ np.where(f == filt) ]):
                                if labelfit: labels, labelfit = 'fits', False
                                else: labels = None
                                if kwargs['ax_xstyle'] == 'rp': xx = (xx-self.t0)/(1+self.z)
                                else:  xx = xx - kwargs['jd_x0']
                                # Note: early models are fitted on normalized LCs
                                if kwargs['ax_ystyle'] == 'norm':                                    
                                    if yshift is not None: _yshift=yshift
                                    elif not filt in ys: continue
                                    else: _yshift = ys[filt]                                    
                                    yy = yy+_yshift                                    
                                elif fm is None:  continue
                                else:  yy = yy/fscale*fm
                                PROP = PROP2[filt].copy()
                                if color is not None: PROP['color'] = color                
                                if ls is not None: PROP['ls'] = ls                
                                self.ax.plot(xx, yy, alpha=kwargs['alphasample'], label=labels, **PROP)
                        else:
                            if labelfit: labels, labelfit = 'fits', False
                            else: labels = None
                            x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                            if kwargs['ax_xstyle'] == 'rp':  xx = (x[np.where(f==filt)]-self.t0)/(1+self.z)
                            else: xx = x[np.where(f==filt)] - kwargs['jd_x0']
                            if kwargs['ax_ystyle'] == 'norm':
                                if fm is None:  continue
                                if yshift is not None: _yshift=yshift
                                elif not filt in ys: continue
                                else: _yshift = ys[filt]                                    
                                yy = yy/fm*fscale+_yshift
                            PROP = PROP2[filt].copy()
                            if color is not None: PROP['color'] = color                
                            if ls is not None: PROP['ls'] = ls     
                            self.ax.plot(xx, yy, alpha=kwargs['alphabest'], label=labels, **PROP)
                            
                if 'multiband_main' in self.fitcls:
                    for model in self.fitcls['multiband_main']:
                        _model = self.fitcls['multiband_main'][model]
                        
                        # samples                        
                        if show_fit_error: 
                            x,y,f=_model.predict_random(
                                limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                                x_pred=None, step=1
                            )
                            for xx, yy in zip(x[ np.where(f == filt) ], y[ np.where(f == filt) ]):
                                if labelfit: labels, labelfit = 'fits', False
                                else: labels = None
                                if kwargs['ax_xstyle'] == 'rp':  xx = (xx-self.t0)/(1+self.z)
                                else:  xx = xx - kwargs['jd_x0']
                                if kwargs['ax_ystyle'] == 'norm':
                                    if fm is None:  continue
                                    if yshift is not None: _yshift=yshift
                                    elif not filt in ys: continue
                                    else: _yshift = ys[filt]                                    
                                    yy = yy/fm*fscale+_yshift
                                PROP = PROP2[filt].copy()
                                if color is not None: PROP['color'] = color                
                                if ls is not None: PROP['ls'] = ls     
                                self.ax.plot(xx, yy, alpha=kwargs['alphasample'], label=labels,**PROP)
                        else:
                            if labelfit: labels, labelfit = 'fits', False
                            else: labels = None                            
                            x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                            if kwargs['ax_xstyle'] == 'rp':  xx = (x[np.where(f==filt)]-self.t0)/(1+self.z)
                            else:  xx = x[np.where(f==filt)] - kwargs['jd_x0']
                            if kwargs['ax_ystyle'] == 'norm':
                                if fm is None:  continue
                                if yshift is not None: _yshift=yshift
                                elif not filt in ys: continue
                                else: _yshift = ys[filt]                            
                                yy = y[np.where(f==filt)]/fm*fscale+_yshift
                            else: yy = y[np.where(f==filt)]
                            PROP = PROP2[filt].copy()
                            if color is not None: PROP['color'] = color                
                            if ls is not None: PROP['ls'] = ls     
                            self.ax.plot(xx, yy, alpha=kwargs['alphabest'], label=labels, **PROP)

            # bol lc fit
            '''
            if 'fitcls' in self.__dict__ and show_bol:
                if 'sed' in self.fitcls:
                    index = 0
                    fl = dict()
                    for source in self.fitcls['sed']:
                        _source = source.split('_')[0]
                        if _source == 'bb': interp = source.split('_')[1]
                        else: continue
                        jd = float(source.split('_')[-1])
                        _phase = (jd-self.t0)/(1+self.z)                
                        model = list(self.fitcls['sed'][source].keys())[index]
                        _model = self.fitcls['sed'][source][model]
                        # get BB SED
                        x, y, w = _model.predict_random(
                            limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                            x_pred=None, step = 1,
                        )

                        # get fitted bolometric mag
                        #for source in self.fitcls['bol_main']:
                        #    for model in self.fitcls['bol_main'][source]:
                        #        _model = self.fitcls['bol_main'][source][model]
                        #        x,y,f=_model.predict(
                        #            limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                        #            x_pred=None, step=1,                            
                        #        )
                                                       
                        for nn, (xx, yy) in enumerate(zip(x, y)): 
                            __ = np.where(xx == central_wavelengths[filt])
                            ff = yy[__]
                            if ff <= 0: continue
                            m, em = self._from_luminosity(ff, filt, el=None, **kwargs)
                            if m < 0 or m > 25: continue
                            f = self.mag_to_flux(m, sigma=kwargs['snrt'], zp=kwargs['zp'])                            
                            if kwargs['ax_xstyle'] == 'rp':  _x = (jd-self.t0)/(1+self.z)
                            else:  _x = jd - kwargs['jd_x0']
                            if kwargs['ax_ystyle'] == 'norm':
                                if fm is None:  continue
                                if yshift is not None: pass
                                elif not filt in ys: continue
                                else: yshift = ys[filt]
                                f = f/fm*fscale+yshift
                            if not nn in fl: fl[nn] = dict()
                            if not 'x' in fl[nn]: fl[nn]['x'] = []
                            if not 'y' in fl[nn]: fl[nn]['y'] = []                            
                            fl[nn]['x'].append( _x )
                            fl[nn]['y'].append( f[0] )                            
                    for n in fl:
                        if nfilt == 0 and n == 0: labels = 'bol fits'
                        else: labels = None
                        xx, yy = np.array(fl[n]['x']), np.array(fl[n]['y'])
                        __ = np.argsort(xx)
                        self.ax.plot(xx[__], yy[__], alpha=kwargs['alphasample'], label=labels, **PROP2[filt])                    
            '''
            # Gaussian process
            if 'gpcls' in self.__dict__ and show_gp:
                if filt in self.gpcls.f_pred:                                        
                    __ = np.where(self.gpcls.f_pred==filt)
                    xx = self.gpcls.x_pred[__]
                    yy = self.gpcls.y_pred[__]
                    yye = self.gpcls.y_prede[__]                                        
                    if labelgp: labels, labelgp = 'GP', False
                    else: labels = None
                    if kwargs['ax_xstyle'] == 'rp':  xx = (xx-self.t0)/(1+self.z)
                    else:  xx = xx - kwargs['jd_x0']
                    if kwargs['ax_ystyle'] == 'norm':
                        if fm is None:  continue
                        if yshift is not None: _yshift=yshift
                        elif not filt in ys: continue
                        else: _yshift = ys[filt]                                    
                        yy = yy/fm*fscale+_yshift
                        yye = yye/fm*fscale
                    PROP = PROP3[filt].copy()
                    if color is not None: PROP['color'] = color                
                    if ls is not None: PROP['ls'] = ls     
                    self.ax.plot(xx, yy, alpha=kwargs['alphabest'], label=labels, **PROP)
                    if show_fit_error: 
                        self.ax.fill_between(
                            xx, yy-yye, yy+yye, alpha=kwargs['alphasample'], **PROP
                        )
            
        #  subplot settings
        if kwargs['ax_xlim'] is not None:
            x1,x2 = min(kwargs['ax_xlim']), max(kwargs['ax_xlim'])
            if kwargs['ax_xstyle'] == 'rp':
                self.ax.set_xlim([x1, x2])
            else:
                self.ax.set_xlim([x1*(1+self.z)+self.t0-kwargs['jd_x0'], x2*(1+self.z)+self.t0-kwargs['jd_x0']])                
        if kwargs['ax_ylim'] is not None: 
            self.ax.set_ylim(kwargs['ax_ylim'])
        if fontsize is None: fontsize = 12
        if show_title:
            self.ax.set_title('%s\n%s z=%.2f'% (self.objid,self.sntype,self.z),fontsize=fontsize)
        
        if self.texp is not None and show_texp:
            if kwargs['ax_xstyle'] == 'rp':
                self.ax.axvline(self.texp[1], ls='--', label='texp', color='k')
                self.ax.axvline(self.texp[0], ls=':', color='k')
                self.ax.axvline(self.texp[2], ls=':', color='k')                                
            else:
                self.ax.axvline(self.texp[1]*(1+self.z)+self.t0-kwargs['jd_x0'], ls='--', label='texp', color='k')
                self.ax.axvline(self.texp[0]*(1+self.z)+self.t0-kwargs['jd_x0'], ls=':', color='k')
                self.ax.axvline(self.texp[2]*(1+self.z)+self.t0-kwargs['jd_x0'], ls=':', color='k')
                
        if kwargs['ax_xstyle'] == 'rp':
            self.ax.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=fontsize)
        else:
            self.ax.set_xlabel('JD - %s (d)'%kwargs['jd_x0'],fontsize=fontsize)
            
        if kwargs['zp'] == 23.9:   yunit = 'Flux ($\mu$Jy)'
        elif kwargs['zp'] == 48.6: yunit = 'Flux (ergs/s/cm2/Hz)'
        else:  yunit = 'Flux (ZP=%s)' % kwargs['zp']                                  
        if ylabel_2right:
            ax=self.ax.twinx()
            if kwargs['ax_ystyle'] == 'norm': ax.set_ylabel('$F_{\mathrm{norm}}$ + $\mathrm{offset}$',fontsize=fontsize)
            else:                             ax.set_ylabel(yunit,fontsize=fontsize)
            ax.set_yticks([])
            self.ax.tick_params(axis="y",direction="in",labelleft="on",pad=-18)
        else:
            if kwargs['ax_ystyle'] == 'norm': self.ax.set_ylabel('$F_{\mathrm{norm}}$ + $\mathrm{offset}$',fontsize=fontsize)
            else:                             self.ax.set_ylabel(yunit,fontsize=fontsize)
        if show_legend: self.ax.legend(fontsize=fontsize, frameon=False)
        
    def _ax1(self, show_title=False, show_telluric=True, show_text=True,
             show_data=True, show_peaks=True, show_fits=True, text_type='date',
             show_fit_error=True, **kwargs):
        """ Spectra plot.
        
        Parameters
        ----------      
        plot_specsources:  `list`
             which sources of spectra to show
        show_stype :     `str`
             spectral data type, options: 'original', 'rest', 'bin', 'continuum', 'flat'
        show_element :    `str`
             1. [full] show full range spectra; 2. [sntype] show characteristic features depends on SN type; 3. else, should be specific element, e.g. 'H~$\alpha$', check and define elements in https://github.com/saberyoung/HAFFET/blob/master/sdapy/constants.py
        v_bounds  :    `list`
             guessed velocity range (unit: 1e3 km/s)
        continuum_method  :  `str`
             The function type for continumm fitting, valid functions are "scalar", "linear", "quadratic", "cubic", "poly", and "exponential"
        continuum_degree  :  `int`
             degree of polynomial when method="poly", for continuum fitting
        pfactor           :  `int`
             threshold used for peak detection
        plot_mcmct  :    `float` between 0 and 1
              a threshold that select good mc samples for plotting, rangiing from 0 to 1, e.g. 0.5 means selecting samplings from the top 50 percent of all samplings relying on the likelihoods
        plot_nsamples  :   `int`
              how many random MC samples to be plotted
        show_telluric :  `bool`
             show telluric lines or not
        show_title :  `bool`
             show subplot title or not       
        show_data :    `bool`
              if show data points
        show_peaks :    `bool`
              if show peaks and valleys of spectra
        show_fit_error :    `bool`
              if False, only show best fit, otherwise, show errors or random samplings
        show_fits :    `bool`
              if show model fittings
        show_text :    `bool`
             show spectra infos to the right of each spectrum        
        alphabest :  `float` between 0 and 1
              matplotlib alpha for best fit fitting
        alphasample  :   `float` between 0 and 1
              matplotlib alpha for random samplings or errors
        text_type :    `str`
             time type in the text: date or phase or jd
        verbose   :   `bool`
               Enable progress report

        Notes
        ---------
        Only working when ``snobject.ax1`` is defined.
        
        See Also
        --------
        snobject.summarize_plot
        """
        kwargs = self.read_kwargs(**kwargs)
        if self.ax1 is None:
            if kwargs['verbose']: print ('Error: no ax1 defined, skipped spectra plots')
            return
        if not 'spec' in self.__dict__:
            if kwargs['verbose']: print ('Error %s: no spec defined, parse spectral data first'%self.objid)
            return
        # which spectra type
        assert kwargs['show_stype'] in ['original', 'rest', 'bin', 'continuum', 'flat']  
        # which elements to show
        element=None
        if kwargs['show_element']=='full':
            pass
        elif kwargs['show_element']=='sntype' and self.sntype in constants.line_forsn:
            element = constants.line_forsn[self.sntype]        
        else:
            element = kwargs['show_element']
            assert element in constants.line_location, 'defined %s in constants.line_location first' % element              
        cw, region = handle_spectrum.parse_element(element, kwargs['v_bounds'])                
        self.spec.sort_spectra()
        for _ in self.spec.data: 
            spec   = self.spec.data[_]['data']
            phase  = float(self.spec.data[_]['phase'])            
            if kwargs['plot_specsources'] is not None and _ not in kwargs['plot_specsources']:
                continue
            _source = '_'.join(_.split())
            ys      = self.spec.data[_]['ys']
            
            # get data, plot them
            wave, flux = spec._norm_spectrum(region=region, stype=kwargs['show_stype'])
            if len(wave) == 0: continue

            if show_data:
                self.ax1.step(wave, flux + ys, **PROP7['data'])
            
            if show_text:
                if text_type == 'date': text = '%s' % _source
                elif text_type == 'phase': text = '%.2f d' % phase
                elif text_type == 'jd': text = '%.2f d' % phase * (1+self.z) + self.t0
                else: text = ''
                if region is not None:
                    self.ax1.text(max(region)+20, ys*1.05, text, color='k')
                else:
                    self.ax1.text(max(wave)+20, ys*1.05, text, color='k')                
                    
            if show_peaks:
                # get peaks
                findp = spec.find_peaks(element=None, region=region, stype=kwargs['show_stype'],
                            pfactor=kwargs['pfactor'], smooth=True, continuum_method=kwargs['continuum_method'],
                            continuum_degree=kwargs['continuum_degree'], line_velocity_bounds=kwargs['v_bounds'])            
                if len(findp) > 0: # if peak exists, show continuum
                    wave, fsmooth = spec.get_local_continuum(region, kwargs['show_stype'],
                                kwargs['continuum_method'], kwargs['continuum_degree'])
                    self.ax1.step(wave, fsmooth + ys, **PROP7['bindata'])
                if 1 in findp:
                    for _w in findp[1]:
                        __ = np.argmin(abs(wave-_w))
                        # en.wikipedia.org/wiki/Template:Unicode_chart_Arrows
                        self.ax1.plot(_w, flux[__] + ys, color='orange', marker=u'$\u2191$', markersize=20)
                if -1 in findp:
                    for _w in findp[-1]:
                        __ = np.argmin(abs(wave-_w))
                        self.ax1.plot(_w, flux[__] + ys, color='cyan', marker=u'$\u2193$', markersize=20)
                        
            # fits 
            if 'fitcls' in self.__dict__ and kwargs['show_stype'] == 'flat' and show_fits:
                if 'specline' in self.fitcls:
                    for sourcename in self.fitcls['specline']:
                        if not _source in sourcename: continue
                        for model in self.fitcls['specline'][sourcename]:                            
                            _model = self.fitcls['specline'][sourcename][model]
                            if _model is None: continue

                            if show_fit_error:
                                x, y, f = _model.predict_random(
                                    limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples']
                                )
                                #if region is not None:
                                #    # match current minima                            
                                #    xmin, ymin = x[0][np.argmin(y[0])], min(y[0])                            
                                #    ys0 = flux[np.argmin(abs(wave-xmin))] - ymin
                                #else:
                                ys0 = 0
                                for xx, yy in zip(x,y):
                                    # there're overlaps?
                                    if region is None: pass
                                    elif min(region) > max(xx) or min(xx) > max(region): continue
                                    self.ax1.plot(xx, yy+ys+ys0, alpha=kwargs['alphasample'], **PROP7['fit'])
                            else:
                                x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                                for xx, yy in zip(x,y):
                                    # there're overlaps?
                                    if region is None: pass
                                    elif min(region) > max(xx) or min(xx) > max(region): continue
                                    self.ax1.plot(xx, yy+ys, alpha=kwargs['alphabest'], **PROP7['fit'])                                    
        if cw is not None: self.ax1.axvline(cw, color='k', ls='--')
        
        # telluric lines?
        if show_telluric:
            for which in telluric_lines:
                tl = telluric_lines[which]
                if region is None: pass
                elif min(region) > max(tl) or min(tl) > max(region): continue            
                self.ax1.axvspan(min(tl), max(tl), alpha=0.2, color='red')                
        #  subplot settings        
        if region is not None: _ticks = np.arange(min(region),max(region)+100,100)
        else: _ticks = np.arange(3000, 9000, 1000)        
        self.ax1.set_xticks( _ticks )
        self.ax1.set_xticklabels(['%.1f'%_ for _ in _ticks], rotation=30)
        if cw is not None and region is not None: # velocity
            ax1 = self.ax1.twiny()
            ax1.set_xlim( region )
            _tickl = ['%.1f'%_ for _ in (_ticks/cw-1)*299792.458/1000.]
            ax1.set_xticks(_ticks)
            ax1.set_xticklabels(_tickl, rotation=30)
            ax1.set_xlabel('velocity ($10^{3}~km~s^{-1}$)', fontsize=12)
        self.ax1.set_xlabel('Rest Wavelength ($\AA$)', fontsize=12)        
        self.ax1.set_yticks([])
        self.ax1.tick_params("both", direction='in', right=True, labelright=True, left=False, labelleft=False)        
        if region is not None:            
            self.ax1.set_xlim([min(region),max(region)])            
        if show_title:
            self.ax1.set_title('%s\n%s z=%.2f'% (self.objid,self.sntype,self.z),fontsize=12) 
        
    def _ax2(self, show_title=False, show_legend=False, ylabel_2right=False,
             show_data=True, show_limits=True, show_fits=True, show_gp=True, show_fit_error=True,
             show_texp=False, color=None, marker=None, markersize=None, label=None,
             ls=None, fillstyle=None, fontsize=None, **kwargs):            
        ''' Magnitude LC plot.
        
        Parameters
        ----------            
        jd_x0      :   `float`
              x axis zeropoint for mag/flux LCs
        ax2_xstyle     :   `str`
              x axis for mag (_ax2) LCs: [rp] rest frame since peak [jd] Junlian date since jd_x0
        ax2_ystyle  :   `str`
              y axis for mag (_ax2) LCs: [app] apparent mag, or [abs] absolute mag
        ax2_xlim    :  `list`
              x limit
        ax2_ylim    :  `list`
              y limit
        plot_sources :   `list`
              which source LCs to show
        plot_bands :   `list`
              which filters to show        
        show_title :   `bool`
              if show title
        show_legend :    `bool`
              if show legend
        ylabel_2right :    `bool`
              if put y label to right 
        corr_mkw   :   `str`
              when calculating absulte mag, if correct milky way extinction if there're any
        corr_host  :   `str`
              when calculating absulte mag, if correct host galaxy extinction if there're any
        show_data :    `bool`
              if show data points
        show_limits :    `bool`
              if show upper limits
        show_fit_error :    `bool`
              if False, only show best fit, otherwise, show errors or random samplings
        show_fits :    `bool`
              if show model fittings
        show_gp  :   `bool`
              if show GP modellings
        show_texp :   `bool`
              if show explosion epochs        
        alphabest :  `float` between 0 and 1
              matplotlib alpha for best fit fitting
        alphasample  :   `float` between 0 and 1
              matplotlib alpha for random samplings or errors
        plot_mcmct  :    `float` between 0 and 1
              a threshold that select good mc samples for plotting, rangiing from 0 to 1, e.g. 0.5 means selecting samplings from the top 50 percent of all samplings relying on the likelihoods
        plot_nsamples  :   `int`
              how many random MC samples to be plotted
        multiband_early_xrangep  :   `list`
              range to reproduce the multiband_early models
        gp_xrangep  :   `list`
              range to reproduce the GP interpolations
        verbose   :   `bool`
               Enable progress report
        
        Notes
        ---------
        Only working when ``snobject.ax2`` is defined.
        
        See Also
        ----------        
        snobject.ax
        '''
        kwargs = self.read_kwargs(**kwargs) 
        if self.ax2 is None:
            if kwargs['verbose']: print ('Error: no ax2 defined, skipped mag plots')
            return
        if not 'lc' in self.__dict__:
            if kwargs['verbose']: print ('Error %s: no lc defined, parse lc data first'%self.objid)
            return
        ebv = 0
        if kwargs['corr_mkw']: ebv += self.mkwebv
        if kwargs['corr_host']: ebv += self.hostebv
        lc = self.lc
        if kwargs['plot_sources'] is not None:
            source = kwargs['plot_sources']
            lc = lc.query('source==@source')
        if kwargs['plot_bands'] is not None:
            plot_bands = kwargs['plot_bands']
        else:
            plot_bands = np.unique(lc['filter'])        
        if kwargs['plot_bands'] is not None:            
            lc = lc.query('filter==@plot_bands')
        labelfit, labelgp = True, True
        for filt in plot_bands:
            ''' for each filter '''
            if filt not in PROP1:
                if kwargs['verbose']: print ('Warning %s: skip %s since it was not registed in the filters'%(self.objid,filt))
                continue
            
            # data points
            __lc = lc.query('filter==@filt and mag<99')
            xx, yy, yye = __lc['jdobs'], __lc['mag'], __lc['emag']

            if kwargs['ax2_ystyle'] == 'abs':
                yy -= self.dm
                yye = np.sqrt(yye**2 + self.dm_error(filt)**2)
            else:
                assert kwargs['ax2_ystyle'] == 'app', 'ax2_ystyle should be either app or abs'
            if kwargs['ax2_xstyle'] == 'rp':
                assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'
                xx = (xx-self.t0)/(1+self.z)
            else:
                assert kwargs['ax2_xstyle'] == 'jd', 'ax2_xstyle should be either rp or jd'
                xx = xx - kwargs['jd_x0']            
            if not filt in Rf:
                print ('Warning: no Rf defined for filter %s, ignore its extinction' % filt)
            else:
                yy -= Rf[filt]*ebv
                
            # show data points
            if show_data:
                PROP = PROP1[filt].copy()
                if color is not None: PROP['color'] = color
                if marker is not None: PROP['marker'] = marker
                if markersize is not None: PROP['markersize'] = markersize
                if label is not None: PROP['label'] = label
                if ls is not None: PROP['ls'] = ls
                if fillstyle is not None: PROP['fillstyle'] = fillstyle                
                self.ax2.errorbar(
                    xx, yy, yerr=yye, alpha=kwargs['alphabest'], **PROP
                )
            
            # upper limits
            # for very early or late epochs
            #if self.texp is not None: tt = self.texp[1]
            #elif self.t0 > 2400000:  tt = self.t0 - 10
            #else: tt = np.inf
            #__lc = lc.query('filter==@filt and mag==99 and jdobs<@tt')
            __lc = lc.query('filter==@filt and mag==99')
            if len(__lc) > 0 and show_limits and 'limmag' in __lc.keys():
                xx, yy = __lc['jdobs'], __lc['limmag']
                if kwargs['ax2_ystyle'] == 'abs': yy -= self.dm            
                if kwargs['ax2_xstyle'] == 'rp':  xx = (xx-self.t0)/(1+self.z)
                else: xx = xx - kwargs['jd_x0']
                if filt in Rf: yy -= Rf[filt]*ebv
                PROP = PROP1l[filt].copy()
                if color is not None: PROP['color'] = color
                if marker is not None: PROP['marker'] = marker
                if markersize is not None: PROP['markersize'] = markersize                
                if ls is not None: PROP['ls'] = ls
                if fillstyle is not None: PROP['fillstyle'] = fillstyle
                self.ax2.plot(xx, yy, alpha=kwargs['alphasample'], **PROP)
                
            # analytic sn lc fit            
            if 'fitcls' in self.__dict__ and show_fits:
                #if 'multiband_early' in self.fitcls:
                if False:
                    for model in self.fitcls['multiband_early']:
                        _model = self.fitcls['multiband_early'][model]
                        
                        # samples
                        if show_fit_error:  
                            x,y,f=_model.predict_random(
                                limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                                x_pred=None, step=1
                            )
                            for xx, yy in zip(x[ np.where(f == filt) ], y[ np.where(f == filt) ]):
                                if labelfit: labels, labelfit = 'fits', False
                                else: labels = None                                                                                
                                m = self.flux_to_mag(yy, dflux=None, sigma=kwargs['snrt'], zp=kwargs['zp'])                                
                                if kwargs['ax2_ystyle'] == 'abs': m -= self.dm            
                                if kwargs['ax2_xstyle'] == 'rp':  xx = (xx-self.t0)/(1+self.z)
                                else: xx = xx - kwargs['jd_x0']            
                                if filt in Rf: m -= Rf[filt]*ebv
                                PROP = PROP2[filt].copy()
                                if color is not None: PROP['color'] = color                
                                if ls is not None: PROP['ls'] = ls     
                                self.ax2.plot(xx, m, alpha=kwargs['alphasample'], label=labels, **PROP)
                        else:
                            if labelfit: labels, labelfit = 'fits', False
                            else: labels = None                            
                            x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                            m = self.flux_to_mag(
                                y[np.where(f==filt)], dflux=None, sigma=kwargs['snrt'], zp=kwargs['zp']
                            )
                            if kwargs['ax2_ystyle'] == 'abs': m -= self.dm            
                            if kwargs['ax2_xstyle'] == 'rp':  xx = (x[np.where(f==filt)]-self.t0)/(1+self.z)
                            else: xx = x[np.where(f==filt)] - kwargs['jd_x0']            
                            if filt in Rf: m -= Rf[filt]*ebv
                            PROP = PROP2[filt].copy()
                            if color is not None: PROP['color'] = color                
                            if ls is not None: PROP['ls'] = ls     
                            self.ax2.plot(xx, m, alpha=kwargs['alphabest'], label=labels, **PROP)
                        
                if 'multiband_main' in self.fitcls:
                    for model in self.fitcls['multiband_main']:
                        _model = self.fitcls['multiband_main'][model]                        
                            
                        # samples
                        if show_fit_error: 
                            x,y,f=_model.predict_random(
                                limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                                x_pred=None, step=1
                            )
                            for xx, yy in zip(x[ np.where(f == filt) ], y[ np.where(f == filt) ]):
                                if labelfit: labels, labelfit = 'fits', False
                                else: labels = None
                                m = self.flux_to_mag(yy, dflux=None, sigma=kwargs['snrt'], zp=kwargs['zp'])
                                if kwargs['ax2_ystyle'] == 'abs': m -= self.dm            
                                if kwargs['ax2_xstyle'] == 'rp':  xx = (xx-self.t0)/(1+self.z)
                                else: xx = xx - kwargs['jd_x0']            
                                if filt in Rf: m -= Rf[filt]*ebv
                                PROP = PROP2[filt].copy()
                                if color is not None: PROP['color'] = color                
                                if ls is not None: PROP['ls'] = ls     
                                self.ax2.plot(xx, m, alpha=kwargs['alphasample'], label=labels, **PROP)
                        else:
                            if labelfit: labels, labelfit = 'fits', False
                            else: labels = None                            
                            x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                            m = self.flux_to_mag(
                                y[np.where(f==filt)], dflux=None, sigma=kwargs['snrt'], zp=kwargs['zp']
                            )
                            if kwargs['ax2_ystyle'] == 'abs': m -= self.dm            
                            if kwargs['ax2_xstyle'] == 'rp':  xx = (x[np.where(f==filt)]-self.t0)/(1+self.z)
                            else: xx = x[np.where(f==filt)] - kwargs['jd_x0']            
                            if filt in Rf: m -= Rf[filt]*ebv
                            PROP = PROP2[filt].copy()
                            if color is not None: PROP['color'] = color                
                            if ls is not None: PROP['ls'] = ls     
                            self.ax2.plot(xx, m, alpha=kwargs['alphabest'], label=labels, **PROP)
                            
            # Gaussian process
            if 'gpcls' in self.__dict__ and show_gp:
                if filt in self.gpcls.f_pred:                                        
                    __ = np.where(self.gpcls.f_pred==filt)
                    xx = self.gpcls.x_pred[__]
                    yy = self.gpcls.y_pred[__]
                    yye = self.gpcls.y_prede[__]
                    if labelgp: labels, labelgp = 'GP', False
                    else: labels = None                                                            
                    m,me,limmag = self.flux_to_mag(yy, dflux=yye, sigma=kwargs['snrt'], zp=kwargs['zp'])
                    __=np.where(m<99)
                    xx, m, me = xx[__], m[__], me[__]
                    if kwargs['ax2_ystyle'] == 'abs': m -= self.dm            
                    if kwargs['ax2_xstyle'] == 'rp':  xx = (xx-self.t0)/(1+self.z)
                    else: xx = xx - kwargs['jd_x0']            
                    if filt in Rf: m -= Rf[filt]*ebv
                    PROP = PROP3[filt].copy()
                    if color is not None: PROP['color'] = color                
                    if ls is not None: PROP['ls'] = ls     
                    self.ax2.plot(xx, m, alpha=kwargs['alphabest'], label=labels, **PROP)
                    if show_fit_error:
                        #self.ax2.fill_between(xx, m-me, m+me, alpha=kwargs['alphasample'], **PROP)
                        self.ax2.errorbar(xx, m, yerr=me, alpha=kwargs['alphasample'], **PROP)
        
        #  subplot settings        
        if kwargs['ax2_xlim'] is not None:
            x1,x2 = min(kwargs['ax2_xlim']), max(kwargs['ax2_xlim'])
            if kwargs['ax2_xstyle'] == 'rp':
                self.ax2.set_xlim([x1, x2])
            else:
                self.ax2.set_xlim([x1*(1+self.z)+self.t0-kwargs['jd_x0'], x2*(1+self.z)+self.t0-kwargs['jd_x0']])
        if kwargs['ax2_ylim'] is not None:
            self.ax2.set_ylim(kwargs['ax2_ylim'])
        else:
            self.ax2.invert_yaxis()
        if fontsize is None: fontsize = 12
        if show_title:
            self.ax2.set_title('%s\n%s z=%.2f'% (self.objid,self.sntype,self.z),fontsize=fontsize)
        if self.texp is not None and show_texp:
            if kwargs['ax2_xstyle'] == 'rp':
                self.ax2.axvline(self.texp[1]/(1+self.z), ls='--', label='texp', color='k')
                self.ax2.axvline(self.texp[0]/(1+self.z), ls=':', color='k')
                self.ax2.axvline(self.texp[2]/(1+self.z), ls=':', color='k')
            else:
                self.ax2.axvline(self.texp[1]*(1+self.z)+self.t0-kwargs['jd_x0'], ls='--', label='texp', color='k')
                self.ax2.axvline(self.texp[0]*(1+self.z)+self.t0-kwargs['jd_x0'], ls=':', color='k')
                self.ax2.axvline(self.texp[2]*(1+self.z)+self.t0-kwargs['jd_x0'], ls=':', color='k')
        if kwargs['ax2_xstyle'] == 'rp':
            self.ax2.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=fontsize)
        else:
            self.ax2.set_xlabel('JD - %s (d)'%kwargs['jd_x0'],fontsize=fontsize)
        if ylabel_2right:
            ax22=self.ax2.twinx()
            if kwargs['ax_ystyle'] == 'abs': ax22.set_ylabel('M$_{abs}$ (mag)',fontsize=fontsize)
            else:   ax22.set_ylabel('M$_{app}$ (mag)',fontsize=fontsize)
            ax22.set_yticks([])
            self.ax2.tick_params(axis="y",direction="in",labelleft="on",pad=-18)
        else:
            if kwargs['ax_ystyle'] == 'abs':self.ax2.set_ylabel('M$_{abs}$ (mag)',fontsize=fontsize)
            else:  self.ax2.set_ylabel('M$_{app}$ (mag)',fontsize=fontsize)
        if show_legend: self.ax2.legend(fontsize=fontsize, frameon=False)
        
    def _ax3(self, xpred=None, index=0, show_title=False, show_legend=False,
             ylabel_2right=True, show_data=True, show_fits=False,
             show_gp=False, show_fit_error=True, show_texp=True,
             color=None, marker=None, markersize=None, label=None,
             ls=None, fillstyle=None, fontsize=None, **kwargs):        
        ''' Colour curve plot.
        
        Parameters
        ----------        
        jd_x0      :   `float`
              x axis zeropoint for mag/flux LCs
        ax3_xstyle     :   `str`
              x axis for mag (_ax2) LCs: [rp] rest frame since peak [jd] Junlian date since jd_x0        
        ax3_xlim    :  `list`
              x limit
        ax3_ylim    :  `list`
              y limit            
        show_title :   `bool`
              if show title
        show_legend :    `bool`
              if show legend
        ylabel_2right :    `bool`
              if put y label to right         
        show_data :    `bool`
              if show data points       
        show_fit_error :    `bool`
              if False, only show best fit, otherwise, show errors or random samplings
        show_fits :    `bool`
              if show model fittings
        show_gp  :   `bool`
              if show GP modellings
        show_texp :   `bool`
              if show explosion epochs        
        alphabest :  `float` between 0 and 1
              matplotlib alpha for best fit fitting
        alphasample  :   `float` between 0 and 1
              matplotlib alpha for random samplings or errors
        plot_mcmct  :    `float` between 0 and 1
              a threshold that select good mc samples for plotting, rangiing from 0 to 1, e.g. 0.5 means selecting samplings from the top 50 percent of all samplings relying on the likelihoods
        plot_nsamples  :   `int`
              how many random MC samples to be plotted        
        verbose   :   `bool`
               Enable progress report
        
        If snobject.colors not available, will run snobject.calc_colors to calculate color first, with parameters below :
        
        xpred     :        `list`
               colour phases, if None will use observing epochs of the color_bands
        color_bands   :        `list`
               two bands of the colour        
        color_interp :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits         
        tdbin :     `float`
               threshold for binning        
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   
                
        Notes
        ---------
        - Only working when ``snobject.ax3`` is defined.

        - If ``snobject.colors`` not available, do ``snobject.calc_colors`` first.
        
        See Also
        ----------        
        snobject.calc_colors, snobject._mag_at_list
        '''
        if self.ax3 is None:
            print ('Error: no ax3 defined, skipped color plots')
            return
        kwargs = self.read_kwargs(**kwargs)       
        
        # calculate colours
        f1, f2 = kwargs['color_bands']
        if not 'colors' in self.__dict__:            
            colors = self.calc_colors(xpred=None, index=index, returnv=True, **kwargs)
        else:
            if 'colorsbands' in self.__dict__: f1, f2 = self.colorsbands
            colors = self.colors
            
        if colors is None: return
        if kwargs['ax3_xstyle'] == 'rp':
            assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'                    
        else:
            assert kwargs['ax3_xstyle'] == 'jd', 'ax3_xstyle should be either rp or jd'
        # show data points
        if show_data:                
            for _ in colors:
                if not _ in kwargs['color_interp']: continue            
                jd,m1,m1e,m2,m2e = colors[_]
                if kwargs['ax3_xstyle'] == 'rp': xx = (jd-self.t0)/(1+self.z)
                else: xx = jd - kwargs['jd_x0']                
                yy = m1-m2
                yye = np.sqrt(m1e**2 + m2e**2)
                PROP = PROP5[_].copy()                
                if color is not None: PROP['color'] = color
                if marker is not None: PROP['marker'] = marker
                if markersize is not None: PROP['markersize'] = markersize
                if label is not None: PROP['label'] = label
                if ls is not None: PROP['ls'] = ls
                if fillstyle is not None: PROP['fillstyle'] = fillstyle
                self.ax3.errorbar(xx, yy, yerr=yye, **PROP)
                
        # fittings        
        xpred = np.arange(min(self.lc.query('mag<99')['jdobs']), max(self.lc.query('mag<99')['jdobs']), 1)
        if kwargs['ax3_xstyle'] == 'rp': xpred = (xpred-self.t0)/(1+self.z)
        else: xpred = xpred - kwargs['jd_x0']
        if show_fits:
            xx, yy, yye = [], [], []
            for phase in xpred:
                c, ce = self._color_at(f1, f2, phase, interpolation='fit', **kwargs)
                if c is None or ce is None: continue
                xx.append(phase)
                yy.append(c)
                yye.append(ce)
            PROP = PROP5['fit_all'].copy()                
            if color is not None: PROP['color'] = color                
            if ls is not None: PROP['ls'] = ls
            if show_fit_error:
                self.ax3.errorbar(xx, yy, yerr=yye, **PROP)
            else:
                self.ax3.plot(xx, yy, **PROP)
            
        # GP
        if show_gp:
            xx, yy, yye = [], [], []
            for phase in xpred:
                c, ce = self._color_at(f1, f2, phase, interpolation='gp', **kwargs)
                if c is None or ce is None: continue
                xx.append(phase)
                yy.append(c)
                yye.append(ce)
            PROP = PROP5['gp_all'].copy()                
            if color is not None: PROP['color'] = color                
            if ls is not None: PROP['ls'] = ls
            if show_fit_error:
                self.ax3.errorbar(xx, yy, yerr=yye, **PROP)
            else:
                self.ax3.plot(xx, yy, **PROP)
        
        #  subplot settings        
        if kwargs['ax3_xlim'] is not None:
            x1,x2 = min(kwargs['ax3_xlim']), max(kwargs['ax3_xlim'])
            if kwargs['ax3_xstyle'] == 'rp':
                self.ax3.set_xlim([x1, x2])
            else:
                self.ax3.set_xlim([x1*(1+self.z)+self.t0-kwargs['jd_x0'], x2*(1+self.z)+self.t0-kwargs['jd_x0']])                
        if kwargs['ax3_ylim'] is not None:
            self.ax3.set_ylim(kwargs['ax3_ylim'])
        if fontsize is None: fontsize = 12
        if show_title:
            self.ax3.set_title('%s\n%s z=%.2f' % (self.objid,self.sntype,self.z),fontsize=fontsize)
        if self.texp is not None and show_texp:
            if kwargs['ax3_xstyle'] == 'rp':
                self.ax3.axvline(self.texp[1], ls='--', label='texp', color='k')
                self.ax3.axvline(self.texp[0], ls=':', color='k')
                self.ax3.axvline(self.texp[2], ls=':', color='k')                                
            else:
                self.ax3.axvline(self.texp[1]*(1+self.z)+self.t0-kwargs['jd_x0'], ls='--', label='texp', color='k')
                self.ax3.axvline(self.texp[0]*(1+self.z)+self.t0-kwargs['jd_x0'], ls=':', color='k')
                self.ax3.axvline(self.texp[2]*(1+self.z)+self.t0-kwargs['jd_x0'], ls=':', color='k')                
        if kwargs['ax3_xstyle'] == 'rp':
            self.ax3.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=fontsize)
        else:
            self.ax3.set_xlabel('JD - %s (d)'%kwargs['jd_x0'],fontsize=fontsize)
        if ylabel_2right:
            ax33=self.ax3.twinx()
            ax33.set_ylabel('%s-%s (mag)'%(f1,f2),fontsize=fontsize)
            ax33.set_yticks([])
            self.ax3.tick_params(axis="y",direction="in",labelleft="on",pad=-18)
        else:
            self.ax3.set_ylabel('%s-%s (mag)'%(f1,f2),fontsize=fontsize)
        if show_legend: self.ax3.legend(fontsize=fontsize, frameon=False)
        
    def _ax4(self, xpred=None, index=0, show_title=False, show_legend=False,
             ylabel_2right=False, show_data=True, show_fits=True,
             show_fit_error=True, show_texp=True, color=None, marker=None,
             markersize=None, label=None, ls=None, fillstyle=None, fontsize=None,
             **kwargs): 
        ''' luminosity LCs plot
        
        Parameters
        ----------        
        make_bol      :  `list`
              how to make bolometric lcs: [lyman] with BC defined by Lyman approach, [bb] blackbody fits on multiband spectra, or [spec] integration on blackbody of absolute calibrated spectra
        jd_x0      :   `float`
              x axis zeropoint for mag/flux LCs
        ax4_xstyle     :   `str`
              x axis for luminosity (_ax4) LCs: [rp] rest frame since peak [jd] Junlian date since jd_x0        
        ax4_ystyle     :   `str`
              y axis for luminosity (_ax4) LCs: [linear] or [log] scale
        ax4_xlim    :  `list`
              x limit
        ax4_ylim    :  `list`
              y limit            
        show_title :   `bool`
              if show title
        show_legend :    `bool`
              if show legend
        ylabel_2right :    `bool`
              if put y label to right         
        show_data :    `bool`
              if show data points        
        show_fit_error :    `bool`
              if False, only show best fit, otherwise, show errors or random samplings
        show_fits :    `bool`
              if show model fittings        
        show_texp :   `bool`
              if show explosion epochs        
        alphabest :  `float` between 0 and 1
              matplotlib alpha for best fit fitting
        alphasample  :   `float` between 0 and 1
              matplotlib alpha for random samplings or errors
        plot_mcmct  :    `float` between 0 and 1
              a threshold that select good mc samples for plotting, rangiing from 0 to 1, 
              e.g. 0.5 means selecting samplings from the top 50 percent of all samplings relying on the likelihoods
        plot_nsamples  :   `int`
              how many random MC samples to be plotted        
        verbose   :   `bool`
               Enable progress report
        
        Notes
        -----------------
        - If snobject.mbol not available and make_bol includes lyman, 
        will run snobject.lyman_bol to calculate Lyman bolometric LCs first, 
        with parameters, e.g. (check snobject.lyman_bol for more):
        
           lyman_bands   :     `list`
               two bands of the colour        
           lyman_interp :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits       
           lyman_bc_phase  :     `str`
               which Lyman BC to be used: [normal] phase, [cool] phase       
        
        If snobject.mbolbb not available and make_bol includes bb, 
        will run snobject.bb_bol to calculate BB bolometric LCs (fitted with multiband photometry) first;
        and the same for spec: if snobject.mbolspec not available and make_bol includes spec, 
        will run snobject.bb_bol to calculate BB bolometric LCs (fited with spectra) first,
        with parameters, e.g. (check snobject.bb_bol for more):
        
        sed_bands   :        `list`
               multiple bands of the BB
        sed_color_interp   :  `list`
               how to make blackbody. 
               1: epochs with all filters within *color_thre* hours. 
               2: epochs with one filter, and the others from analytic fits.
               3. epochs with one filter, and the other from Gaussian process
        sed_xrangep  :    `str`
               after fitting, range to reproduce the SED       
        
        See Also
        ----------        
        snobject.calc_colors, snobject._mag_at_list, snobject.lyman_bol, snobject.bb_bol
        '''
        if self.ax4 is None:
            print ('Error: no ax4 defined, skipped luminosity plots')
            return        
        kwargs = self.read_kwargs(**kwargs)
        if kwargs['ax4_xstyle'] == 'rp':
            assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP'                    
        else:
            assert kwargs['ax4_xstyle'] == 'jd', 'ax3_xstyle should be either rp or jd'        
        if 'lyman' in kwargs['make_bol']:
            if 'mbol' in self.__dict__:
                mbol = self.mbol
            else:
                _par = kwargs.copy()
                if not 'lyman' in _par['make_bol']: _par['make_bol'].append('lyman')
                mbol = self.lyman_bol(xpred=xpred, index=index, returnv=True, **_par)            
            if show_data and mbol is not None:
                for _ in mbol:
                    if len(mbol[_][0]) == 0: continue
                    if kwargs['ax4_xstyle'] == 'rp':
                        xx = (get_numpy(mbol[_][0])-self.t0)/(1+self.z)
                    else:
                        xx = get_numpy(mbol[_][0]) - kwargs['jd_x0']
                    yy = get_numpy(mbol[_][1])
                    yye = get_numpy(mbol[_][2])
                    PROP = PROP6['lyman_%s'%_].copy()                    
                    if color is not None: PROP['color'] = color
                    if marker is not None: PROP['marker'] = marker
                    if markersize is not None: PROP['markersize'] = markersize
                    if label is not None: PROP['label'] = label
                    if ls is not None: PROP['ls'] = ls
                    if fillstyle is not None: PROP['fillstyle'] = fillstyle
                    self.ax4.errorbar(xx,yy,yerr=yye,**PROP)
            elif kwargs['verbose']: print ('skip lyman bol at ax4')
        if 'bb' in kwargs['make_bol']:            
            if 'mbolbb' in self.__dict__:
                mbolbb = self.mbolbb
            else:
                _par = kwargs.copy()
                if not 'bb' in _par['make_bol']: _par['make_bol'].append('bb')                
                mbolbb, mbolspec = self.bb_bol(index=index, returnv=True, fastsedfitting=True, **kwargs)
                
            if show_data and mbolbb is not None:  
                for _ in mbolbb:
                    if len(mbolbb[_][0]) == 0: continue
                    if kwargs['ax4_xstyle'] == 'rp':
                        xx = (get_numpy(mbolbb[_][0])-self.t0)/(1+self.z)
                    else:
                        xx = get_numpy(mbolbb[_][0]) - kwargs['jd_x0']
                    yy = get_numpy(mbolbb[_][1])
                    yye = get_numpy(mbolbb[_][2])
                    PROP = PROP6['bb_%s'%_].copy()
                    if color is not None: PROP['color'] = color
                    if marker is not None: PROP['marker'] = marker
                    if markersize is not None: PROP['markersize'] = markersize
                    if label is not None: PROP['label'] = label
                    if ls is not None: PROP['ls'] = ls
                    if fillstyle is not None: PROP['fillstyle'] = fillstyle
                    self.ax4.errorbar(xx,yy,yerr=yye,**PROP)
            elif kwargs['verbose']: print ('skip bb bol at ax4')
        if 'spec' in kwargs['make_bol']:
            if 'mbolspec' in self.__dict__:
                mbolspec = self.mbolspec
            else:
                kwargs['make_bol'] = ['spec']
                mbolbb, mbolspec = self.bb_bol(index=index, returnv=True, fastsedfitting=True, **kwargs)                
            if show_data and mbolspec is not None:
                if kwargs['ax4_xstyle'] == 'rp':
                    xx = (get_numpy(mbolspec[0])-self.t0)/(1+self.z)
                else:
                    xx = get_numpy(mbolspec[0]) - kwargs['jd_x0']
                yy = get_numpy(mbolspec[1])
                yye = get_numpy(mbolspec[2])
                PROP = PROP6['spec'].copy()
                if color is not None: PROP['color'] = color
                if marker is not None: PROP['marker'] = marker
                if markersize is not None: PROP['markersize'] = markersize
                if label is not None: PROP['label'] = label
                if ls is not None: PROP['ls'] = ls
                if fillstyle is not None: PROP['fillstyle'] = fillstyle
                self.ax4.errorbar(xx,yy,yerr=yye,**PROP)
            elif kwargs['verbose']: print ('skip spec bol at ax4')
            
        # fitting models
        if 'fitcls' in self.__dict__ and show_fits:
            if 'bol_early' in self.fitcls:
                for source in self.fitcls['bol_early']:
                    #if not source in kwargs['make_bol']: continue
                    for model in self.fitcls['bol_early'][source]:
                        _model = self.fitcls['bol_early'][source][model]
                        
                        # samples
                        PROP = PROP6['early_fit'].copy()
                        if color is not None: PROP['color'] = color                
                        if ls is not None: PROP['ls'] = ls
                        if show_fit_error: 
                            x,y,f=_model.predict_random(
                                limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                                x_pred=None, step=1, 
                            )
                            for nn, (xx, yy) in enumerate(zip(x, y)):                                
                                if nn != 0: PROP['label']=None
                                if kwargs['ax4_xstyle'] == 'rp':
                                    xx = (xx-self.t0)/(1+self.z)
                                else:
                                    xx = xx - kwargs['jd_x0']
                                self.ax4.plot(xx, yy, alpha=kwargs['alphasample'], **PROP)
                        else:
                            x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                            if kwargs['ax4_xstyle'] == 'rp':
                                xx = (x-self.t0)/(1+self.z)
                            else:
                                xx = x - kwargs['jd_x0']
                            self.ax4.plot(xx, y, alpha=kwargs['alphabest'], **PROP)
            if 'bol_main' in self.fitcls:
                for source in self.fitcls['bol_main']:
                    #if not source in kwargs['make_bol']: continue
                    for model in self.fitcls['bol_main'][source]:
                        _model = self.fitcls['bol_main'][source][model]                        
                        
                        # samples
                        PROP = PROP6['main_fit'].copy()
                        if color is not None: PROP['color'] = color                
                        if ls is not None: PROP['ls'] = ls
                        if show_fit_error:                            
                            x,y,f=_model.predict_random(
                                limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                                x_pred=None, step=1,                            
                            )
                            for nn, (xx, yy) in enumerate(zip(x, y)):                                
                                if nn != 0: PROP['label']=None
                                if kwargs['ax4_xstyle'] == 'rp':
                                    xx = (xx-self.t0)/(1+self.z)
                                else:
                                    xx = xx - kwargs['jd_x0']
                                self.ax4.plot(xx, yy, alpha=kwargs['alphasample'], **PROP)
                        else:
                            x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                            if kwargs['ax4_xstyle'] == 'rp':
                                xx = (x-self.t0)/(1+self.z)
                            else:
                                xx = x - kwargs['jd_x0']
                            self.ax4.plot(xx, y, alpha=kwargs['alphabest'], **PROP)
            if 'bol_tail' in self.fitcls:
                for source in self.fitcls['bol_tail']:
                    #if not source in kwargs['make_bol']: continue
                    for model in self.fitcls['bol_tail'][source]:
                        _model = self.fitcls['bol_tail'][source][model]                        
                        
                        # samples
                        PROP = PROP6['tail_fit'].copy()
                        if color is not None: PROP['color'] = color                
                        if ls is not None: PROP['ls'] = ls
                        if show_fit_error:
                            x,y,f=_model.predict_random(
                                limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                                x_pred=None, step=1,
                            )
                            for nn, (xx, yy) in enumerate(zip(x, y)):                                
                                if nn != 0: PROP['label']=None
                                if kwargs['ax4_xstyle'] == 'rp':
                                    xx = (xx-self.t0)/(1+self.z)
                                else:
                                    xx = xx - kwargs['jd_x0']
                                self.ax4.plot(xx, yy, alpha=kwargs['alphasample'], **PROP)
                        else:
                            x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                            if kwargs['ax4_xstyle'] == 'rp':
                                xx = (x-self.t0)/(1+self.z)
                            else:
                                xx = x - kwargs['jd_x0']
                            self.ax4.plot(xx, y, alpha=kwargs['alphabest'], **PROP)
            if 'bol_full' in self.fitcls:
                for source in self.fitcls['bol_full']:
                    #if not source in kwargs['make_bol']: continue
                    for model in self.fitcls['bol_full'][source]:
                        _model = self.fitcls['bol_full'][source][model]
                                                
                        # samples
                        PROP = PROP6['joint_fit'].copy()
                        if color is not None: PROP['color'] = color                
                        if ls is not None: PROP['ls'] = ls
                        if show_fit_error:
                            x,y,f=_model.predict_random(
                                limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                                x_pred=None, step=1,
                            )
                            for nn, (xx, yy) in enumerate(zip(x, y)):                                
                                if nn != 0: PROP['label']=None
                                if kwargs['ax4_xstyle'] == 'rp':
                                    xx = (xx-self.t0)/(1+self.z)
                                else:
                                    xx = xx - kwargs['jd_x0']
                                self.ax4.plot(xx, yy, alpha=kwargs['alphasample'], **PROP)
                        else:
                            x, y, y1, y2, f = _model.predict(x_pred=None, step = 1, returnv=True, quant=[.5,.5,.5])
                            if kwargs['ax4_xstyle'] == 'rp':
                                xx = (x-self.t0)/(1+self.z)
                            else:
                                xx = x - kwargs['jd_x0']
                            self.ax4.plot(xx, y, alpha=kwargs['alphabest'], **PROP)
        
        #  subplot settings
        if kwargs['ax4_xlim'] is not None:
            x1,x2 = min(kwargs['ax4_xlim']), max(kwargs['ax4_xlim'])
            if kwargs['ax4_xstyle'] == 'rp':
                self.ax4.set_xlim([x1, x2])
            else:
                self.ax4.set_xlim([x1*(1+self.z)+self.t0-kwargs['jd_x0'], x2*(1+self.z)+self.t0-kwargs['jd_x0']])                        
        if kwargs['ax4_ylim'] is not None:
            self.ax4.set_ylim(kwargs['ax4_ylim'])
        if kwargs['ax4_ystyle'] == 'log': self.ax4.set_yscale('log')
        else: assert kwargs['ax4_ystyle'] == 'linear'
        if fontsize is None: fontsize = 12
        if show_title:
            self.ax4.set_title('%s\n%s z=%.2f'% (self.objid,self.sntype,self.z),fontsize=fontsize)
        if self.texp is not None and show_texp:
            if kwargs['ax4_xstyle'] == 'rp':
                self.ax4.axvline(self.texp[1]/(1+self.z), ls='--', label='texp', color='k')
                self.ax4.axvline(self.texp[0]/(1+self.z), ls=':', color='k')
                self.ax4.axvline(self.texp[2]/(1+self.z), ls=':', color='k')
            else:
                self.ax4.axvline(self.texp[1]*(1+self.z)+self.t0-kwargs['jd_x0'], ls='--', label='texp', color='k')
                self.ax4.axvline(self.texp[0]*(1+self.z)+self.t0-kwargs['jd_x0'], ls=':', color='k')
                self.ax4.axvline(self.texp[2]*(1+self.z)+self.t0-kwargs['jd_x0'], ls=':', color='k')
        if kwargs['ax4_xstyle'] == 'rp':
            self.ax4.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=fontsize)
        else:
            self.ax4.set_xlabel('JD - %s (d)'%kwargs['jd_x0'],fontsize=fontsize)
        if ylabel_2right:
            ax44=self.ax4.twinx()
            ax44.set_ylabel('L$_{bol}$ (erg $s^{-1}$)',fontsize=fontsize)
            ax44.set_yticks([])
            self.ax4.tick_params(axis="y",direction="in",labelleft="on",pad=-18)
        else:
            self.ax4.set_ylabel('L$_{bol}$ (erg $s^{-1}$)',fontsize=fontsize) 
        if show_legend: self.ax4.legend(fontsize=fontsize, frameon=False)        
        
    def _ax5(self, phase, index=0, interpolation=None, show_title=False, show_legend=False,
             show_data=True, show_fits=True, show_fit_error=True, **kwargs): 
        ''' show constructed SED of specific epochs 
        
        Parameters
        ----------         
        ax5_ystyle  :   `str`
              y axis for sed (_ax5) LCs: [linear] or [log] scale
        ax5_xlim    :  `list`
              x limit
        ax5_ylim    :  `list`
              y limit
        phase :       `float`
               rest frame phase (days) relative to t0        
        tdbin   :        `float`
               threshold for binning
        interpolation :     `str`
               estimate flux with data epoch less than than **tdbin**, or interpolation from GP/fits
        corr_mkw    :     `bool`
               if correct milky way extinction
        corr_host    :     `bool`
               if host galaxy extinction
        index      :  `int`
               if multiple models available, which of them to be used
        quantile   :  `list`
               use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors
        clobber   :   `bool`
               Redo analysis
        zp :     `float`
              zeropoint to convert magnitude to flux.
              zp = 23.9 for micro Jansky to AB mag
              zp = 48.6 for ergs/s/cm2/Hz to AB mag
              e.g. mab = -2.5 * log10(fv[Jy]/3631) = -2.5 * log10(fv[mJy]) + 2.5*log10(3631*1e6)
        snrt  :    `float`
              SNR threshold to distinguish detection/limit   
        do_kcorr    :  `bool`
              K corrected the photometry or spectra
        ab2vega     :  `bool`
              set it True if input UV/IR data in AB, if they're vega, set to False
        sed_type   :        `list`
               which data to fit the SED: [bb] photometric points from varies of interpolations, or [spec] the spectrum
        sed_color_interp   : `list`
               for SED fitting with multiband photometry, how to estimate the missing bands, 1: epochs with all filters within *color_thre* hours. 2: epochs with one filter, and the others from analytic fits.  3. epochs with one filter, and the other from Gaussian process
        sed_bands          : `list`
               filters to make blackbody, give a list. e.g.['g','r','i'], or None will use all possible bands        
        sed_abscal_bands   : `list`
               which photometry bands to be used to absolute calibrate the selected spectra
        sed_abscal_method  : `list`
               how to use photometry to calibrate the spectra: [cal]ibrate or [mangle] the spetra
        sed_mangle_func    : `list`
               mangling function
        spec_snr          :  `float`
               if no error from spectra, add noise level of snr        
        
        See Also
        ----------  
        snobject._mag_at, snobject._to_luminosity, self._cal_spectrum
        '''
        kwargs = self.read_kwargs(**kwargs) 
        if self.ax5 is None:
            if kwargs['verbose']: print ('Error: no ax5 defined, skipped sed plots')
            return                
        assert self.t0 > 2400000, '!!!either input jdpeak or do GP first and set jdpeak with GP' 
        if 'bb' in kwargs['sed_type'] and show_data:
            filters = kwargs['sed_bands']        
            if filters is None: filters = np.unique(self.lc['filter'])            
            for interpolation in kwargs['sed_color_interp']:                
                for n, f in enumerate(filters):                    
                    m, em = self._mag_at(f, phase, interpolation=interpolation, index=index, **kwargs)
                    if m is None: continue
                    _flux, _dflux, wlref, bandwidths = self._to_luminosity(m, f, em=em, **kwargs)
                    if em is None: # limits
                        self.ax5.plot(wlref, _flux,
                            color = PROP6['bb_%s'%interpolation]['color'],
                            marker='v', fillstyle='none', ls='none')
                    else: # detections
                        PROP = PROP6['bb_%s'%interpolation].copy()
                        if n>0: PROP['label'] = None
                        if show_fit_error:
                            self.ax5.errorbar(wlref, _flux, yerr = _dflux, xerr=bandwidths, **PROP)
                        else:
                            self.ax5.plot(wlref, _flux, **PROP6['bb_%s'%interpolation])
        if 'spec' in kwargs['sed_type'] and 'spec' in self.__dict__ and show_data:
            for _ in self.spec.data:
                spec   = self.spec.data[_]['data']
                _phase  = float(self.spec.data[_]['phase'])
                if abs(_phase - phase) > kwargs['tdbin']: continue
                
                # get original data
                ws, fluxes = spec.wave, spec.flux
                efluxes = spec._add_noise(fluxes, kwargs['spec_snr'])
                if len(ws) == 0: continue                
                fluxes, efluxes = self._cal_spectrum(ws, fluxes, efluxes, phase, **kwargs)
                if show_fit_error:
                    self.ax5.errorbar(ws, fluxes, yerr = efluxes, **PROP6['spec'])
                else:
                    self.ax5.plot(ws, fluxes, **PROP6['spec'])
        if show_fits and 'fitcls' in self.__dict__ and 'sed' in self.fitcls:
            for source in self.fitcls['sed']:
                _source = source.split('_')[0]
                if _source == 'bb': interpolation = source.split('_')[1]
                jd = float(source.split('_')[-1])
                _phase = (jd-self.t0)/(1+self.z)
                if abs(_phase - phase) > kwargs['tdbin']: continue
                
                model = list(self.fitcls['sed'][source].keys())[index]            
                _model = self.fitcls['sed'][source][model]
                if kwargs['sed_xrangep'] is None:
                    x, y, w = _model.predict_random(
                        limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                        x_pred=None, step = 1,
                    )
                else:
                    x, y, w = _model.predict_random(
                        limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                        x_pred=np.arange(min(kwargs['sed_xrangep']), max(kwargs['sed_xrangep']), 1), 
                    )
                if 'bb' in kwargs['sed_type'] and _source == 'bb':
                    pass
                elif 'spec' in kwargs['sed_type'] and _source == 'spec':
                    pass
                else:
                    continue
                for nn, (xx, yy) in enumerate(zip(x, y)):                    
                    if _source == 'spec': PROP = PROP6['spec_fit']
                    else: PROP = PROP6['cbb_%s_fit'%interpolation]                    
                    self.ax5.plot(xx, yy, alpha=kwargs['alphasample'], **PROP)
                    if show_fit_error: break
                    
        #  subplot settings
        if kwargs['ax5_xlim'] is not None:
            x1,x2 = min(kwargs['ax5_xlim']), max(kwargs['ax5_xlim'])            
            self.ax5.set_xlim([x1, x2])
        elif kwargs['sed_xrangep'] is not None:            
            self.ax5.set_xlim([min(kwargs['sed_xrangep']), max(kwargs['sed_xrangep'])])
        else:
            self.ax5.set_xlim([500, 20000])            
        if kwargs['ax5_ylim'] is not None:
            self.ax5.set_ylim(kwargs['ax5_ylim'])            
        if show_title:
            self.ax5.set_title('%s\n%s z=%.2f' % (self.objid,self.sntype,self.z),fontsize=12)
        if kwargs['do_kcorr']: self.ax5.set_xlabel(r'Rest Frame Wavelength ($\AA$)')
        else:  self.ax5.set_xlabel(r'Wavelength ($\AA$)')                        
        if kwargs['ax5_ystyle'] == 'log':
            self.ax5.set_yscale('log')
            self.ax5.set_ylabel('log Luminosity (erg s-1 A-1)')
        else:
            assert kwargs['ax5_ystyle'] == 'linear'
            self.ax5.set_ylabel('Luminosity (erg s-1 A-1)')
        if show_legend: self.ax5.legend(fontsize=10, frameon=False)
        
    def _ax6(self, show_title=False, show_legend=True, show_data=True,
             show_fits=True, show_fit_error=True, color=None, marker=None,
             markersize=None, label=None, ls=None, fillstyle=None, fontsize=None,
             **kwargs):
        ''' show spectral velocity evolution.
        
        Parameters
        ----------         
        sn_line   :   `str`
              which line to fit, e.g. 'H~$\alpha$', 'He~5876$\AA$', 'O~7774$\AA$'        
        ax6_xlim    :  `list`
              x limit
        ax6_ylim    :  `list`
              y limit        
        
        See Also
        ----------  
        snobject.set_vexp
        '''
        kwargs = self.read_kwargs(**kwargs) 
        if self.ax6 is None:
            if kwargs['verbose']: print ('Error: no ax5 defined, skipped sed plots')
            return
        if show_data:
            if 'fitcls' not in self.__dict__: return
            if not 'specline' in self.fitcls.keys(): return        
            assert self.t0 > 2400000, 'set t0 first'
            
            # spectral line name
            _sline=None
            if kwargs['sn_line']=='full':
                print ('!!! Warning: define a line first')
                return
            elif kwargs['sn_line']=='sntype' and self.sntype in line_forsn:
                _sline = line_forsn[self.sntype]
            elif len(kwargs['sn_line']) == 0 or kwargs['sn_line'] is None:
                if self.sntype in line_forsn:
                    _sline = line_forsn[self.sntype]
                else:
                    print ('!!! Warning: define a line first')
                    return
            else:
                _sline = kwargs['sn_line']
                assert _sline in line_location, 'define line location for %s, otherwise use %s'%\
                    (_sline, line_location.keys())
            intrinstic_wave = line_location[_sline]

            # get spectral fittings
            xx, yy, yye = [], [], []
            def func(a): return -a/intrinstic_wave*299792.458/1000.
            quant=kwargs['quantile']
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
            PROP = PROP7['velocity'].copy()
            if color is not None: PROP['color'] = color
            if marker is not None: PROP['marker'] = marker
            if markersize is not None: PROP['markersize'] = markersize
            if label is not None: PROP['label'] = label
            if ls is not None: PROP['ls'] = ls
            if fillstyle is not None: PROP['fillstyle'] = fillstyle                
            self.ax6.errorbar(xx,yy,yerr=yye,**PROP)
            
        if show_fits and 'fitcls' in self.__dict__ and 'specv_evolution' in self.fitcls:
            for model in self.fitcls['specv_evolution']:
                _model = self.fitcls['specv_evolution'][model]
                if kwargs['specv_evolution_xrangep'] is None:
                    x, y, w = _model.predict_random(
                        limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                        x_pred=None, step = 1,
                    )
                else:
                    x, y, w = _model.predict_random(
                        limit=kwargs['plot_mcmct'], plotnsamples=kwargs['plot_nsamples'],
                        x_pred=np.arange(min(kwargs['specv_evolution_xrangep']),
                                         max(kwargs['specv_evolution_xrangep']), 1), 
                    ) 
                for nn, (xx, yy) in enumerate(zip(x, y)):
                    PROP = PROP7['velocity_model'].copy()
                    if color is not None: PROP['color'] = color                    
                    if ls is not None: PROP['ls'] = ls                    
                    self.ax6.plot(xx, yy, alpha=kwargs['alphasample'], **PROP)
                    
                # peak velocity
                v = self.set_vexp(model_name=model, returnv=True)
                self.ax6.errorbar(0, v[1], yerr=(abs(v[0]-v[1])+abs(v[2]-v[1]))/2.,
                                  marker='o', fillstyle='none', color='red')
                
        #  subplot settings        
        if kwargs['ax6_xlim'] is not None:
            x1,x2 = min(kwargs['ax6_xlim']), max(kwargs['ax6_xlim'])            
            self.ax6.set_xlim([x1, x2])                  
        if kwargs['ax6_ylim'] is not None:
            self.ax6.set_ylim(kwargs['ax5_ylim'])
        if fontsize is None: fontsize = 12
        if show_title:
            self.ax6.set_title('%s\n%s z=%.2f' % (self.objid,self.sntype,self.z),fontsize=fontsize)
        self.ax6.set_xlabel('$t - T_{r,\mathrm{max}} \; (\mathrm{restframe \; d})$',fontsize=fontsize)
        self.ax6.set_ylabel('$v \; (1000 \mathrm{km/s})$',fontsize=fontsize)
        if show_legend: self.ax6.legend(fontsize=fontsize, frameon=False)
        
    def _to_luminosity(self, m, f, em=None, **kwargs):
        ''' Convert magnitude to luminosity.

        Parameters
        ----------    
        m      :  `float`
             apparent (AB) magnitude.
             For UV/NIR, if input is AB, set ``ab2vega`` True, otherwise leave it False.
        em     :  `float`
             apparent magnitude error
        f     :  `str`
             filter
        ab2vega     :  `bool`
             set it True if input UV/IR data in AB, if they're vega, set to False
        do_kcorr    :  `bool`
             K corrected the photometry or spectra
        
        Returns
        -----------
        flux   :  `float`
             flux, unit in erg/s/A
        dflux  :  `float`
             flux error, unit in erg/s/A
        w      :  `float`
             central wavelength, unit in A
        bandwidth  :  `float`
             band width, unit in A
        '''
        kwargs = self.read_kwargs(**kwargs)
        SN_distance = 10**((self.dm-25)/5)*3.086e24 # in cm                               
        if m is None: return        
        # If UVOT bands are in AB, need to convert to Vega
        if f in ab_to_vega.keys() and kwargs['ab2vega']:
            m -= ab_to_vega[f]
        wlref = central_wavelengths[f]
        fref = zp[f] * 1e-11
        bandwidths = filter_width[f]                        
        if kwargs['do_kcorr']:
            wlref /= (1+self.z)
            fref *= (1+self.z)
            bandwidths /= (1+self.z)
        _flux = 4*np.pi*SN_distance**2*fref*10**(-0.4*m)
        if em is None: _dflux = None
        else: _dflux = 2.5/np.log(10) * _flux * em
        return _flux, _dflux, wlref, bandwidths
    
    def _from_luminosity(self, l, f, el=None, **kwargs):
        ''' Convert luminosity to (AB) magnitude.
        
        Parameters
        ----------    
        l      :  `float`
             luminosity (erg/s)
        el     :  `float`
             luminosity error
        f     :  `str`
             filter        
        ab2vega     :  `bool`
             set it True if want UV/IR data in AB, othersie will return vega mag
        do_kcorr    :  `bool`
             K corrected the photometry or spectra
        
        Returns
        -----------
        m   :  `float`
             AB magnitude, unit in mag
        em  :  `float`
             AB magnitude error, unit in mag        
        '''
        kwargs = self.read_kwargs(**kwargs)
        SN_distance = 10**((self.dm-25)/5)*3.086e24 # in cm                               
        if l is None: return
        wlref = central_wavelengths[f]
        fref = zp[f] * 1e-11
        bandwidths = filter_width[f]                        
        if kwargs['do_kcorr']:
            wlref /= (1+self.z)
            fref *= (1+self.z)
            bandwidths /= (1+self.z)
        m = -2.5 * np.log10(l/(4*np.pi*SN_distance**2*fref))        
        # If UVOT bands are in AB, need to convert to Vega
        if f in ab_to_vega.keys() and kwargs['ab2vega']:
            m += ab_to_vega[f]                
        if el is None: em = None
        else: em = el * np.log(10) / 2.5 / l
        return m, em
    
    def _cal_spectrum(self, ws, fluxes, efluxes, phase, interpolation=None, index=0, mangling_function='linear', **kwargs):
        ''' Absolute calibrate or mangle the spetrum
        
         Parameters
        ----------    
        ws :  `list`
             spectrum wavelengths, unit in A
        fluxes :  `list`
             spectrum fluxes, arbitrary unit
        efluxes :  `list`
             spectrum flux errors, arbitrary unit
        phase  :  `float`
             spectrum phase, used to get absolute mag from interpolations
        sed_abscal_method :  `str`
             how to use photometry to calibrate the spectra: 

            - [cal]ibrate

            - [mangle] the spetra
                     
        sed_abscal_bands :  `str`             
             which photometry bands to be used to absolute calibrate the selected spectra
        mangling_function :  `str`
             mangling function   
        maxfev  :  `float`
             for scipy. The maximum number of calls to the function. 

        Notes
        ---------------
        More parameters can be input for snobject._mag_at, snobject.sym_mag, snobject._to_luminosity.
        
        Returns
        ------------
        fluxes  :  `list`
             spectrum fluxes, unit in erg/s/A
        efluxes :  `list`
             spectrum flux errors, unit in erg/s/A
        
        See Also
        -------------
        snobject._mag_at, snobject.sym_mag, snobject._to_luminosity
        '''
        kwargs = self.read_kwargs(**kwargs)
        filters = kwargs['sed_abscal_bands']
        if filters is None: filters = np.unique(self.lc['filter'])
        
        assert len(ws) == len(fluxes) and len(ws) > 0        
        wls, ratios = [], []
        for filt in filters:
            m, em = self._mag_at(filt, phase, interpolation=interpolation, index=index, **kwargs)
            if m is None or em is None: continue
            _flux, _dflux, wlref, bandwidths = self._to_luminosity(m, filt, em=em, **kwargs)
            ff = self.sym_mag(ws, fluxes, filt) 
            if ff is None: continue
            wls.append(wlref)
            ratios.append ( _flux/ff )                
        if len(ratios) == 0: return None, None
        if kwargs['sed_abscal_method'] == 'cal':
            ratio = np.mean(ratios)
            fluxes *= ratio
            efluxes *= ratio
            return fluxes, efluxes
        else:
            assert mangling_function in ['linear', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6']
            func = eval('models.polynomial.%s' % mangling_function)
            params, covar = curve_fit(func, wls, ratios, method='trf', maxfev=kwargs['maxfev'])
            perr = np.sqrt(np.diag(covar))
            fit = func(ws, *params)
            fluxes *= fit
        return fluxes, efluxes
    
    def show_corner(self, ax, index=0, gp=False, engine=None,
                    model=None, source=None, filts=None, **kwargs):
        ''' Make contour plots for model samplings.
        
         Parameters
        ----------    
        ax :  `matplotlib.subplot`
             matplotlib.subplot         
        index   :   `int`
              if there're multple models can provide, e.g. the nickel mass mni, which one to use
        gp    :   `bool`
              if show samples from GP or from fittings
        engine  :   `str`
              which engine models to be used
        model   :   `str`
              which models to be used
        source  :   `str`
              which source to be used
        filts   :   `str`
              for `multiband` engines, which filter to show
        clobber  :   `bool`
              Redo analysis
        verbose  :   `bool`
              Enable progress report
        figsize  :   `tuple`
              figure size         
        dpi     :   `int`
              matplotlib resolution
        
        See Also
        -------------
        snobject.get_model
        '''        
        kwargs = self.read_kwargs(**kwargs)
        if gp:
            assert 'gpcls' in self.__dict__, 'Error: no GP had been done...'
        else:
            assert 'fitcls' in self.__dict__, 'Error: no fittings had been done...'            
        kwargs = self.read_kwargs(**kwargs)
        
        if gp: # gp
            figpath = '{}_gpcorner'.format(self.objid)
            fignames = self.gpcls.save_corner(figpath,
                    datadir='%s/plots/'%LOCALSOURCE, clobber=kwargs['clobber'])
        else:  # fit                        
            models, modelnames, enginenames, sourcenames = self.get_model(
                engine=engine, model=model, source=source
            )            
            if len(models) == 0:
                if kwargs['verbose']: print ('No models found')
                return
            elif len(models) > 1:
                if kwargs['verbose']: print ('!!! More than one model found: %s'%modelnames)             
            _model, _modelname, _enginename, _sourcename = models[index],\
                modelnames[index], enginenames[index], sourcenames[index]
            if _sourcename is not None:
                figpath = '{}_fitcorner_{}_{}'.format(self.objid, _modelname, _sourcename)
            else:
                figpath = '{}_fitcorner_{}'.format(self.objid, _modelname)            
            fignames = _model.save_corner(figpath,
                datadir='%s/plots/'%LOCALSOURCE, clobber=kwargs['clobber'], filts=filts)
        if fignames is None or len(fignames) == 0:
            if kwargs['verbose']: print ('No Contours found')
            return
        plt.figure(constrained_layout=True, figsize=kwargs['figsize'], dpi=kwargs['dpi'])
        if len(fignames) > 1: divider = make_axes_locatable(ax)
        for nn, figname in enumerate(fignames):
            if kwargs['verbose']: print ('show corner plot %i- %s'%(nn,figname))
            img = mpimg.imread(figname)
            if nn == 0:
                ax.imshow(img)
                ax.axis("off") 
            else:              
                _ax = divider.append_axes("bottom", size="100%", pad="2%")
                _ax.imshow(img)
                _ax.axis("off")            
                
    def all_fittings(self):
        ''' Show all fittings done so far.
        
        Parameters
        ----------    
        N/A
        '''
        for engine_name in self.fitcls:
            print ('->', engine_name, ':')
            if engine_name in ['specline', 'bol_early', 'bol_main', 'bol_tail', 'bol_full', 'sed']:
                for source_name in self.fitcls[engine_name]:                    
                    for model_name in self.fitcls[engine_name][source_name]:  
                        print(' '*10, source_name, ': ', model_name )
            else:
                for model_name in self.fitcls[engine_name]:  
                    print( ' '*10, model_name )                
        
    def get_model(self, engine=None, model=None, source=None):
        ''' Get model with specific model name, and/or engine name, and/or source name.
        
        Parameters
        ----------    
        engine  :   `str`
              which engine models to be used
        model   :   `str`
              which models to be used
        source  :   `str`
              which source to be used

        Returns
        ----------
        models   :   `list`
              model list
        modelnames :   `list`
              model name list
        enginenames :   `list`
              engine name list
        sourcenames  :   `list`
              source name list
        '''
        assert 'fitcls' in self.__dict__
        models, modelnames, enginenames, sourcenames = [], [], [], []
        
        for engine_name in self.fitcls:            
            if engine is not None and engine != engine_name: continue
            if engine_name in ['specline', 'bol_early', 'bol_main',
                               'bol_tail', 'bol_full', 'sed']:
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
            
    def get_par(self, par, index=0, engine=None, model=None,
                source=None, filt='r', k_opt=None, **kwargs):
        ''' Get specific parameter value and error, of a fitted model.
        
        Parameters
        ----------    
        par  :   `str`
              parameter name
        index   :   `int`
              if there're multple models can provide, e.g. the nickel mass mni, which one to use         
        engine  :   `str`
              which engine models to be used
        model   :   `str`
              which models to be used
        source  :   `str`
              which source to be used
        filts   :   `str`
              for `multiband` engines, which filter
        quantile   :  `list`
              use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors        
        k_opt  :  `float`
              for Arnett alike models, which opacity to be used:

              - None: will use sntype instead to set opacity, if not defined, use 0.1
        
              - float: opacity number
        
        Returns
        ----------
        pars   :   `dict`        
        
        See Also
        ----------
        snobject._get_par
        '''
        kwargs = self.read_kwargs(**kwargs)
        quant=kwargs['quantile']
        _dict = self._get_par(par, None, par2=None, engine=engine, model=model,
                              source=source, filt=filt, quant=quant)                     
        if par == 'mni':
            def func(a,b): return a*b
            _dict1 = self._get_par('fni', func, 'mej', engine, model, source, filt, quant)
        elif par == 'taum':
            def func(a,b): return models.arnett_tail.Mej_Ek_to_taum(a, b, k_opt=k_opt, sntype=self.sntype)
            _dict1 = self._get_par('mej', func, 'ekin', engine, model, source, filt, quant)
        elif par in ['ekin','mej'] and self.vexp is not None:
            def func(a):
                v_ej = self.vexp[1]
                vej_err = (abs(self.vexp[0]-self.vexp[1])+abs(self.vexp[0]-self.vexp[1]))/2.
                vm, vme = dessart_vej_to_vm(v_ej, vej_err, self.sntype, False)
                if vm is not None:                    
                    mej, ek = models.arnett_tail.taum_to_Mej_Ek(
                        a, vm, taum_err=None, vej_err=vme, k_opt=k_opt, sntype=self.sntype)
                else:                    
                    mej, ek = models.arnett_tail.taum_to_Mej_Ek_1(
                        a, v_ej, taum_err=None, vej_err=vej_err, k_opt=k_opt, sntype=self.sntype)
                if par == 'ekin': return ek
                else: return mej
            _dict1 = self._get_par('taum', func, None, engine, model, source, filt, quant)
        else:
            _dict1 = dict()
        if _dict is None: return _dict1
        else:  return {**_dict, **_dict1}
        
    def _get_par(self, par, operator, par2=None, engine=None, model=None,
                source=None, filt='r', quant=[.05,.5,.95]):
        ''' Get specific parameter value and error, of a fitted model.
        
        Parameters
        ----------    
        par  :   `str`
              parameter name
        operator   :   `func`
              used when the inferred parameter is a function of par and/or par2
        par2    :  `str`
              parameter 2 name
        engine  :   `str`
              which engine models to be used
        model   :   `str`
              which models to be used
        source  :   `str`
              which source to be used
        filts   :   `str`
              for `multiband` engines, which filter
        quant   :  `list`
              use 50 percentile as mean, and 1 sigma (68%% -> 16%% - 84%%) as errors        
        
        Returns
        ----------
        pars   :   `dict`        
        
        See Also
        ----------
        snobject.get_par
        '''
        _dict = dict()
        if not 'fitcls' in self.__dict__:
            #print ('Warning: do fitting first')
            return        
        models, modelnames, enginenames, sourcenames = self.get_model(
            engine=engine, model=model, source=source
        ) 
        for _model, _modelname, _enginename, _sourcename in zip(models,\
            modelnames, enginenames, sourcenames):
            pars = get_pars(_modelname, with_alias=True)['par']            
            if not par in pars: continue
            if _model is None: continue
            v1 = _model.get_par(filt=filt, quant=quant, parname=par)            
            if operator is None:
                _dict['%s %s' % (_modelname, _sourcename)] = v1
            elif par2 is None:                                
                _dict['%s %s' % (_modelname, _sourcename)] = \
                    (operator(v1[0]), operator(v1[1]), operator(v1[2]))
            elif par2 in pars:                 
                v2 = _model.get_par(filt=filt, quant=quant, parname=par2)            
                _dict['%s %s' % (_modelname, _sourcename)] = \
                    (operator(v1[0],v2[0]), operator(v1[1],v2[1]), operator(v1[2],v2[2]))
        return _dict
    
    def savefig(self, **kwargs):
        ''' Save snobject.fig to data directory.
        
        Parameters
        ----------    
        figsize  :   `tuple`
              figure size         
        figpath :   `str`
              figure file name  
        dpi     :   `int`
              matplotlib resolution
        verbose     :   `bool`
              Enable progress report

        See Also
        ------------
        snobject.showfig
        '''
        kwargs = self.read_kwargs(**kwargs)
        if self.fig is None:
            if kwargs['verbose']: print ('Error: no fig to be used to savefig')
            return
        self.fig.set_size_inches(kwargs['figsize'][0], kwargs['figsize'][1])                
        if '%s' in kwargs['figpath']:
            figpath = '{}/plots/{}'.format(LOCALSOURCE, kwargs['figpath']%self.objid)
        else:
            figpath = '{}/plots/{}'.format(LOCALSOURCE, kwargs['figpath'])
        self.fig.savefig(figpath, dpi=kwargs['dpi'], bbox_inches='tight')
        if kwargs['verbose']: print ('saved fig to %s'%figpath)
        
    def showfig(self, ax=None, **kwargs):
        ''' Show the stored figure, that was saved by the snobject.savefig.
        
        Parameters
        ----------    
        ax   :   `matplotlib.subplot`
              matplotlib.subplot, if None, will create one instead.
        figsize  :   `tuple`
              figure size         
        figpath :   `str`
              figure file name  
        dpi     :   `int`
              matplotlib resolution
        verbose     :   `bool`
              Enable progress report

        See Also
        ------------
        snobject.savefig
        '''
        kwargs = self.read_kwargs(**kwargs)        
        if '%s' in kwargs['figpath']:                    
            figpath = '{}/plots/{}'.format(LOCALSOURCE, kwargs['figpath']%self.objid)
        else:
            figpath = '{}/plots/{}'.format(LOCALSOURCE, kwargs['figpath'])
        if not os.path.exists(figpath):
            if kwargs['verbose']: print ('savefig first')
            return
        plt.figure(constrained_layout=True, figsize=kwargs['figsize'], dpi=kwargs['dpi'])
        img = mpimg.imread(figpath)
        if ax is None:
            imgplot = plt.imshow(img)
            ax = plt.gca()
        else:
            ax.imshow(img)
        ax.axis("off")
    
    def dm_error(self, filt, verbose=False):
        ''' Estimate the distance module errors, including errors from milky ebv, redshift and hubble constants.
        
        Parameters
        ----------    
        filt   :   `str`
              filter        
        verbose     :   `bool`
              Enable progress report

        Returns
        ------------
        dme  :  `float`
              distance module error, unit in mag   
        '''        
        if filt not in Rf:
            if verbose: print('Warning: no R defined for %s, ignore its extinction then'%filt)
            R = 0
        else:
            R=Rf[filt]
        
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
    def read_c10(filename='c10_template.txt', verbose=True):
        ''' Read colour templates.
        
        Parameters
        ----------    
        filename   :   `str`
              template file name.
        verbose     :   `bool`
              Enable progress report

        Returns
        ---------
        color   :   `dict`
              template colors at 10d post peak
         '''
        if os.path.isfile( filename ):
            if verbose: print (">> WARNING: using local color file (%s)"%filename)            
        else:
            filename = '%s/%s' % (LOCALSOURCE, filename)
        assert os.path.exists(filename)
        c10 = dict()
        for ll in open(filename, encoding='utf-8').readlines():
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
        ''' Convert from magnitude to fluxes.

        Parameters
        ----------    
        mag   :   `list`
              magnitudes
        magerr     :   `list`
              magnitude errors
        limmag     :   `list`
              limiting magnitude
        sigma  :   `float`
              signal noise ratio threshold for detections 
        units  :   `str`
              unit system:

              - [zp] zero point
        
              - [phys] counts
        
        zp   :   `float`
              zero points from mag to flux
        wavelength  :   `float`
              central wavelength
        
        Returns
        ---------
        flux   :   `list`
              fluxes
        dflux   :   `list`
              flux errors

        See Also
        ---------
        snobject.flux_to_mag
        '''
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
        ''' Converts from fluxes into magnitudes.

        Parameters
        ----------    
        flux   :   `list`
              fluxes
        dflux     :   `list`
              flux errors        
        sigma  :   `float`
              signal noise ratio threshold for detections 
        units  :   `str`
              unit system:

              - [zp] zero point
        
              - [phys] counts
        
        zp   :   `float`
              zero points from mag to flux
        wavelength  :   `float`
              central wavelength
        
        Returns
        ---------
        mag   :   `list`
              magnitudes
        dmag   :   `list`
              magnitude errors
        limmag     :   `list`
              limiting magnitude
        
        See Also
        ---------
        snobject.mag_to_flux
        '''        
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
    def sym_mag(w, f, filt, verbose=False):
        ''' Synthetic magnitudes with filter transmission curves.
        
        Parameters
        ---------- 
        w   :        `list`
               spectral wavelengths
        f :       `list`
               spectral fluxes
        filt   :        `float`
               filter
        verbose     :   `bool`
              Enable progress report
        
        Notes
        ---------
        Currently we had only ztf g/r/i filter transmission file avaliable at 
        https://github.com/saberyoung/HAFFET/tree/master/sdapy/data/filter_transmission,
        One can add their own filters into the directory.
        
        Return
        ---------- 
        mag    :        `float`
               synthetic magnitudes
        '''         
        tfile = '%s/data/filter_transmission/ztf_%s_band.csv' % (srcpath, filt)        
        if os.path.exists(tfile):
            df = pd.read_csv(tfile, names=['w','T'])
            wnm = w/10. # AA to nm
    
            # interpolation
            tfunc = interp1d(df['w'], df['T'], kind="zero")
            _idin = np.logical_and(wnm > min(df['w']), wnm < max(df['w']))
            __a = np.trapz(tfunc(wnm[_idin])*wnm[_idin], wnm[_idin])
            __b = np.trapz(f[_idin]*tfunc(wnm[_idin])*wnm[_idin], wnm[_idin])
            __flux =  __b / __a
            return __flux
        elif verbose: print ('Error: %s is missing' % tfile)
        return
        
    @staticmethod
    def bin_df(df, deltah = 1., xkey='jdobs'):        
        ''' Bin dataframe. 
        Will add one column, ``jdbin``, into the input dataframe.
        
        Parameters
        ---------- 
        df   :        `pandas dataframe`
               snobject.lc
        deltah :       `float`
               binning factor, unit in days
        xkey   :        `str`
               jd key                
        
        Return
        ---------- 
        df    :        `pandas dataframe`
               snobject.lc
        '''
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
        ''' Merge dataframes
        
        Parameters
        ---------- 
        df   :        `pandas dataframe`               

        Returns
        ---------- 
        df   :        `pandas dataframe`
        '''
        rows = __df.shape[0]
        __df1, __n = __df.values[0], 1
        for __nk, _kk in enumerate(range(rows)):
            if __nk==0:continue
            __df1 += np.array(__df.values[_kk])
            __n += 1
        _arr = [__df2[0:int(len(__df2)/__n)] if type(__df2) is str else float(__df2)/__n for __df2 in __df1]
        return _arr

def main():    
    print_logo()
    description="run analysis on a list of SNe"    
    parser = argparse.ArgumentParser(description=description,\
         formatter_class=argparse.ArgumentDefaultsHelpFormatter)   
    parser.add_argument('-s', '--source', help='meta source', default='BTS')
    parser.add_argument('--syntax', help='meta syntax', default='type in ["Ib","Ic"]')
    parser.add_argument('-z', '--ztfid', help='specify ztfid of sn')
    parser.add_argument('-i', '--iauid', help='specify iauid of sn')
    parser.add_argument('-t', '--sntype', help='specify type of sn')
    parser.add_argument("-v", "--verbose",dest="verbose",default=False,\
                        action="store_true",help='Enable progress report')
    parser.add_argument("-c", "--clobber",dest="clobber",default=False,\
                        action="store_true",help='Redo analysis')
    parser.add_argument("-d", "--debug",dest="debug",default=False,\
                        action="store_true",help='If failed fittings, ignore or crashed')
    parser.add_argument("--checkdir",dest="checkdir",default=False,\
                        action="store_true",help='Check data directory')    
    args = parser.parse_args()  
    pars = vars(args)    
    
    # run    
    if args.checkdir:
        check_dir(check=False, askifdirmissed=True)        
    else:
        if pars['ztfid'] is not None:  pars.update(syntax = 'objid =="%s"'%pars['ztfid'])
        if pars['iauid'] is not None:  pars.update(syntax = 'alias =="%s"'%pars['iauid'])
        if pars['sntype'] is not None: pars.update(syntax = 'type in ["%s"]' %pars['sntype'])    
        cls = snelist( **pars )
        cls.run(debug=args.debug)
