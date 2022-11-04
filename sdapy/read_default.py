#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/read_default.py
# Author            : syang <sheng.yang@astro.su.se>
# Date              : 12.11.2019
# Last Modified Date: 11.02.2022
# Last Modified By  : syang <sheng.yang@astro.su.se>

import os, configparser, sdapy
filepath = os.getenv('ZTFDATA',"./Data/")

def read_default(filename, wdir=None, verbose=False,
                 return_comments=False, return_type=False):
    ''' Read default parameters from a file.
    
    Parameters
    ----------                
    filename :        `str`
              file name
    wdir :       `str`
              data directory, if None, will be set as HAFFET default data directory, i.e. $ZTFDATA
    verbose   :       `bool`
              show report
    return_comments   :       `bool`
              if True, return comments of each parameters
    return_type   :       `bool`
              if True, return types of each parameters
    
    Returns
    ---------- 
    default_infos   :   `dict`
    '''
    if wdir is None: wdir = filepath
    if os.path.isfile( filename ):
        print (">> WARNING: using local default file (%s)"%filename)
        wdir = '.'
    config = configparser.ConfigParser()
    config.read('%s/%s'%(wdir, filename))
    if verbose: print ('read %s from %s'%(filename, wdir))
    optlist = {}
    for s in config.sections():        
        optlist[s] = {}
        for o in config.options(s):            
            ll = config.get(s,o).split('#')
            assert len(ll) == 3, 'Error: check %s %s' % (s, o)
            val, typ, comment = ll           
            if return_comments:
                optlist[s][o] = comment.strip()
            else:                
                assert typ.strip() in ['str', 'eval']
                if return_type:
                    optlist[s][o] = typ.strip()
                else:                   
                    if typ.strip() == 'str':  optlist[s][o] = str(val.strip())
                    else:  optlist[s][o] = eval(val.strip())
    return optlist

def get_keypairs(filename='auth.txt', wdir=None,
                 verbose=False, return_comments=False, return_type=False):
    ''' get keypairs that would be used by online services, e.g. forced photometry
    
    Parameters
    ----------                
    filename :        `str`
              file name
    wdir :       `str`
              data directory, if None, will be set as HAFFET default data directory, i.e. $ZTFDATA
    verbose   :       `bool`
              show report
    return_comments   :       `bool`
              if True, return comments of each parameters
    return_type   :       `bool`
              if True, return types of each parameters
    
    Returns
    ---------- 
    default_infos   :   `dict`
    '''
    optlist = read_default(filename, wdir, verbose, return_comments, return_type)
    if len(optlist) == 0:
        print ("!!! no content for %s"%filename)
    return optlist

def get_parameters(filename='default_par.txt', wdir=None,
        verbose=False, return_comments=False, return_type=False, keylist=None):
    ''' get default parameters for each SN
    
    Parameters
    ----------                
    filename :        `str`
              file name
    wdir :       `str`
              data directory, if None, will be set as HAFFET default data directory, i.e. $ZTFDATA
    verbose   :       `bool`
              show report
    return_comments   :       `bool`
              if True, return comments of each parameters
    return_type   :       `bool`
              if True, return types of each parameters
    keylist : `list`
              if specified, return those demanded keywords
    
    Returns
    ---------- 
    default_infos   :   `dict`
    '''
    optlist = read_default(filename, wdir, verbose, return_comments, return_type)
    if len(optlist) == 0:
        print ("!!! no content for %s"%filename)
    if keylist is not None:
        _list = dict()
        for k in optlist:
            if k in keylist: _list[k] = optlist[k]
        optlist = _list    
    return optlist
