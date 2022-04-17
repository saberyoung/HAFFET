#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/ztfanalysis.py
# Author            : syang <sheng.yang@astro.su.se>
# Date              : 12.11.2019
# Last Modified Date: 11.02.2022
# Last Modified By  : syang <sheng.yang@astro.su.se>

import os, configparser
filepath = os.path.join(os.path.dirname(__file__))
filename = 'auth.txt'

def read_default(wdir=None):

    if wdir is None: wdir = filepath
    if  os.path.isfile( filename ):
        print (">> WARNING: using local default file")
        wdir = ''        
    config = configparser.ConfigParser()
    config.read('%s/%s'%(wdir, filename))
    optlist = {}
    for s in config.sections():
        optlist[s] = {}
        for o in config.options(s):
            optlist[s][o] = config.get(s,o)
    return optlist

def get_keypairs(wdir=None):
    
    optlist = read_default(wdir)
    if len(optlist) == 0:
        print ("!!! no content from default file")
    return optlist
