# !/usr/bin/env python3
# -*- coding: UTF-8 -*-
from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os,string,re,sys,subprocess,shlex,math,requests,wget
from scipy.spatial import distance
import sewpy
try:    import pyfits
except: from astropy.io import fits as pyfits
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.wcs import WCS
from astroquery.sdss import SDSS
from PIL import Image
from io import StringIO

class handle_trigger:    
    def __init__(self, ra, dec, img=None, ax=None, snname=None,
                 z1=None, z2=None, fov=20, source='ps1', maglim1=10, maglim2=24,
                 filt='r', clobber=False):
        """Handle finder
        """            
        # ----- read meta ----- #           
        self.ra   = ra    # ra
        self.dec  = dec   # dec
        self.img  = img    # file name of a .fits image to use as the finding chart
        self.name = snname # name of the supernova (default is the OBJECT keyword in the header)
        self.z1   = z1     # z-scale in main image (default same as ds9 zscale)
        self.z2   = z2     # z-scale in main image (default same as ds9 zscale)            
        self.fov  = fov  # fov size of main image in arcsec
        self.ax   = ax
        assert source in ['ps1', 'dss']
        self.source = source
        assert filt in 'ugriz'
        self.filt = filt
        self.clobber = clobber
        self.maglim1 = maglim1
        self.maglim2 = maglim2
        
    def run(self):
        
        self.readradec()
        if self.img is None:           
            if self.source == 'ps1':  self.readpsimg()        
        if self.img is None: return
        self.plot()
        
    def readradec(self):
        # read ra dec
        raunit, decunit = u.deg, u.deg
        try:    float(self.ra)            
        except: raunit = u.hourangle
        _p = SkyCoord(self.ra, self.dec, unit=(raunit, decunit))           
        rasec = '{0:0{width}}'.format(float("%.2f"%_p.ra.hms[2]), width=5)   
        self.rahms = '%.2i:%.2i:%s'%(_p.ra.hms[0],abs(_p.ra.hms[1]),rasec)
        self.decdms = '%.2i:%.2i:%.2f'%(_p.dec.dms[0],abs(_p.dec.dms[1]),abs(_p.dec.dms[2]))
        self.radeg = float('%.5f'%_p.ra.deg)
        self.decdeg = float('%.5f'%_p.dec.deg)         

    def readpsimg(self):
        size = self.fov / 0.25
        self.img = download_ps1(self.radeg, self.decdeg, size=size, filt=self.filt)
        if self.img is None or not os.path.exists(self.img):
            print ('Failed download PS1 image...')        
        
    def plot(self):   

        if self.img is None:return
        
        # read data and header of image
        ext, resolution = 0, 1
        if self.source == 'ps1': ext, resolution = 0, .25
        
        image_data = pyfits.getdata(self.img, ext=ext)
        hdr = pyfits.getheader(self.img, ext=ext)
        
        # ds9 zscale
        z1, z2 = zscale(image_data)
        
        # read fits wcs from img
        w = WCS(hdr)        
    
        # image vertices
        v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = \
            vertices(self.radeg, self.decdeg, self.fov/3600., self.fov/3600.)
        ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra],\
                                     [v1_dec, v2_dec, v4_dec, v3_dec])
        left = (ra_vertices[0] + ra_vertices[3])/2.
        right = (ra_vertices[1] + ra_vertices[2])/2.
        bottom = (dec_vertices[2] + dec_vertices[3])/2.
        top = (dec_vertices[0] + dec_vertices[1])/2.
        
        # show image          
        self.ax.imshow(image_data, cmap='gray', vmin=z1, vmax=z2, aspect='equal',
                       interpolation='nearest', origin='lower', extent=(left, right, bottom, top))
        self.ax.set_xlabel('R.A. (degrees)')
        self.ax.set_ylabel('Dec. (degrees)')
        self.ax.get_xaxis().get_major_formatter().set_useOffset(False)
        self.ax.get_yaxis().get_major_formatter().set_useOffset(False)
        
        # mark SN        
        self.ax.plot(self.radeg, self.decdeg, color='r', marker='o',
            fillstyle='none', markersize=12, ls='', label='SN %s %s'%(self.radeg, self.decdeg))
        
        # query SDSS field stars
        cata = retrieve_sloan(self.radeg, self.decdeg, self.fov, self.maglim1, self.maglim2, False)        
        if not cata is None:            
            for n, (_ra, _dec) in enumerate(zip(cata['ra'], cata['dec'])):            
                if n == 0: label='SDSS field stars'
                else: label = None
                self.ax.plot(_ra, _dec, color='g', marker='o',
                             fillstyle='none', markersize=12, ls='', label=label)
        
        self.ax.set_xlim([self.radeg-self.fov/3600./2., self.radeg+self.fov/3600./2.])
        self.ax.set_ylim([self.decdeg-self.fov/3600./2., self.decdeg+self.fov/3600./2.])
        self.ax.legend()  
        
if True:
    #fig,ax=plt.subplots(1,1)
    a=handle_trigger(20, 30, fov = 40, ax=ax, img=None)
    a.run()
    #plt.show()
