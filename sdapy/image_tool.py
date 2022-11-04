# !/usr/bin/env python3
# -*- coding: UTF-8 -*-
from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os,string,re,sys,subprocess,shlex,math,requests,wget
from scipy.spatial import distance
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

class handle_image:
    """Make finder by downloading images from online archive, e.g. PS1
    """
    def __init__(self, ra, dec, img=None, ax=None,
                 z1=None, z2=None, fov=20, source='ps1', maglim1=10, maglim2=24,
                 filt='r', clobber=False):                 
        # ----- read meta ----- #           
        self.ra   = ra    # ra
        self.dec  = dec   # dec
        self.img  = img    # file name of a .fits image to use as the finding chart        
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
        if self.img is None or not os.path.exists(self.img):           
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
        download_ps1(self.radeg, self.decdeg, size=size, filt=self.filt, fname=self.img)
        if not os.path.exists(self.img):
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
        
def vertices(ra,dec,fovw,fovh):
    """finding the vertices of a FoV by giving
    the central location and the FoV size
    """
    fovw,fovh = fovw/2.,fovh/2.
    vert_ra,vert_dec=[],[]
    ra_rad,dec_rad,fovw_rad,fovh_rad = np.deg2rad(ra), np.deg2rad(dec),\
        np.deg2rad(fovw), np.deg2rad(fovh)
    for i,j in zip([-fovw_rad, fovw_rad, fovw_rad, -fovw_rad],\
                   [fovh_rad, fovh_rad, -fovh_rad, -fovh_rad]):
        arg = -i/(np.cos(dec_rad)-j*np.sin(dec_rad))
        v_ra = np.rad2deg(ra_rad+np.arctan(arg))       
        v_dec = np.rad2deg(np.arcsin((np.sin(dec_rad)+\
                                j*np.cos(dec_rad))/(1+i**2+j**2)**0.5))
        vert_ra.append(v_ra)
        vert_dec.append(v_dec)
    return vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], \
        vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]        

def zscale(image, nsamples=1000, contrast=0.25, max_reject=0.5, min_npixels=5, krej=2.5, max_iterations=5):

    # Sample the image
    image = np.asarray(image)
    image = image[np.isfinite(image)]
    stride = int(max(1.0, image.size / nsamples))
    samples = image[::stride][:nsamples]
    samples.sort()
    
    npix = len(samples)
    zmin = samples[0]
    zmax = samples[-1]

    # Fit a line to the sorted array of samples
    minpix = max(min_npixels, int(npix * max_reject))
    x = np.arange(npix)
    ngoodpix = npix
    last_ngoodpix = npix + 1

    # Bad pixels mask used in k-sigma clipping
    badpix = np.zeros(npix, dtype=bool)

    # Kernel used to dilate the bad pixels mask
    ngrow = max(1, int(npix * 0.01))
    kernel = np.ones(ngrow, dtype=bool)

    for niter in range(max_iterations):
        if ngoodpix >= last_ngoodpix or ngoodpix < minpix:
            break

        fit = np.polyfit(x, samples, deg=1, w=(~badpix).astype(int))
        fitted = np.poly1d(fit)(x)
        
        # Subtract fitted line from the data array
        flat = samples - fitted

        # Compute the k-sigma rejection threshold
        threshold = krej * flat[~badpix].std()

        # Detect and reject pixels further than k*sigma from the fitted line
        badpix[(flat < - threshold) | (flat > threshold)] = True

        # Convolve with a kernel of length ngrow
        badpix = np.convolve(badpix, kernel, mode='same')

        last_ngoodpix = ngoodpix
        ngoodpix = np.sum(~badpix)

    slope, intercept = fit

    if ngoodpix >= minpix:
        if contrast > 0:
            slope = slope / contrast
        center_pixel = (npix - 1) // 2
        median = np.median(samples)
        zmin = max(zmin, median - (center_pixel - 1) * slope)
        zmax = min(zmax, median + (npix - center_pixel) * slope)
    return zmin,zmax

def download_sdss(img,ra,dec,size,filtro,sdssimg='tmpsdss.fits.bz2'):     
    if os.path.exists(sdssimg):os.remove(sdssimg)
    t = Table.read(sch_dir+'/sdss_fields.fits')  
    _dist = np.sqrt((t['RA']-ra)**2+(t['DEC']-dec)**2)
    _csize = True
    _chc = 0
    while _csize:       
        _tab = t[np.where(_dist==sorted(_dist)[_chc])]
        rerun, run, camcol, field = _tab['RERUN'], _tab['RUN'], _tab['CAMCOL'], _tab['FIELD']
        fitsurl = 'http://data.sdss3.org/sas/dr12/boss/photoObj/frames/'+\
                  '%s/%s/%s/frame-%s-00%s-%s-00%s.fits.bz2'%\
                  (rerun[0], run[0], camcol[0],filtro,run[0],camcol[0],field[0])    
        wget.download(fitsurl, out=imgnew)
        fsize = os.path.getsize(imgnew)
        if fsize < 10000:
            os.remove(imgnew)
            _chc+=1            
        else:_csize=False       
    return imgnew

def download_ps1(ra,dec,size=240,filt="r",fname=None,clobber=False):
    table = getpsimages([ra], [dec], size=int(size), filters=filt,
                        format="fits", imagetypes="stack")
    for row in table:
        ra = row['ra']
        dec = row['dec']
        projcell = row['projcell']
        subcell = row['subcell']
        filter = row['filter']
        
        # create a name for the image -- could also include the projection cell or other info
        if fname is None:
            fname = "t{:08.4f}{:+07.4f}{:.1f}.{}.fits".format(ra,dec,size,filter)
        if os.path.exists(fname):
            if clobber:  os.remove(fname)
            else: return fname
        url = row["url"]
        r = requests.get(url)
        open(fname,"wb").write(r.content)
        
def getpsimages(tra, tdec, size=240, filters="grizy", format="fits", imagetypes="stack"):
     
    """Query ps1filenames.py service for multiple positions to get a list of images
    This adds a url column to the table to retrieve the cutout.
     
    tra, tdec = list of positions in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    format = data format (options are "fits", "jpg", or "png")
    imagetypes = list of any of the acceptable image types.  Default is stack;
        other common choices include warp (single-epoch images), stack.wt (weight image),
        stack.mask, stack.exp (exposure time), stack.num (number of exposures),
        warp.wt, and warp.mask.  This parameter can be a list of strings or a
        comma-separated string.
 
    Returns an astropy table with the results
    """
    ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    fitscut = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    # if imagetypes is a list, convert to a comma-separated string
    if not isinstance(imagetypes,str):
        imagetypes = ",".join(imagetypes)
    # put the positions in an in-memory file object
    cbuf = StringIO()
    cbuf.write('\n'.join(["{} {}".format(ra, dec) for (ra, dec) in zip(tra,tdec)]))
    cbuf.seek(0)
    # use requests.post to pass in positions as a file
    r = requests.post(ps1filename, data=dict(filters=filters, type=imagetypes),
        files=dict(file=cbuf))
    r.raise_for_status()
    tab = Table.read(r.text, format="ascii")
 
    urlbase = "{}?size={}&format={}".format(fitscut,size,format)
    tab["url"] = ["{}&ra={}&dec={}&red={}".format(urlbase,ra,dec,filename)
            for (filename,ra,dec) in zip(tab["filename"],tab["ra"],tab["dec"])]
    return tab

def retrieve_sloan(racen, deccen, rwidth, maglim1, maglim2, _write):

    dew = rwidth/3600./2.
    raw = dew/np.cos(deccen*np.pi/180.)
    rara = [racen-raw, racen+raw]
    decra = [deccen-dew, deccen+dew]

    sloan = SDSS.query_sql('SELECT ra,dec,u,err_u,g,err_g,r,err_r,i,err_i,z,err_z FROM Star WHERE ra BETWEEN %f AND %f AND dec BETWEEN %f and %f AND g BETWEEN %f AND %f ' % (
        rara[0], rara[1], decra[0], decra[1], maglim1, maglim2), data_release=14)

    if not sloan:
        print("!!! WARNING: there are no Sloan stars in field")
        return None

    for b in 'ugriz':
        sloan.rename_column('err_'+b, b+'err')
    sloan = Table(sloan, masked=True)
    if _write:
        write_catalog(sloan)

    return sloan

if False:
    fig,ax=plt.subplots(1,1)
    a=handle_image(20, 30, fov = 40, ax=ax, img='tmp.fits')
    a.run()
    plt.show()
