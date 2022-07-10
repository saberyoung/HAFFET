# !/usr/bin/env python3
# -*- coding: UTF-8 -*-

import extinction
import numpy as np
from scipy.ndimage.filters import gaussian_filter, generic_filter, median_filter
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
#import logging, random
#from uncertainties.unumpy import nominal_values, std_devs
#from scipy.optimize import curve_fit
#from astropy import units
#from astropy.constants import c
from sdapy.functions import *
from sdapy.filters import *
#from sdapy.corner_hack import quantile
#from sdapy.constants import speclines
from astropy.time import Time

class handle_spectra:    
    def __init__(self, z, ebv, rv=3.1, fnorm=1, t0=None, **kwargs):
        """Handle spectra
        """            
        # ----- read meta ----- #
        self.kwargs = dict()
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        self.data = dict()        
        self.ebv = ebv
        self.z = z        
        self.rv = rv       
        self.fnorm = fnorm
        self.ys = 0
        # set t0 as mjd
        if t0 is None: self.t0 = 0
        elif t0 > 2400000: self.t0 = t0-2400000.5
        else: self.t0 = t0
        
    def define_phase(self, epoch):
        t = '%s-%s-%s'%(epoch[:4],epoch[4:6],epoch[6:])
        at = Time(t, format='isot', scale='utc')
        phase = (at.mjd-self.t0)/(1+self.z)
        return '%.2f' % phase
    
    def add_spectrum(self, wave, flux, instru, epoch, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        sid = '%s %s'%(epoch, instru)
        if sid in self.data and not kwargs['clobber']:
            print ('%s %s already added, skip'%(epoch,instru))
            return
        specdata = handle_spectrum(
            wave, flux, self.z, self.ebv, rv=self.rv, fnorm=self.fnorm, **kwargs
        )        
        self.data[sid] = {
            'epoch' : epoch,
            'phase' : self.define_phase(epoch),
            'instru': instru,
            'data'  : specdata,
            'ys'    : self.ys,
        }
        # y shift
        self.ys += self.fnorm/kwargs['spec_shift']
    
class handle_spectrum:    
    def __init__(self, wave, flux, z, ebv, flux_err=None, rv=3.1, fnorm=1, **kwargs):
        """Handle spectrum
        """            
        # ----- read meta ----- #
        self.kwargs = dict()
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        self.wave = wave
        self.flux = flux        
        self.ebv = ebv
        self.z = z
        self.rv = rv       
        self.fnorm = fnorm
        
        # after extinction
        self.rest_wave = None
        self.rest_flux = None
        # after binning
        self.bin_wave = None
        self.bin_flux = None        
        # continuum
        self.continuum_wave = None
        self.continuum_flux = None
        # after continuum correction
        self.flat_wave = None
        self.flat_flux = None  
        
        # correct extinction and time dilation
        self._correct_extinction(**kwargs)
        self._bin_spectrum(**kwargs)
        self._estimate_continuum(**kwargs)
        self._correct_continuum(**kwargs)   
        
    def _correct_extinction(self, **kwargs):
        """Rest frame spectra and correct for MW extinction
        """
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        # Determine extinction        
        mag_ext = extinction.fitzpatrick99(self.wave, self.rv * self.ebv, self.rv)
        
        # Correct flux to rest-frame
        self.rest_wave = self.wave / (1 + self.z)
        self.rest_flux = self.flux * 10 ** (0.4 * mag_ext)        
        
    def _bin_spectrum(self, **kwargs):
        """Bin a spectrum to a given resolution
        """
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        assert self.rest_wave is not None, 'do extinction correction first'
        wave, flux = self.rest_wave, self.rest_flux        
        assert kwargs['bin_method'] in ['sum', 'average', 'gauss', 'median', 'savgol']
        self.bin_wave = wave
        
        if kwargs['bin_method'] == 'sum':
            self.bin_flux = generic_filter(flux, sum, kwargs['bin_size'])

        elif kwargs['bin_method'] == 'average':
            self.bin_flux = generic_filter(flux, np.average, kwargs['bin_size'])

        elif kwargs['bin_method'] == 'gauss':
            self.bin_flux = gaussian_filter(flux, kwargs['bin_size'])

        elif kwargs['bin_method'] == 'median':
            self.bin_flux = median_filter(flux, kwargs['bin_size'])

        elif kwargs['bin_method'] == 'savgol':           
            self.bin_flux = savgol_filter(flux, kwargs['bin_size'],
                                kwargs['savgol_order'], mode='nearest')

        else:
            raise ValueError(f'Unknown method')

    def _estimate_continuum(self,**kwargs):
        """fit continuum
        """
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert self.rest_wave is not None, 'do extinction correction first'
        
        method = kwargs['continuum_method']
        degree = kwargs['continuum_degree']
        
        # use binning spectra for continuum and test on unbinned spectral wavelength               
        continuum_fit = handle_spectrum._calc_continuum(self.bin_wave,
                        self.bin_flux, method, degree, waveon=self.rest_wave)        
        self.continuum_wave = self.rest_wave
        self.continuum_flux = continuum_fit

    def _correct_continuum(self,**kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        assert self.bin_wave is not None, 'do extinction correction first'
        assert self.continuum_wave is not None, 'do extinction correction first'

        self.flat_wave = self.bin_wave
        self.flat_flux = self.bin_flux - self.continuum_flux
        
    def _norm_spectrum(self, region=None, stype='bin', return_func=False, **kwargs):
        ''' 
        normalize spectrum
        if region available, cut region and then normlalization
        '''
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        wave, flux = self.get_data(stype)        
        if region is not None:
            # Get indices for beginning and end of the feature
            feat_start = min(region)
            feat_end = max(region)                             
            __ = np.logical_and(self.rest_wave >= feat_start, self.rest_wave <= feat_end)
            wave, flux = wave[__], flux[__]            
            assert len(wave) > 10, 'Range too small. Please select a wider range'
            
        # normalize spectrum
        if return_func:
            def func(x): return (x-min(flux))/(max(flux)-min(flux))*self.fnorm
            return func
        else:
            return  wave, (flux-min(flux))/(max(flux)-min(flux))*self.fnorm

    def get_data(self, stype):
        assert stype in ['original', 'rest', 'bin', 'continuum', 'flat']        
        if stype == 'original':
            wave, flux = self.wave, self.flux
        elif stype == 'bin':            
            wave, flux = self.bin_wave, self.bin_flux
        elif stype == 'rest':            
            wave, flux = self.rest_wave, self.rest_flux
        elif stype == 'continuum':            
            wave, flux = self.continuum_wave, self.continuum_flux
        elif stype == 'flat':            
            wave, flux = self.flat_wave, self.flat_flux
        return wave, flux
    
    def _find_peaks(self, region=None, stype='bin', **kwargs):
        # find peaks in binning spectrum
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        wave, flux = self.get_data(stype)
        if region is not None:
            # Get indices for beginning and end of the feature
            feat_start = min(region)
            feat_end = max(region)                             
            __ = np.logical_and(self.rest_wave >= feat_start, self.rest_wave <= feat_end)
            wave, flux = wave[__], flux[__]            
            assert len(wave) > 10, 'Range too small. Please select a wider range'
        pfactor = kwargs['pfactor']
        return handle_spectrum._calc_peaks(wave, flux, pfactor)        
    
    @staticmethod
    def _add_noise(flux, snr):
        """Add noise level of snr to the flux of the spectrum."""            
        sigma = flux / snr
        # Add normal distributed noise at the SNR level.        
        return np.random.normal(0, sigma)

    @staticmethod
    def _calc_continuum(wave, flux, method, degree, waveon=None):
        """fit continuum
        """        
        if method not in ("scalar", "linear", "quadratic", "cubic", "poly", "exponential"):
            raise ValueError("Incorrect method for polynomial fit.")
        
        if method == "poly" and degree is None:
            raise ValueError("No degree specified for continuum method 'poly'.")
        
        poly_degree = {"scalar": 0, "linear": 1, "quadratic": 2, "cubic": 3, "poly": degree}
        if waveon is None: waveon = wave
        if method == "exponential":
            z = np.polyfit(wave, np.log(flux), deg=1, w=np.sqrt(flux))
            p = np.poly1d(z)
            continuum_fit = np.exp(p(waveon))
        else:
            z = np.polyfit(wave, flux, deg=poly_degree[method])
            p = np.poly1d(z)
            continuum_fit = p(waveon)
        return continuum_fit

    @staticmethod
    def _calc_peaks(wave, flux, pfactor):
        """find peaks
        """ 
        xl = dict()
        for symbol in [1,-1]:
            peaks, properties = find_peaks(
                flux*symbol, prominence=(max(flux) - min(flux)) / pfactor,
            )
            if len(peaks) == 0: continue        
            _xl = []
            for pp in peaks:  _xl.append( wave[pp] )
            xl[symbol] = np.array(_xl)
        return xl
