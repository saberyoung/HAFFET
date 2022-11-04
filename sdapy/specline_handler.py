# !/usr/bin/env python3
# -*- coding: UTF-8 -*-

import extinction
import numpy as np
from scipy.ndimage.filters import gaussian_filter, generic_filter, median_filter
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from astropy.time import Time
from sdapy.functions import *
from sdapy.filters import *
from sdapy.constants import line_location

class handle_spectra:
    """Handle spectra for one SN
    """
    
    def __init__(self, z, ebv, t0=None, **kwargs):                 
        # ----- read meta ----- #
        self.kwargs = dict()
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        self.data = dict()        
        self.ebv = ebv
        self.z = z        
        self.rv = kwargs.get('rv', 3.1)
        self.fnorm = 1        # flux normalizarion
        self.spec_shift = 2.  # spectra shift as fnorm / *shift*
        self.ys = 0
        # set t0 as mjd
        if t0 is None: self.t0 = 0
        elif t0 > 2400000: self.t0 = t0-2400000.5
        else: self.t0 = t0
        
    def define_phase(self, epoch):
        if is_date(epoch):
            if '-' in epoch: t = epoch
            else: t = '%s-%s-%s'%(epoch[:4],epoch[4:6],epoch[6:])
            at = Time(t, format='isot', scale='utc')
            phase = (at.mjd-self.t0)/(1+self.z)
        elif is_number(epoch):
            if float(epoch) < 2400000.5:
                phase = (float(epoch)-self.t0)/(1+self.z)
            else:
                phase = (float(epoch)-2400000.5-self.t0)/(1+self.z)
        else:
            print ('!!!', epoch)
        return '%.2f' % phase
    
    def add_spectrum(self, wave, flux, instru, epoch, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        verbose = kwargs['verbose']        
        sid = '%s %s'%(epoch, instru)
        if sid in self.data and not kwargs['clobber']:
            if verbose: print ('%s %s already added, skip'%(epoch,instru))
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
        self.ys += self.fnorm/self.spec_shift
        
    def sort_spectra(self, returnepochs=False):
        # rearrange ys based on phase
        epochs, phases, yss = [], [], []
        for _ in self.data:
            phases.append( float(self.data[_]['phase']) )
            epochs.append( _ )
            yss.append( self.data[_]['ys'] )        
        __ = np.argsort(phases)[::-1]
        if returnepochs: return np.array(epochs)[__][::-1]
        for _e, _ys in zip(np.array(epochs)[__], sorted(yss)):
            self.data[_e]['ys'] = _ys
        
class handle_spectrum:
    """Handle one single spectrum
    """  
    def __init__(self, wave, flux, z, ebv, flux_err=None, rv=3.1, fnorm=1, **kwargs):                  
        # ----- read meta ----- #
        self.kwargs = dict()
        for _key in self.kwargs:
            kwargs.setdefault(_key, self.kwargs[_key])
        
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
        self._correct_extinction()
        self._bin_spectrum(**kwargs)
        self._estimate_continuum(**kwargs)
        self._correct_continuum()   
        
    def _correct_extinction(self):
        """Rest frame spectra and correct for MW extinction
        """                
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

    def _estimate_continuum(self, **kwargs):
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

    def _correct_continuum(self):
        ''' correct continuum
        ''' 
        assert self.bin_wave is not None, 'do extinction correction first'
        assert self.continuum_wave is not None, 'do extinction correction first'

        self.flat_wave = self.bin_wave
        self.flat_flux = self.bin_flux - self.continuum_flux
        
    def _norm_spectrum(self, region=None, stype='bin', return_func=False):
        ''' 
        normalize spectrum
        if region available, cut region and then normlalization
        '''        
        wave, flux = self.get_data(stype)        
        if region is not None:
            # Get indices for beginning and end of the feature
            feat_start = min(region)
            feat_end = max(region)                             
            __ = np.logical_and(self.rest_wave >= feat_start, self.rest_wave <= feat_end)
            wave, flux = wave[__], flux[__]            
                        
        # normalize spectrum
        if return_func:
            def func(x): return (x-min(flux))/(max(flux)-min(flux))*self.fnorm
            return func
        else:
            if len(wave) < 10:
                print ('!!! Range too small. Please select a wider range')
            if len(wave) == 0: return [], []
            else: return  wave, (flux-min(flux))/(max(flux)-min(flux))*self.fnorm
            
    def get_data(self, stype):
        ''' get different type of spectra
        '''
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

    def get_local_continuum(self, region, stype, continuum_method, continuum_degree):
        ''' smooth data locally '''
        assert region is not None        
        wave, flux = self._norm_spectrum(
            region=region, stype=stype, return_func=False
        )
        
        fsmooth = self._calc_continuum(
            wave, flux, continuum_method, continuum_degree,
        )
        return wave, fsmooth        
    
    def find_peaks(self, element=None, region=None, line_velocity_bounds=None,
                   stype='bin', pfactor=20, smooth=True, continuum_method=None,
                   continuum_degree=None):        
        ''' find peaks for one element or in a spectral range'''
        if region is not None: pass
        elif element is not None and line_velocity_bounds is not None:
            cw, region = self.parse_element(element, line_velocity_bounds)
        else: return dict()
        
        if smooth and continuum_method is not None and continuum_degree is not None:                    
            # find peaks on smoothed spectra
            wave, fsmooth = self.get_local_continuum(region, stype, continuum_method, continuum_degree)
            if len(wave) < 5: return
            return  self._calc_peaks(wave, fsmooth, pfactor)
        else:
            wave, flux = self._norm_spectrum(region=region, stype=stype, return_func=False)
            if len(wave) < 5: return
            return handle_spectrum._calc_peaks(wave, flux, pfactor)
        
    @staticmethod
    def parse_element(element, line_velocity_bounds):
        ''' obtain central wavelength and fitting range for an element '''
        if element is None: return None, None
        assert element in line_location, '%s not included in %s' % (element, line_location)
        cw = line_location[element]
        # blueside: twice of the allowed maximum velocity
        # redside : line intrinstic wavelength
        blueside, redside = calc_wav(max(line_velocity_bounds)*2,cw), cw    
        region = [blueside, cw + abs(blueside-redside)/10.]
        return cw, region
    
    @staticmethod
    def _add_noise(flux, snr):
        """Add noise level of snr to the flux of the spectrum."""        
        sigma = flux / snr
        
        # clear data            
        _ = np.where(sigma <= 0)
        __ = np.where(sigma > 0)
        sigma[_] = min(sigma[__])
        sigma[np.isnan(sigma)] = min(sigma[__])                
        sigma[np.isinf(sigma)] = min(sigma[__])        
        # Add normal distributed noise at the SNR level.        
        return np.random.normal(0, abs(sigma))

    @staticmethod
    def _calc_continuum(wave, flux, method, degree, waveon=None):
        """fit continuum
        """
        if len(wave) < 5: return
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
