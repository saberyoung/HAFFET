# !/usr/bin/env python3
# -*- coding: UTF-8 -*-

import extinction
import numpy as np
from scipy.ndimage.filters import gaussian_filter, generic_filter, median_filter
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import logging, random
from uncertainties.unumpy import nominal_values, std_devs
from scipy.optimize import curve_fit
from astropy import units
from astropy.constants import c
from .functions import *
from .models import *
from .filters import *
from .corner_hack import quantile

class handle_spectrum:    
    def __init__(self, wave, flux, z, ebv, cw, region=None,
                 ax=None, instru=None, phase=None,
                 epoch=None, source=None, rv=3.1, fnorm=1,
                 ys=0, mcmc_h5_file=None, clobber=False, **kwargs):
        """Handle spectrum, plot, fit and measures gaussian minima and pEW of spectral features
        """            
        # ----- read meta ----- #
        self.kwargs = dict()
        for _key in kwargs: self.kwargs[_key] = kwargs[_key]        
        
        self.wave = wave
        self.flux = flux
        self.ebv = ebv
        if self.ebv is None: self.ebv = 0
        self.z = z
        self.ax = ax
        self.instru = instru
        self.epoch = epoch
        self.source = source
        self.rv = rv
        self.cw = cw
        self.region = region
        self.fnorm = fnorm
        self.ys = ys
        self.phase = phase
        if mcmc_h5_file is not None and kwargs['datadir'] is not None:
            self.sampfile = '%s/%s.h5' % (kwargs['datadir'], mcmc_h5_file)
        else:
            self.sampfile = None
        self.clobber = clobber
        
        # Place holders for results of intermediate analyses
        self.rest_wave = None
        self.rest_flux = None
        self.bin_flux = None
        self.norm_flux = None
        self.feature_fit = None
        self.velocity = None
        
    def run(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        self._correct_extinction(**kwargs)                
        if self._cut_spectrum(**kwargs):
            self._norm_spectrum(**kwargs)
            self._bin_spectrum(**kwargs)        
            self._find_peaks(**kwargs)
            self._feature_properties(**kwargs)
            if self.feature_fit is not None:
                q, q1, q2 = self.feature_fit.get_par(quant=kwargs['quantile'])
                if q is not None:
                    self.velocity = ( calc_vel(q[1],self.cw), calc_vel(q1[1],self.cw), calc_vel(q2[1],self.cw) )
            else:
                self.velocity = None
                
    def _cut_spectrum(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if self.rest_wave is None:
            raise RuntimeError('Spectrum must be corrected with extinction before going to features')
        
        if self.region is not None:
            # Get indices for beginning and end of the feature
            feat_start = min(self.region)
            feat_end = max(self.region)                             
            __ = np.logical_and(self.rest_wave >= feat_start,
                                self.rest_wave <= feat_end)
            self.wave, self.flux = self.wave[__], self.flux[__]
            self.rest_wave, self.rest_flux = self.rest_wave[__], self.rest_flux[__]            
            if len(self.rest_wave) <= 10:
                print ('Range too small. Please select a wider range')
                return False
        return True
    
    def _norm_spectrum(self, **kwargs):
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        
        if self.rest_flux is None:
            raise RuntimeError('Spectrum must be corrected with extinction before normalization')        
        self.norm_flux = (self.rest_flux-min(self.rest_flux))/(max(self.rest_flux)-min(self.rest_flux))*self.fnorm
    
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
        if self.norm_flux is None:
            raise RuntimeError('Spectrum must be normalized before binning')        
        assert kwargs['bin_method'] in ['sum', 'average', 'gauss', 'median', 'savgol']
        
        if self.rest_wave is None or self.norm_flux is None:
            raise RuntimeError('Spectrum must be corrected for extinction before binning')

        if kwargs['bin_method'] == 'sum':
            self.bin_flux = generic_filter(self.norm_flux, sum, kwargs['bin_size'])

        elif kwargs['bin_method'] == 'average':
            self.bin_flux = generic_filter(self.norm_flux, np.average, kwargs['bin_size'])

        elif kwargs['bin_method'] == 'gauss':
            self.bin_flux = gaussian_filter(self.norm_flux, kwargs['bin_size'])

        elif kwargs['bin_method'] == 'median':
            self.bin_flux = median_filter(self.norm_flux, kwargs['bin_size'])

        elif kwargs['bin_method'] == 'savgol':
            polyorder = 2
            self.bin_flux = savgol_filter(self.norm_flux, kwargs['bin_size'], polyorder, mode='nearest')

        else:
            raise ValueError(f'Unknown method')
    
    def _find_peaks(self,**kwargs):
        # find peaks in binning spectrum
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])                
        x, y = self.rest_wave, self.bin_flux
        peaks, properties = find_peaks(
            y, prominence=(max(y) - min(y)) / kwargs['pfactor'],
        )
        if len(peaks) == 0: return
        
        _xl, _yl = [], []
        for pp in peaks: 
            _xl.append( x[pp] )
            _yl.append( y[pp] )
        _xl, _yl = np.array(_xl), np.array(_yl)
        
        # the peak closest to intrinstic line location
        if kwargs['spec_guess_red']:
            __=np.argmin(abs(_xl - self.cw))
            self.peaks = ( _xl[__],_yl[__] )
        else:
            __=np.argmin(abs(x - self.cw))
            self.peaks = ( self.cw,y[__] )
    
    def _feature_properties(self, **kwargs):
        """Calculate the properties of a single feature in a spectrum
        """
        for _key in self.kwargs: kwargs.setdefault(_key, self.kwargs[_key])
        if not 'peaks' in self.__dict__: return
        
        # feature red side endpoint
        xln,yln = self.peaks
        xl1 = min(self.region)

        # feature blue side endpoint limit
        __ = np.logical_and(self.rest_wave>=xl1, self.rest_wave<=xln)        
        xx, yy = self.rest_wave[__], self.norm_flux[__]
        
        vmin, vmax = min(kwargs['spec_fitr']), max(kwargs['spec_fitr'])
        self.feature_fit = FeatureVelocities(xx, yy, self.cw, nstep=kwargs['scipysamples'],
            opt_routine=kwargs['spec_routine'], sampfile=self.sampfile,
            maxfev=kwargs['maxfev'], clobber=self.clobber,vmin=vmin, vmax=vmax
        )
        self.feature_fit.train()
        
class FeatureVelocities:
    """Represents the velocity calculation for a spectroscopic feature"""

    def __init__(self, wave, flux, rest_frame, opt_routine='trf',
                 sampfile='tmp.h5', maxfev=20000, nstep=10, clobber=False,
                 vmin=0, vmax=50000):
        """Calculates the pEW of a spectroscopic feature

        Args:
            wave (ndarray): The wavelength values of the feature
            flux (ndarray): The flux values for each feature
        """
        assert opt_routine in ['mcmc', 'trf', 'dogbox']
        
        self.rest_frame = rest_frame
        self.wave = wave
        self.flux = flux
        self.opt_routine = opt_routine
        self.maxfev = maxfev
        self.sampfile = sampfile
        self.nstep = nstep
        self.clobber = clobber
        self.vmin = vmin
        self.vmax = vmax
        
        self._gauss_params = None
        self._cov = None
        self._velocity = None
        self.samples = []
        self.lnprob = []
        
    def train(self):
        if os.path.exists( self.sampfile ) and not self.clobber:
            self.samples, self.lnprob = get_samples_scipy( self.sampfile )
        else:
            for i in range(self.nstep):                           
                # random select feature blue side endpoint            
                _x = (max(self.wave)-min(self.wave))/1000.*random.randint(0,1000)+min(self.wave)
                #__ = np.where(self.wave >= _x)
                #xx, yy = self.wave[__], self.flux[__]
                xx, yy = self.wave, self.flux
                if len(xx) <=5: continue
                
                feature = FeatureVelocity(xx, yy, self.rest_frame,
                        opt_routine=self.opt_routine, maxfev=self.maxfev, minimum=_x)

                if feature.velocity is None: continue
                if feature.velocity > self.vmax or feature.velocity < self.vmin: continue                
                self.samples.append( feature._fit_gauss_params() )
                self.lnprob.append( redchisqg(feature.flux, feature.gaussian_fit(), deg=2, sd=None) )
                
            self.samples, self.lnprob = np.array(self.samples), np.array(self.lnprob)                        
            if os.path.exists( self.sampfile ): os.remove( self.sampfile )
            # save sample files            
            hf = h5py.File(self.sampfile, 'w')
            hf.create_dataset('samples', data=self.samples)
            hf.create_dataset('lnprob', data=self.lnprob)
            hf.close()

    def get_par(self, quant=[.05,.5,.95]):
        assert 'samples' in self.__dict__
        if len(self.samples) == 0: return None, None, None
        mean, p1, p2 = [], [], []
        for _ in [0, 1, 2, 3]:           
            xs = self.samples[:,_]
            ql,q,qu = quantile(np.atleast_1d(xs), quant, weights=None)
            mean.append(q)
            p1.append(ql)
            p2.append(qu)        
        return mean, p1, p2

    def predict(self, w_pred=None, step = 1e-1, quant=[.05,.5,.95]):
        ''' output fitting products '''
        if w_pred is None:
            w_pred = np.arange(self.wave.min(), self.wave.max(), step)
        if len(self.samples) == 0: return None, None, None, None
        q1,q,q2 = self.get_par(quant=quant)                
        y = gaussian( w_pred, *q )
        y1 = gaussian( w_pred, *q1 )
        y2 = gaussian( w_pred, *q2 )        
        return np.array(w_pred), np.array(y), np.array(y1), np.array(y2)
    
class FeatureVelocity:
    """Represents the velocity calculation for a spectroscopic feature"""

    def __init__(self, wave, flux, rest_frame,
                 opt_routine='trf', maxfev=20000, minimum=None):
        """Calculates the pEW of a spectroscopic feature

        Args:
            wave (ndarray): The wavelength values of the feature
            flux (ndarray): The flux values for each feature
        """
        assert opt_routine in ['mcmc', 'trf', 'dogbox']
        
        self.rest_frame = rest_frame
        self.wave = wave
        self.flux = flux
        self.opt_routine = opt_routine
        self.maxfev = maxfev
        self.minimum = minimum
        if self.minimum is None:
            self.minimum = np.median(self.wave)
        
        self._gauss_params = None
        self._cov = None
        self._velocity = None
        
    def _fit_gauss_params(self):
        """Fitted an negative gaussian to the binned flux

        Returns:
            A list of fitted parameters
        """

        if self._gauss_params is not None:
            return self._gauss_params

        try:
            self._gauss_params, self._cov = curve_fit(
                f=gaussian,
                xdata=self.wave,
                ydata=self.flux,
                p0=[max(abs(self.flux)), self.minimum, 50., 0],
                method=self.opt_routine,
                maxfev=self.maxfev,
                bounds=((max(abs(self.flux))/5., min(self.wave), 0, -np.inf),
                        (max(abs(self.flux))*5, max(self.wave), 500, np.inf) ))
        except:
            self._gauss_params = None, None, None, None

        return self._gauss_params
            
    @property
    def gauss_amplitude(self):
        """The fitted gaussian amplitude"""

        return self._fit_gauss_params()[0]

    @property
    def gauss_avg(self):
        """The fitted gaussian average"""

        return self._fit_gauss_params()[1]

    @property
    def gauss_stddev(self):
        """The fitted gaussian standard deviation"""

        return self._fit_gauss_params()[2]

    @property
    def gauss_offset(self):
        """The fitted gaussian offset"""

        return self._fit_gauss_params()[3]

    def gaussian_fit(self):
        """Evaluate the gaussian fit of the normalized flux for feature wavelengths

        Returns:
            An array of normalized flux values
        """

        return gaussian(self.wave, *self._fit_gauss_params())

    @property
    def velocity(self):
        """Calculate the velocity of a feature

        Fit a feature with a negative gaussian and determine the feature's
        velocity. Returned value is ``np.nan`` if the fit fails.

        Returns:
            The velocity of the feature in km / s
        """

        if self._velocity is None:
            gauss_avg = self.gauss_avg
            if gauss_avg is None: return
            
            # Calculate velocity
            unit = units.km / units.s
            speed_of_light = c.to(unit).value
            
            self._velocity = speed_of_light * (
                    ((((self.rest_frame - gauss_avg) / self.rest_frame) + 1) ** 2 - 1) /
                    ((((self.rest_frame - gauss_avg) / self.rest_frame) + 1) ** 2 + 1)
            )

        return self._velocity
