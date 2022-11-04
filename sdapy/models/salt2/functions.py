# borrow from https://github.com/sncosmo/sncosmo/blob/master/sncosmo/models.py
# for inteporlations, instead Cpython, here I use pickle to cache the inteporlations

import os
import numpy as np
from sdapy import __path__, constants, models, filters
from joblib import dump, load
srcpath = __path__[0]

class interp1d(object):
    """
    Utility class for providing functionality like
    scipy.interpolate.interp1d while also providing left/right fill values
    like numpy.interp.  Sigh.
    """
    def __init__(self, x, y):
        """
        Return a callable f such that y = f(x)
        """
        self.x = x
        self.y = y
        
    def __call__(self, x):
        return np.interp(x, self.x, self.y)
    
class TimeSeries:
    """
    A series of values associated with a phase and a wavelength,
    e.g. a time series of spectra or a time series of their errors
    """    
    def __init__(self, filename, cachefile='./tmp.clf'):
        """
        Initialize with ASCII file with grid data in the form
          phase wavelength value
        """
        if os.path.exists(cachefile):
            self._wavelengths, self._phases, self._spectra = load(cachefile)['w'], load(cachefile)['p'], load(cachefile)['s']
        else:        
            self._wavelengths = None  #- wavelength of first day, assume others are the same
            self._phases = []         #- Phases in the model file
            self._spectra = []        #- One spectrum interpolation function for each phase
            
            currentday = None
            w = []
            flux = []
            for line in open(filename):
                day, wavelength, value = map(float, line.split())
                if currentday == None:
                    currentday = day

                if day != currentday:
                    self._phases.append(currentday)
                    self._spectra.append(interp1d(w, flux))
                    if self._wavelengths is None:
                        self._wavelengths = np.array(w)
                    currentday = day
                    w = []
                    flux = []
            
                w.append(wavelength)
                flux.append(value)
            
            #- Get the last day of information in there 
            self._phases.append(currentday)
            self._spectra.append(interp1d(w, flux))

            # dump them
            dump({'w': self._wavelengths, 'p': self._phases, 's' : self._spectra}, cachefile)
            
    def spectrum(self, phase, wavelengths=None, extend=True):
        """
        Return spectrum at requested phase and wavelengths.
        Raise ValueError if phase is out of range of model unless extend=True.
        """
        #- Bounds check first
        if phase < self._phases[0] and not extend:
            raise ValueError("phase %.2f before first model phase %.2f" % (phase, self._phases[0]))
        if phase > self._phases[-1] and not extend:
            raise ValueError("phase %.2f after last model phase %.2f" % (phase, self._phases[-1]))

        #- Use default wavelengths if none are specified
        if wavelengths is None:
            wavelengths = self._wavelengths

        #- Check if requested phase is out of bounds or exactly in the list
        if phase in self._phases:
            iphase = self._phases.index(phase)
            return self._spectra[iphase](wavelengths)
        elif phase < self._phases[0]:
            return self._spectra[0](wavelengths)
        elif phase > self._phases[-1]:
            return self._spectra[-1](wavelengths)
            
        #- If we got this far, we need to interpolate phases
        i = np.searchsorted(self._phases, phase)
        speclate = self._spectra[i](wavelengths)
        specearly = self._spectra[i-1](wavelengths)
        dphase = (phase - self._phases[i-1]) / (self._phases[i] - self._phases[i-1] )
        dspec = speclate - specearly
        spec = specearly + dphase*dspec
        
        return spec

    def wavelengths(self):
        """Return array of wavelengths sampled in the model"""
        return self._wavelengths.copy()
    
    def phases(self):
        """Return array of phases sampled in the model"""
        return np.array(self._phases)
        
    def grid(self):
        """Return a 2D array of spectrum[phase, wavelength]"""
        nspec = len(self._phases)
        nwave = len(self._wavelengths)
        z = np.zeros((nspec, nwave), dtype='float64')
        for i, spec in enumerate(self._spectra):
            z[i, :] = spec(self._wavelengths)
        return z

class Salt2Model:
    """The SALT2 Type Ia supernova spectral timeseries model.

    The spectral flux density of this model is given by

    .. math::

       F(t, \\lambda) = x_0 (M_0(t, \\lambda) + x_1 M_1(t, \\lambda))
                       \\times 10^{-0.4 CL(\\lambda) c}

    where ``x0``, ``x1`` and ``c`` are the free parameters of the model,
    ``M_0``, ``M_1`` are the zeroth and first components of the model, and
    ``CL`` is the colorlaw, which gives the extinction in magnitudes for
    ``c=1``.

    Parameters
    ----------
    saltpath : str, optional
        Directory path containing model component files. Default is `None`,
        which means that no directory is prepended to filenames when
        determining their path.   

    Notes
    -----
    The "2-d grid" files have the format ``<phase> <wavelength>
    <value>`` on each line.

    The phase and wavelength values of the various components don't
    necessarily need to match. (In the most recent salt2 model data,
    they do not all match.) The phase and wavelength values of the
    first model component (in ``m0file``) are taken as the "native"
    sampling of the model, even though these values might require
    interpolation of the other model components.

    """        
    def __init__(self, saltpath='data/salt2-4/'):
        """
        Create a SALT2 model, using the model files found under modeldir,       
        """
        modeldir = '%s/%s' % (srcpath, saltpath)
        if not os.path.isdir(modeldir):
            raise OSError("Error: can't find salt models %s!" % (modeldir))
        
        self._model = dict()
        self._model['M0'] = TimeSeries(modeldir+'/salt2_template_0.dat', modeldir+'/salt2_template_0.clf')
        self._model['M1'] = TimeSeries(modeldir+'/salt2_template_1.dat', modeldir+'/salt2_template_1.clf')

        self._model['V00'] = TimeSeries(modeldir+'/salt2_spec_variance_0.dat', modeldir+'/salt2_spec_variance_0.clf')
        self._model['V11'] = TimeSeries(modeldir+'/salt2_spec_variance_1.dat', modeldir+'/salt2_spec_variance_1.clf')
        self._model['V01'] = TimeSeries(modeldir+'/salt2_spec_covariance_01.dat', modeldir+'/salt2_spec_covariance_01.clf')
        errorscale = modeldir+'/salt2_spec_dispersion_scaling.dat'
        if os.path.exists(errorscale):
            errorscalecache = modeldir+'/salt2_spec_dispersion_scaling.clf'
            self._model['errorscale'] = TimeSeries(errorscale, errorscalecache)
        
        self._wavelengths = self._model['M0'].wavelengths()
        
    def wavelengths(self):
        """Return wavelength coverage array of the model"""
        return self._wavelengths.copy()

    def flux(self, phase, wavelengths=None, x0=1.0, x1=0.0, c=0.0):
        """Return the model (wavelength, flux) spectrum for these parameters"""
        if wavelengths is None:
            wavelengths = self._wavelengths
        f0 = self._model['M0'].spectrum(phase, wavelengths=wavelengths)
        f1 = self._model['M1'].spectrum(phase, wavelengths=wavelengths)
        flux = x0*(f0 + x1*f1)
        return flux

    def error(self, phase, wavelengths=None, x0=1.0, x1=0.0, c=0.0):
        """
        Return flux error spectrum for given parameters
        """
        if wavelengths is None:
            wavelengths = self._wavelengths
            
        v00 = self._model['V00'].spectrum(phase, wavelengths=wavelengths)
        v11 = self._model['V11'].spectrum(phase, wavelengths=wavelengths)
        v01 = self._model['V01'].spectrum(phase, wavelengths=wavelengths)

        sigma = x0 * numpy.sqrt(v00 + x1*x1*v11 + 2*x1*v01)
        sigma *= self._extinction(wavelengths, c)
        ### sigma *= 1e-12   #- Magic constant from SALT2 code
        
        if 'errorscale' in self._model:
            sigma *= self._model['errorscale'].spectrum(phase, wavelengths=wavelengths)
        
        #- Hack adjustment to error (from SALT2 code)
        if phase < -7 or phase > +12:
            xx = numpy.nonzero(wavelengths < 3400)
            sigma[xx] *= 1000
        
        return sigma

    def variance(self, phase, wavelengths=None, x0=1.0, x1=0.0, c=0.0):
        """
        Return model flux variance for given parameters.
        """
        return self.error(phase, wavelengths=wavelengths, x0=x0, x1=x1, c=c)**2

def salt2_fit(time, t0=0, x0=1.0, x1=0.0, c=0.0, saltpath='data/salt2-4/', filt='r'):
    """
    :param time: time in days in source frame
    :param t0: peak time in days
    :param filt: filter
    :param x0/x1/c: salt2 parameters    
    :param saltpath: salt2 model folder
    :return: flux
    """    
    saltcls = Salt2Model(saltpath=saltpath)
    wavelengths = filters.central_wavelengths[filt]
    flux = []
    for t in time:
        f = saltcls.flux(t-t0, wavelengths, x0, x1, c)
        flux.append(f)
    return flux
