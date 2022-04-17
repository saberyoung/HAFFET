import numpy as np
import emcee, os
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, NullLocator
from matplotlib.colors import LinearSegmentedColormap, colorConverter
from matplotlib.ticker import ScalarFormatter
        
def Lbol_to_Mbol(L):
    # bolometric luminosity to bolometric mag
    return -2.5*( np.log10(L) - (71.21+17.5)*0.4 )

def Mbol_to_Lbol(M):
    # bolometric mag to bolometric luminosity
    return 10**(-0.4*M + (71.21+17.5)*0.4)

def dilute(T, a0, a1, a2):
    # dilute facotr
    return a0+a1*(10**4/T)+a2*(10**4/T)**2

def bbody(x,T,R,lambda_cutoff=3000,alpha=1,
                   a0=0.711,a1=-0.476,a2=0.308):
    '''
    Calculate the blackbody radiance for a set
    of wavelengths given a temperature and radiance.
    Modified in the UV
    Parameters
    ---------------
    lam: Reference wavelengths in Angstroms
    T:   Temperature in Kelvin
    R:   Radius in cm
    Output
    ---------------
    Spectral radiance in units of erg/s/Angstrom
    (calculation and constants checked by Sebastian Gomez)
    '''

    T *= 1000
    R *= 1e15

    # Planck Constant in cm^2 * g / s
    h = 6.62607E-27
    # Speed of light in cm/s
    c = 2.99792458E10

    # Convert wavelength to cm
    lam_cm = x * 1E-8

    # Boltzmann Constant in cm^2 * g / s^2 / K
    k_B = 1.38064852E-16

    # Calculate Radiance B_lam, in units of (erg / s) / cm ^ 2 / cm
    exponential = (h * c) / (lam_cm * k_B * T)
    B_lam = ((2 * np.pi * h * c ** 2) / (lam_cm ** 5)) / (np.exp(exponential) - 1)
    B_lam[x <= lambda_cutoff] *= (x[x <= lambda_cutoff]/lambda_cutoff)**alpha

    # Multiply by the surface area
    A = 4*np.pi*R**2

    # Output radiance in units of (erg / s) / Angstrom
    Radiance = B_lam * A / 1E8

    return Radiance / 1e40 * dilute(T, a0, a1, a2)**2

def redchisqg(ydata,ymod,deg=2,sd=None):  
    if sd is None:  chisq=np.sum((ydata-ymod)**2)  
    else:  chisq=np.sum( ((ydata-ymod)/sd)**2 ) 
    nu=ydata.size-1-deg  
    return chisq/nu

# pseudo equivalent width
def calc_pew(w,f,fc): return scipy.integrate.trapz(1-f/fc,w)

# minima to velocity
def calc_vel(x,x0): return -(x-x0)/x0*299792.458/1000.

# velocity to wavelength
def calc_wav(v,x0): return -v*1000./299792.458*x0+x0

def get_samples_mc(h5_file):    
    if not os.path.exists(h5_file):  return    
    reader = emcee.backends.HDFBackend(h5_file)
    tau = reader.get_autocorr_time(tol=0)   
    burnin = int(5*np.max(tau))
    samples = reader.get_chain(discard=burnin, thin=np.max(int(np.max(tau)), 0), flat=True)
    lnpost = reader.get_log_prob(discard=burnin, thin=np.max([int(np.max(tau)), 1]), flat=True)    
    try: return samples.value, lnpost.value
    except: return samples, lnpost
    
def get_samples_scipy(h5_file):    
    if not os.path.exists(h5_file):  return    
    reader = h5py.File(h5_file, 'r')    
    samples = reader['samples']
    lnpost = reader['lnprob']    
    try: return samples.value, lnpost.value
    except: return samples, lnpost

def get_numpy(v):
    try:    return v.to_numpy() # pandas dataframe
    except: return v            # numpy array

def is_seq(o):
    """Check if the object is a sequence

    Parameters
    ----------
    o : any object
           The object to check
        
    Returns
    -------
        is_seq : bool, scalar
           True if *o* is a sequence, False otherwise
    """
    return hasattr(o, "__len__")
