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

    return Radiance / 1e40 #* dilute(T, a0, a1, a2)**2

def redchisqg(ydata,ymod,deg=2,sd=None):  
    if sd is None:  chisq=np.sum((ydata-ymod)**2)  
    else:  chisq=np.sum( ((ydata-ymod)/sd)**2 ) 
    nu=ydata.size-1-deg  
    return chisq/nu

def BC_Lyman(color, colortype='g-r', phase='normal', sntype='Ic'):
    '''
    Lyman analytic bolometric correction for SE/II SNe
    https://ui.adsabs.harvard.edu/abs/2014MNRAS.437.3848L/abstract
    table 2, 3, 4
    '''
    if not sntype in ['Ib','Ic','IIb','Ic-BL','SLSN','II']:  return None, None
    if not phase in ['normal','cool']:  return None, None
    if not colortype in ['g-r','g-i','B-V','B-R','B-I','V-R','V-I']:  return None, None
    if phase == 'normal' and sntype == ['Ib','Ic','IIb','Ic-BL','SLSN']:
        if colortype == 'g-r':
            crange = [-.3,1.]
            mbol = .054-.195*color-.719*color**2
        elif colortype == 'g-i':
            crange = [-.8,1.1]            
            mbol = -.029-.404*color-.230*color**2
        elif colortype == 'B-V':
            crange = [0.,1.3]
            mbol = -.083-.139*color-.691*color**2
        elif colortype == 'B-R':
            crange = [.1,2.]
            mbol = -.029-.302*color-.224*color**2
        elif colortype == 'B-I':
            crange = [-.4,2.3]
            mbol = -.055-.240*color-.154*color**2
        elif colortype == 'V-R':
            crange = [-.2,.7]
            mbol = .197-.183*color-.419*color**2
        elif colortype == 'V-I':
            crange = [-.7,1.1]
            mbol = .213-.203*color-.079*color**2       
    elif phase == 'normal' and sntype == 'II':
        if colortype == 'g-r':
            crange = [-.2,1.3]
            mbol = .053-.089*color-.736*color**2
        elif colortype == 'g-i':
            crange = [-.5,1.4]
            mbol = -0.007-.359*color-.336*color**2
        elif colortype == 'B-V':
            crange = [0,1.6]
            mbol = -.139-.013*color-.649*color**2
        elif colortype == 'B-R':
            crange = [.1,2.5]
            mbol = .004-.303*color-.213*color**2
        elif colortype == 'B-I':
            crange = [0,2.8]
            mbol = .004-.297*color-.149*color**2
        elif colortype == 'V-R':
            crange = [0,.9]
            mbol = .073-.902*color-1.796*color**2
        elif colortype == 'V-I':
            crange = [0,1.2]
            mbol = .057+.708*color-.912*color**2            
    else:
        if colortype == 'g-r':
            crange = [-.3,.3]
            mbol = -.146+.479*color-2.257*color**2
        elif colortype == 'g-i':
            crange = [-.7,.1]
            mbol = -.158-.459*color-1.599*color**2
        elif colortype == 'B-B':
            crange = [-.2,.5]
            mbol = -.393+.786*color-2.124*color**2
        elif colortype == 'B-R':
            crange = [-.2,.8]
            mbol = -.463+.790*color-1.034*color**2
        elif colortype == 'B-I':
            crange = [-.2,.8]
            mbol = -.473+.830*color-1.064*color**2
        elif colortype == 'V-R':
            crange = [0,.4]
            mbol = -.719+4.093*color-6.419*color**2
        elif colortype == 'V-I':
            crange = [0,.4]
            mbol = -.610+2.244*color-2.107*color**2 
    _ = np.logical_and(color >= min(crange), color <= max(crange))
    return mbol, _

def dessart_vej_to_vm(vpeak, vpeake, sntype, verbose):
    '''
    convert vej at peak to photospheirc vm, unit: 10*3 km/s
     following Dessart 16, convert photospheric velocity (vej) to characteristic velicity (vm)
        vej is O 7772 line velocity for SNe Ic, and He 5886 for Ib, at the peak epoch.
        https://academic.oup.com/mnras/article/458/2/1618/2589109
    '''
    if sntype in ['SN Ib', 'SN IIb']:
        vm, dvm = (vpeak-2.64) / 0.765, vpeake / 0.765
    elif sntype in ['SN Ic', 'SN Ic-BL']:
        vm, dvm = (vpeak-2.99) / 0.443, vpeake / 0.443
    else:
        if verbose: print ('check if Dessart et al suitable for your sn type')
        vm, dvm = None, None
    return vm, dvm

# pseudo equivalent width
def calc_pew(w,f,fc): return scipy.integrate.trapz(1-f/fc,w)

# wavelength (A) to velocity (1e3 km/s)
def calc_vel(x,x0): return -(x-x0)/x0*299792.458/1000.

# velocity (1e3 km/s) to wavelength (A)
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

def get_samples_nstep(h5_file, thin_by=1):    
    if not os.path.exists(h5_file):  return    
    reader = emcee.backends.HDFBackend(h5_file)
    nsteps = thin_by*np.shape(reader.get_chain())[0]
    return nsteps

def get_samples_scipy(h5_file):    
    if not os.path.exists(h5_file):  return    
    reader = h5py.File(h5_file, 'r')    
    samples = reader['samples']
    lnpost = reader['lnprob']    
    try: return samples.value, lnpost.value
    except: return samples, lnpost

def get_numpy(v):
    try:    return v.to_numpy() # pandas dataframe
    except: return np.array(v)  # numpy array

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

def penc_to_errors(mean, perc_low, perc_up):        
    return abs(max(mean-perc_low, perc_up-mean))
