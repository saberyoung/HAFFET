import numpy as np

def dilute(T, a0, a1, a2):
    # dilute facotr
    return a0+a1*(10**4/T)+a2*(10**4/T)**2

def bbody(x,T,R,lambda_cutoff=3000,alpha=1):
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

    return Radiance

def bbody_dilute(x,T,R,lambda_cutoff=3000,alpha=1,a0=0.711,a1=-0.476,a2=0.308):
    ''' diluted blackbody
    '''
    return bbody(x,T1,R1,lambda_cutoff,alpha) * dilute(T, a0, a1, a2)**2

def bbody_double(x,T1,R1,T2,R2,lambda_cutoff=3000,alpha=1):
    ''' two components, each od them is BB
    '''
    return bbody(x,T1,R1,lambda_cutoff,alpha) + \
        bbody(x,T2,R2,lambda_cutoff,alpha)
