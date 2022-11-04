import numpy as np
from sdapy import constants, models

def shock_fit(time, Me, Re, Ee, texp=None):
    ''' shock cooling fit with Piro et al 2020 model (https://arxiv.org/pdf/2007.08543.pdf)
    
    Parameters
    ----------       
    time :  `array`
        Independent values.
    Me    : `float`
         envolop mass (unit: solar mass)
    Re    : `float`
         envolop radius (unit: solar radius)
    Ee     : `float`
         envolop energy (unit: foe) 
    texp    : `float`
        explosion time, time between first light to the peak epoch.
    '''
    times = time
    if texp is not None: times = time - texp
    t = times * constants.day_to_s
    E = Ee * 1e51    
    M = Me * constants.M_sun
    R = Re * constants.R_sun
    
    n = 10
    delta = 1.1
    K = (n-3) * (3-delta) / (4 * np.pi * (n-delta)) # K = 0.119
    kappa = 0.2 # Motivated by the helium-rich composition of this event
                # normally for Ibc we use 0.07
    vt = ((n-5) * (5-delta) / ((n-3) * (3-delta)))**0.5 * (2 * E / M)**0.5
    td = (3 * kappa * K * M / ((n-1) * vt * constants.c))**0.5 # in second    
    
    prefactor = np.pi*(n-1)/(3*(n-5)) * constants.c * R * vt**2 / kappa 
    L1 = prefactor * (td/t)**(4/(n-2))
    L2 = prefactor * np.exp(-0.5 * ((t/td)**2 - 1))
    Ls = np.zeros(len(t))
    ix1 = t < td
    Ls[ix1] = L1[ix1]
    Ls[~ix1] = L2[~ix1]
    return Ls

def shock_arnett_fit(time, Me, Re, Ee, mni, taum, texp=None):
    ''' shock cooling + Arnett (mni, taum) fit
    
    Parameters
    ----------       
    time :  `array`
        Independent values.
    Me    : `float`
         envolop mass (unit: solar mass)
    Re    : `float`
         envolop radius (unit: solar radius)
    Ee     : `float`
         envolop energy (unit: foe) 
    mni    : `float`
        Arnett model parameter, Nickel mass, unit in solar mass.
    taum    : `float`
        Arnett model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.
    texp    : `float`
        explosion time, time between first light to the peak epoch.
    '''
    times = time
    if texp is not None: times = time - texp
    La = models.arnett_tail.Arnett_fit_taum(times, mni, taum)
    Ls = shock_fit(times, Me, Re, Ee)
    return La + Ls

def shock_arnett_mejek_fit(time, Me, Re, Ee, f_ni, mej, ek, texp=None, k_opt=None):
    ''' shock cooling + Arnett (mni, mej, ek) fit
    
    Parameters
    ----------       
    time :  `array`
        Independent values.
    Me    : `float`
         envolop mass (unit: solar mass)
    Re    : `float`
         envolop radius (unit: solar radius)
    Ee     : `float`
         envolop energy (unit: foe) 
    f_ni    : `float`
        fraction of nickel mass
    Mej    : `float`
         ejecta mass (unit: solar mass)
    Ek     : `float`
         kinetic energy (unit: foe)
    texp    : `float`
        explosion time, time between first light to the peak epoch.
    k_opt   : `float`
        diffusion opacity
    '''
    times = time
    if texp is not None: times = time - texp
    if k_opt is None: k_opt = constants.k_opt
    La = models.arnett_tail.Arnett_fit_Mej_Ek(times, f_ni, ek, mej, texp=texp, k_opt=k_opt)
    Ls = shock_fit(times, Me, Re, Ee)
    return La + Ls
