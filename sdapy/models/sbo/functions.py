import numpy as np
from sdapy import constants

def shock_fit(times, Me, Re, Ee):
    '''
    shock cooling fit with Piro et al 2020 model
    https://arxiv.org/pdf/2007.08543.pdf
      times: day
      Me   : solar mass
      Re   : solar radius
      Ee   : foe
    '''
    t = times * 86400.
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

def shock_fit_texp(times, Me, Re, Ee, texp):
    '''
    shock cooling fit with Piro et al 2020 model
      with texp as start of shock cooling phase
    ''' 
    return shock_fit(times-texp, Me, Re, Ee)

def shock_fit_texp_tc(times, Me, Re, Ee, texp, tc):
    '''
    shock cooling fit with Piro et al 2020 model
      with tc as end of shock cooling phase
    '''    
    ix = times < tc
    times = times[ix]
    return shock_fit_texp(times, Me, Re, Ee, texp)
