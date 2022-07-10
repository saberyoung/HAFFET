import numpy as np
from sdapy import constants

def tail_fit_t0(t, mni, t0):
    '''
    fit radioactive tail with Wygoda 2019 Eq 10, 11 and 12
    https://arxiv.org/pdf/1711.00969.pdf       
    '''
    t = np.array([t]) if (isinstance(t, np.float64) or isinstance(t, float)) else t    
    decay_factor = 1-np.e**(-(t0/t)**2)
    tau_Ni = constants.tau_Ni
    tau_Co = constants.tau_Co
    Lgamma = mni*( 6.45*np.e**(-t/tau_Ni) + 1.38*np.e**(-t/tau_Co) )
    Lpos = 4.64*mni*( -np.e**(-t/tau_Ni) + np.e**(-t/tau_Co) )    
    return Lgamma * decay_factor *1e43 + Lpos * 1e41
    
def t0_to_Mej_Ek(t0, vej, t0err=None, vejerr=None):
    '''
    from t0 to kinetic energy and ejecta mass seperately
    t0 = sqrt(c*kgamma*Mej**2/Ekin)
    vej = sqrt(2*Ekin/Mej)
    '''
    vej *= 1e8
    if vejerr is not None: vejerr *= 1e8    
    t0 *= 86400.
    if t0err is not None: t0err *= 86400.

    Mej = 0.5*vej**2*t0**2/constants.C/constants.k_gamma
    Ekin = 0.25*vej**4*t0**2/constants.C/constants.k_gamma    
    if t0err is None or vejerr is None: return Mej/constants.M_sun, Ekin/1e51
    
    # error    
    Mejerr = Mej*((2*t0err/t0)**2 + (2*vejerr/vej)**2)**0.5
    Ekerr = Ekin*((2*t0err/t0)**2 + (4*vejerr/vej)**2)**0.5    
    return Mej/constants.M_sun, Mejerr/constants.M_sun, Ekin/1e51, Ekerr/1e51

def Mej_Ek_to_t0(Mej, Ek):
    '''    
    from kinetic energy and ejecta mass to t0
    '''
    Mej *= constants.M_sun
    Ek *= 1e51
    t0 = np.sqrt(constants.C*constants.k_gamma*Mej**2/Ek)
    return t0 / 86400.

def tail_fit_Mej_Ek(times, m_ni, Ek, Mej):
    '''
    fit tail model on Ek, Mej, with vm as free parameter
    ''' 
    t0 = Mej_Ek_to_t0(Mej, Ek)
    return tail_fit_t0(times, m_ni, t0)

def tail_fit_t0_ts(times, m_ni, t0, ts):
    '''
    fit tail model with one more free parameter, ts,
    which is the changing time between Arnett model and tail model
    '''    
    return tail_fit_t0(times[np.where(times>ts)], m_ni, t0)

def tail_fit_Mej_Ek_ts(times, m_ni, Ek, Mej, ts):
    '''
    fit tail model Ek,Mej with ts as free parameter
    '''
    t0 = Mej_Ek_to_t0(Mej, Ek)
    return tail_fit_t0(times[np.where(times>ts)], m_ni, t0)
