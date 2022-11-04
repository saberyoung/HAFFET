import numpy as np
from sdapy import constants
from scipy.integrate import simps
import math

'''
1- Arnett fits on main peak
'''
def Arnett_fit_taum(times, m_ni, taum, texp=None):
    ''' output Arnett bolometric luminosities
    
    Parameters
    ----------
    time :  `array`
        Independent values.
    m_ni    : `float`
        Arnett model parameter, Nickel mass, unit in solar mass.
    taum    : `float`
        Arnett model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.
    texp    : `float`
        explosion time, time between first light to the peak epoch.
    '''
    time = times
    if texp is not None: time = times - texp    
    mni = m_ni * constants.M_sun
    t = time * constants.day_to_s
    t = np.array([t]) if (isinstance(t, np.float64) or isinstance(t, float)) else t    
    tau_Ni = constants.tau_Ni * constants.day_to_s
    tau_Co = constants.tau_Co * constants.day_to_s
    tau_m  = taum * constants.day_to_s
    int_A=np.zeros(len(t)) 
    int_B=np.zeros(len(t)) 
    L_ph=np.zeros(len(t))
    
    x = t/tau_m
    y = tau_m/(2.*tau_Ni)
    s = tau_m*(tau_Co-tau_Ni)/(2.*tau_Co*tau_Ni)
    
    for i in range(len(t)):
        z  = np.arange(100)*x[i]/100.
        Az = 2.*z*np.exp(-2.*z*y+np.square(z))
        Bz = 2.*z*np.exp(-2.*z*y+2.*z*s+np.square(z))
        int_A[i] = simps(Az,z)
        int_B[i] = simps(Bz,z)
        _L_ph = (mni*np.exp(-1.*np.square(x[i])))*((constants.e_Ni-constants.e_Co)*int_A[i]+constants.e_Co*int_B[i])
        if math.isnan(_L_ph): continue
        L_ph[i]=_L_ph
    return L_ph

def taum_to_MejEk(tau_m, taum_err=None, k_opt=None, sntype='Ic'):
    ''' from characteristic time, taum (unit: day) to the product of kinetic energy (unit: foe) and ejecta mass (unit: solar mass),
    i.e. M$_{ej}^{3/4}$ E$_{kin}^{-1/4}$
    
    Parameters
    ----------   
    tau_m    : `float`
        Arnett model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.
    taum_err   : `float`
        Error of taum, if not setted, will return only value, otherwise will return value as well as error.
    k_opt   : `float`
        diffusion opacity
    '''
    if k_opt is None:
        assert sntype in constants.k_opt
        k_opt = constants.k_opt[sntype]    
    taum = tau_m * constants.day_to_s
    M_ejE_K = taum/((k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    unit = (constants.M_sun**3/1e51)**(1./4.)
    
    if taum_err is not None:
        taumerr = taum_err * constants.day_to_s
        M_ejE_KE = taumerr/((k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
        return M_ejE_K/unit, M_ejE_KE/unit
    else:        
        return M_ejE_K/unit

def taum_to_Mej_Ek(tau_m, v_ej, taum_err=None, vej_err=None, k_opt=None, sntype='Ic'):
    ''' break the degenracy of kinetic energy and ejecta mass, with the help of velocity.
        for Thin Shell:  vej = sqrt(2*Ekin/Mej)

    Parameters
    ----------   
    tau_m    : `float`
        Arnett model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.
    v_ej     : `float`
        here vej is the photospheric velocity (unit: 10**3 km/s), which can be related to
        line velocities (He 5876 and O 7772 line velocity at peak epoch for SN Ib and Ic correspondingly)
        via Dessart et al 2014.
    taum_err   : `float`
        Error of taum.
    vej_err    : `float`
        Error of vej.
    k_opt   : `float`
        diffusion opacity
    '''
    if k_opt is None:
        assert sntype in constants.k_opt
        k_opt = constants.k_opt[sntype]   
    vej = v_ej * 1e8 # to cgs units
    if vej_err is not None: vejerr = vej_err * 1e8    
    taum = tau_m * constants.day_to_s
    if taum_err is not None: taumerr = taum_err * constants.day_to_s
    
    Mej = vej*taum**2*constants.beta*constants.c/k_opt/(2.*constants.gamma)**(.5)
    Ek = vej**3*taum**2*constants.beta*constants.c/k_opt/(8.*constants.gamma)**(.5)
    if taum_err is None or vej_err is None: return Mej/constants.M_sun, Ek/1e51
    
    # error
    M_ejE_K = taum/((k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    M_ejE_KE = taumerr/((k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    Mejerr = Mej*((2*M_ejE_KE/M_ejE_K)**2 + (vejerr/vej)**2)**0.5
    Ekerr = Ek*((2*vejerr/vej)**2 + (Mejerr/Mej)**2)**0.5    
    return Mej/constants.M_sun, Mejerr/constants.M_sun, Ek/1e51, Ekerr/1e51

def taum_to_Mej_Ek_1(tau_m, v_ej, taum_err=None, vej_err=None, k_opt=None, sntype='Ic'):
    ''' break the degenracy of kinetic energy and ejecta mass, with the help of velocity.
        for Homologous Expansion: vm = vej = sqrt(Ekin/Mej*10./3.)
    
    Parameters
    ----------   
    tau_m    : `float`
        Arnett model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.
    v_ej     : `float`
        here vej is the photospheric velocity (unit: 10**3 km/s), which can be assumed as line velocities.
    taum_err   : `float`
        Error of taum.
    vej_err    : `float`
        Error of vej.
    k_opt   : `float`
        diffusion opacity
    '''
    if k_opt is None:
        assert sntype in constants.k_opt
        k_opt = constants.k_opt[sntype]   
    vej = v_ej * 1e8
    if vej_err is not None: vejerr = vej_err * 1e8    
    taum = tau_m * constants.day_to_s
    if taum_err is not None: taumerr = taum_err * constants.day_to_s
    
    Mej = vej*taum**2*constants.beta*constants.c/k_opt/(10./3.*constants.gamma)**(.5)
    Ek = vej**3*taum**2*constants.beta*constants.c/k_opt/((10./3.)**3*constants.gamma)**(.5)
    if taum_err is None or vej_err is None: return Mej/constants.M_sun, Ek/1e51
    
    # error
    M_ejE_K = taum/((k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    M_ejE_KE = taumerr/((k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    Mejerr = Mej*((2*M_ejE_KE/M_ejE_K)**2 + (vejerr/vej)**2)**0.5
    Ekerr = Ek*((2*vejerr/vej)**2 + (Mejerr/Mej)**2)**0.5    
    return Mej/constants.M_sun, Mejerr/constants.M_sun, Ek/1e51, Ekerr/1e51

def Mej_Ek_to_taum(Mej, Ek, k_opt=None, sntype='Ic'):
    '''  from kinetic energy (unit: foe) and ejecta mass (unit: solar mass) to taum (unit: day)
    
    Parameters
    ----------   
    Mej    : `float`
         ejecta mass (unit: solar mass)
    Ek     : `float`
         kinetic energy (unit: foe)
    k_opt   : `float`
        diffusion opacity
    '''
    if k_opt is None:
        assert sntype in constants.k_opt
        k_opt = constants.k_opt[sntype]
    M_ejE_K = ( (Mej*constants.M_sun)**3 / (Ek*1e51) )**(1./4.)
    taum = M_ejE_K * ((k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))           
    return taum / constants.day_to_s

def Mej_Ek_to_vej(mej, ek):
    '''    from kinetic energy (unit: foe) and ejecta mass (unit: solar mass) to photospheric velocity (unit: 10**3 km/s)

    Parameters
    ----------   
    mej    : `float`
         ejecta mass (unit: solar mass)
    ek     : `float`
         kinetic energy (unit: foe)
    '''
    Mej = mej * constants.M_sun
    Ek = ek * 1e51
    vej = np.sqrt(2*Ek/Mej)         
    return vej * 1e-8

def Mej_vej_to_Ek(mej, vej):
    '''    from ejecta mass (unit: solar mass) and photospheric velocity (unit: 10**3 km/s) to kinteic energy (unit: foe)
    
    Parameters
    ----------   
    ej    : `float`
         ejecta mass (unit: solar mass)
    vej     : `float`
        here vej is the photospheric velocity (unit: 10**3 km/s), which can be assumed as line velocities.
    '''
    Vej = vej * 1e8
    Mej = mej * constants.M_sun
    return Vej**2*Mej/2/1e51

def Arnett_fit_Mej_Ek(times, f_ni, Ek, Mej, texp=None, k_opt=None, sntype='Ic'):
    '''   fit Arnett model with Ek (unit: foe), Mej (unit: solar mass),
          as free parameter
    
    Parameters
    ----------   
    time :  `array`
        Independent values.
    f_ni    : `float`
        fraction of nickel mass
    Mej    : `float`
        ejecta mass (unit: solar mass)
    Ek     : `float`
         kinetic energy (unit: foe)
    k_opt   : `float`
        diffusion opacity
    texp    : `float`
        explosion time, time between first light to the peak epoch.
    '''
    time = times
    if texp is not None: time = times - texp    
    m_ni = f_ni * Mej
    tau_m = Mej_Ek_to_taum(Mej, Ek, k_opt=k_opt, sntype=sntype)
    return Arnett_fit_taum(time, m_ni, tau_m) 

'''
2- Colbat decay fits on radioactive rail
'''
def tail_fit_t0(t, mni, t0):
    '''   fit radioactive tail with Wygoda 2019 Eq 10, 11 and 12
    (https://arxiv.org/pdf/1711.00969.pdf)

    Parameters
    ----------   
    t :  `array`
        Independent values.
    m_ni    : `float`
        Tail model parameter, Nickel mass, unit in solar mass.
    t0    : `float`
        Tail model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.    
    '''    
    t = np.array([t]) if (isinstance(t, np.float64) or isinstance(t, float)) else t    
    decay_factor = 1-np.e**(-(t0/t)**2)
    tau_Ni = constants.tau_Ni
    tau_Co = constants.tau_Co
    Lgamma = mni*( 6.45*np.e**(-t/tau_Ni) + 1.38*np.e**(-t/tau_Co) )
    Lpos = 4.64*mni*( -np.e**(-t/tau_Ni) + np.e**(-t/tau_Co) )    
    return Lgamma * decay_factor *1e43 + Lpos * 1e41

def t0_to_Mej_Ek(t_0, v_ej, t0_err=None, vej_err=None, k_gamma=None):
    '''   from tail characteristic time, t0 to kinetic energy and ejecta mass seperately
    t0 = sqrt(c*kgamma*Mej**2/Ekin)
    vej = sqrt(2*Ekin/Mej)

    Parameters
    ----------          
    t_0    : `float`
        Tail model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.
    v_ej     : `float`
        here vej is the photospheric velocity (unit: 10**3 km/s), which can be assumed as line velocities.
    t0_err   : `float`
        Error of t0.
    vej_err    : `float`
        Error of vej.
    k_gamma   : `float`
        gamma ray opacity
    '''
    if k_gamma is None: k_gamma = constants.k_gamma
    vej = v_ej * 1e8
    if vej_err is not None: vejerr = vej_err * 1e8
    t0 = t_0 * constants.day_to_s
    if t0_err is not None: t0err = t0_err * constants.day_to_s
    
    Mej = 0.5*vej**2*t0**2/constants.C/k_gamma
    Ekin = 0.25*vej**4*t0**2/constants.C/k_gamma    
    if t0_err is None or vejerr is None: return Mej/constants.M_sun, Ekin/1e51
    
    # error    
    Mejerr = Mej*((2*t0err/t0)**2 + (2*vejerr/vej)**2)**0.5
    Ekerr = Ekin*((2*t0err/t0)**2 + (4*vejerr/vej)**2)**0.5    
    return Mej/constants.M_sun, Mejerr/constants.M_sun, Ekin/1e51, Ekerr/1e51

def Mej_Ek_to_t0(Mej, Ek, k_gamma=None):
    '''    from kinetic energy and ejecta mass to t0
    
    Parameters
    ----------       
    Mej    : `float`
         ejecta mass (unit: solar mass)
    Ek     : `float`
         kinetic energy (unit: foe)  
    k_gamma   : `float`
        gamma ray opacity 
    '''
    if k_gamma is None: k_gamma = constants.k_gamma
    M_ej = Mej * constants.M_sun
    E_k = Ek * 1e51
    t0 = np.sqrt(constants.C*k_gamma*M_ej**2/E_k)
    return t0 / constants.day_to_s

def tail_fit_Mej_Ek(times, f_ni, Ek, Mej, k_gamma=None):
    '''  fit tail model on Ek, Mej, with photospheric velocity, vm as free parameter as well

    Parameters
    ----------   
    time :  `array`
        Independent values.
    f_ni    : `float`
        fraction of nickel mass
    Mej    : `float`
         ejecta mass (unit: solar mass)
    Ek     : `float`
         kinetic energy (unit: foe)       
    k_gamma   : `float`
        gamma ray opacity
    '''
    m_ni = f_ni * Mej    
    t0 = Mej_Ek_to_t0(Mej, Ek, k_gamma=k_gamma)
    return tail_fit_t0(times, m_ni, t0)

'''
3- Arnett + tail fits on full range
'''
def joint_fit_taum_t0(times, mni, taum, t0, ts, texp=None):
    '''  for the full range, *Arnett_fit_taum* fit + *tail_fit_t0* fit, with ts as free parameter
    
    Parameters
    ----------   
    time :  `array`
        Independent values.
    m_ni    : `float`
        Arnett model parameter, Nickel mass, unit in solar mass.
    taum    : `float`
        Arnett model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.
    t0    : `float`
        Tail model parameter, characteristic time, unit in days, decided by ejecta mass and kinetic enegies.
    ts    : `float`
        boundary between Arnett model and tail model
    texp    : `float`
        explosion time, time between first light to the peak epoch.
    '''    
    t = times
    if texp is not None: t = times - texp
    L1 = Arnett_fit_taum(t, mni, taum)
    L2 = tail_fit_t0(t, mni, t0)
    Ls = np.zeros(len(t))
    ix1 = times < ts
    Ls[ix1] = L1[ix1]
    Ls[~ix1] = L2[~ix1]
    return Ls

def joint_fit_Mej_Ek(times, mni, Mej, Ek, ts, texp=None, k_opt=None, sntype='Ic'):
    '''  for the full range, *Arnett_fit_Mej_Ek* fit + *tail_fit_Mej_Ek* fit, with ts as free parameter
    
    Parameters
    ----------   
    time :  `array`
        Independent values.
    m_ni    : `float`
        Arnett model parameter, Nickel mass, unit in solar mass.
    Mej    : `float`
         ejecta mass (unit: solar mass)
    Ek     : `float`
         kinetic energy (unit: foe)   
    ts    : `float`
        boundary between Arnett model and tail model   
    texp    : `float`
        explosion time, time between first light to the peak epoch.
    '''
    t = times
    if texp is not None: t = times - texp
    t0 = Mej_Ek_to_t0(Mej, Ek)
    taum = Mej_Ek_to_taum(Mej, Ek, k_opt=k_opt, sntype=sntype) 
    L1 = Arnett_fit_taum(t, mni, taum)
    L2 = tail_fit_t0(t, mni, t0)
    Ls = np.zeros(len(t))
    ix1 = times < ts
    Ls[ix1] = L1[ix1]
    Ls[~ix1] = L2[~ix1]
    return Ls
