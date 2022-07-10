import numpy as np
from sdapy import constants
from scipy.integrate import simps
import math

def Arnett_fit_taum(times, m_ni, tau_m):
    '''
    decide Arnett bolometric output
    with nickel mass (unit: msun) and taum (unit: day)
    '''
    mni = m_ni * constants.M_sun
    t = times * 86400.    
    t = np.array([t]) if (isinstance(t, np.float64) or isinstance(t, float)) else t    
    tau_Ni = constants.tau_Ni * 86400.
    tau_Co = constants.tau_Co * 86400    
    tau_m  *= 86400.
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

def taum_to_MejEk(taum, taumerr=None):
    '''    
    from taum (unit: day) to the product of kinetic energy (unit: foe) and ejecta mass (unit: solar mass),
    i.e. M$_{ej}^{3/4}$ E$_{kin}^{-1/4}$
    '''
    taum *= 86400.
    M_ejE_K = taum/((constants.k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    unit = (constants.M_sun**3/1e51)**(1./4.)
    
    if taumerr is not None:
        taumerr *= 86400.    
        M_ejE_KE = taumerr/((constants.k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
        return M_ejE_K/unit, M_ejE_KE/unit
    else:        
        return M_ejE_K/unit

def taum_to_Mej_Ek(taum, vej, taumerr=None, vejerr=None):
    '''
    break degenracy with velocity
    from taum (unit: day) to kinetic energy (unit: foe) and ejecta mass (unit: solar mass) seperately
       here vej is the photospheric vm (unit: 10**3 km/s)
        can be obtained via Dessart et al with peak velocity in He 5876 (Ib) or O 7772 (Ic)
        vej: 1e3 km/s       
        vej = sqrt(2*Ekin/Mej)
    '''
    vej *= 1e8 # to cgs units
    if vejerr is not None: vejerr *= 1e8    
    taum *= 86400.
    if taumerr is not None: taumerr *= 86400.
    
    Mej = vej*taum**2*constants.beta*constants.c/constants.k_opt/(2.*constants.gamma)**(.5)
    Ek = vej**3*taum**2*constants.beta*constants.c/constants.k_opt/(8.*constants.gamma)**(.5)
    if taumerr is None or vejerr is None: return Mej/constants.M_sun, Ek/1e51
    
    # error
    M_ejE_K = taum/((constants.k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    M_ejE_KE = taumerr/((constants.k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    Mejerr = Mej*((2*M_ejE_KE/M_ejE_K)**2 + (vejerr/vej)**2)**0.5
    Ekerr = Ek*((2*vejerr/vej)**2 + (Mejerr/Mej)**2)**0.5    
    return Mej/constants.M_sun, Mejerr/constants.M_sun, Ek/1e51, Ekerr/1e51

def taum_to_Mej_Ek_1(taum, vej, taumerr=None, vejerr=None):
    '''
    break degenracy with velocity
    from taum (unit: day) to kinetic energy (unit: foe) and ejecta mass (unit: solar mass) seperately
       here Ekin/Mej = 3/10*vej**2
       ->  vm = vej = sqrt(Ekin/Mej*10./3.)
    '''
    vej *= 1e8
    if vejerr is not None: vejerr *= 1e8    
    taum *= 86400.
    if taumerr is not None: taumerr *= 86400.
    
    Mej = vej*taum**2*constants.beta*constants.c/constants.k_opt/(10./3.*constants.gamma)**(.5)
    Ek = vej**3*taum**2*constants.beta*constants.c/constants.k_opt/((10./3.)**3*constants.gamma)**(.5)
    if taumerr is None or vejerr is None: return Mej/constants.M_sun, Ek/1e51
    
    # error
    M_ejE_K = taum/((constants.k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    M_ejE_KE = taumerr/((constants.k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))
    Mejerr = Mej*((2*M_ejE_KE/M_ejE_K)**2 + (vejerr/vej)**2)**0.5
    Ekerr = Ek*((2*vejerr/vej)**2 + (Mejerr/Mej)**2)**0.5    
    return Mej/constants.M_sun, Mejerr/constants.M_sun, Ek/1e51, Ekerr/1e51

def Mej_Ek_to_taum(Mej, Ek):
    '''    
    from kinetic energy (unit: foe) and ejecta mass (unit: solar mass) to taum (unit: day)
    '''    
    M_ejE_K = ( (Mej*constants.M_sun)**3 / (Ek*1e51) )**(1./4.)
    taum = M_ejE_K * ((constants.k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))           
    return taum / 86400.

def Mej_Ek_to_vej(Mej, Ek):
    '''    
    from kinetic energy (unit: foe) and ejecta mass (unit: solar mass) to photospheric velocity (unit: 10**3 km/s)
    '''
    Mej *= constants.M_sun
    Ek *= 1e51
    vej = np.sqrt(2*Ek/Mej)         
    return vej * 1e-8

def Arnett_fit_Mej_Ek(times, m_ni, Ek, Mej):
    '''
    fit Arnett model on Ek (unit: foe), Mej (unit: solar mass), with vm (unit: 1e3 km/s) as free parameter
    ''' 
    tau_m = Mej_Ek_to_taum(Mej, Ek)
    return Arnett_fit_taum(times, m_ni, tau_m)    

def Arnett_fit_taum_texp(times, m_ni, tau_m, t_exp):
    '''
    fit Arnett model taum with explosion epoch as free parameter
    '''
    t = times - t_exp
    return Arnett_fit_taum(t, m_ni, tau_m)

def Arnett_fit_Mej_Ek_texp(times, m_ni, Ek, Mej, t_exp):
    '''
    fit Arnett model Ek,Mej with explosion epoch as free parameter
    '''
    t = times - t_exp
    tau_m = Mej_Ek_to_taum(Mej, Ek)
    return Arnett_fit_taum(t, m_ni, tau_m)
