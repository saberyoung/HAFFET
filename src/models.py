from scipy.integrate import simps
import math
import numpy as np
from filters import central_wavelengths

class constants:
    # covert things to cgs units
    M_sun=1.989e33
    R_sun=696340e5
    c=3.e10
    tau_Ni=8.8
    tau_Co=111.3
    e_Ni=3.9e10
    e_Co=6.8e9
    k_opt=0.07 # k_opt=0.07 g/cm^2, fixed optical opacity (corresponds to electron scattering) for SE SNe
    beta=13.8  # beta=13.8, constant of integration (Arnett 1982)
    gamma=6./5.
                        
def BC_Lyman(color, phase=1, sntype=1):
    '''
    Lyman analytic bolometric correction for SE SNe: BC_g vs g-r
        sntype = 1: SE SN
        sntype = 2: II
        phase = 1: normal phase
        phase = 2: cooling phase
    '''
    assert sntype in [1,2]
    assert phase in [1,2]
    if phase == 1 and sntype == 1:
        crange = [-.3,1.]
        mbol = 0.054-.195*color-.719*color**2
    elif phase == 1 and sntype == 2:
        crange = [-.2,1.3]
        mbol = 0.053-.089*color-.736*color**2
    else:
        crange = [-.3,.3]
        mbol = -.146+.479*color-2.257*color**2  
    _ = np.logical_and(color >= min(crange), color <= max(crange))
    return mbol, _

def Arnett_fit(times, m_ni, tau_m):
    '''
    decide Arnett bolometric output
    with nickel mass (msun) and taum (day)
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

def Arnett_taum(taum, taumerr=None):
    '''    
    from taum to to the product of kinetic energy and ejecta mass      
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

def Arnett_taumr(Mej, Ek):
    '''    
    from kinetic energy and ejecta mass to taum       
    '''    
    M_ejE_K = ( (Mej*constants.M_sun)**3 / (Ek*1e51) )**(1./4.)
    taum = M_ejE_K * ((constants.k_opt/(constants.beta*constants.c))**0.5)*(constants.gamma**(0.25))           
    return taum # secs

def Arnett_mej_ek(taum, vej, taumerr=None, vejerr=None):
    '''
    break degenracy with velocity 
    so that decide kinetic energy and ejecta mass seperately
       here vej is vm, obtained via Dessart et al with peak velocity in He 5876 (Ib) or O 7772 (Ic)
       vej: 1e3 km/s       
       vej = sqrt(2*Ekin/Mej)
    '''
    vej *= 1e8
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

def Arnett_mej_ek_1(taum, vej, taumerr=None, vejerr=None):
    '''
       Ekin/Mej = 3/10*vej**2
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

def Arnett_fit_withoutv(times, m_ni, Ek, Mej):
    '''
    fit Arnett model together with v
    '''    
    mni = m_ni * constants.M_sun
    t = times * 86400.    
    t = np.array([t]) if (isinstance(t, np.float64) or isinstance(t, float)) else t    
    tau_Ni = constants.tau_Ni * 86400.
    tau_Co = constants.tau_Co * 86400.
    
    int_A=np.zeros(len(t)) 
    int_B=np.zeros(len(t)) 
    L_ph=np.zeros(len(t))
    tau_m = Arnett_taumr(Mej, Ek)
    
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

def Arnett_fit_withoutt(times, m_ni, tau_m, t_exp):
    '''
    fit Arnett model together with explosion epoch
    '''
    t = times - t_exp    
    return Arnett_fit(t, m_ni, tau_m)

def tail_nickel(t, mni, t0):
    # Wygoda 2019 Eq 10, 11 and 12
    # https://arxiv.org/pdf/1711.00969.pdf       
    t = np.array([t]) if (isinstance(t, np.float64) or isinstance(t, float)) else t    
    decay_factor = 1-np.e**(-(t0/t)**2)
    tau_Ni = constants.tau_Ni
    tau_Co = constants.tau_Co
    Lgamma = mni*( 6.45*np.e**(-t/tau_Ni) + 1.38*np.e**(-t/tau_Co) )
    Lpos = 4.64*mni*( -np.e**(-t/tau_Ni) + np.e**(-t/tau_Co) )    
    return Lgamma * decay_factor *1e43 + Lpos * 1e41

def tail_nickel_fitt(t, mni, t0, ts):    
    return tail_nickel(t[np.where(t>ts)], mni, t0)

def arnett_tail(times, mni, taum, t0, texp, td):
    t = times - texp    
    L1 = Arnett_fit(t, mni, taum)
    L2 = tail_nickel(t, mni, t0)
    Ls = np.zeros(len(t))
    ix1 = times < td
    Ls[ix1] = L1[ix1]
    Ls[~ix1] = L2[~ix1]
    return Ls

def Piro_sbofit(times, Me, Re, Ee):
    '''
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

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gaussian(x, depth, avg, std, offset):
    """Evaluate a negative gaussian

    f = -depth * e^(-((x - avg)^2) / (2 * std ** 2)) + offset

    Args:
        x    (ndarray): Values to evaluate the gaussian at
        depth  (float): Amplitude of the gaussian
        avg    (float): Average of the gaussian
        std    (float): Standard deviation of the gaussian
        offset (float): Vertical offset

    Returns:
        The evaluated gaussian
    """

    return -depth * np.exp(-((x - avg) ** 2) / (2 * std ** 2)) + offset

def double_gauss(x, H, A1, x01, sigma1, A2, x02, sigma2): 
    return H + A1 * np.exp(-(x - x01) ** 2 / (2 * sigma1 ** 2)) \
        + A2 * np.exp(-(x - x02) ** 2 / (2 * sigma2 ** 2))

def poly6(x, a, b, c, d, e, f, g):
    return a*x**6+b*x**5+c*x**4+d*x**3+e*x**2+f*x+g

def poly5(x, a, b, c, d, e, f):
    return a*x**5+b*x**4+c*x**3+d*x**2+e*x+f

def poly4(x, a, b, c, d, e):
    return a*x**4+b*x**3+c*x**2+d*x+e

def poly3(x, a, b, c, d):
    return a*x**3+b*x**2+c*x+d

def poly2(x, a, b, c):
    return a*x**2+b*x+c

def linear(x, a, b):
    return a*x+b

def exp(t, a, t0, b, c):
    return a*np.power(t-t0,b) + c  

# analytic model flux form from Bazin et al 2009
def bazin1(time, A, t0, tfall):
    X = np.exp(-(time - t0) / tfall)
    bazinfunc = A**.5 * X
    return bazinfunc

def bazin2(time, A, t0, trise):
    X = 1 / (1 + np.exp(-(time - t0) / trise))
    bazinfunc = A**.5 * X
    return bazinfunc

def bazin(time, A, t0, tfall, trise, C):
    X = bazin1(time, A, t0, tfall)*bazin2(time, A, t0, trise)
    bazinfunc = X + C
    return bazinfunc
    
def villar(time, a, b, t0, t1, tfall, trise):
    """
    Villar et al 2019
    """
    _flux = []
    for _time in time:
        if _time < t1:
            X = 1 / (1 + np.exp(-(_time - t0) / trise))
            _flux.append( (a + b*(_time-t0)) * X )
        else:
            X = np.exp(-(_time - t1) / tfall) / (1 + np.exp(-(_time - t0) / trise))
            _flux.append( (a + b*(t1-t0)) * X )
    return _flux

def powerlaw_post_baseline(times, amplitude=25, t_0=0, alpha_r=2):    
    return amplitude * (times - t_0)**alpha_r

def powerlaw_full(times, amplitude=25, t_0=0, alpha_r=2, c=0):
    f = []
    for t in times:
        if t<t_0:
            f.append( c )
        else:
            f.append( powerlaw_post_baseline(t, amplitude, t_0, alpha_r) + c )
    return f
    
'''
Define likelihood for power law fits
'''
def nll_pls(theta, f, t, f_err):
    return -1*lnlikelihood_pls(theta, f, t, f_err)

def lnlikelihood_pls(theta, f, t, f_err):
    t_0, a, a_prime, alpha_r, f_sigma = theta
    
    pre_exp = np.logical_not(t > t_0)
    model = -np.inf*np.ones_like(f)
    model[pre_exp] = a
    
    time_term = (t[~pre_exp] - t_0)    
    model[~pre_exp] = a + a_prime * (time_term)**alpha_r * 10.**(-alpha_r)
    assert np.all(model > -np.inf),"fewer model values than flux values\n{}\n{}\na{}A'{}alpha{}f_sigma{}".format(model, t,a,a_prime,alpha_r,f_sigma)    
    ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
    return ln_l
    
def lnprior_pls(theta):
    t_0, a, a_prime, alpha_r, f_sigma = theta
    if (a_prime < 0 or 
        f_sigma < 0 or 
        t_0 < -100 or
        t_0 > 0 or 
        alpha_r < 0 or
        alpha_r > 1e8 or
        a < -1e8 or
        a > 1e8):
        return -np.inf
    else:
        return -np.log(a_prime) - np.log(f_sigma) - alpha_r*np.log(10)

def lnposterior_pls(theta, f, t, f_err):
    lnp = lnprior_pls(theta)    
    lnl = lnlikelihood_pls(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp

def nll_pl(theta, f, t, f_err, filters):
    return -1*lnlikelihood_pl(theta, f, t, f_err, filters)

def lnlikelihood_pl(theta, f, t, f_err, filters):
    n_filt = len(np.unique(filters))
    ln_l = 0
    f_num = 0    
    for filt in central_wavelengths:
        # !!! arange filter order depending on central_wavelengths
        if not filt in filters: continue                
        __theta = np.array([theta[0],                            
                            theta[1 + 4*f_num],                            
                            theta[2 + 4*f_num],                            
                            theta[3 + 4*f_num],                            
                            theta[4 + 4*f_num]])
        f_obs = np.where(filters == filt)        
        ln_l += lnlikelihood_pls(__theta, f[f_obs], t[f_obs], f_err[f_obs])
        f_num += 1
    return ln_l

def lnprior_pl(theta, filters):            
    n_filt = len(np.unique(filters))
    ln_p = 0
    f_num = 0
    for filt in central_wavelengths:
        if not filt in filters: continue                     
        __theta = np.array([theta[0],                            
                            theta[1 + 4*f_num],                            
                            theta[2 + 4*f_num],                            
                            theta[3 + 4*f_num],                            
                            theta[4 + 4*f_num]])
        ln_p += lnprior_pls(__theta)
        f_num += 1
    return ln_p

def lnposterior_pl(theta, f, t, f_err, filters):    
    lnp = lnprior_pl(theta, filters)
    if not np.isfinite(lnp): return -np.inf
    lnl = lnlikelihood_pl(theta, f, t, f_err, filters)
    if not np.isfinite(lnl): return -np.inf
    return lnl + lnp
    
'''
Define likelihood for Bazin fit
'''                       
def nll_bzs(theta, f, t, f_err):
    return -1*lnlikelihood_bzs(theta, f, t, f_err)

def lnlikelihood_bzs(theta, f, t, f_err):
    A, t0, tfall, trise, C, f_sigma = theta
        
    model = bazin(t, A, t0, tfall, trise, C)    
    assert np.all(model > -np.inf),"fewer model values than flux values\n{}\n{}\nA{}\nt0{}\ntfall{}\ntrise{}\nC{}\nf_sigma{}".format(model, t, A, t0, tfall, trise, C, f_sigma)      
    ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
    return ln_l

def lnprior_bzs(theta):
    A, t0, tfall, trise, C, f_sigma = theta
    if (A < 0 or
        A < C or
        C < 0 or
        A > 1e8 or
        f_sigma < 0 or 
        t0 < -100 or
        t0 > 100 or 
        trise > 60 or
        trise < 0 or
        tfall < 0 or
        tfall > 120):
        return -np.inf
    else:
        return 0

def lnposterior_bzs(theta, f, t, f_err):
    lnp = lnprior_bzs(theta)
    lnl = lnlikelihood_bzs(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp

def nll_bz(theta, f, t, f_err, filters):
    return -1*lnlikelihood_bz(theta, f, t, f_err, filters)

def lnlikelihood_bz(theta, f, t, f_err, filters):
    n_filt = len(np.unique(filters))
    ln_l = 0
    f_num = 0    
    for filt in central_wavelengths:
        # !!! arange filter order depending on central_wavelengths
        if not filt in filters: continue                
        __theta = np.array([theta[0 + 5*f_num],
                            theta[1 + 5*f_num],                            
                            theta[2 + 5*f_num],                            
                            theta[3 + 5*f_num],                            
                            theta[4 + 5*f_num],
                            theta[5 + 5*f_num]])
        f_obs = np.where(filters == filt)        
        ln_l += lnlikelihood_bzs(__theta, f[f_obs], t[f_obs], f_err[f_obs])
        f_num += 1
    return ln_l

def lnprior_bz(theta, filters):            
    n_filt = len(np.unique(filters))
    ln_p = 0
    f_num = 0
    for filt in central_wavelengths:
        if not filt in filters: continue   
        __theta = np.array([theta[0 + 5*f_num],
                            theta[1 + 5*f_num],                            
                            theta[2 + 5*f_num],                            
                            theta[3 + 5*f_num],                            
                            theta[4 + 5*f_num],
                            theta[5 + 5*f_num]])        
        ln_p += lnprior_bzs(__theta)
        f_num += 1
    return ln_p

def lnposterior_bz(theta, f, t, f_err, filters):        
    lnp = lnprior_bz(theta, filters)
    if not np.isfinite(lnp): return -np.inf
    lnl = lnlikelihood_bz(theta, f, t, f_err, filters)
    if not np.isfinite(lnl): return -np.inf
    return lnl + lnp

'''
Define likelihood for Arnett fit
'''
def nll_arnett(theta, f, t, f_err):
    return -1*lnlikelihood_arnett(theta, f, t, f_err)

def lnlikelihood_arnett(theta, f, t, f_err):
    mni, taum, f_sigma = theta    
    model = Arnett_fit(t, mni, taum)
    assert np.all(model > -np.inf),"fewer model values than flux values\n{}\n{}\nMni{}taum{}f_sigma{}".format(model, t, mni, taum, f_sigma)    
    ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
    return ln_l

def lnprior_arnett(theta):
    mni, taum, f_sigma = theta
    if (mni < 0 or
        mni > 10 or 
        f_sigma < 0 or 
        taum < 0 or
        taum > 100):
        return -np.inf
    else:
        return 0

def lnposterior_arnett(theta, f, t, f_err):
    lnp = lnprior_arnett(theta)    
    lnl = lnlikelihood_arnett(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp

'''
Define likelihood for Arnett fit without v
'''
def nll_arnett_fitv(theta, f, t, f_err):
    return -1*lnlikelihood_arnett_fitv(theta, f, t, f_err)

def lnlikelihood_arnett_fitv(theta, f, t, f_err):
    mni, Ek, Mej, f_sigma = theta    
    model = Arnett_fit_withoutv(t, mni, Ek, Mej)
    assert np.all(model > -np.inf),"fewer model values than flux values\n{}\n{}\nMni{}taum{}f_sigma{}".format(model, t, mni, Ek, Mej, f_sigma)    
    ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
    return ln_l

def lnprior_arnett_fitv(theta):
    mni, Ek, Mej, f_sigma = theta
    vm = np.sqrt(2 * Ek / Mej)
    if (mni < 0 or
        mni > 5 or 
        f_sigma < 0 or 
        Ek < 0.1 or
        Ek > 30 or
        Mej < 3 or
        Mej > 60 or
        vm < 0 or
        vm > 20):
        return -np.inf
    else:
        return 0

def lnposterior_arnett_fitv(theta, f, t, f_err):
    lnp = lnprior_arnett_fitv(theta)    
    lnl = lnlikelihood_arnett_fitv(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp

'''
Define likelihood for Arnett fit without t exp
'''
def nll_arnett_fitt(theta, f, t, f_err):
    return -1*lnlikelihood_arnett_fitt(theta, f, t, f_err)

def lnlikelihood_arnett_fitt(theta, f, t, f_err):
    mni, taum, texp, f_sigma = theta
    model = Arnett_fit_withoutt(t, mni, taum, texp)
    assert np.all(model > -np.inf),"fewer model values than flux values\n{}\n{}\nMni{}taum{}texp{}f_sigma{}".format(model, t, mni, taum, texp, f_sigma)
    ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
    return ln_l

def lnprior_arnett_fitt(theta):
    mni, taum, texp, f_sigma = theta
    if (mni < 0 or
        mni > 5 or 
        f_sigma < 0 or 
        taum < 0 or
        taum > 100 or
        texp < -100 or
        texp > 0):
        return -np.inf
    else:
        return 0

def lnposterior_arnett_fitt(theta, f, t, f_err):
    lnp = lnprior_arnett_fitt(theta)    
    lnl = lnlikelihood_arnett_fitt(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp

'''
Define likelihood for Tail fit
'''
def nll_tail(theta, f, t, f_err):
    return -1*lnlikelihood_tail(theta, f, t, f_err)

def lnlikelihood_tail(theta, f, t, f_err):
    mni, taum, f_sigma = theta    
    model = tail_nickel(t, mni, taum)
    assert np.all(model > -np.inf),"fewer model values than flux values\n{}\n{}\nMni{}taum{}f_sigma{}".format(model, t, mni, taum, f_sigma)
    ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
    return ln_l

def lnprior_tail(theta):
    mni, taum, f_sigma = theta
    if (mni < 0 or
        mni > 10 or 
        f_sigma < 0 or
        taum > 1e3 or 
        taum < 0):
        return -np.inf
    else:
        return 0

def lnposterior_tail(theta, f, t, f_err):
    lnp = lnprior_tail(theta)    
    lnl = lnlikelihood_tail(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp

'''
Define likelihood for Tail fit with t start as well
'''
def nll_tail_fitt(theta, f, t, f_err):
    return -1*lnlikelihood_tail_fitt(theta, f, t, f_err)

def lnlikelihood_tail_fitt(theta, f, t, f_err):
    mni, taum, ts, f_sigma = theta
    __ = np.where( t>=ts )    
    model = tail_nickel(t[__], mni, taum)
    f, f_err = f[__], f_err[__]
    assert np.all(model > -np.inf),"fewer model values than flux values\n{}\n{}\nMni{}taum{}ts{}f_sigma{}".format(model, t, mni, taum, ts, f_sigma)
    ln_l = -0.5*np.sum((f - model)**2 / ((f_sigma*f_err)**2)) - np.sum(np.log(f_sigma*f_err)) - 0.5*len(model)*np.log(2*np.pi)
    return ln_l

def lnprior_tail_fitt(theta):
    mni, taum, ts, f_sigma = theta
    if (mni < 0 or
        mni > 10 or
        ts < 40 or
        ts > 100 or 
        f_sigma < 0 or
        taum > 1e3 or 
        taum < 0):
        return -np.inf
    else:
        return 0

def lnposterior_tail_fitt(theta, f, t, f_err):
    lnp = lnprior_tail_fitt(theta)    
    lnl = lnlikelihood_tail_fitt(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp

'''
Define likelihood for Arnett + Tail fit
'''
def nll_arnett_tail(theta, f, t, f_err):
    return -1*lnlikelihood_arnett_tail(theta, f, t, f_err)

def lnlikelihood_arnett_tail(theta, f, t, f_err):
    mni, taum, t0, texp, td = theta    
    model = arnett_tail(t, mni, taum, t0, texp, td)
    chi2_term = -1/2*np.sum((f - model)**2/f_err**2)
    error_term = np.sum(np.log(1/np.sqrt(2*np.pi*f_err**2)))    
    ln_l = chi2_term + error_term
    return ln_l

def lnprior_arnett_tail(theta):
    mni, taum, t0, texp, td = theta
    if ((0 < mni < 10) and (0 < taum < 1e3) and (0 < t0 < 1e3) and (-100 < texp < 0) and (20 < td < 100)):
        return 0.0
    return -np.inf

def lnposterior_arnett_tail(theta, f, t, f_err):
    lnp = lnprior_arnett_tail(theta)    
    lnl = lnlikelihood_arnett_tail(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp

'''
Define likelihood for SBO fit using the Piro model
https://arxiv.org/pdf/2007.08543.pdf
'''
def nll_arnett_sbo(theta, f, t, f_err):
    return -1*lnlikelihood_arnett_sbo(theta, f, t, f_err)

def lnlikelihood_arnett_sbo(theta, f, t, f_err):
    Me, Re, Ee, texp = theta
    t_ = t - texp
    model = Piro_sbofit(t_, Me, Re, Ee)    
    chi2_term = -1/2*np.sum((f - model)**2/f_err**2)
    error_term = np.sum(np.log(1/np.sqrt(2*np.pi*f_err**2)))    
    ln_l = chi2_term + error_term
    return ln_l

def lnprior_arnett_sbo(theta):
    Me, Re, Ee, texp = theta
    if (Me < 0 or
        Me > 1 or
        Re < 100 or
        Re > 10000 or 
        Ee < 1e-3 or
        Ee > 1 or 
        texp < -100 or
        texp > 0):
        return -np.inf
    else:
        return 0

def lnposterior_arnett_sbo(theta, f, t, f_err):
    lnp = lnprior_arnett_sbo(theta)    
    lnl = lnlikelihood_arnett_sbo(theta, f, t, f_err)
    if not np.isfinite(lnl):
        return -np.inf
    if not np.isfinite(lnp):
        return -np.inf
    return lnl + lnp
