import numpy as np
import scipy.special as ss
from collections import namedtuple
from scipy.interpolate import interp1d
from scipy.integrate import quad, cumtrapz
from inspect import isfunction
from sdapy import constants, models

# assume a constant optical opacity for SLSN?
k_opt_slsn = 0.2

def Diffusion(time, luminosity, kappa, kappa_gamma, mej, vej):
    """ 
    inspired from https://github.com/guillochon/MOSFiT/blob/master/mosfit/modules/transforms/diffusion.py        
    """
    N_INT_TIMES = 100
    MIN_LOG_SPACING = -3
    DIFF_CONST = 2.0 * constants.M_sun / (13.7 * constants.c * constants.km_cgs)
    TRAP_CONST = 3.0 * constants.M_sun / (4*np.pi * constants.km_cgs ** 2)    
    
    _tau_diff = np.sqrt(DIFF_CONST * kappa * mej / vej) / constants.day_to_s
    _trap_coeff = (
        TRAP_CONST * kappa_gamma * mej / vej ** 2) / constants.day_to_s ** 2
    td2, A = _tau_diff ** 2, _trap_coeff
    
    new_lums = np.zeros_like(time)
    min_te = min(time)
    tb = max(0.0, min_te)
    linterp = interp1d(time, luminosity, copy=False, assume_sorted=True)
    
    uniq_times = np.unique(time[time >= tb])    
    lu = len(uniq_times)
    
    num = int(round(N_INT_TIMES / 2.0))
    lsp = np.logspace(
        np.log10(_tau_diff / time[-1]) +
        MIN_LOG_SPACING, 0, num)
    xm = np.unique(np.concatenate((lsp, 1 - lsp)))
    
    int_times = np.clip(
        tb + (uniq_times.reshape(lu, 1) - tb) * xm, tb, time[-1])

    int_te2s = int_times[:, -1] ** 2
    int_lums = linterp(int_times)
    int_args = int_lums * int_times * np.exp(
        (int_times ** 2 - int_te2s.reshape(lu, 1)) / td2)
    int_args[np.isnan(int_args)] = 0.0
    
    uniq_lums = np.trapz(int_args, int_times)
    uniq_lums *= -2.0 * np.expm1(-A / int_te2s) / td2
    
    new_lums = uniq_lums[np.searchsorted(uniq_times, time)]
    
    return new_lums
    
def basic_magnetar(time, p0, bp, mass_ns, theta_pb):
    """
    https://github.com/guillochon/MOSFiT/blob/master/mosfit/modules/engines/magnetar.py    
    """    
    Ep = 2.6e52 * (mass_ns / 1.4) ** (3. / 2.) * p0 ** (-2)
    # ^ E_rot = 1/2 I (2pi/P)^2, unit = erg

    tp = 1.3e5 * bp ** (-2) * p0 ** 2 * (
        mass_ns / 1.4) ** (3. / 2.) * (np.sin(theta_pb)) ** (-2)
    # ^ tau_spindown = P/(dP/dt), unit = s
    # Magnetic dipole: power = 2/(3c^3)*(R^3 Bsin(theta))^2 * (2pi/P)^4
    # Set equal to -d/dt(E_rot) to derive tau

    # ^ From Ostriker and Gunn 1971 eq 4    
    return 2 * Ep / tp / (1. + 2 * time / tp) ** 2

def basic_magnetar_powered_bolometric(times, p0, bp, mass_ns,
                theta_pb, mej, ek, texp=None, k_opt=None, k_gamma=None):    
    time = times
    if texp is not None: time = times - texp
    if k_opt is None: k_opt = k_opt_slsn
    if k_gamma is None: k_gamma = constants.k_gamma    
    lbol = basic_magnetar(time=time * constants.day_to_s, p0=p0, bp=bp, mass_ns=mass_ns, theta_pb=theta_pb)
    vej = models.arnett_tail.Mej_Ek_to_vej(mej, ek)
    return Diffusion(time, lbol, k_opt, k_gamma, mej, vej*1e3)

'''
import matplotlib.pyplot as plt
times = np.arange(-30, 80, 1)
p0 = 5
bp = 1
mass_ns = 1.4
theta_pb = 1
mej = 1
ek = 1
l = basic_magnetar_powered_bolometric(times, p0, bp, mass_ns,
            theta_pb, mej, ek, texp=-20, k_opt=None, k_gamma=None)
plt.plot(times, l, 'ko')
plt.savefig('test.png')
'''
