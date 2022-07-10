import numpy as np
from sdapy import constants

def joint_fit_taum_t0(times, mni, taum, t0, ts):
    t = times
    L1 = Arnett_fit_taum(t, mni, taum)
    L2 = tail_fit_t0(t, mni, t0)
    Ls = np.zeros(len(t))
    ix1 = times < ts
    Ls[ix1] = L1[ix1]
    Ls[~ix1] = L2[~ix1]
    return Ls

def joint_fit_taum_t0_texp(times, mni, taum, t0, texp, ts):
    t = times - texp
    L1 = Arnett_fit_taum(t, mni, taum)
    L2 = tail_fit_t0(t, mni, t0)
    Ls = np.zeros(len(t))
    ix1 = times < ts
    Ls[ix1] = L1[ix1]
    Ls[~ix1] = L2[~ix1]
    return Ls

def joint_fit_Mej_Ek(times, mni, Mej, Ek, ts):
    t0 = Mej_Ek_to_t0(Mej, Ek)
    taum = Mej_Ek_to_taum(Mej, Ek)
    t = times
    L1 = Arnett_fit_taum(t, mni, taum)
    L2 = tail_fit_t0(t, mni, t0)
    Ls = np.zeros(len(t))
    ix1 = times < ts
    Ls[ix1] = L1[ix1]
    Ls[~ix1] = L2[~ix1]
    return Ls

def joint_fit_Mej_Ek_texp(times, mni, Mej, Ek, texp, ts):
    t0 = Mej_Ek_to_t0(Mej, Ek)
    taum = Mej_Ek_to_taum(Mej, Ek)
    t = times - texp
    L1 = Arnett_fit_taum(t, mni, taum)
    L2 = tail_fit_t0(t, mni, t0)
    Ls = np.zeros(len(t))
    ix1 = times < ts
    Ls[ix1] = L1[ix1]
    Ls[~ix1] = L2[~ix1]
    return Ls    
