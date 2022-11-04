'''
Follow superbol (https://github.com/mnicholl/superbol/blob/master/superbol.py)

SDSS filters and AB mags:
These effective wavelengths for SDSS filters are from Fukugita et al. (1996, AJ, 111, 1748) and are
the wavelength weighted averages (effective wavelengths in their Table 2a, first row)

For Swift UVOT: S=UVW2, D=UVM2, A=UVW1
For GALEX: F=FUV, N=NUV
For NEOWISE: W=W1, Q=W2
For Spitzer: C=Ch1, E=Ch2

 The below zeropoints are needed to convert magnitudes to fluxes
For AB mags,
    m(AB) = -2.5 log(f_nu) - 48.60. f_nu is in units of ergs/s/cm2/Hz
    m(AB) = -2.5 log(f_nu) - 23.90. f_nu is in units of micor Jansky
Therefore, AB magnitudes are directly related to a physical flux.
 Working through the conversion to ergs/s/cm2/Angs, gives
 f_lam = 0.1089/(lambda_eff^2)  where lambda_eff is the effective wavelength of the filter in angstroms
 Note then that the AB flux zeropoint is defined ONLY by the choice of effective wavelength of the bandpass

 However, not all bands here are AB mag, so for consistency across all filters the zeropoints are stored in the following dictionary

 Matt originally listed the following from  Paul Martini's page : http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
 That is not an original source, for AB mags it simply uses the f_lam =0.1089/(lambda_eff^2) relation, and the effective wavelengths from Fukugita et al.

 ugriz and GALEX NUV/FUV are in AB mag system, UBVRI are Johnson-Cousins in Vega mag, JHK are Glass system Vega mags, Swift UVOT SDA and WISE WQ are in Vega mag system

The values for UBVRIJHK are for the Johnson-Cousins-Glass system and are taken directly from Bessell et al. 1998, A&A, 333, 231 (Paul Martini's page lists these verbatim)
Note that these Bessell et al. (1998) values were calculated not from the spectrum of Vega itself, but from a Kurucz model atmosphere of an AOV star.
GALEX effective wavelengths from here: http://galex.stsci.edu/gr6/?page=faq
WISE data from http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=WISE&asttype=
ATLAS values taken from Tonry et al 2018
'''

# order: FSDNAuUBgcVwrRoiIzyJHKWQ
# central wavelengths (in Angs)
central_wavelengths = {
    'u': 3560,  'g': 4830, 'r': 6260, 'i': 7670,
    'z': 8890, 'y': 9600, 'w':5985, 
    'U': 3600,  'B': 4380, 'V': 5450, 'R': 6410,
    'G': 6730, 'I': 7980, 'J': 12200, 'H': 16300,
    'K': 21900, 'S': 2030, 'D': 2231, 'A': 2634,
    'F': 1516, 'N': 2267, 'o': 6790, 'c': 5330,
    'W': 33526, 'Q': 46028, 'C': 35075, 'E': 44366,
}

# If UVOT bands are in AB, need to convert to Vega
ab_to_vega={    
    'V':-0.01,'B':-0.13,'U':1.02,'A':1.51,'D':1.69,'S':1.73,'White':0.80,    
    'J':0.929,'H':1.394,'K':1.859,'W':2.699,'Q':3.339,'C':2.779,'E':3.264,
}

#All values in 1e-11 erg/s/cm2/Angs
zp = {'u': 859.5, 'g': 466.9, 'r': 278.0, 'i': 185.2, 'z': 137.8, 'y': 118.2, 'w': 245.7, 
      'U': 417.5, 'B': 632.0, 'V': 363.1, 'R': 217.7, 'G': 240.0, 'I': 112.6, 'J': 31.47, 'H': 11.38,
      'K': 3.961, 'S': 536.2, 'D': 463.7, 'A': 412.3, 'F': 4801., 'N': 2119., 'o': 236.2, 'c': 383.3,
      'W': 0.818, 'Q': 0.242, 'C': 0.676, 'E': 0.273}

#Filter widths (in Angs)
filter_width = {'u': 458,  'g': 928, 'r': 812, 'i': 894,  'z': 1183, 'y': 628, 'w': 2560,
                'U': 485,  'B': 831, 'V': 827, 'R': 1389, 'G': 4203, 'I': 899, 'J': 1759, 'H': 2041,
                'K': 2800, 'S': 671, 'D': 446, 'A': 821,  'F': 268,  'N': 732, 'o': 2580, 'c': 2280,
                'W': 6626, 'Q': 10422, 'C': 7432, 'E': 10097}

#Extinction coefficients in A_lam / E(B-V). Uses York Extinction Solver (http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/YorkExtinctionSolver/coefficients.cgi)
Rf = {'u': 4.786,  'g': 3.587, 'r': 2.471, 'i': 1.798,  'z': 1.403, 'y': 1.228, 'w':2.762,
      'U': 4.744,  'B': 4.016, 'V': 3.011, 'R': 2.386, 'G': 2.216, 'I': 1.684, 'J': 0.813, 'H': 0.516,
      'K': 0.337, 'S': 8.795, 'D': 9.270, 'A': 6.432,  'F': 8.054,  'N': 8.969, 'o': 2.185, 'c': 3.111,
      'W': 0.190, 'Q': 0.127, 'C': 0.237, 'E': 0.137}

# flux shift, FSDNAuUBgcVwrRoiIzyJHKWQ, SDAuUBgcVrRoiI
#ys = {'g':0, 'c':20, 'V':30, 'r':40, 'R':50, 'o': 60, 'i':80, 'I':70,
#      'u':-40, 'U': -30, 'B':-20,'S': -100, 'D': -80, 'A': -60, 'z':90}
#ys = {'g':0, 'c':40, 'r':80, 'o': 120, 'i':160}
ys = {'F':-50,'S':-40,'D':-30,'N':-20,'A':-10,'u':0,'U':10,'B':20,
      'g':30,'c':40,'V':50,'w':60,'r':70,'R':80,'o':90,'i':100,
      'I':110,'z':120,'y':130,'J':140,'H':150,'K':160,'W':170,'Q':180,'C':190,'E':200}

# E(X-Y) to E(B-V)
# https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2188Y/abstract table 1
toebv = {'F-N': 0.154, 'N-u':2.460, 'u-g': 0.984,  'g-r': 0.901, 'r-i': 0.557,
         'i-z': 0.496, 'z-J':0.554, 'J-H': 0.279,  'H-K': 0.188, 'K-W': 0.149, 'W-Q': 0.075}

# matplotlib pyplot style
## !!! don't put alpha for them since it was already defined in the def_params.txt file

## for data points
PROP1 = {
    "r":dict(color='red', marker='o', ls='none', fillstyle='none', label='r'),
    "g":dict(color='green', marker='o', ls='none', fillstyle='none', label='g'),
    "R":dict(color='gold', marker='o', ls='none', fillstyle='none', label='R'),
    "B":dict(color='skyblue', marker='o', ls='none', fillstyle='none', label='B'),
    "i":dict(color='k', marker='o', ls='none', fillstyle='none', label='i'),
    "z":dict(color='blue', marker='o', ls='none', fillstyle='none', label='z'),
    "u":dict(color='brown', marker='o', ls='none', fillstyle='none', label='u'),
    "c":dict(color='cyan', marker='o', ls='none', fillstyle='none', label='c'),
    "o":dict(color='orange', marker='o', ls='none', fillstyle='none', label='o'),
    "S":dict(color='purple', marker='o', ls='none', fillstyle='none', label='UVW2'),
    "D":dict(color='salmon', marker='o', ls='none', fillstyle='none', label='UVM2'),
    "A":dict(color='pink', marker='o', ls='none', fillstyle='none', label='UVM1'),
    "U":dict(color='grey', marker='o', ls='none', fillstyle='none', label='U'),
    "V":dict(color='magenta', marker='o', ls='none', fillstyle='none', label='V'),
    "H":dict(color='lime', marker='o', ls='none', fillstyle='none', label='H'),
    "J":dict(color='navy', marker='o', ls='none', fillstyle='none', label='J'),
    "K":dict(color='peru', marker='o', ls='none', fillstyle='none', label='K'),
    "Q":dict(color='violet', marker='o', ls='none', fillstyle='none', label='W2'),
    "W":dict(color='gold', marker='o', ls='none', fillstyle='none', label='W1'),
    "C":dict(color='indigo', marker='o', ls='none', fillstyle='none', label='CH1'),
    "E":dict(color='tan', marker='o', ls='none', fillstyle='none', label='CH2'),
}

## for data upper limits
PROP1l = {
    "r":dict(color='red', marker='v', ls='none', fillstyle='none',),
    "g":dict(color='green', marker='v', ls='none', fillstyle='none',),
    "R":dict(color='gold', marker='v', ls='none', fillstyle='none',),
    "B":dict(color='skyblue', marker='v', ls='none', fillstyle='none',),
    "i":dict(color='k', marker='v', ls='none', fillstyle='none',),
    "z":dict(color='blue', marker='v', ls='none', fillstyle='none',),
    "u":dict(color='brown', marker='v', ls='none', fillstyle='none',),
    "c":dict(color='cyan', marker='v', ls='none', fillstyle='none',),
    "o":dict(color='orange', marker='v', ls='none', fillstyle='none',),
    "S":dict(color='purple', marker='v', ls='none', fillstyle='none'),
    "D":dict(color='salmon', marker='v', ls='none', fillstyle='none'),
    "A":dict(color='pink', marker='v', ls='none', fillstyle='none'),
    "U":dict(color='grey', marker='v', ls='none', fillstyle='none'),
    "V":dict(color='magenta', marker='v', ls='none', fillstyle='none'),
    "H":dict(color='lime', marker='v', ls='none', fillstyle='none'),
    "J":dict(color='navy', marker='v', ls='none', fillstyle='none'),
    "K":dict(color='peru', marker='v', ls='none', fillstyle='none'),
    "Q":dict(color='violet', marker='v', ls='none', fillstyle='none'),
    "W":dict(color='gold', marker='v', ls='none', fillstyle='none'),
    "C":dict(color='indigo', marker='v', ls='none', fillstyle='none'),
    "E":dict(color='tan', marker='v', ls='none', fillstyle='none'),
}

## for multiband fittings sampling curves
PROP2 = {
    "r":dict(color='red', ls='-', zorder=9),
    "g":dict(color='green', ls='-', zorder=9),
    "R":dict(color='gold', ls='-', zorder=9),
    "B":dict(color='skyblue', ls='-', zorder=9),
    "i":dict(color='k', ls='-', zorder=9),
    "z":dict(color='blue', ls='-', zorder=9),
    "u":dict(color='brown', ls='-', zorder=9),
    "c":dict(color='cyan', ls='-', zorder=9),
    "o":dict(color='orange', ls='-', zorder=9),
    "S":dict(color='purple', ls='-', zorder=9),
    "D":dict(color='salmon', ls='-', zorder=9),
    "A":dict(color='pink', ls='-', zorder=9),
    "U":dict(color='grey', ls='-', zorder=9),
    "V":dict(color='magenta', ls='-', zorder=9),
    "H":dict(color='lime', ls='-', zorder=9),
    "J":dict(color='navy', ls='-', zorder=9),
    "K":dict(color='peru', ls='-', zorder=9),
    "Q":dict(color='violet', ls='-', zorder=9),
    "W":dict(color='gold', ls='-', zorder=9),
    "C":dict(color='indigo', ls='-', zorder=9),
    "E":dict(color='tan', ls='-', zorder=9),
}

## for Gaussian process sampling curve
## !!! don't put label and alpha here
PROP3 = {
    "r":dict(color='red', ls=':', zorder=7),
    "g":dict(color='green', ls=':', zorder=7),
    "R":dict(color='gold', ls=':', zorder=7),
    "B":dict(color='skyblue', ls=':', zorder=7),
    "i":dict(color='k', ls=':', zorder=7),
    "z":dict(color='blue', ls=':', zorder=7),
    "u":dict(color='brown', ls=':', zorder=7),
    "c":dict(color='cyan', ls=':', zorder=7),
    "o":dict(color='orange', ls=':', zorder=7),
    "S":dict(color='purple', ls=':', zorder=7),
    "D":dict(color='salmon', ls=':', zorder=7),
    "A":dict(color='pink', ls=':', zorder=7),
    "U":dict(color='grey', ls=':', zorder=7),
    "V":dict(color='magenta', ls=':', zorder=7),
    "H":dict(color='lime', ls=':', zorder=7),
    "J":dict(color='navy', ls=':', zorder=7),
    "K":dict(color='peru', ls=':', zorder=7),
    "Q":dict(color='violet', ls=':', zorder=7),
    "W":dict(color='gold', ls=':', zorder=7),
    "C":dict(color='indigo', ls=':', zorder=7),
    "E":dict(color='tan', ls=':', zorder=7),
}

## for fitting peak axvline
PROP4 = {
    "r":dict(color='red', ls='--', zorder=8, label='r peak'),
    "g":dict(color='green', ls='--', zorder=8, label=None),
    "R":dict(color='gold', ls='--', zorder=8, label='R peak'),
    "B":dict(color='skyblue', ls='--', zorder=8, label=None),
    "i":dict(color='k', ls='--', zorder=8, label=None),
    "z":dict(color='blue', ls='--', zorder=8, label=None),
    "u":dict(color='brown', ls='--', zorder=8, label=None),
    "c":dict(color='cyan', ls='--', zorder=8, label=None),
    "o":dict(color='orange', ls='--', zorder=8, label=None),
    "S":dict(color='purple', ls='--', zorder=8, label=None),
    "D":dict(color='salmon', ls='--', zorder=8, label=None),
    "A":dict(color='pink', ls='--', zorder=8, label=None),
    "U":dict(color='grey', ls='--', zorder=8, label=None),
    "V":dict(color='magenta', ls='--', zorder=8, label=None),
    "H":dict(color='lime', ls='--', zorder=8, label=None),
    "J":dict(color='navy', ls='--', zorder=8, label=None),
    "K":dict(color='peru', ls='--', zorder=8, label=None),
    "Q":dict(color='violet', ls='--', zorder=8, label=None),
    "W":dict(color='gold', ls='--', zorder=8, label=None),
    "C":dict(color='indigo', ls='--', zorder=8, label=None),
    "E":dict(color='tan', ls='--', zorder=8, label=None),
}

## for colour points
PROP5 = {
    "bin":dict(color='k', marker='o', ls='none', fillstyle='none', zorder=3, label='bin'),
    "fit":dict(color='g', marker='o', ls='none', fillstyle='none', zorder=2, label='fit'),
    "gp" :dict(color='b', marker='o', ls='none', fillstyle='none', zorder=1, label='GP'),
    "bin_all":dict(color='k', ls='--', zorder=3),
    "fit_all":dict(color='g', ls='--', zorder=2),
    "gp_all" :dict(color='b', ls='--', zorder=1),
}

## for bol fits
PROP6 = {
    # data
    "lyman_bin":dict(color='k', fillstyle='none', marker = 'o', ls='none', zorder=3, label='BC bin'),
    "lyman_fit":dict(color='g', fillstyle='none', marker = 'o', ls='none', zorder=2, label='BC fit'),
    "lyman_gp": dict(color='b', fillstyle='none', marker = 'o', ls='none', zorder=1, label='BC GP'),
    "bb_bin":   dict(color='k', fillstyle='none', marker = 's', ls='none', zorder=3, label='BB bin'),
    "bb_fit":   dict(color='g', fillstyle='none', marker = 's', ls='none', zorder=2, label='BB fit'),
    "bb_gp":    dict(color='b', fillstyle='none', marker = 's', ls='none', zorder=1, label='BB GP'),
    "bb_mix":   dict(color='r', fillstyle='none', marker = 's', ls='none', zorder=1, label='BB mix'),
    "spec":     dict(color='orange', fillstyle='none', marker = 'v', ls='none', zorder=1, label='spec'),
    # sed fit
    "cbb_bin_fit":  dict(color='k', ls='--', zorder=1),
    "cbb_gp_fit":   dict(color='b', ls='--', zorder=1),
    "cbb_fit_fit":  dict(color='g', ls='--', zorder=1),
    "cbb_mix_fit":  dict(color='r', ls='--', zorder=1),
    "spec_fit": dict(color='orange', ls='--', zorder=1),
    # bol fit
    "early_fit":dict(color='orange',ls='--',zorder=5,label='bol early fit'),
    "main_fit":dict(color='cyan', ls='--', zorder=4,label='bol main fit'),
    "tail_fit":dict(color='r',ls='--',zorder=4,label='bol tail fit'),
    "joint_fit":dict(color='grey',ls='--',zorder=6,label='bol full fit'),
}

## for spectra fits
PROP7 = {
    "fit":dict(color='blue',ls='-',label=None),
    "data":dict(color='grey',ls='-',label=None),
    "bindata":dict(color='brown',ls='--',label=None),
    "velocity":dict(color='k',ls='none',marker='o',label=None,fillstyle='none'),
    "velocity_model":dict(color='r',ls='-',label=None),
}
