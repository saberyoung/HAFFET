# central wavelengths (in Angs)
central_wavelengths = {'u': 3560,  'g': 4830, 'r': 6260, 'i': 7670,
                       'z': 8890, 'y': 9600, 'w':5985, 'Y': 9600,
                       'U': 3600,  'B': 4380, 'V': 5450, 'R': 6410,
                       'G': 6730, 'I': 7980, 'J': 12200, 'H': 16300,
                       'K': 21900, 'S': 2030, 'D': 2231, 'A': 2634,
                       'F': 1516, 'N': 2267, 'o': 6790, 'c': 5330,
                       'W': 33526, 'Q': 46028
}

#All values in 1e-11 erg/s/cm2/Angs
zp = {'u': 859.5, 'g': 466.9, 'r': 278.0, 'i': 185.2, 'z': 137.8, 'y': 118.2, 'w': 245.7, 'Y': 118.2,
      'U': 417.5, 'B': 632.0, 'V': 363.1, 'R': 217.7, 'G': 240.0, 'I': 112.6, 'J': 31.47, 'H': 11.38,
      'K': 3.961, 'S': 536.2, 'D': 463.7, 'A': 412.3, 'F': 4801., 'N': 2119., 'o': 236.2, 'c': 383.3,
      'W': 0.818, 'Q': 0.242}

#Filter widths (in Angs)
filter_width = {'u': 458,  'g': 928, 'r': 812, 'i': 894,  'z': 1183, 'y': 628, 'w': 2560, 'Y': 628,
                'U': 485,  'B': 831, 'V': 827, 'R': 1389, 'G': 4203, 'I': 899, 'J': 1759, 'H': 2041,
                'K': 2800, 'S': 671, 'D': 446, 'A': 821,  'F': 268,  'N': 732, 'o': 2580, 'c': 2280,
                'W': 6626, 'Q': 10422}

# flux shift
ys = {'u':-40, 'g':0, 'B':0, 'r':40, 'R':40, 'i':80, 'z':200, 'c':120, 'o':160, }

# matplotlib pyplot style
## for data points
PROP1 = {
    "r":dict(color='red', marker='o', ls='none', fillstyle='none', label='r'),
    "g":dict(color='green', marker='o', ls='none', fillstyle='none', label='g'),
    "R":dict(color='red', marker='o', ls='none', fillstyle='none', label='R'),
    "B":dict(color='green', marker='o', ls='none', fillstyle='none', label='B'),
    "i":dict(color='k', marker='o', ls='none', fillstyle='none', label='i'),
    "z":dict(color='blue', marker='o', ls='none', fillstyle='none', label='z'),
    "u":dict(color='brown', marker='o', ls='none', fillstyle='none', label='u'),
    "c":dict(color='cyan', marker='o', ls='none', fillstyle='none', label='c'),
    "o":dict(color='orange', marker='o', ls='none', fillstyle='none', label='o'),
}

## for data points in flux plot
PROP1f = {
    "r":dict(color='red', marker='o', ls='none', fillstyle='none'),
    "g":dict(color='green', marker='o', ls='none', fillstyle='none'),
    "R":dict(color='red', marker='o', ls='none', fillstyle='none'),
    "B":dict(color='green', marker='o', ls='none', fillstyle='none'),
    "i":dict(color='k', marker='o', ls='none', fillstyle='none'),
    "z":dict(color='blue', marker='o', ls='none', fillstyle='none'),
    "u":dict(color='brown', marker='o', ls='none', fillstyle='none'),
    "c":dict(color='cyan', marker='o', ls='none', fillstyle='none'),
    "o":dict(color='orange', marker='o', ls='none', fillstyle='none'),
}

## for data upper limits
PROP1l = {
    "r":dict(color='red', marker='v', ls='none', fillstyle='none',),
    "g":dict(color='green', marker='v', ls='none', fillstyle='none',),
    "R":dict(color='red', marker='v', ls='none', fillstyle='none',),
    "B":dict(color='green', marker='v', ls='none', fillstyle='none',),
    "i":dict(color='k', marker='v', ls='none', fillstyle='none',),
    "z":dict(color='blue', marker='v', ls='none', fillstyle='none',),
    "u":dict(color='brown', marker='v', ls='none', fillstyle='none',),
    "c":dict(color='cyan', marker='v', ls='none', fillstyle='none',),
    "o":dict(color='orange', marker='v', ls='none', fillstyle='none',),    
}

## for fittings sampling curves
PROP2 = {
    "r":dict(color='red', ls='-', zorder=9),
    "g":dict(color='green', ls='-', zorder=9),
    "R":dict(color='red', ls='-', zorder=9),
    "B":dict(color='green', ls='-', zorder=9),
    "i":dict(color='k', ls='-', zorder=9),
    "z":dict(color='blue', ls='-', zorder=9),
    "u":dict(color='brown', ls='-', zorder=9),
    "c":dict(color='cyan', ls='-', zorder=9),
    "o":dict(color='orange', ls='-', zorder=9),
}

## for fitting peak axvline
PROP2p = {
    "r":dict(color='red', ls='--', zorder=8, label='r peak'),
    "g":dict(color='green', ls='--', zorder=8, label=None),
    "R":dict(color='red', ls='--', zorder=8, label='R peak'),
    "B":dict(color='green', ls='--', zorder=8, label=None),
    "i":dict(color='k', ls='--', zorder=8, label=None),
    "z":dict(color='blue', ls='--', zorder=8, label=None),
    "u":dict(color='brown', ls='--', zorder=8, label=None),
    "c":dict(color='cyan', ls='--', zorder=8, label=None),
    "o":dict(color='orange', ls='--', zorder=8, label=None),
}

## for Gaussian process sampling curve
PROP3 = {
    "r":dict(color='red', ls=':', zorder=7),
    "g":dict(color='green', ls=':', zorder=7),
    "R":dict(color='red', ls=':', zorder=7),
    "B":dict(color='green', ls=':', zorder=7),
    "i":dict(color='k', ls=':', zorder=7),
    "z":dict(color='blue', ls=':', zorder=7),
    "u":dict(color='brown', ls=':', zorder=7),
    "c":dict(color='cyan', ls=':', zorder=7),
    "o":dict(color='orange', ls=':', zorder=7),
}

## for early power law fits
PROP4 = {
    "r":dict(color='red', ls='-', zorder=1),
    "g":dict(color='green', ls='-', zorder=1),
    "R":dict(color='red', ls='-', zorder=1),
    "B":dict(color='green', ls='-', zorder=1),
    "i":dict(color='k', ls='-', zorder=1),
    "z":dict(color='blue', ls='-', zorder=1),
    "u":dict(color='brown', ls='-', zorder=1),
    "c":dict(color='cyan', ls='-', zorder=1),
    "o":dict(color='orange', ls='-', zorder=1),
}

## for colour points
PROP5 = {
    "gr1":dict(color='k', marker='o', ls='none', fillstyle='none', zorder=3, label='<1d'),
    "gr2":dict(color='g', marker='o', ls='none', fillstyle='none', zorder=2, label='+fit'),
    "gr3":dict(color='b', marker='o', ls='none', fillstyle='none', zorder=1, label='+GP'),
}

## for bol Arnett fits
PROP6 = {
    "fit":dict(color='cyan', ls='--', zorder=4),
    "lyman_data_1":dict(color='k', fillstyle='none', marker = 'o',
                        ls='none', alpha=.5, label='Lyman', zorder=3),
    "lyman_data_2":dict(color='g', fillstyle='none', marker = 'o',
                        ls='none', alpha=.5, label=None, zorder=2),
    "lyman_data_3":dict(color='b', fillstyle='none', marker = 'o',
                        ls='none', alpha=.5, label=None, zorder=1),   
    "bb_data_1":   dict(color='k', fillstyle='none', marker = 's',
                        ls='none', alpha=.5, label='BB pseudo', zorder=3),
    "bb_data_2":   dict(color='g', fillstyle='none', marker = 's',
                        ls='none', alpha=.5, label=None, zorder=2),
    "bb_data_3":   dict(color='b', fillstyle='none', marker = 's',
                        ls='none', alpha=.5, label=None, zorder=1),
    "bb_data_1f":  dict(color='k', fillstyle='none', marker = 'v',
                        ls='none', alpha=.5, label='BB SB', zorder=3),
    "bb_data_2f":  dict(color='g', fillstyle='none', marker = 'v',
                        ls='none', alpha=.5, label=None, zorder=2),
    "bb_data_3f":  dict(color='b', fillstyle='none', marker = 'v',
                        ls='none', alpha=.5, label=None, zorder=1),
}

## for Tail fits
PROP7 = {
    "fit":dict(color='r',ls='--',zorder=4),
}

## for spectra fits
PROP8 = {
    "fit":dict(color='purple',ls='--',alpha=.2),
    "data":dict(color='grey',ls='-',alpha=.8),
    "bindata":dict(color='k',ls='-',alpha=1),  
}


''' R from literature 
from Table 2 in Yuan, Liu, Xiang, MNRAS 2013
Rf = {
    "g": 3.294,
    "r": 2.393,
    "i": 1.836,    
}
'''

#Extinction coefficients in A_lam / E(B-V). Uses York Extinction Solver (http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/YorkExtinctionSolver/coefficients.cgi)
Rf = {'u': 4.786,  'g': 3.587, 'r': 2.471, 'i': 1.798,  'z': 1.403, 'y': 1.228, 'w':2.762, 'Y': 1.228,
      'U': 4.744,  'B': 4.016, 'V': 3.011, 'R': 2.386, 'G': 2.216, 'I': 1.684, 'J': 0.813, 'H': 0.516,
      'K': 0.337, 'S': 8.795, 'D': 9.270, 'A': 6.432,  'F': 8.054,  'N': 8.969, 'o': 2.185, 'c': 3.111,
      'W': 0.190, 'Q': 0.127}
