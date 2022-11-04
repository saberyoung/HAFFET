''' ##################
define constants and prior boundries for models
''' ##################
import numpy as np

# covert things to cgs units
day_to_s = 86400.
M_sun=1.989e33
R_sun=696340e5
c=3.e10
km_cgs=100000
au_cgs = 14959787070000
tau_Ni=8.8
tau_Co=111.3
e_Ni=3.9e10
e_Co=6.8e9
# unit in g/cm^2, fixed optical opacity (corresponds to electron scattering) for SE SNe, Nagy2018, table 6
k_opt = {
    'II': 0.21,
    'IIP': 0.21,
    'IIb': 0.20,
    'Ib':  0.18,
    'Ic':  0.10,
    'Ic-BL': 0.07,    
}
C=.05
# gamma opacity
k_gamma = 0.027
# beta=13.8, constant of integration (Arnett 1982)
beta=13.8
gamma=6./5.

# define best fit and boundaries for parameters
# fsigma is the factor to enlarge errors
fsigma_p0 = 1
fsigma_bounds = [0, 10]

# spectral lines
line_location = {
    r'H~$\alpha$': 6563,
    r'H~$\beta$': 4861,
    r'H~$\gamma$': 4340,
    r'H~$\delta$': 4102,
    r'H~$\epsilon$': 3970,
    r'He~I~3889$\AA$': 3889,
    r'He~I~5876$\AA$': 5876,
    r'He~I~6678$\AA$': 6678,
    r'He~I~7065$\AA$': 7065,
    r'He~II~4686$\AA$': 4686,
    r'O~I~5577$\AA$': 5577,
    r'O~I~7774$\AA$': 7774,
    r'O~I~doublet(6300, 6364 $\AA$)': (6300 + 6364)/2,
    r'O~I~8446$\AA$': 8446,
    r'Na~ID': (5896 + 5890)/2.,
    r'Ca~II~doublet(7291,7324$\AA$)': (7291 + 7324)/2,
    r'Ca~II~triplet(8498,8542,8662$\AA$)': (8498 + 8542 + 8662)/3,
}

telluric_lines = {
    r'telluric~$I$': [6860, 6890],
    r'telluric~$II$': [7600, 7630],
    r'telluric~$II$': [7170, 7350],
}

line_forsn = {
    'Ib': r'He~I~5876$\AA$',
    'IIb': r'He~I~5876$\AA$',
    'Ic': r'O~I~7774$\AA$',
    'Ic-BL': r'O~I~7774$\AA$',
}
