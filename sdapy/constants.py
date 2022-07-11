''' ##################
define constants and prior boundries for models
''' ##################
import numpy as np

# covert things to cgs units
M_sun=1.989e33
R_sun=696340e5
c=3.e10
tau_Ni=8.8
tau_Co=111.3
e_Ni=3.9e10
e_Co=6.8e9
k_opt=0.07 # k_opt=0.07 g/cm^2, fixed optical opacity (corresponds to electron scattering) for SE SNe
C=.05
k_gamma = 0.027
beta=13.8  # beta=13.8, constant of integration (Arnett 1982)
gamma=6./5.

# define best fit and boundaries for parameters
# fsigma is the factor to enlarge errors
fsigma_p0 = 1
fsigma_bounds = [0, 10]
# power law fits
a_p0 = 100.
a_bounds = [0,10000]
c_p0 = 60.
c_bounds = [0, 10000]
alpha_p0 = 2
alpha_bounds = [0, 100]
# Bazin fits     
dt_p0 = 0
dt_bounds = [-100, 100]
tfall_p0 = 10
tfall_bounds = [0, 120]
trise_p0 = 10
trise_bounds = [0, 60]
#Villar fits
b_p0 = -100
b_bounds = [-10000, 0]
ts_p0 = -20
ts_bounds = [-100, 100]
te_p0 = 20
te_bounds = [-100, 100]
# nickel fits
mni_p0 = .2
mni_bounds = [0, 5]
taum_p0 = 10
taum_bounds = [0, 100]
t0_p0 = 100
t0_bounds = [0, 1e3]
mej_p0 = 1
mej_bounds = [.1, 50]
ek_p0 = 1
ek_bounds = [.1, 100]    
vm_p0 = 10
vm_bounds = [1, 40]
texp_p0 = -18
texp_bounds = [-100, 0]
ts_p0 = 60
ts_bounds = [20, 120]
tc_p0 = -18
tc_bounds = [-100, 0]
# SBO fits
me_p0 = .1
me_bounds = [0, 1]
re_p0 = 100
re_bounds = [0, 1e3]
ee_p0 = 1e-2
ee_bounds = [1e-3, 10]
# spectral lines
line_location = {
    'H~$\alpha$': 6563,
    'H~$\beta$': 4861,
    'H~$\gamma$': 4340,
    'H~$\delta$': 4102,
    'H~$\epsilon$': 3970,
    'He~5876$\AA$': 5876,
    'O~7774$\AA$': 7774,
    'Na~I~D1': 5896,
    'Na~I~D2': 5890,
}
line_forsn = {
    'SN Ib': 'He~5876$\AA$',       
    'SN IIb': 'He~5876$\AA$',      
    'SN Ic': 'O~7774$\AA$',          
    'SN Ic-BL': 'O~7774$\AA$',       
}

line_velocity_p0 = 10 # *1e3 km/s
line_velocity_bounds = [2, 30]
line_width_p0 = 50
line_width_bounds = [0,500]
# calculate gaussian parameters
gh_p0 = 1 # normalized offset
gh_bounds = [0, 2] 
ga_p0 = -1 # normalized amplititude, - for absorption
ga_bounds = [-2, 0]
gx_p0 = -200  # line location relative to intrinstic wavelength in AA
gx_bounds = [-1000, 0]
gs_p0 = line_width_p0  # line width in AA
gs_bounds = line_width_bounds
