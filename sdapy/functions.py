import numpy as np
import emcee, os, re, math
import h5py
from six import string_types
from dateutil.parser import parse
from sdapy import filters
import pandas as pd
# ztfquery source path
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

def Lbol_to_Mbol(L):
    # bolometric luminosity to bolometric mag
    return -2.5*( np.log10(L) - (71.21+17.5)*0.4 )

def Mbol_to_Lbol(M):
    # bolometric mag to bolometric luminosity
    return 10**(-0.4*M + (71.21+17.5)*0.4)

def redchisqg(ydata,ymod,deg=2,sd=None):  
    if sd is None:  chisq=np.sum((ydata-ymod)**2)  
    else:  chisq=np.sum( ((ydata-ymod)/sd)**2 ) 
    nu=ydata.size-1-deg  
    return chisq/nu

def BC_Lyman(color, colortype='g-r', phase='normal', sntype='Ic'):
    '''
    Lyman analytic bolometric correction for SE/II SNe
    https://ui.adsabs.harvard.edu/abs/2014MNRAS.437.3848L/abstract
    table 2, 3, 4
    '''
    bcfile = '%s/bc_table.txt' % LOCALSOURCE
    if not os.path.exists:
        print ('Error: BC file, %s, not found'%bcfile)
        return None, None
    typemap = dict()
    typemap['Ib'] = 'SESNe'
    typemap['Ibn'] = 'SESNe'
    typemap['Ic'] = 'SESNe'
    typemap['IIb'] = 'SESNe'
    typemap['Ic-BL'] = 'SESNe'
    typemap['II'] = 'II'
    if not sntype in typemap.keys():
        print ('Error: no BC found for type %s'%sntype)
        return None, None
    _type = typemap[sntype]
    if not phase in ['normal','cool']:
        print ('Error: wrong phase as %s'%phase)
        return None, None
    if phase == 'cool': _type = 'cool'
    if not colortype in ['g-r','g-i','B-V','B-R','B-I','V-R','V-I']:
        print ('Error: colortype %s not registered'%colortype)
        return None, None
    c1,c2 = colortype.split('-')
    _ = []
    for nn,ll in enumerate(open(bcfile).readlines()):
        if ll[0]=='#' or len(ll)==0:_.append(nn)
    bctab = pd.read_csv(bcfile, skiprows=_, delim_whitespace=True).query('type==@_type and x==@c1 and y==@c2')
    crange = [float(bctab['range'].to_list()[0].split('/')[0]),
              float(bctab['range'].to_list()[0].split('/')[1])]
    mbol = float(bctab['bc_c0'].to_list()[0]) + \
        float(bctab['bc_c1'].to_list()[0]) * color + \
        float(bctab['bc_c2'].to_list()[0]) * color**2    
    _ = np.logical_and(color >= min(crange), color <= max(crange))
    return mbol, _

def dessart_vej_to_vm(vpeak, vpeake, sntype, verbose):
    '''
    convert vej at peak to photospheirc vm, unit: 10*3 km/s
     following Dessart 16, convert photospheric velocity (vej) to characteristic velicity (vm)
        vej is O 7772 line velocity for SNe Ic, and He 5886 for Ib, at the peak epoch.
        https://academic.oup.com/mnras/article/458/2/1618/2589109
    '''
    if sntype in ['Ib', 'IIb']:
        vm, dvm = (vpeak-2.64) / 0.765, vpeake / 0.765
    elif sntype in ['Ic', 'Ic-BL']:
        vm, dvm = (vpeak-2.99) / 0.443, vpeake / 0.443
    else:
        if verbose: print ('check if Dessart et al suitable for your sn type')
        vm, dvm = None, None
    return vm, dvm

# pseudo equivalent width
def calc_pew(w,f,fc): return scipy.integrate.trapz(1-f/fc,w)

# wavelength (A) to velocity (1e3 km/s)
def calc_vel(x,x0): return -(x-x0)/x0*299792.458/1000.

# velocity (1e3 km/s) to wavelength (A)
def calc_wav(v,x0): return -v*1000./299792.458*x0+x0

def get_samples_mc(h5_file):    
    if not os.path.exists(h5_file):  return None, None
    reader = emcee.backends.HDFBackend(h5_file)
    tau = reader.get_autocorr_time(tol=0)   
    try: burnin = int(5*np.max(tau))
    except: return None, None
    samples = reader.get_chain(discard=burnin, thin=np.max(int(np.max(tau)), 0), flat=True)
    lnpost = reader.get_log_prob(discard=burnin, thin=np.max([int(np.max(tau)), 1]), flat=True)    
    try: return samples.value, lnpost.value
    except: return samples, lnpost

def get_samples_nstep(h5_file, thin_by=1):    
    if not os.path.exists(h5_file):  return    
    reader = emcee.backends.HDFBackend(h5_file)
    nsteps = thin_by*np.shape(reader.get_chain())[0]
    return nsteps

def get_samples_scipy(h5_file):    
    if not os.path.exists(h5_file):  return    
    reader = h5py.File(h5_file, 'r')    
    samples = reader['samples']
    lnpost = reader['lnprob']    
    try: return samples.value, lnpost.value
    except: return samples, lnpost

def get_numpy(v):
    try:    return v.to_numpy() # pandas dataframe
    except: return np.array(v)  # numpy array

def is_seq(o):
    """Check if the object is a sequence

    Parameters
    ----------
    o : any object
           The object to check
        
    Returns
    -------
        is_seq : bool, scalar
           True if *o* is a sequence, False otherwise
    """
    return hasattr(o, "__len__")

def is_coordinate(s):
    """Check if input is a coordinate."""
    if type(s) is str: 
        if len(s.split()) != 1: return False
    else:
        return False        
    matches = re.search(
        r'([+-]?([0-9]{1,2}):([0-9]{2})(:[0-9]{2})?\.?([0-9]+)?)', s)

    return False if matches is None else True

def is_date(s):
    """Check if input is a valid date."""
    try:
        parse(s)
        return True
    except ValueError:
        return False
    
def penc_to_errors(mean, perc_low, perc_up):        
    return abs(max(mean-perc_low, perc_up-mean))

def detect_symbol_1(ll, pattern='\"'): 
    quote_locations = []
    done, startIndex = False, 0
    while not done:
        startIndex = ll.find(pattern, startIndex)
        endIndex = ll.find(pattern, startIndex + 1)        
        if startIndex != -1 and endIndex != -1:
            quote_locations.append([startIndex, endIndex])
        else:
            done = True
        startIndex = endIndex + 1        
    return quote_locations

def detect_symbol_2(ll, pattern='\"'): 
    quote_locations = []
    done, startIndex = False, 0
    while not done:
        startIndex = ll.find(pattern, startIndex)      
        if startIndex != -1:
            quote_locations.append(startIndex)
        else:
            done = True
        startIndex += 1        
    return quote_locations

def replace_str_index(text,index=0,replacement=''):
    return '%s%s%s'%(text[:index],replacement,text[index+1:])

def remove_quote(ll): # detect quote
    nq = detect_symbol_1(ll, pattern='\"')
    if len(nq) == 0: return ll
    nc = detect_symbol_2(ll, pattern=',')    
    for ncid in nc:        
        for q in nq:
            if ncid > q[0] and ncid < q[1]:
                ll = replace_str_index(ll,index=ncid,replacement=' ')
    return ll.replace('\"','')

def is_number(s):
    """Check if input is a number."""
    if s is None: return False
    if isinstance(s, list) and not isinstance(s, string_types):
        try:
            for x in s:
                if isinstance(x, string_types) and ' ' in x:
                    raise ValueError
            _ = [float(x) for x in s]
            return True
        except ValueError:
            return False
    else:
        try:
            if isinstance(s, string_types) and ' ' in s:
                raise ValueError
            float(s)
            return True
        except ValueError:
            return False

def is_integer(s):
    """Check if input is an integer."""
    if isinstance(s, list) and not isinstance(s, string_types):
        try:
            _ = [int(x) for x in s]
            return True
        except ValueError:
            return False
    else:
        try:
            int(s)
            return True
        except ValueError:
            return False

def str_clean(s):
    # for oac meta, at some point there're multiple objects in one column
    if type(s) is str:  return s.split()[0]
    else: return s
    
def filter_clean(filt):
    if type(filt) != str: return '-'
    filt = filt.strip().lower().replace('_','').replace(':','')
    if filt.startswith('sdss'): filt = filt.replace('sdss', '')
    if filt.startswith('ztf'): filt = filt.replace('ztf', '')
    if filt.startswith('atlas'): filt = filt.replace('atlas', '')
    if filt.startswith('uvot'): filt = filt.replace('uvot', '')
    if filt == 'uvm2': filt = 'D'
    if filt == 'uvw1': filt = 'A'
    if filt == 'uvw2': filt = 'S'
    if filt == 'fuv': filt = 'F'
    if filt == 'nuv': filt = 'N'
    if filt == 'w1': filt = 'W'
    if filt == 'w2': filt = 'Q'
    filt = filt[0]
    if not filt.lower() in filters.central_wavelengths and not filt.upper() in filters.central_wavelengths:
        print ('!!! Error: check filter: %s/%s not registered correctly' % (filt.upper(), filt.lower()))
        return '-'
    if filt.lower() in filters.central_wavelengths:
        return filt.lower()
    else:        
        return filt.upper()
    
def type_clean(sntype):
    if type(sntype) == str:
        return sntype.replace('SLSN ','').replace('SN ','').replace('AT ','').strip()
    elif type(sntype) == float:
        if math.isnan(sntype): return '-'
        else: return str(sntype)
        
# Borrowed from astrocats' supernova catalog module.
def name_clean(name):
    """Clean a transient's name."""
    if name.strip() in ['-', 'None', 'nan']: return
    
    newname = name.strip(' ;,*.')
    if newname.startswith('NAME '):
        newname = newname.replace('NAME ', '', 1)
    if newname.endswith(' SN'):
        newname = newname.replace(' SN', '')
    if newname.endswith(':SN'):
        newname = newname.replace(':SN', '')
    if newname.startswith('MASJ'):
        newname = newname.replace('MASJ', 'MASTER OT J', 1)
    if (newname.startswith('MASTER') and len(newname) > 7 and
            is_number(newname[7])):
        newname = newname.replace('MASTER', 'MASTER OT J', 1)
    if (newname.startswith('MASTER OT') and len(newname) > 10 and
            is_number(newname[10])):
        newname = newname.replace('MASTER OT', 'MASTER OT J', 1)
    if newname.startswith('MASTER OT J '):
        newname = newname.replace('MASTER OT J ', 'MASTER OT J', 1)
    if newname.startswith('OGLE '):
        newname = newname.replace('OGLE ', 'OGLE-', 1)
    if newname.startswith('OGLE-') and len(newname) != 16:
        namesp = newname.split('-')
        if (len(namesp) == 4 and len(namesp[1]) == 4 and
                is_number(namesp[1]) and is_number(namesp[3])):
            newname = 'OGLE-' + namesp[1] + '-SN-' + namesp[3].zfill(3)
        elif (len(namesp) == 2 and is_number(namesp[1][:2]) and
              not is_number(namesp[1][2:])):
            newname = 'OGLE' + namesp[1]
    if newname.startswith('SN SDSS'):
        newname = newname.replace('SN SDSS ', 'SDSS', 1)
    if newname.startswith('SDSS '):
        newname = newname.replace('SDSS ', 'SDSS', 1)
    if newname.startswith('SDSS'):
        namesp = newname.split('-')
        if (len(namesp) == 3 and is_number(namesp[0][4:]) and
                is_number(namesp[1]) and is_number(namesp[2])):
            newname = namesp[0] + '-' + namesp[1] + '-' + namesp[2].zfill(3)
    if newname.startswith('SDSS-II SN'):
        namesp = newname.split()
        if len(namesp) == 3 and is_number(namesp[2]):
            newname = 'SDSS-II SN ' + namesp[2].lstrip('0')
    if newname.startswith('SN CL'):
        newname = newname.replace('SN CL', 'CL', 1)
    if newname.startswith('SN HiTS'):
        newname = newname.replace('SN HiTS', 'SNHiTS', 1)
    if newname.startswith('SNHiTS '):
        newname = newname.replace('SNHiTS ', 'SNHiTS', 1)
    if newname.startswith('GAIA'):
        newname = newname.replace('GAIA', 'Gaia', 1)
    if newname.startswith('KSN-'):
        newname = newname.replace('KSN-', 'KSN', 1)
    if newname.startswith('KSN'):
        newname = 'KSN' + newname[3:].lower()
    if newname.startswith('Gaia '):
        newname = newname.replace('Gaia ', 'Gaia', 1)
    if newname.startswith('Gaia'):
        newname = 'Gaia' + newname[4:].lower()
    if newname.startswith('GRB'):
        newname = newname.replace('GRB', 'GRB ', 1)
    if newname.startswith('GRB ') and is_number(newname[4:].strip()):
        newname = 'GRB ' + newname[4:].strip() + 'A'
    if newname.startswith('ESSENCE '):
        newname = newname.replace('ESSENCE ', 'ESSENCE', 1)
    if newname.startswith('LSQ '):
        newname = newname.replace('LSQ ', 'LSQ', 1)
    if newname.startswith('LSQ') and is_number(newname[3]):
        newname = newname[:3] + newname[3:].lower()
    if newname.startswith('DES') and is_number(newname[3]):
        newname = newname[:7] + newname[7:].lower()
    if newname.startswith('SNSDF '):
        newname = newname.replace(' ', '')
    if newname.startswith('SNSDF'):
        namesp = newname.split('.')
        if len(namesp[0]) == 9:
            newname = namesp[0] + '-' + namesp[1].zfill(2)
    if newname.startswith('HFF '):
        newname = newname.replace(' ', '')
    if newname.startswith('SN HST'):
        newname = newname.replace('SN HST', 'HST', 1)
    if newname.startswith('HST ') and newname[4] != 'J':
        newname = newname.replace('HST ', 'HST J', 1)
    if newname.startswith('SNLS') and newname[4] != '-':
        newname = newname.replace('SNLS', 'SNLS-', 1)
    if newname.startswith('SNLS- '):
        newname = newname.replace('SNLS- ', 'SNLS-', 1)
    if newname.startswith('CRTS CSS'):
        newname = newname.replace('CRTS CSS', 'CSS', 1)
    if newname.startswith('CRTS MLS'):
        newname = newname.replace('CRTS MLS', 'MLS', 1)
    if newname.startswith('CRTS SSS'):
        newname = newname.replace('CRTS SSS', 'SSS', 1)
    if newname.startswith(('CSS', 'MLS', 'SSS')):
        newname = newname.replace(' ', ':').replace('J', '')
    if newname.startswith('SN HFF'):
        newname = newname.replace('SN HFF', 'HFF', 1)
    if newname.startswith('SN GND'):
        newname = newname.replace('SN GND', 'GND', 1)
    if newname.startswith('SN SCP'):
        newname = newname.replace('SN SCP', 'SCP', 1)
    if newname.startswith('SN UDS'):
        newname = newname.replace('SN UDS', 'UDS', 1)
    if newname.startswith('SCP') and newname[3] != '-':
        newname = newname.replace('SCP', 'SCP-', 1)
    if newname.startswith('SCP- '):
        newname = newname.replace('SCP- ', 'SCP-', 1)
    if newname.startswith('SCP-') and is_integer(newname[7:]):
        newname = 'SCP-' + newname[4:7] + str(int(newname[7:]))
    if newname.startswith('PS 1'):
        newname = newname.replace('PS 1', 'PS1', 1)
    if newname.startswith('PS1 SN PS'):
        newname = newname.replace('PS1 SN PS', 'PS', 1)
    if newname.startswith('PS1 SN'):
        newname = newname.replace('PS1 SN', 'PS1', 1)
    if newname.startswith('PS1') and is_number(newname[3]):
        newname = newname[:3] + newname[3:].lower()
    elif newname.startswith('PS1-') and is_number(newname[4]):
        newname = newname[:4] + newname[4:].lower()
    if newname.startswith('PSN K'):
        newname = newname.replace('PSN K', 'K', 1)
    if newname.startswith('K') and is_number(newname[1:5]):
        namesp = newname.split('-')
        if len(namesp[0]) == 5:
            newname = namesp[0] + '-' + namesp[1].zfill(3)
    if newname.startswith('Psn'):
        newname = newname.replace('Psn', 'PSN', 1)
    if newname.startswith('PSNJ'):
        newname = newname.replace('PSNJ', 'PSN J', 1)
    if newname.startswith('TCPJ'):
        newname = newname.replace('TCPJ', 'TCP J', 1)
    if newname.startswith('SMTJ'):
        newname = newname.replace('SMTJ', 'SMT J', 1)
    if newname.startswith('PSN20J'):
        newname = newname.replace('PSN20J', 'PSN J', 1)
    if newname.startswith('SN ASASSN'):
        newname = newname.replace('SN ASASSN', 'ASASSN', 1)
    if newname.startswith('ASASSN-20') and is_number(newname[9]):
        newname = newname.replace('ASASSN-20', 'ASASSN-', 1)
    if newname.startswith('ASASSN '):
        newname = newname.replace('ASASSN ', 'ASASSN-', 1).replace('--', '-')
    if newname.startswith('ASASSN') and newname[6] != '-':
        newname = newname.replace('ASASSN', 'ASASSN-', 1)
    if newname.startswith('ASASSN-') and is_number(newname[7]):
        newname = newname[:7] + newname[7:].lower()
    if newname.startswith('ROTSE3J'):
        newname = newname.replace('ROTSE3J', 'ROTSE3 J', 1)
    if newname.startswith('MACSJ'):
        newname = newname.replace('MACSJ', 'MACS J', 1)
    if newname.startswith('MWSNR'):
        newname = newname.replace('MWSNR', 'MWSNR ', 1)
    if newname.startswith('SN HUNT'):
        newname = newname.replace('SN HUNT', 'SNhunt', 1)
    if newname.startswith('SN Hunt'):
        newname = newname.replace(' ', '')
    if newname.startswith('SNHunt'):
        newname = newname.replace('SNHunt', 'SNhunt', 1)
    if newname.startswith('SNhunt '):
        newname = newname.replace('SNhunt ', 'SNhunt', 1)
    if newname.startswith('ptf'):
        newname = newname.replace('ptf', 'PTF', 1)
    if newname.startswith('SN PTF'):
        newname = newname.replace('SN PTF', 'PTF', 1)
    if newname.startswith('PTF '):
        newname = newname.replace('PTF ', 'PTF', 1)
    if newname.startswith('PTF') and is_number(newname[3]):
        newname = newname[:3] + newname[3:].lower()
    if newname.startswith('IPTF'):
        newname = newname.replace('IPTF', 'iPTF', 1)
    if newname.startswith('iPTF '):
        newname = newname.replace('iPTF ', 'iPTF', 1)
    if newname.startswith('iPTF') and is_number(newname[4]):
        newname = newname[:4] + newname[4:].lower()
    if newname.startswith('PESSTOESO'):
        newname = newname.replace('PESSTOESO', 'PESSTO ESO ', 1)
    if newname.startswith('snf'):
        newname = newname.replace('snf', 'SNF', 1)
    if newname.startswith('SNF '):
        newname = newname.replace('SNF ', 'SNF', 1)
    if (newname.startswith('SNF') and is_number(newname[3:]) and
            len(newname) >= 12):
        newname = 'SNF' + newname[3:11] + '-' + newname[11:]
    if newname.startswith(('MASTER OT J', 'ROTSE3 J')):
        prefix = newname.split('J')[0]
        coords = newname.split('J')[-1].strip()
        decsign = '+' if '+' in coords else '-'
        coordsplit = coords.replace('+', '-').split('-')
        if ('.' not in coordsplit[0] and len(coordsplit[0]) > 6 and
                '.' not in coordsplit[1] and len(coordsplit[1]) > 6):
            newname = (
                prefix + 'J' + coordsplit[0][:6] + '.' + coordsplit[0][6:] +
                decsign + coordsplit[1][:6] + '.' + coordsplit[1][6:])
    if (newname.startswith('Gaia ') and is_number(newname[3:4]) and
            len(newname) > 5):
        newname = newname.replace('Gaia ', 'Gaia', 1)
    if (newname.startswith('AT ') and is_number(newname[3:7]) and
            len(newname) > 7):
        newname = newname.replace('AT ', 'AT', 1)
    if len(newname) <= 4 and is_number(newname):
        newname = 'SN' + newname + 'A'
    if (len(newname) > 4 and is_number(newname[:4]) and
            not is_number(newname[4:])):
        newname = 'SN' + newname
    if (newname.startswith('Sn ') and is_number(newname[3:7]) and
            len(newname) > 7):
        newname = newname.replace('Sn ', 'SN', 1)
    if (newname.startswith('sn') and is_number(newname[2:6]) and
            len(newname) > 6):
        newname = newname.replace('sn', 'SN', 1)
    if (newname.startswith('SN ') and is_number(newname[3:7]) and
            len(newname) > 7):
        newname = newname.replace('SN ', 'SN', 1)
    if (newname.startswith('SN') and is_number(newname[2:6]) and
            len(newname) == 7 and newname[6].islower()):
        newname = 'SN' + newname[2:6] + newname[6].upper()
    elif (newname.startswith('SN') and is_number(newname[2:6]) and
          (len(newname) == 8 or len(newname) == 9) and newname[6:].isupper()):
        newname = 'SN' + newname[2:6] + newname[6:].lower()
    if (newname.startswith('AT') and is_number(newname[2:6]) and
            len(newname) == 7 and newname[6].islower()):
        newname = 'AT' + newname[2:6] + newname[6].upper()
    elif (newname.startswith('AT') and is_number(newname[2:6]) and
          (len(newname) == 8 or len(newname) == 9) and newname[6:].isupper()):
        newname = 'AT' + newname[2:6] + newname[6:].lower()

    newname = (' '.join(newname.split())).strip()
    return newname
