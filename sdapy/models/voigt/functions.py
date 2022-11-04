import numpy as np
from scipy.special import wofz

def voigt(x, H, A, x0, sigma, gamma):
    ''' voigt function
    '''
    return H + A * np.real(wofz((x - x0 + 1j*gamma)/sigma/np.sqrt(2))) / sigma /np.sqrt(2*np.pi)
