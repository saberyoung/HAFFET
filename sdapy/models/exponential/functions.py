import numpy as np

def exp(t, a, b, c):
    ''' exponential function
    '''
    return a * np.exp(-b * t) + c
