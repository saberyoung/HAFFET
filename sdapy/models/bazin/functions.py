import numpy as np

def bazin1(time, A, t0, tfall):
    '''bazin model part I
    '''
    X = np.exp(-(time - t0) / tfall)
    bazinfunc = A**.5 * X    
    return bazinfunc

def bazin2(time, A, t0, trise):
    '''bazin model part II
    '''
    X = 1 / (1 + np.exp(-(time - t0) / trise))
    bazinfunc = A**.5 * X
    return bazinfunc

def bazin(time, A, t0, tfall, trise, C):
    '''bazin 2009 et al model
    '''
    X = bazin1(time, A, t0, tfall)*bazin2(time, A, t0, trise)
    bazinfunc = X + C
    return bazinfunc
