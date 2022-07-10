import numpy as np

# analytic model flux form from Bazin et al 2009
def bazin1(time, A, t0, tfall):
    X = np.exp(-(time - t0) / tfall)
    bazinfunc = A**.5 * X    
    return bazinfunc

def bazin2(time, A, t0, trise):
    X = 1 / (1 + np.exp(-(time - t0) / trise))
    bazinfunc = A**.5 * X
    return bazinfunc

def bazin(time, A, t0, tfall, trise, C):
    X = bazin1(time, A, t0, tfall)*bazin2(time, A, t0, trise)
    bazinfunc = X + C
    return bazinfunc
