import numpy as np

def exp(t, a, t0, b, c):
    return a*np.power(t-t0,b) + c  
