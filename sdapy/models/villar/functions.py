import numpy as np

def villar(time, a, b, t0, t1, tfall, trise, c):
    """
    Villar et al 2019 (https://iopscience.iop.org/article/10.3847/1538-4357/ab418c/pdf), function 1    

    Parameters
    ----------
    time : array
        Independent values.    
    a    : float
        A, Amplitude
    b    : float
        beta (flux/day), Plateau slope
    t0   : float
        “Start” time
    t1   : float
        plateau onset of LC
    tfall: float
        rise time of LC
    trise: float
        fall time of LC
    c    : float
        baseline flux
    """
    _flux = []
    for _time in time:
        if _time < t1:
            X = 1 / (1 + np.exp(-(_time - t0) / trise))
            _flux.append( (a + b*(_time-t0)) * X + c )
        else:
            X = np.exp(-(_time - t1) / tfall) / (1 + np.exp(-(_time - t0) / trise))
            _flux.append( (a + b*(t1-t0)) * X + c )
    return _flux
