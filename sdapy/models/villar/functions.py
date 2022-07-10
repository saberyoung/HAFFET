def villar(time, a, b, t0, t1, tfall, trise):
    """
    Villar et al 2019
    """
    _flux = []
    for _time in time:
        if _time < t1:
            X = 1 / (1 + np.exp(-(_time - t0) / trise))
            _flux.append( (a + b*(_time-t0)) * X )
        else:
            X = np.exp(-(_time - t1) / tfall) / (1 + np.exp(-(_time - t0) / trise))
            _flux.append( (a + b*(t1-t0)) * X )
    return _flux
