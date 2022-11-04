def powerlaw_post_baseline(times, t_0=0, amplitude=25, alpha_r=2):
    ''' power law function
    '''
    return amplitude * (times - t_0)**alpha_r

def powerlaw_full(times, t_0=0, amplitude=25, alpha_r=2, c=0):
    ''' power law fits (based on Miller et al, https://ui.adsabs.harvard.edu/abs/2020ApJ...902...47M/abstract), constant (background) before explosion, and power law post explosion.
    '''
    f = []
    for t in times:
        if t<t_0:
            f.append( c )
        else:
            f.append( powerlaw_post_baseline(t, t_0, amplitude, alpha_r) + c )
    return f
