def powerlaw_post_baseline(times, t_0=0, amplitude=25, alpha_r=2):    
    return amplitude * (times - t_0)**alpha_r

def powerlaw_full(times, t_0=0, amplitude=25, alpha_r=2, c=0):
    f = []
    for t in times:
        if t<t_0:
            f.append( c )
        else:
            f.append( powerlaw_post_baseline(t, t_0, amplitude, alpha_r) + c )
    return f
