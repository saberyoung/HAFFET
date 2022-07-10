modelmeta = {
    'powerlaw_multiple':
    {
        'engine': 'multiband_early',
        'alias': [
            'pl', 'powerlaw', 'Powerlaw',
        ],
        'description': 'Power law fits',
        'func': 'functions.powerlaw_full',       
        'parname' : [
            r'$t_\mathrm{fl}$', r"$A$", r"$\alpha$", r"$C$",
        ],
        'par' : [
            'texp', 'a', 'alpha', 'c',
        ],
        'same_parameter': ['texp'],
        'fit_error' : True,
    },
    'powerlaw_single':
    {
        'engine': 'multiband_early',
        'alias': [
            'pl2', 'powerlaw2', 'Powerlaw2',
            'pl_single', 'powerlaw_single', 'Powerlaw_single',
        ],
        'description': 'Power law fits',
        'func': 'functions.powerlaw_full',        
        'parname' : [
            r'$t_\mathrm{fl}$', r"$A$", r"$\alpha$", r"$C$",
        ],
        'par' : [
            'texp', 'a', 'alpha', 'c',
        ],        
    },   
}
