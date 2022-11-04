modelmeta = {
    'powerlaw_multiple':
    {
        'engine': 'multiband_early',
        'alias':
        [
            'pl',
            'powerlaw',
            'Powerlaw',
        ],
        'description': 'Power law fits',
        'func': 'functions.powerlaw_full',       
        'parname' :
        [
            r'$t_\mathrm{fl}$', r"$A$", r"$\alpha$", r"$C$",
        ],
        'par' :
        [
            'texp', 'a', 'alpha', 'c',
        ],
        'bestv':
        {
            'texp' : -18,
            'a' : 100,
            'alpha' : 2,           
            'c' : 60,
        },
        'bounds':
        {
            'texp' : [-100, 0],
            'a' : [0, 10000],
            'alpha' : [0, 10],
            'c' : [0, 10000],
        },
        'same_parameter': ['texp'],
        'fit_error' : True,
    },
    'powerlaw_single':
    {
        'engine': 'multiband_early',
        'alias':
        [
            'pl2',
            'powerlaw2',
            'Powerlaw2',
            'pl_single',
            'Powerlaw_single',
        ],
        'description': 'Power law fits',
        'func': 'functions.powerlaw_full',        
        'parname' :
        [
            r'$t_\mathrm{fl}$', r"$A$", r"$\alpha$", r"$C$",
        ],
        'par' :
        [
            'texp', 'a', 'alpha', 'c',
        ],
        'bestv':
        {
            'texp' : -18,
            'a' : 100,
            'alpha' : 2,           
            'c' : 60,
        },
        'bounds':
        {
            'texp' : [-100, 0],
            'a' : [0, 10000],
            'alpha' : [0, 10],
            'c' : [0, 10000],
        },
    },
}
