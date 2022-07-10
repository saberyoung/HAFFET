modelmeta = {
    'bazin1':
    {
        'engine': 'multiband_main',
        'alias': [
            'bazin',
            'bazin1',
            'bz',
            'bz1',
        ],
        'description': 'Bazin fits',
        'func': 'functions.bazin',        
        'parname' : [
            r'$A$', r'$t_\mathrm{0}$', r'$\tau_\mathrm{fall}$',
            r'$\tau_{rise}$', r'$C$',
        ],
        'par' : [
            'a', 'dt', 'tfall', 'trise', 'c',
        ],        
    },
    'bazin2':
    {
        'engine': 'multiband_main',
        'alias': [
            'bazin2',
            'bz2',
        ],
        'description': 'Bazin fit with f_sigma',
        'func': 'functions.bazin',        
        'parname' : [
            r'$A$', r'$t_\mathrm{0}$', r'$\tau_\mathrm{fall}$',
            r'$\tau_{rise}$', r'$C$',
        ],
        'par' : [
            'a', 'dt', 'tfall', 'trise', 'c',
        ],
        'same_parameter': ['dt'],
        'fit_error' : True,
    }
}
