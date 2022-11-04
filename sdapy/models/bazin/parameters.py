modelmeta = {
    'bazin1':
    {
        'engine': 'multiband_main',
        'alias':
        [
            'bazin',            
            'bz',
            'bz1',
        ],
        'description': 'Bazin fits',
        'func': 'functions.bazin',        
        'parname' :
        [
            r'$A$', r'$t_\mathrm{0}$',
            r'$\tau_\mathrm{fall}$',
            r'$\tau_{rise}$', r'$C$',
        ],
        'par' :
        [
            'a', 'dt', 'tfall', 'trise', 'c',
        ],
        'bestv':
        {
            'a' : 100.,
            'dt' : 0,
            'tfall' : 20,
            'trise' : 10,
            'c' : 60,
        },
        'bounds':
        {
            'a' : [0,10000],
            'dt' : [-100, 100],
            'tfall' : [0, 120],
            'trise' : [0, 60],
            'c' : [0, 10000],
        },
    },
    'bazin2':
    {
        'engine': 'multiband_main',
        'alias':
        [           
            'bz2',
        ],
        'description': 'Bazin fit with f_sigma',
        'func': 'functions.bazin',        
        'parname' :
        [
            r'$A$', r'$t_\mathrm{0}$', r'$\tau_\mathrm{fall}$',
            r'$\tau_{rise}$', r'$C$',
        ],
        'par' :
        [
            'a', 'dt', 'tfall', 'trise', 'c',
        ],
        'bestv':
        {
            'a' : 100.,
            'dt' : 0,
            'tfall' : 20,
            'trise' : 10,
            'c' : 60,
        },
        'bounds':
        {
            'a' : [0,10000],
            'dt' : [-100, 100],
            'tfall' : [0, 120],
            'trise' : [0, 60],
            'c' : [0, 10000],
        },
        'same_parameter': ['dt'],
        'fit_error' : True,
    }
}
