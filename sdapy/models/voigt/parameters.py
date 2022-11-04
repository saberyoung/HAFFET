modelmeta = {
    'voigt':
    {
        'engine': 'specline',
        'alias':
        [
            'voigt',
        ],
        'description': 'voigt fits',
        'func': 'functions.voigt',
        'parname' :
        [
            r'$H$', r'$A$', r'$x_\mathrm{0}$', r'$\sigma$', r'$\gamma$',
        ],
        'par' :
        [
            'h', 'a', 'x0', 'sigma', 'gamma'
        ],
        'bestv':
        {
            'h' : 1,
            'a' : -100.,
            'x0' : -200,
            'sigma' : 50,
            'gamma' : 1,
        },
        'bounds':
        {
            'h' : [0, 2],
            'a' : [-1000, 0],
            'x0' : [-1000, 0],
            'sigma' : [0, 500],
            'gamma' : [0, 500],
        },
    },    
}
