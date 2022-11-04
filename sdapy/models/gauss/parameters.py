modelmeta = {
    'gauss':
    {
        'engine': 'specline',
        'alias':
        [
            'gauss', 'singlegauss', 'gaussian'
        ],
        'description': 'Gaussain',
        'func': 'functions.gauss',
        'parname' :
        [
            r'$H$', r'$A$', r'$x_\mathrm{0}$', r'$\sigma$',
        ],
        'par' :
        [
            'h', 'a', 'x0', 'sigma',
        ],
        'bestv':
        {
            'h' : 1,
            'a' : -1.,
            'x0' : -200,
            'sigma' : 50,                    
        },
        'bounds':
        {
            'h' : [0, 2],
            'a' : [-2, 0],
            'x0' : [-1000, 0],
            'sigma' : [0, 1000],                
        },
    },
    'doublegauss':
    {
        'engine': 'specline',
        'alias':
        [
            'doublegauss', 'doublegaussian'
        ],
        'description': 'Double Gaussian',
        'func': 'functions.double_gauss',
        'parname' :
        [
            r'$\mathrm{H}$',
            r'$A_\mathrm{1}$', r'$x_\mathrm{0,1}$', r'$\sigma_\mathrm{1}$',
            r'$A_\mathrm{2}$', r'$x_\mathrm{0,2}$', r'$\sigma_\mathrm{2}$',
        ],
        'par' :
        [
            'h', 'a', 'x0', 'sigma', 'a2', 'x02', 'sigma2',
        ],
        'bestv':
        {
            'h'  : 1,
            'a' : -1.,
            'x0' : -200,
            'sigma' : 50,
            'a2' : -1.,
            'x02': -200,
            'sigma2' : 50,   
        },
        'bounds':
        {
            'h'  : [0, 2],
            'a' : [-2, 0],
            'x0': [-1000, 0],
            'sigma' : [0, 1000],
            'a2' : [-2, 0],
            'x02': [-1000, 0],
            'sigma2' : [0, 1000],  
        },
    },    
}
