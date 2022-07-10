modelmeta = {
    'gauss':
    {
        'engine': 'specline',
        'alias': [
            'gauss', 'singlegauss', 'gaussian'
        ],
        'description': 'Gaussain',
        'func': 'functions.gauss',
        'parname' : [
            r'$H$', r'$A$', r'$x_\mathrm{0}$', r'$\sigma$',
        ],
        'par' : [
            'gh', 'ga', 'gx', 'gs',
        ],       
    },
    'doublegauss':
    {
        'engine': 'specline',
        'alias': [
            'doublegauss', 'doublegaussian'
        ],
        'description': 'Double Gaussian',
        'func': 'functions.double_gauss',
        'parname' : [
            r'$H_\mathrm{1}$', r'$A_\mathrm{1}$', r'$x_\mathrm{0,1}$', r'$\sigma_\mathrm{1}$',
            r'$H_\mathrm{2}$', r'$A_\mathrm{2}$', r'$x_\mathrm{0,2}$', r'$\sigma_\mathrm{2}$',
        ],
        'par' : [
            'a', 'dt', 'tfall', 'trise', 'c',
            'a', 'dt', 'tfall', 'trise', 'c',
        ],       
    },    
}
