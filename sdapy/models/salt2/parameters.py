modelmeta = {
    'salt':
    {
        'engine': 'multiband_main',
        'alias':
        [
            'salt2',
            'salt2_fit',
        ],
        'description': 'Salt2 fit',
        'func': 'functions.salt2_fit',
        'parname':
        [
            r'$t_\mathrm{0}$', r'$x_\mathrm{0}$',
            r"$x_\mathrm{1}$", r'$\mathrm{C}$',
        ],
        'par' :
        [
            't0', 'x0', 'x1', 'c',
        ],
        'bestv':
        {
            't0' : 0,
            'x0' : 1000.,
            'x1' : 0,
            'c' : 10,                    
        },
        'bounds':
        {
            't0' : [-20, 20],
            'x0' : [0, 10000],
            'x1' : [-10000, 10000],
            'c'  : [0, 10000],                
        },
    },    
}
