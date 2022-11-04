modelmeta = {   
    'villar':
    {
        'engine': 'multiband_main',
        'alias':
        [
            'villar',            
        ],
        'description': 'Villar fits',
        'func': 'functions.villar',        
        'parname' :
        [
            r'$A$', r'$B$', r'$t_\mathrm{s}$', r'$t_\mathrm{e}$',
            r'$\tau_\mathrm{fall}$', r'$\tau_{rise}$', r'$\tau_{C}$',
        ],
        'par' : [            
            'a', 'b', 'ts', 'te', 'tfall', 'trise', 'c',
        ],
        'bestv':
        {
            'a' : 100.,
            'b' : -10.,
            'ts' : 30,
            'te' : 30,
            'tfall' : 40,
            'trise' : 20,
            'c' : 0, 
        },
        'bounds':
        {
            'a' : [0,5000],
            'b' : [-1000, 0],
            'ts' : [-100, 100],
            'te' : [-100, 100],
            'tfall' : [0, 120],
            'trise' : [0, 60],
            'c' :  [0,5000],
        },
    },
}
