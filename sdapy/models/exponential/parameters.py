modelmeta = {
    'exponential':
    {
        'engine': 'specv_evolution',
        'alias':
        [
            'exp', 'exponential',
        ],
        'description': 'exponential fits',
        'func': 'functions.exp',
        'parname':
        [
            r'$\mathrm{A}$', r"$\mathrm{B}$", r'$\mathrm{C}$',
        ],
        'par' :
        [
            'a', 'b', 'c',
        ],
        'bestv':
        {
            'a' : 10,           
            'b' : .1,
            'c' : 1,                    
        },
        'bounds':
        {
            'a' : [0, 100],            
            'b' : [0, 10],
            'c' : [0, 10],                
        },
    },    
}
