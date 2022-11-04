import numpy as np

modelmeta = {
    'poly1':
    {
        'engine': 'multiband_tail',
        'alias':
        [
            'linear', 'polynomial1',
        ],
        'description': '1 order polynomial',
        'func': 'functions.linear',
        'parname' :
        [
            r'$a$', r'$b$',
        ],
        'par' :
        [
            'a', 'b',
        ],
        'bestv':
        {
            'a' : 1,
            'b' : 0.,                             
        },
        'bounds':
        {
            'a' : [-np.inf, np.inf],
            'b' : [-np.inf, np.inf],                         
        },
    },
}
