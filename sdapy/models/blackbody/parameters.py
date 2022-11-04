modelmeta = {
    'blackbody':
    {
        'engine': 'sed',
        'alias':
        [
            'bb', 'bb_single',
        ],
        'description': 'Blackbody fits',
        'func': 'functions.bbody',        
        'parname' :
        [
            r'$T$', r'$R$',
        ],
        'par' :
        [
            't', 'r',
        ],
        'bestv':
        {
            't' : 10,
            'r' : 2,
        },
        'bounds':
        {
            't' : [0, 100],
            'r' : [0, 100],
        },
    },
    'dilute_blackbody':
    {
        'engine': 'sed',
        'alias':
        [
            'bb_dilute',                        
        ],
        'description': 'diluted Blackbody fits',
        'func': 'functions.bbody_dilute',        
        'parname' :
        [
            r'$T$', r'$R$',
        ],
        'par' :
        [
            't', 'r',
        ],
        'bestv':
        {
            't' : 10,
            'r' : 2,
        },
        'bounds':
        {
            't' : [0,100],
            'r' : [0, 100],
        },
    },
    'double_blackbody':
    {
        'engine': 'sed',
        'alias':
        [           
            'bb2',
        ],
        'description': 'blackbody with 2 conponents',
        'func': 'functions.bbody_double',        
        'parname' :
        [
            r'$T1$', r'$R1$',
            r'$T2$', r'$R2$',
        ],
        'par' :
        [
            't', 'r',
        ],
        'bestv':
        {
            't' : 10,
            'r' : 2,
            't1' : 10,
            'r1' : 2,
        },
        'bounds':
        {
            't' : [0,100],
            'r' : [0, 100],
            't1' : [0,100],
            'r1' : [0, 100],
        },
    }
}
