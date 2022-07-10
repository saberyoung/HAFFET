modelmeta = {   
    'villar':
    {
        'engine': 'multiband_main',
        'alias': [
            'villar',            
        ],
        'description': 'Villar fits',
        'func': 'functions.villar',        
        'parname' : [
            r'$A$', r'$t_\mathrm{0}$', r'$\tau_\mathrm{fall}$',
            r'$\tau_{rise}$', r'$C$',
        ],
        'par' : [            
            'a', 'b', 'ts', 'te', 'tfall', 'trise',
        ],        
    },
}
