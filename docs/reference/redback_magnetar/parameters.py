modelmeta = {   
    'magnetar':
    {
        'engine': 'bol_main',
        'alias':
        [
            'redback_magnetar',
            'basic_magnetar_powered_bolometric'
        ],
        'description': 'magnetar model fit',
        'func': 'functions.basic_magnetar_powered_bolometric',  
        'parname':
        [            
            r'$P_{0}\ [ms]$',
            r'$B_{p}\ [10^{14}G$]',
            r'$M_{\mathrm{NS}} [M_{\odot}]$',
            r'$\theta_{P-B}$',            
        ],
        'par' :
        [
            'p0', 'bp', 'mass_ns', 'theta_pb',
        ],
        'bestv':
        {
            'p0': 5,
            'bp': 1,
            'mass_ns': 1.4,
            'theta_pb': 1,                     
        },
        'bounds':
        {
            'p0': [1,10],
            'bp': [.1,10],
            'mass_ns':  [1.1,2.2],
            'theta_pb': [0, 3.14/2],                    
        },
    },
}
