modelmeta = {
    'magnetar_spindown':
    {
        'engine': 'bol_main',
        'alias': [],
        'description': 'magnetar spin donw',
        'func': 'functions.basic_magnetar',  
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
    'magnetar1':
    {
        'engine': 'bol_main',
        'alias':
        [
            'magnetar',
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
            r'$M_{\mathrm{ej}}\ [M_{\odot}]$',
            r'$E_{\mathrm{kin}}\ [foe]$',
        ],
        'par' :
        [
            'p0', 'bp', 'mass_ns', 'theta_pb', 'mej', 'ek',
        ],
        'bestv':
        {
            'p0': 5,
            'bp': 1,
            'mass_ns': 1.4,
            'theta_pb': 1,
            'mej': 1,
            'ek': 1,           
        },
        'bounds':
        {
            'p0': [1,10],
            'bp': [.1,10],
            'mass_ns':  [1.1,2.2],
            'theta_pb': [0, 3.14/2],
            'mej': [.1, 100],
            'ek': [.1, 100],           
        },
    },
    'magnetar1_texp':
    {
        'engine': 'bol_main',
        'alias':
        [
            'magnetar',
            'basic_magnetar_powered_bolometric_texp'
        ],
        'description': 'magnetar model fit',
        'func': 'functions.basic_magnetar_powered_bolometric',  
        'parname':
        [            
            r'$P_{0}\ [ms]$',
            r'$B_{p}\ [10^{14}G$]',
            r'$M_{\mathrm{NS}} [M_{\odot}]$',
            r'$\theta_{P-B}$',
            r'$M_{\mathrm{ej}}\ [M_{\odot}]$',
            r'$E_{\mathrm{kin}}\ [foe]$',
            r'$t_\mathrm{fl} [d]$',
        ],
        'par' :
        [
            'p0', 'bp', 'mass_ns', 'theta_pb', 'mej', 'ek', 'texp',
        ],
        'bestv':
        {
            'p0': 5,
            'bp': 1,
            'mass_ns': 1.4,
            'theta_pb': 1,
            'mej': 1,
            'ek': 1,
            'texp': -18,
        },
        'bounds':
        {
            'p0': [1,10],
            'bp': [.1,10],
            'mass_ns':  [1.1,2.2],
            'theta_pb': [0, 3.14/2],
            'mej': [.1, 100],
            'ek': [.1, 100],
            'texp' : [-100, 0],
        },
    },
    'magnetar1_texp_opacity':
    {
        'engine': 'bol_main',
        'alias':
        [
            'magnetar',
            'basic_magnetar_powered_bolometric_texp_opacity'
        ],
        'description': 'magnetar model fit',
        'func': 'functions.basic_magnetar_powered_bolometric',  
        'parname':
        [            
            r'$P_{0}\ [ms]$',
            r'$B_{p}\ [10$^{14}G]$',
            r'$M_{\mathrm{NS}} [M_{\odot}]$',
            r'$\theta_{P-B}$',
            r'$M_{\mathrm{ej}}\ [M_{\odot}]$',
            r'$E_{\mathrm{kin}}\ [foe]$',
            r'$t_\mathrm{fl} [d]$',
            r'$\kappa$',
            r'$\kappa_{\gamma}$'
        ],
        'par' :
        [
            'p0', 'bp', 'mass_ns', 'theta_pb', 'mej', 'ek', 'texp', 'k_opt', 'k_gamma',
        ],
        'bestv':
        {
            'p0': 5,
            'bp': 1,
            'mass_ns': 1.4,
            'theta_pb': 1,
            'mej': 1,
            'ek': 1,
            'texp' :-18,
            'k_opt': 0.07,
            'k_gamma': 0.027,
        },
        'bounds':
        {
            'p0': [1,10],
            'bp': [.1,10],
            'mass_ns':  [1.1,2.2],
            'theta_pb': [0, 3.14/2],
            'mej': [.1, 100],
            'ek': [.1, 100],
            'texp' : [-100, 0],
            'k_opt': [0.05, 2],
            'k_gamma' : [1e-4, 10],
        },
    },    
}
