modelmeta = {
    'pure_shock':
    {
        'engine': 'bol_early',
        'alias': ['pure_sbo'],
        'description': 'Pure Shock break out fit',
        'func': 'functions.shock_fit',  
        'parname':
        [
            r'$M_\mathrm{ext}$',
            r'$R_\mathrm{ext}$',
            r"$E_\mathrm{ext}$",
        ],
        'par' :
        [
            'me', 're', 'ee',
        ],
        'bestv':
        {
            'me' : .1,
            're' : 200.,
            'ee' : .1,                           
        },
        'bounds':
        {
            'me' : [0, 1],
            're' : [0, 1000],
            'ee' : [1e-3, 10],                         
        },
    },
    'pure_shock_texp':
    {
        'engine': 'bol_early',
        'alias': ['pure_sbo_texp'],
        'description': 'Pure Shock break out fit',
        'func': 'functions.shock_fit',  
        'parname':
        [
            r'$M_\mathrm{ext}$',
            r'$R_\mathrm{ext}$',
            r"$E_\mathrm{ext}$",
            r'$t_\mathrm{fl}$',
        ],
        'par' :
        [
            'me', 're', 'ee', 'texp'
        ],
        'bestv':
        {
            'me' : .1,
            're' : 200.,
            'ee' : .1,
            'texp' : -18,
        },
        'bounds':
        {
            'me' : [0, 1],
            're' : [0, 1000],
            'ee' : [1e-3, 10],
            'texp' : [-100, 0], 
        },
    },
    'shock_fit_taum':
    {
        'engine': 'bol_early',
        'alias':
        [
            'shock_fit1',
            'shock_fit',
            'shock',
            'sbo',
            'shock1',
            'sbo1',
        ],
        'description': 'Shock break out fit',
        'func': 'functions.shock_arnett_fit',  
        'parname':
        [
            r'$M_\mathrm{ext}$', r'$R_\mathrm{ext}$',
            r"$E_\mathrm{ext}$", r'$f_\mathrm{Ni}$',
            r'$\tau_{m}$', r'$t_\mathrm{fl}$',
        ],
        'par' :
        [
            'me', 're', 'ee', 'fni', 'taum', 
        ],
        'bestv':
        {
            'me' : .1,
            're' : 200.,
            'ee' : .1,     
            'fni': .2,
            'taum': 10,                              
        },
        'bounds':
        {
            'me' : [0, 1],
            're' : [0, 1000],
            'ee' : [1e-3, 10],
            'fni' : [0, 1],
            'taum' : [0, 100],              
        },
    },
    'shock_fit_mejek':
    {
        'engine': 'bol_early',
        'alias':
        [
            'shock_fit2',
            'shock2',
            'sbo2',      
        ],
        'description': 'Shock break out fit',
        'func': 'functions.shock_arnett_mejek_fit',  
        'parname':
        [
            r'$M_\mathrm{ext}$', r'$R_\mathrm{ext}$',
            r"$E_\mathrm{ext}$", r"$f_\mathrm{Ni}$",
            r"$M_\mathrm{ej}$", r'$E_\mathrm{kin}$',            
        ],
        'par' :
        [
            'me', 're', 'ee', 'fni', 'mej', 'ek',
        ],
        'bestv':
        {
            'me' : .1,
            're' : 200.,
            'ee' : .1,     
            'fni': .2,
            'ek' : 1,
            'mej' : 1,                
        },
        'bounds':
        {
            'me' : [0, 1],
            're' : [0, 1000],
            'ee' : [1e-3, 10],
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],           
        },
    },
    'shock_fit_taum_texp':
    {
        'engine': 'bol_early',
        'alias':
        [
            'shock_fit3'            
            'shock3',
            'sbo3',
        ],
        'description': 'Shock break out fit',
        'func': 'functions.shock_arnett_fit',  
        'parname':
        [
            r'$M_\mathrm{ext}$', r'$R_\mathrm{ext}$',
            r"$E_\mathrm{ext}$", r'$f_\mathrm{Ni}$',
            r'$\tau_{m}$', r'$t_\mathrm{fl}$',
        ],
        'par' :
        [
            'me', 're', 'ee', 'fni', 'taum', 'texp',
        ],
        'bestv':
        {
            'me' : .1,
            're' : 200.,
            'ee' : .1,     
            'fni': .2,
            'taum': 10,
            'texp' : -18,                    
        },
        'bounds':
        {
            'me' : [0, 1],
            're' : [0, 1000],
            'ee' : [1e-3, 10],
            'fni' : [0, 1],
            'taum' : [0, 100],
            'texp' : [-100, 0],                
        },
    },
    'shock_fit_mejek_texp':
    {
        'engine': 'bol_early',
        'alias':
        [
            'shock_fit4',
            'shock4',
            'sbo4',      
        ],
        'description': 'Shock break out fit',
        'func': 'functions.shock_arnett_mejek_fit',  
        'parname':
        [
            r'$M_\mathrm{ext}$', r'$R_\mathrm{ext}$',
            r"$E_\mathrm{ext}$", r"$f_\mathrm{Ni}$",
            r"$M_\mathrm{ej}$", r'$E_\mathrm{kin}$',
            r'$t_\mathrm{fl}$',
        ],
        'par' :
        [
            'me', 're', 'ee', 'fni', 'mej', 'ek', 'texp',
        ],
        'bestv':
        {
            'me' : .1,
            're' : 200.,
            'ee' : .1,     
            'fni': .2,
            'ek' : 1,
            'mej' : 1,
            'texp' : -18,
        },
        'bounds':
        {
            'me' : [0, 1],
            're' : [0, 1000],
            'ee' : [1e-3, 10],
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],
            'texp' : [-100, 0],  
        },
    },
    'shock_fit_mejek_texp_opacity':
    {
        'engine': 'bol_early',
        'alias':
        [
            'shock_fit5',
            'shock5',
            'sbo5',      
        ],
        'description': 'Shock break out fit',
        'func': 'functions.shock_arnett_mejek_fit',  
        'parname':
        [
            r'$M_\mathrm{ext}$', r'$R_\mathrm{ext}$',
            r"$E_\mathrm{ext}$", r"$f_\mathrm{Ni}$",
            r"$M_\mathrm{ej}$", r'$E_\mathrm{kin}$',
            r'$t_\mathrm{fl}$', r'$\kappa$',
        ],
        'par' :
        [
            'me', 're', 'ee', 'fni', 'mej', 'ek', 'texp', 'k_opt',
        ],
        'bestv':
        {
            'me' : .1,
            're' : 200.,
            'ee' : .1,     
            'fni': .2,
            'ek' : 1,
            'mej' : 1,
            'texp' : -18,
            'k_opt': 0.07,
        },
        'bounds':
        {
            'me' : [0, 1],
            're' : [0, 1000],
            'ee' : [1e-3, 10],
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],
            'texp' : [-100, 0],
            'k_opt': [0.05, 2],
        },
    },
}
