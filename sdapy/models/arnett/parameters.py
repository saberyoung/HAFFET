modelmeta = {
    'arnett_fit_taum':
    {
        'engine': 'bol_main',
        'alias': [
            'arnett_fit_taum',
            'arnett1',
            'arnett',
        ],
        'description': 'Arnett fit, free parameters: mni, taum',
        'func': 'functions.Arnett_fit_taum',        
        'parname' : [
            r'$M_\mathrm{Ni}$', r'$\tau_{m}$',
        ],
        'par' : [
            'mni', 'taum',
        ],        
    },
    'arnett_fit_taum_texp':
    {
        'engine': 'bol_main',
        'alias': [
            'arnett_fit_taum_texp',
            'arnett2',
        ],
        'description': 'Arnett fit, free parameters: mni, taum, texp',
        'func': 'functions.Arnett_fit_taum_texp',        
        'parname' : [
            r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'$t_\mathrm{fl}$',
        ],
        'par' : [
            'mni', 'taum', 'texp',
        ],
    }, 
    'arnett_fit_Mej_Ek':
    {
        'engine': 'bol_main',
        'alias': [
            'arnett_fit_mej_ek',
            'arnett3',
        ],
        'description': 'Arnett fit, free parameters: mni, mej, ek',
        'func': 'functions.Arnett_fit_Mej_Ek',        
        'parname' : [
            r'$M_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$',
        ],
        'par' : [
            'mni', 'ek', 'mej',
        ],
    },
    'arnett_fit_Mej_Ek_texp':
    {
        'engine': 'bol_main',
        'alias': [
            'arnett_fit_mej_ek_texp',
            'arnett4',
        ],
        'description': 'Arnett fit, free parameters: mni, mej, ek, texp',
        'func': 'functions.Arnett_fit_Mej_Ek_texp',        
        'parname': [
            r'$M_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'$t_\mathrm{fl}$',
        ],
        'par' : [
            'mni', 'ek', 'mej', 'texp',
        ],
    },    
}
