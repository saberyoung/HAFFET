modelmeta = {
    'tail_fit_t0':
    {
        'engine': 'bol_tail',
        'alias': [
            'tail_fit_t0',
            'tail1',
            'tail',
        ],
        'description': 'Tail fit, free parameters: mni, t0',
        'func': 'functions.tail_fit_t0',        
        'parname': [
            r'$M_\mathrm{Ni}$', r'$t_{0}$',
        ],
        'par' : [
            'mni', 't0',
        ],
    },    
    'tail_fit_t0_ts':
    {
        'engine': 'bol_tail',
        'alias': [
            'tail_fit_t0_ts',
            'tail2',
        ],
        'description': 'Tail fit, free parameters: mni, t0, ts',
        'func': 'functions.tail_fit_t0_ts',        
        'parname': [
            r'$M_\mathrm{Ni}$', r'$t_{0}$', r'$t_{s}$',
        ],
        'par' : [
            'mni', 't0', 'ts',
        ],
    },      
    'tail_fit_Mej_Ek':
    {
        'engine': 'bol_tail',
        'alias': [
            'tail_fit_mek_ek',
            'tail3',
        ],
        'description': 'Tail fit, free parameters: mni, mej, ek',
        'func': 'functions.tail_fit_Mej_Ek',       
        'parname': [
            r'$M_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$',
        ],
       'par' : [
            'mni', 'ek', 'mej',
        ],
    },    
    'tail_fit_Mej_Ek_ts':
    {
        'engine': 'bol_tail',
        'alias': [
            'tail_fit_mek_ek_ts',
            'tail4',
        ],
        'description': 'Tail fit, free parameters: mni, mej, ek, ts',
        'func': 'functions.tail_fit_Mej_Ek_ts',        
        'parnames': [
            r'$M_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'$t_{s}$',
        ],
        'par' : [
            'mni', 'ek', 'mej', 'ts',
        ],
    },    
}
