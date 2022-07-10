modelmeta = {
    'joint_fit_taum_t0':
    {
        'engine': 'bol_full',
        'alias': [
            'joint_fit_taum_t0',
            'arnett_tail',
            'joint',
            'arnett_tail1',
            'joint1',
        ],
        'description': 'Arnett+Tail fit',
        'func': 'functions.joint_fit_taum_t0',        
        'parname': [
            r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'$t_{0}$', r'$t_{turn}$',
        ],
        'par' : [
            'mni', 'taum', 't0', 'ts',
        ],
    },    
    'joint_fit_taum_t0_texp':
    {
        'engine': 'bol_full',
        'alias': [
            'joint_fit_taum_t0_texp',
            'arnett_tail2',
            'joint2',
        ],
        'description': 'Arnett+Tail fit',
        'func': 'functions.joint_fit_taum_t0_texp',
        'likelihood':
        {
            'nll': 'likelihood2.nll',
            'lnposterior': 'likelihood2.lnposterior',
        },
        'parname': [
            r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'$t_{0}$', r'$t_\mathrm{fl}$', r'$t_{turn}$',
        ],
        'par' : [
            'mni', 'taum', 't0', 'texp', 'ts',
        ],
    },    
    'joint_fit_Mej_Ek':
    {
        'engine': 'bol_full',
        'alias': [
            'joint_fit_Mej_Ek',         
            'arnett_tail3',
            'joint3',
        ],
        'description': 'Arnett+Tail fit',
        'func': 'functions.joint_fit_Mej_Ek',        
        'parname': [
            r'$M_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'$t_{turn}$',
        ],
        'par' : [
            'mni', 'ek', 'mej', 'ts'
        ],
    }, 
    'joint_fit_Mej_Ek_texp':
    {
        'engine': 'bol_full',
        'alias': [
            'joint_fit_Mej_Ek_texp',         
            'arnett_tail4',
            'joint4',
        ],
        'description': 'Arnett+Tail fit',
        'func': 'functions.joint_fit_Mej_Ek_texp',       
        'parname': [
            r'$M_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'$t_\mathrm{fl}$', r'$t_{turn}$',
        ],
        'par' : [
            'mni', 'ek', 'mej', 'texp', 'ts',
        ],
    },    
}
