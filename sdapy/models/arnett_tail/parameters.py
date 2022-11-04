modelmeta = {    
    'arnett_fit_taum':
    {
        'engine': 'bol_main',
        'alias':
        [           
            'arnett1', 'arnett',
        ],
        'description': 'Arnett fit, free parameters: mni, taum',
        'func': 'functions.Arnett_fit_taum',        
        'parname' :
        [
            r'$M_\mathrm{Ni}$', r'$\tau_{m}$',
        ],
        'par' :
        [
            'mni', 'taum',
        ],
        'bestv':
        {
            'mni' : .2,
            'taum' : 10,
        },
        'bounds':
        {
            'mni' : [0, 5],
            'taum' : [0, 100],
        },
    },
    'arnett_fit_taum_texp':
    {
        'engine': 'bol_main',
        'alias':
        [
            'arnett2',
        ],
        'description': 'Arnett fit, free parameters: mni, taum, texp',
        'func': 'functions.Arnett_fit_taum',        
        'parname' :
        [
            r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'$t_\mathrm{fl}$',
        ],
        'par' :
        [
            'mni', 'taum', 'texp',
        ],
        'bestv':
        {
            'mni' : .2,
            'taum' : 10,
            'texp' : -18,
        },
        'bounds':
        {
            'mni' : [0, 5],
            'taum' : [0, 100],
            'texp' : [-100, 0],
        },
    },
    'arnett_fit_Mej_Ek':
    {
        'engine': 'bol_main',
        'alias':
        [           
            'arnett3',
        ],
        'description': 'Arnett fit, free parameters: fni, mej, ek',
        'func': 'functions.Arnett_fit_Mej_Ek',        
        'parname' :
        [
            r'$f_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$',
        ],
        'par' :
        [
            'fni', 'ek', 'mej',
        ],
        'bestv':
        {
            'fni' : .2,
            'ek' : 1,
            'mej' : 1,
        },
        'bounds':
        {
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],
        },
    },
    'arnett_fit_Mej_Ek_texp':
    {
        'engine': 'bol_main',
        'alias':
        [           
            'arnett4',
        ],
        'description': 'Arnett fit, free parameters: fni, mej, ek, texp',
        'func': 'functions.Arnett_fit_Mej_Ek',        
        'parname':
        [
            r'$f_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'$t_\mathrm{fl}$',
        ],
        'par' :
        [
            'fni', 'ek', 'mej', 'texp',
        ],
        'bestv':
        {
            'fni' : .2,
            'ek' : 1,
            'mej' : 1,
            'texp' : -18,
        },
        'bounds':
        {
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],
            'texp' : [-100, 0],
        },
    },
    'arnett_fit_Mej_Ek_texp_kopt':
    {
        'engine': 'bol_main',
        'alias':
        [           
            'arnett5',
        ],
        'description': 'Arnett fit, free parameters: fni, mej, ek, texp, kopt',
        'func': 'functions.Arnett_fit_Mej_Ek',        
        'parname':
        [
            r'$f_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'$t_\mathrm{fl}$', r'$\kappa$',
        ],
        'par' :
        [
            'fni', 'ek', 'mej', 'texp', 'k_opt',
        ],
        'bestv':
        {
            'fni' : .2,
            'ek' : 1,
            'mej' : 1,
            'texp' : -18,
            'k_opt': 0.07,
        },
        'bounds':
        {
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],
            'texp' : [-100, 0],
            'k_opt': [0.05, 2],
        },
    }, 
    'tail_fit_t0':
    {
        'engine': 'bol_tail',
        'alias':
        [            
            'tail1',
            'tail',
        ],
        'description': 'Tail fit, free parameters: mni, t0',
        'func': 'functions.tail_fit_t0',        
        'parname':
        [
            r'$M_\mathrm{Ni}$', r'$t_{0}$',
        ],
        'par' :
        [
            'mni', 't0',
        ],
        'bestv':
        {
            'mni' : .2,
            't0' : 100,            
        },
        'bounds':
        {
            'mni' : [0, 5],
            't0' : [0, 10000],
        },
    },        
    'tail_fit_Mej_Ek':
    {
        'engine': 'bol_tail',
        'alias':
        [           
            'tail3',
        ],
        'description': 'Tail fit, free parameters: fni, mej, ek',
        'func': 'functions.tail_fit_Mej_Ek',       
        'parname':
        [
            r'$f_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$',
        ],
       'par' :
        [
            'fni', 'ek', 'mej',
        ],
        'bestv':
        {
            'fni' : .2,
            'ek' : 1,
            'mej' : 1,
        },
        'bounds':
        {
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],
        },
    },        
    'joint_fit_taum_t0':
    {
        'engine': 'bol_full',
        'alias':
        [            
            'arnett_tail',
            'joint',
            'arnett_tail1',
            'joint1',
        ],
        'description': 'Arnett+Tail fit',
        'func': 'functions.joint_fit_taum_t0',        
        'parname':
        [
            r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'$t_{0}$', r'$t_{turn}$',
        ],
        'par' :
        [
            'mni', 'taum', 't0', 'ts',
        ],
        'bestv':
        {
            'mni' : .2,
            'taum' : 10,
            't0' : 100,
            'ts' : 60,
        },
        'bounds':
        {
            'mni' : [0, 5],
            'taum' : [0, 100],
            't0' : [0, 10000],
            'ts' : [20, 120],
        },        
    },    
    'joint_fit_taum_t0_texp':
    {
        'engine': 'bol_full',
        'alias':
        [           
            'arnett_tail2',
            'joint2',
        ],
        'description': 'Arnett+Tail fit',
        'func': 'functions.joint_fit_taum_t0',
        'likelihood':
        {
            'nll': 'likelihood2.nll',
            'lnposterior': 'likelihood2.lnposterior',
        },
        'parname':
        [
            r'$M_\mathrm{Ni}$', r'$\tau_{m}$', r'$t_{0}$', r'$t_{turn}$', r'$t_\mathrm{fl}$',
        ],
        'par' : [
            'mni', 'taum', 't0', 'ts', 'texp', 
        ],
        'bestv':
        {
            'mni' : .2,
            'taum' : 10,
            't0' : 100,            
            'ts' : 60,
            'texp' : -18,
        },
        'bounds':
        {
            'mni' : [0, 5],
            'taum' : [0, 100],
            't0' : [0, 10000],            
            'ts' : [20, 120],
            'texp' : [-100, 0],
        },
    },    
    'joint_fit_Mej_Ek':
    {
        'engine': 'bol_full',
        'alias':
        [                
            'arnett_tail3',
            'joint3',
        ],
        'description': 'Arnett+Tail fit',
        'func': 'functions.joint_fit_Mej_Ek',        
        'parname':
        [
            r'$f_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'$t_{turn}$',
        ],
        'par' :
        [
            'fni', 'ek', 'mej', 'ts'
        ],
        'bestv':
        {
            'fni' : .2,
            'ek' : 1,
            'mej' : 1,
            'ts' : 60,
        },
        'bounds':
        {
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],
            'ts' : [20, 120],
        },
    }, 
    'joint_fit_Mej_Ek_texp':
    {
        'engine': 'bol_full',
        'alias':
        [               
            'arnett_tail4',
            'joint4',
        ],
        'description': 'Arnett+Tail fit',
        'func': 'functions.joint_fit_Mej_Ek', 
        'parname':
        [
            r'$f_\mathrm{Ni}$', r'E$_{kin}$', r'M$_{ej}$', r'$t_{turn}$', r'$t_\mathrm{fl}$',
        ],
        'par' :
        [
            'fni', 'ek', 'mej', 'ts', 'texp',
        ],
        'bestv':
        {
            'fni' : .2,
            'ek' : 1,
            'mej' : 1,
            'ts' : 60,
            'texp' : -18,
        },
        'bounds':
        {
            'fni' : [0, 1],
            'ek' : [.1, 100],
            'mej' : [.1, 50],
            'ts' : [20, 120],
            'texp' : [-100, 0],
        },
    },    
}
