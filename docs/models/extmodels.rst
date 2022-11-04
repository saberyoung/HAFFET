.. _esternalmodels:

Add your own models
----------------------------------------

Besides those build-in models, users can use their own models as well. We show an example below, that is a linear model fittting for the tail of SNe multi-band lightcurves. One can add models directly to the exsisting models, e.g. `ZTFDATA`/models/polynomial/ directory. Another option is to make a directory/model, e.g. `test`, under `ZTFDATA`/models/, and registered it in the `ZTFDATA`/models/__init__.py with:

.. code-block:: bash

   from .test import *

Then inside the directory, one can define analytic models (in `ZTFDATA`/models/test/functions.py), e.g.

.. code-block:: bash

   def linear(x, a, b):
      return a*x+b

and claim model name, :ref:`engine <engines>`, parameter descriptions, etc in the `ZTFDATA`/models/test/parameters.py, e.g.

.. code-block:: bash
		
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
	 'same_parameter': ['a'],
	 'fit_error' : True,
      }
   }

Here, the **same_parameter** item decides whether to fit multiple band LCs simultaneously.  If the **same_parameter** is left blank, HAFFET will fit the LC tails in different bands one by one. If some parameters are specified, such as **a** for example, HAFFET fit all bands with a same **a**, which is the slope of tail.

Demo codes can be downloaded :download:`here <../reference/poly.zip>`

Import external models
----------------------------------------

Another option is to import models from other open public python packages. As an example, we introduce magnetar model from `redback`, which can be then used to fit on bolometric lightcurves of SLSNe. After installing `redback` and registering it to `HAFFET` correctly, one can import functions at `ZTFDATA`/models/test2/functions.py as:

.. code-block:: bash

   from redback.transient_models.supernova_models import basic_magnetar_powered_bolometric

and claim model settings as:

.. code-block:: bash

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

Demo codes can be downloaded :download:`here <../reference/redback_magnetar.zip>`
