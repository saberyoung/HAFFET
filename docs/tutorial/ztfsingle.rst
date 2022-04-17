.. _single:
   
:mod:`ztfsingle: for specific object`
===========================================

| :ref:`Next <multiple>`
| :ref:`1. input photometric and spectral data <single1>`
| :ref:`2. interpolation, define t0 as peak epoch <single2>`
| :ref:`3. pre peak fit, define t0 as explosion epoch <single3>`
| :ref:`4. build bolometric LCs <single4>`
| :ref:`5. fit bolometric LCs with analytic models <single5>`
| :ref:`6. photospheric velocities <single6>`
	   
.. _single1:

1- input data
----------------------------------------

The first step is to load observational data. In `sdapy`, we use `pandas.DataFrame <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html>`_, and asked at least 3 columns ('mag', 'emag' and 'jdobs') to be setted.

An easiest way to input photometric data is to manually add it:

.. code-block:: bash

   >>> # import module
   >>> from sdapy.ztfanalysis import *

   >>> # kwargs description can be found at
   >>> # https://github.com/saberyoung/sn_data_analysis/blob/a4db8a221eb46bc3141f481665f08b71c4d70f23/sdapy/ztfanalysis.py#L237
   >>> # and we will detailly describe them later in the tutorial
   >>> ztfp = ztfsingle(ztfid, iauid=None, z=None, ra=None, dec=None,
   >>> ...             mkwebv=None, sntype=None, dm=None, jdpeak=None,
   >>> ...             logger=None, axes=None, **kwargs**kwargs)
   
   >>> # make a pandas dataframe from a csv file
   >>> import pandas as pd
   >>> df = pd.readcsv(...)
   
   >>> # make sure there're 'mag', 'emag' and 'jdobs' columns in the dataframe
   >>> ztfp.lc = df
   
   >>> # or use add_lc, that will check the column names
   >>> ztfp.add_lc(df)
   
For specific sources, there're built-in functions that can help download and load local data. `ztfquery` is used to deal with ZTF Growth marshal/fritz data, thus we decide to follow their data structure for other data, e.g. ATLAS, TNS, as well. As shown below, the `targetdir` is a defined parameter, and one should define it as an existing folder, construct it with the following strcuture. You can also find my data as example at https://stockholmuniversity.box.com/s/2c3z8yrgvd9zumm4c35u8jtvapyva3e1

.. code-block:: bash

   `targetdir`/
      marshal/
         lightcurves/
	    `ztfid`/
	       marshal_plot_lc_lightcurve_`ztfid`.csv 
	 spectra/
	    `ztfid`/
	       `ztfid`_***.ascii
      fritz/
         lightcurves/
	    `ztfid`/
	       fritz_`ztfid`_lc.csv
	 sample/
	    fritz_groups.json
	 spectra/
	    `ztfid`/
	       fritz_spectrum_`instrument`_spec_auto_***.ascii
      ForcePhot/
         `ztfid`/
	    forcedphotometry_`ztfid`_lc.csv	    
      ForcePhot_atlas/
         `ztfid`/
	    forcedphotometry_`ztfid`_lc.csv
	    forcedphotometry_`ztfid`_lc_atlas_fp_stacked_`bindays`_days.txt
      TNS/
         `ztfid`/
	    `ztfid`_***.txt

After define and build up the `targetdir` correctly, the `sdapy` is possible to be used to download and reload ZTF/ATLAS/TNS data. We show examples below (they're also avaliable in a jupyter-notebook at https://github.com/saberyoung/sn_data_analysis/blob/master/notebook/loaddata.ipynb):

.. code-block:: bash

   >>> # import module
   >>> from sdapy.ztfanalysis import *

   >>> # initialize a ztfsingle class
   >>> # ra, dec are used if one need to query ZTF/ATLAS forced photometry
   >>> # ztfid is needed by marshal/fritz pjotometry
   >>> # one can set them, e.g. ra, to None, and might have assertion errors
   >>> #    when doing some specific tasks
   >>> ztfp = ztfsingle(ztfid='ZTF20aajcdad', iauid='SN2020bcq', 
   >>> ...     z=0.0186, dm=34.6, mkwebv=0.01387, sntype='SN Ib', 
   >>> ...     ra='13:26:29.65', dec='+36:00:31.1', jdpeak=None, 
   >>> ...     logger=None, axes=None)

   >>> # make a pandas dataframe
   >>> import pandas as pd
   >>> df = pd.readcsv(...)
   
   >>> # make sure there're 'mag', 'emag' and 'jdobs' columns in the dataframe
   >>> # and then give it to ztfsingle
   >>> ztfp.add_lc(df)

.. _single2:

2- interpolation
----------------------------------------


.. _single3:

3- power law fit
----------------------------------------


.. _single4:

4- build bolometric LCs
----------------------------------------


.. _single5:

5- compare bolometric LCs to analytic models
------------------------------------------


.. _single6:

6- photospheric velocities
----------------------------------------

| :ref:`Next <multiple>`
| :ref:`Top <single>`
