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

ztfquery is used to deal with Growth marshal/fritz data, and we follow their data structure for other data, e.g. ATLAS, TNS, as well. 

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
