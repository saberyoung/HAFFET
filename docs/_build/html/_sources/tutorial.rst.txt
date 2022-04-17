Tutorial
===================================

The aim of `sdapy` is to handle observational data for a set of targets, to estimate their physical parameters, and visualize the population of inferred parameters. Therefore, there're two classes defined, i.e. `ztfsingle` is to deal with one specific object, and `ztfmultiple` is to organise the overall runnings for multiple objects, store their fitting results, which can be visualize with the help of `plotter`. The inheritance scheme of `sdapy` is shown as directed flowchart as followed:

.. image:: static/sdapy.png
   :width: 800

:ref:`1. ztfsingle  <single>`
--------------------------------------------------
| :ref:`1.1 input photometric and spectral data <single1>`
| :ref:`1.2 interpolation, define t0 as peak epoch <single2>`
| :ref:`1.3 pre peak fit, define t0 as explosion epoch <single3>`
| :ref:`1.4 build bolometric LCs <single4>`
| :ref:`1.5 fit bolometric LCs with analytic models <single5>`
| :ref:`1.6 photospheric velocities <single6>`

:ref:`2. ztfmultiple and plotter <multiple>`
----------------------------------------------


Check more detials in the corresponding API page.
