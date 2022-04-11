*****************************************
Welcome to the sdapy documentation
*****************************************

What is sdapy?
==================
	   
`SN Data Analysis python package (sdapy) <https://sngyang.com/sdapy>`_ is an open source Python package to help analyze SN photometric and spectroscopic data.

`sdapy` provides utilities to:

* parse meta data from `BTS page <https://sites.astro.caltech.edu/ztf/bts/bts.php>`_.
* download/restore ZTF alert/forced photometry and spectra (via `ztfquery <https://github.com/MickaelRigault/ztfquery/tree/master/ztfquery>`_), `ATLAS forced phtometry <https://fallingstar-data.com/forcedphot/>`_, and the `TNS <https://www.wis-tns.org/>`_ spectra.
* intepolated multi band photometry via Gaussian Process (with `george <https://george.readthedocs.io/en/latest//>`_).
* intepolated each band with SN alike analytic models, e.g. Bazin et al, Villar et al.
* power law fits on early LCs in multi bands simultaneously, to determine first light (developed based on `<https://github.com/adamamiller/ztf_early_Ia_2018>`_).
* build colour curves, and bolometric LCs, via bolometric corrections defined in Lyman et al 2014, or diluted black body fits.
* estimate host galaxy extinction by comparing to intrinsic colours (e.g. Drout et al), or from the pEW of Na ID line doublets.
* Bolometric LC fitting with the Arnett models, and gamma ray leakage tail.
* Gaussian/viglot fitting for spectral lines for photospheric velocity.
* shock cooling fitting with Piro 2020.

Documentation
==============================

.. toctree::
   :maxdepth: 1

   install
   tutorial
   issues
   todo
   faq/index
   license
	
Reference/API
====================
     
.. toctree::
   :maxdepth: 1
   
   ztfanalysis   
   
Indices and tables
====================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
  
Links
======

* `Source code <https://github.com/saberyoung/sn_data_analysis>`_
* `Docs <https://sdapy.readthedocs.io/>`_
* `Issues <https://github.com/saberyoung/sn_data_analysis/issues>`_

Author
====================

`Sheng Yang <http://www.sngyang.com>`_: saberyoung@gmail.com
