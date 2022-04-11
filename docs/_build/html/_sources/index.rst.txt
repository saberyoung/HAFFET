*****************************************
Welcome to the sdapy documentation
*****************************************
Documentation of `sdapy` is hosted by `Read the Docs <http://sdapy.readthedocs.org/en/stable/>`_.

Background
==================
The progenitor scenarios for supernovae (SNe) are still open questions, and one approach to diagnose their physical origins is to investigate the shape and luminosity of bolometric light curves of a large set of SNe, fit to semianalytic or hydrodynamic models to estimate possible physical parameters, and populate these together with what model allowed. There're various of survey groups doing these, however different analysis will bring systematic errors that would make their conclusion suspicious. A generic code is needed to handle fittings for SNe with different type, from different survey, in different casence, etc, and provide reliable results for comparison. As part of the `ZTF collaboration <https://www.ztf.caltech.edu/>`_, I develop `sdapy` to analysis a subset of Ib/Ic SNe classified by `BTS survey <https://sites.astro.caltech.edu/ztf/bts/bts.php>`_, and hope it could be afterwards used by more groups so that we could compare their common properties to understand their real physics behind.


What is sdapy?
==================
	   
`SN Data Analysis python package (sdapy) <https://sngyang.com/sdapy>`_ is an open source Python package to help analyze SN photometric and spectroscopic data. 

`sdapy` provides utilities to:

* parse SN meta data from `BTS page <https://sites.astro.caltech.edu/ztf/bts/bts.php>`_.
* download/restore ZTF alert photometry/spectra from Growth marshal/fritz (via `ztfquery <https://github.com/MickaelRigault/ztfquery/tree/master/ztfquery>`_), `ZTF forced photometry <https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi/>`_ , `ATLAS forced phtometry <https://fallingstar-data.com/forcedphot/>`_, and the `TNS spectra <https://www.wis-tns.org/>`_ (and also other LCs/spectra by manual input them).
* intepolated multi band photometry via Gaussian Process (with `george <https://george.readthedocs.io/en/latest//>`_).
* intepolated each band with SN alike analytic models, e.g. `Bazin et al 2009 <https://ui.adsabs.harvard.edu/abs/2009A%26A...499..653B/abstract>`_, `Villar et al 2019 <https://iopscience.iop.org/article/10.3847/1538-4357/ab418c>`_.
* power law fits on early LCs in multi bands simultaneously, to determine first light (developed based on `<https://github.com/adamamiller/ztf_early_Ia_2018>`_).
* build colour curves, and bolometric LCs, via bolometric corrections defined in `Lyman et al 2014 <https://academic.oup.com/mnras/article/437/4/3848/1011706>`_, or the `diluted black body fits <https://en.wikipedia.org/wiki/Black_body>`_.
* estimate host galaxy extinction by comparing colours to intrinsic colours (e.g. `Stritzinger et al 2018 <https://ui.adsabs.harvard.edu/abs/2018A&A...609A.135S>`_), or from the `pEW of Na ID line doublets <https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.1465P>`_.
* Bolometric LC fitting with the `Arnett models <https://ui.adsabs.harvard.edu/abs/1982ApJ...253..785A>`_, and/or the `gamma ray leakage tail <https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.3941W>`_.
* Gaussian/viglot fitting for absorption minima of He 5876 and O 7772 (for SN Ib and Ic correspondingly), which were used to estimate photospheric velocities, using an analytic function defined by `Dessart et al 2016 <https://academic.oup.com/mnras/article/458/2/1618/2589109>`_.
* shock cooling phase fitting with `Piro et al 2020 <https://arxiv.org/pdf/2007.08543.pdf>`_.

For fittings, there're two approaches implemented as hyperparameter optimization routine:

* scipy.optimize (https://docs.scipy.org/doc/scipy/reference/optimize.html)
* emcee (https://emcee.readthedocs.io/en/stable/)

Documentation
==============================

.. toctree::
   :maxdepth: 1

   install
   tutorial
   gallery
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
