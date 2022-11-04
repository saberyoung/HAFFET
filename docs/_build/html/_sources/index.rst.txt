*****************************************
Welcome to the HAFFET documentation
*****************************************

Background
==================
The progenitor scenarios of supernovae (SNe) are still open questions, and one approach to diagnose their physical origins is to investigate the bolometric light curves of a large set of SNe, and fit them to theoretical models to estimate their physical parameter distributions.
Such analysis from different studies often use different approaches and codes which makes the comparisons more difficult. A generic code-package to handle light curve fitting for transients with different types, from different surveys, in different cadences, is therefore useful to provide reliable results for comparison. For this purpose, we present `HAFFET`, a data-driven model fitter for transients.


What is HAFFET?
==================
	   
`HAFFET: Hybrid Analytic Flux FittEr for Transients <https://github.com/saberyoung/HAFFET>`_ is an open source Python package to help analyze SN photometric and spectroscopic data. 

The aim of `HAFFET` is to handle observational data for a set of targets, to estimate their physical parameters, and visualize the population of inferred parameters. Therefore, there are two classes defined, i.e. `snobject` is to deal with data and fittings for one specific object, and `snelist` is to organise the overall running for a list of objects. The inheritance scheme of `HAFFET` is shown as a directed flowchart below:

.. image:: static/sdapy.png
   :width: 800

As shown, `HAFFET` provides utilities to:

* download SN data from online sources:

  #. ZTF alert photometry/spectra from Growth marshal/fritz via `ztfquery <https://github.com/MickaelRigault/ztfquery/tree/master/ztfquery>`_ (an account is needed; for ZTF internal collaborators).
  
  #. `ZTF forced photometry services <https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi/>`_ (an account is needed; open public).
  
  #. `ATLAS forced phtometry services <https://fallingstar-data.com/forcedphot/>`_ (an account is needed; open public).
     
  #. `Open Astronomy Catalog  <https://github.com/astrocatalogs/OACAPI>`_ (open public).
     
  #. Private lightcurves/spectra from users (take care of data format).
  
* intepolated multi band lightcurves with:
  
  #. Gaussian Process regression (via `george <https://george.readthedocs.io/en/latest//>`_).
  
  #. fittings to analytic SNe lightcurve models, e.g. `Bazin et al. (2009) <https://ui.adsabs.harvard.edu/abs/2009A%26A...499..653B/abstract>`_, `Villar et al. (2019) <https://iopscience.iop.org/article/10.3847/1538-4357/ab418c>`_, and the `SALT <http://supernovae.in2p3.fr/salt/doku.php>`_ model for SNe Ia.
  
* characterise the first light and rising of SNe with power law fits:

  #. on multi band photometry simultaneously, developed based on `Miller et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...902...47M/abstract>`_ (`Github <https://github.com/adamamiller/ztf_early_Ia_2018>`_).

  #. on different bands seperately.
  
* match epochs of different bands, via: 

  #. binning

  #. GP interpolation

  #. model fittings

  calculate colours, and/or construct the spectral energy distribution (SED) by assuming e.g. a blackbody distribution.
  
* estimate bolometric LCs via:

  #. bolometric corrections defined in `Lyman et al. (2014) <https://academic.oup.com/mnras/article/437/4/3848/1011706>`_ for stripped envelope SNe or SNe II, and `Chen et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022arXiv220202059C>`_ for SLSNe.    
  
  #. `diluted black body fits <https://en.wikipedia.org/wiki/Black_body>`_ on the estimated SEDs constructed from multi-band photometry. The blackbody fits are developed following `superbol <https://github.com/mnicholl/superbol>`_, and the dilution factor are based on `Dessart & Hillier (2005) <https://ui.adsabs.harvard.edu/abs/2005A%26A...439..671D/abstract>`_.

  #. integration of the absolute calibrated spectra.
  
* estimate host galaxy extinction by:

  #. comparing colours to intrinsic colours, e.g. `Stritzinger et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018A&A...609A.135S>`_, `Taddia et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015A%26A...574A..60T/abstract>`_, etc.

  #. the `pEW of Na ID line doublets <https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.1465P>`_.
  
* fit the constructed bolometric lightcurves to different models, e.g.:

  #. the `Arnett models <https://ui.adsabs.harvard.edu/abs/1982ApJ...253..785A>`_ for SNe Ia and core collapse during their main peaks.

  #. the `gamma ray leakage tail model <https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.3941W>`_ for SNe Ia and core collapse at tail phases.

  #. the `shock cooling emission model <https://arxiv.org/pdf/2007.08543.pdf>`_ for SNe IIb or some Ibc that have early shock cooling tails.

* identify and fit the absorption minima of spectral lines with, e.g.:
  
  #. `Gaussian <https://en.wikipedia.org/wiki/Normal_distribution>`_ function.

  #. `Viglot <https://en.wikipedia.org/wiki/Voigt_profile>`_ function.

* fit the spectral line velocity evolution with, e.g.:

  #. `Exponential <https://en.wikipedia.org/wiki/Exponential_function>`_ function.
  
* scatter all above features into parameter spaces for sample exploration.

It should be noted that besides above build-in models, we provide possibility for users to add their own models, or import models from other python packages (e.g. `MOSFIT <https://mosfit.readthedocs.io/en/latest/>`_, `redback <https://redback.readthedocs.io/en/latest/>`_, etc) into `HAFFET`, to fit on data prepared by our build-in engines. For instance, for those transients that are unlikely to be explained by the Arnett model, one could try fit their bolometric lightcurves with the magnetar/CSM models for the SLSNe, kilonovae/afterglow models for the fast transients, etc. Some demo codes show how to implement these can be found :ref:`here <esternalmodels>`.

In `HAFFET`, there are two approaches implemented as hyperparameter optimization routine:

* scipy.optimize (https://docs.scipy.org/doc/scipy/reference/optimize.html), which is fast.
  
* emcee (https://emcee.readthedocs.io/en/stable/), which is reliable.

Documentation
==============================

.. toctree::
   :maxdepth: 1

   install
   tutorial   
   reference
   version
   faq
   license

API
====================

.. toctree::
   :maxdepth: 1
   
   snelist
   snobject   
   models   
   
Indices and tables
====================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Contact
===================================

If you need additional help, you're welcome to join `our slack channel <https://join.slack.com/t/haffets/shared_invite/zt-1f9qkkapx-DDKdv3vhAZe9fgdY8evWnw>`_, or report any issues `here <https://github.com/saberyoung/sn_data_analysis/issues>`_, or drop us an `email <saberyoung@gmail.com>`_.
