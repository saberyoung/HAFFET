# sn_data_analysis: codes to play with SN alike LCs and spectra

**tutorial (not finished yet)**
https://ztfdataanalysis.readthedocs.io/en/latest/

functions
=============
* Parse local ZTF forced/alert photometry, spectra, etc via ztfquery (https://github.com/MickaelRigault/ztfquery/tree/master/ztfquery)
* Take meta data from BTS catalog (https://sites.astro.caltech.edu/ztf/bts/bts.php)
* Power law fits on early LCs in multi bands simultaneously (developed based on https://github.com/adamamiller/ztf_early_Ia_2018)
* Gaussian Process intepolation (with https://george.readthedocs.io/en/latest/)
* SN alike LC fittings with analytic models, e.g. Bazin et al, Villar et al, etc
* Build colour curve and bolometric LCs, via bolometric corrections defined in Lyman et al 2014, or diluted black body fits if there're at least 3 band fluxes
* Host galaxy extinction estimation by comparing to intrinsic colours, or from the pEW of Na ID line doublets
* Bolometric LC fitting with Arnett models around peak phase, and gamma ray leakage at tail
* Gaussian/viglot fitting for spectral lines for photospheric velocity, which is used to break the Arnett degenracy

hyperparameter optimization routine
=============
Two approaches are implemented:
* minimize: scipy.optimize (https://docs.scipy.org/doc/scipy/reference/optimize.html)
* monte carlo: emcee (https://emcee.readthedocs.io/en/stable/)

structure
=============

data
--
directory for data and cached files

notebook
--
provide a few examples

src
--
* ztfanalysis.py:  
   main program with 2 classes defined: 
   - *ztfsingle* deal with one ZTF SN. Parameters can be redefined in data/priors.txt
   - *ztfmultiple* read meta, deal with a sample of SNe in order
* model_fitters.py: define fitting codes
* gaussian_process.py: define GP codes
* filters.py: define filter information and matplotlib options for different bands
* functions.py: define generic functions
* models.py: define functions, likelihoods for various of models

usage
=============
- data firectory is too large for github, find them here: https://stockholmuniversity.box.com/s/2c3z8yrgvd9zumm4c35u8jtvapyva3e1

- dowload data for your SN, and organize the structure as mine. Carefully define the *datadir* parameter.

- you can find examples in notebook/

- or do: python ztfanalysis.py -h, there're a few infos in helper.

Q&A
=============
Let me (sheng.yang@astro.su.se) know if you had any questions
