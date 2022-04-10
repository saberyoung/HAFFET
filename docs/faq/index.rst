Frequently Asked Questions
============================

``ImportError: /lib64/libstdc++.so.6``
---------------------------------------

- reason: the Gcc dynamic library version is too old.

- how to solve::
  #edit bash
  LD_LIBRARY_PATH=/home/feng/anaconda3/lib:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH

``ERROR: setuptools 1.0 or later is required by astropy-helpers``
-------------------------------------------------------------------

- reason: when installing astroquery, setuptools version is low

- how to solve::

   >>> import setuptools
   >>> setuptools.__path__

   see where setuptools is called, and then upgrade it: conda install setuptools/pip install setuptools --upgrade/etc
