Installation
===================================

Requirements
-------------

`HAFFET` depends on several python libraries, e.g.
`matplotlib <https://matplotlib.org/>`_, `astropy <https://www.astropy.org/>`_, etc
(see `full list <https://github.com/saberyoung/sn_data_analysis/blob/master/requirements.txt>`_),
and all of them could be installed via `pip`, e.g.::

  pip install -r requirements.txt

.. note::
   * **requirements.txt** can be found at https://github.com/saberyoung/HAFFET/blob/master/requirements.txt.
   * Linux and MAC OS systems have been tested, and for Windows users, there's a `docker version <https://github.com/saberyoung/sn_data_analysis/tree/Docker>`_ in preparation.
   * `HAFFET` is tested under Python 3.

Source installation with Pypi (RECOMMENDED)
------------------------------------------------------------------

It is possible to build the latest ``HAFFET`` with `pip <http://www.pip-installer.org>`_ ::

    pip install --user haffet

If you have installed with ``pip``, you can keep your installation up to date
by upgrading from time to time::

    pip install --user --upgrade haffet

Make sure your ``pip`` is up to date, if it cannot find ``HAFFET``, upgrate your ``pip``::
  
    python -m pip install -U pip

Almost-as-quick installation from official source release
----------------------------------------------------------

``HAFFET`` is also available in the
`Github <https://github.com/saberyoung/HAFFET>`_. You can
download and build it with::

    python setup.py install --user    

Check
-----

If everything goes well, you will have it::

    python
    >>> import sdapy
    >>> sdapy.__version__
    '0.1.0'    

Clean
-----

When you run "python setup.py", temporary build products are placed in the
"build" directory. If you want to clean out and remove the ``build`` directory,
then run::

    python setup.py clean --all

Uninstall
-----------

For uninstallation, one can easily delete the files directly.
In order to know the file path, you should start python in correct environment and do::

>>> import sdapy
>>> sdapy.__path__

Another approach is to remove it via ``pip``::
  
    pip uninstall haffet
