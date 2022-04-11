Installation
===================================

Requirements
-------------

`sdapy` depends on several python libraries, e.g.
`matplotlib <https://matplotlib.org/>`_, `astropy <https://www.astropy.org/>`_, etc
(see `full list <https://github.com/saberyoung/sn_data_analysis/blob/master/requirements.txt>`_),
and all of them could be installed via `pip`, e.g.::

  pip install -r requirements.txt

.. note:: 
   * Only Linux and MAC OS have been tested, for Windows everything reamins unknown.
   * Considering the Python version, the tool is working for Python 3.

Source installation with Pypi (RECOMMENDED but not available yet)
------------------------------------------------------------------

It is possible to build the latest ``sdapy`` with `pip <http://www.pip-installer.org>`_ ::

    pip install --user sdapy

If you have installed with ``pip``, you can keep your installation up to date
by upgrading from time to time::

    pip install --user --upgrade sdapy

Almost-as-quick installation from official source release
----------------------------------------------------------

``sdapy`` is also available in the
`Github <https://github.com/saberyoung/sn_data_analysis>`_. You can
download and build it with::
  
    python setup.py install --user
    

Check
-----

If everything goes fine, you can test it::

    python
    >>> import sdapy
    >>> sdapy.__version__
    '0.0.0.1'
    

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

Another approach is to remove it via pip::

    pip uninstall sdapy
