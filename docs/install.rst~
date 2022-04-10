Installation
===================================

Requirements
-------------

`KOBE` depends on several python libraries, e.g. `healpy <https://healpy.readthedocs.io/en/latest/>`_, 
`matplotlib <https://matplotlib.org/>`_, `astropy <https://www.astropy.org/>`_, etc
(see `full list <https://github.com/saberyoung/kobe/blob/master/requirements.txt>`_),
and all of them could be installed via `pip`, e.g.::

  pip install -r requirements.txt

.. note:: 
   * Only Linux and MAC OS have been tested, for Windows everything reamins unknown.
   * Considering the Python version, the tool is working for Python 2/3.

Source installation with Pypi (RECOMMENDED)
--------------------------------------------

It is possible to build the latest ``kobe`` with `pip <http://www.pip-installer.org>`_ ::

    pip install --user kobe

If you have installed with ``pip``, you can keep your installation up to date
by upgrading from time to time::

    pip install --user --upgrade kobe

Almost-as-quick installation from official source release
----------------------------------------------------------

KOBE is also available in the
`Python Package Index (PyPI) <https://pypi.org/project/kobe/>`_. You can
download it with::

    curl -O https://files.pythonhosted.org/packages/1c/d5/42cb34cd80b1049b2f4352f17e9277344a25ee07edc8557d9abf9e963147/kobe-0.0.1.tar.gz

and build it with::

    tar -xzf kobe-0.0.1.tar.gz    
    python setup.py install --user

Bash installation 
---------------------------------------------

You can download the source files from `KOBE repository <https://github.com/saberyoung/kobe.git>`_ with 
the ``git`` command::

    git clone https://github.com/saberyoung/kobe.git
    
which would create one directory in the current path, e.g. ``/home/xxx/kobe/``.
Then, you need to source the shell file::

    source /home/xxx/kobe/kobe.bash

This approach would not install the pipeline, but instead call python via envirnmental defination,
which is light and easy.

Check
-----

If everything goes fine, you can test it::

    python
    >>> import kobe
    >>> kobe.__version__
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

>>> import kobe
>>> kobe.__path__

Another approach is to remove it via pip::

    pip uninstall kobe
