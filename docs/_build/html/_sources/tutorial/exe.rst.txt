.. _exe:
   
:mod:`The executable file`
===========================================

| :ref:`Previous: The Graphical user interface <gui>`
| :ref:`Next: Python package, snerun <package>`
| :ref:`Back: Getting started <getstart>`

There's an executable file available in the bin directory, which can be called with:

.. code-block:: bash

   >>> sdapy_run -h
   welcome to ->
   ██╗  ██╗ █████╗ ███████╗███████╗███████╗████████╗
   ██║  ██║██╔══██╗██╔════╝██╔════╝██╔════╝╚══██╔══╝
   ███████║███████║█████╗  █████╗  █████╗     ██║
   ██╔══██║██╔══██║██╔══╝  ██╔══╝  ██╔══╝     ██║
   ██║  ██║██║  ██║██║     ██║     ███████╗   ██║
   ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝     ╚══════╝   ╚═╝

   usage: sdapy_run [-h] [-s SOURCE] [--syntax SYNTAX] [-z ZTFID] [-i IAUID]
                 [-t SNTYPE] [-v] [-c]

   run analysis on a list of SNe

   optional arguments:
   -h, --help            show this help message and exit
   -s SOURCE, --source SOURCE
                        meta source (default: BTS)
   --syntax SYNTAX       meta syntax (default: type in ["Ib","Ic"])
   -z ZTFID, --ztfid ZTFID
                        specify ztfid of sn (default: None)
   -i IAUID, --iauid IAUID
                        specify iauid of sn (default: None)
   -t SNTYPE, --sntype SNTYPE
                        specify type of sn (default: None)
   -v, --verbose         Enable progress report (default: False)
   -c, --clobber         Redo analysis (default: False)     
   -d, --debug           If failed fittings, ignore or crashed (default: False)
   --checkdir            Check data directory (default: False)
   
The above command will prompt the user to define a list of candidates to run.
More detailed options can be set via the **default_par.txt** (for general) and **individual_par.txt** (for individual objects), in the data directory.

| :ref:`Top <exe>`
