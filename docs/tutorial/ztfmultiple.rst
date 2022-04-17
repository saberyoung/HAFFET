.. _multiple:
   
:mod:`ztfmultiple and plotter`
===========================================

| :ref:`Previous <single>`
| :ref:`1. parameters <multiple1>`
| :ref:`2. run <multiple2>`
| :ref:`3. plotter <multiple3>`

On the basis of telescope pointings setup (which is described in the :ref:`previous chapter <kbpoint>`), here in this chapter we will describe the construction and usage of `KOBE telescope` and `observatory`.
     
.. _multiple1:

parameters
----------------------------------------

`KOBE telescope` object contains 2 attributes:

* `pointings` - telescope pointings
* `telname` - telescope name

.. code-block:: bash
		
   >>> from kobe import telescope
   >>> a=telescope()
   >>> a.telname
   >>> a.pointings

.. _multiple2:

run
----------------------------------------

.. _multiple3:

plotter
----------------------------------------

   

| :ref:`Previous <single>`
| :ref:`Top <multiple>`
