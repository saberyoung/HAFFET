.. _single:
   
:mod:`For specific object`
===========================================

| :ref:`Next <multiple>`
| :ref:`1. input photometric and spectral data <single1>`
| :ref:`2. interpolation, define t0 as peak epoch <single2>`
| :ref:`3. pre peak fit, define t0 as explosion epoch <single3>`
| :ref:`4. build bolometric LCs <single4>`
| :ref:`5. fit bolometric LCs with analytic models <single5>`
| :ref:`6. photospheric velocities <single6>`
| :ref:`7. driven physical mechanisms <single7>`

.. _single1:

input data
----------------------------------------

The aim of `KOBE` is to define and optimize observing strategies for telescopes.
We mimic different scenarios, while in the process, we realized there're lots of tasks are used multiple times.
In order to make the code clear and simple, we define python classes for major concepts, adopting Python inheritance to make them reusability.

The inheritance scheme of `KOBE` is shown as directed flowchart as followed:

.. image:: ../static/kobesheme.png
   :width: 800

.. note::
   
   It should be awared that functions from the parent classes (e.g. utils), could be used by their son classes (e.g. schedule). In the tutorial, we show examples for the current classes, but they're also applicable to its successors. For example:

   We parsed one `trigger` object from a source:
   
   .. code-block:: bash

      >>> from kobe import trigger
      >>> a=trigger()
      >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')

   Since `trigger` is thus inheriented by `schedule`, `url` function is also available in `schedule`:
   
   .. code-block:: bash

      >>> from kobe import schedule
      >>> a = schedule()       
      >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')

As shown, we define generic functions in `utils`, which would be then inherited by `visualization` and `circulate`.
These 3 classes would provide basic functions for `KOBE` project, assisting with visualization, circulate and other tasks.

Afterwards, we define 3 major classes, i.e. `trigger`, `candidates`, and `pointings`, to mimic the followup search process.

**Search without trigger**:

Before receiving an alert, `KOBE` provide `pointings` class that could help generate pointing lists for specific telescopes, depending on different strategy, i.e. tiling search or galaxy search.
On this basis, a `telescope` object is defined after setting its name.
Finally, user could define several `observatory` objects, while each of them is composed by a location, observatory name, and a list of `telescope` objects.

At this phase, `schedule` could generate observing strategies for telescopes, depending only on the visibility of each pointings (i.e. airmass, galactic plane and bright object constrains).

Meanwhile, user could define a series of candidates via `candidates` and submit them to `schedule`, asking `KOBE` to assess them.

**Search with trigger**:

The contourpart search starts with receiving a trigger from one interesting source.
In `trigger` class, we adopt `Healpy/HEALPix <https://healpy.readthedocs.io/en/latest/index.html>`_ to help parse trigger information from alert files.
These extracted time and space informations were then transferred to `schedule`.
Together with pre-submitted pointing lists and user options, `schedule` would generate OBs for a series of telescopes, try coordinate their observations in order to increase the overall efficiency of EM counterpart searches.

**Evaluation for trigger search**:

In order to assess the schedule, `KOBE` provide routes to calculate the overall detection efficiency, by simulating a series of targets randomly according to trigger, see how many of injections could be finally detected, and the detection efficiency is defines as the ratio.

.. _single2:

Interpolation
----------------------------------------

Parameter is one basic element of programming that allow interaction of the program with external entities.
Before use of `KOBE`, we would like to describe how `KOBE` would read and handle with parameters:

1. In most cases, `KOBE` put all parameters into optional parameters and provide them default values.
   If the result is not as what you expected, read carefully its documentation, set as request, and redo the tasks.

2. Some parameters would be used by multiple classes and functions.
   In order to reduce duplicated inputs, `KOBE` defines an attribute named `defkwargs`, which defines a series of default parameter settings.
   One need to modify only once at the beginning of importing a class, and would be hence works for the rest calls.
  
As an example, we try building a tiling network with `KOBE`:
     
.. code-block:: bash
		
   >>> import kobe as kb
   >>> b=kb.tilings()
   >>> b.defkwargs
   {'wdir': './', 'clobber': False, 'rot_theta': 0, 'rot_phi': 0, 'rot_psi': 0, 'num': 1, 'emailpass': '', 'emailsmtp': '', 'subject': 'GW Alerts', 'fromaddr': '', 'toaddrs': '', 'slacktoken': None, 'slackto': None}
   >>> b.generatep(limdec=[-90,90])
   >>> b.generatep(limdec=[0,90])
   INFO:kobe.KBpointings:Warning: tiling data already parsed
  
As shown, the tile generation would be ignored if tiling is already exists.
Then one can set `clobber` parameter either at the call or at the beginning when initialize the class:

.. code-block:: bash

   # call clobber   
   >>> b.generatep(limdec=[0,90],clobber=True)
   >>> b.data

   # An alternative approach is to set it when initialize the `tilings` class		
   >>> b=kb.tilings(class=True)
   >>> b.defkwargs['clobber']
   True   
   >>> b.generatep(limdec=[-90,90])
   >>> b.generatep(limdec=[0,90])
   
| :ref:`Next <multiple>`
| :ref:`Top <single>`
