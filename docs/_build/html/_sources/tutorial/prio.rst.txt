.. _kbprio:
   
:mod:`Priorization and Optimization`
===========================================

:ref:`Previous <kbtrigger>`
    
Priorization
----------------------------------------

`Yang et al, 2019 <https://iopscience.iop.org/article/10.3847/1538-4357/ab0e06/meta>`_

.. code-block:: bash
   >>> from kobe import schedule
   >>> a.schedule()

1. generate OBs
   
   In `KOBE`, we adopt the concept, `OB`, to describe a (set of) observation(s) in detials.
   Apart from pointing, `KOBE OB` include also the filter, observing time, and/or the corresponding visibilities, trigger probabilities, etc.   

   .. code-block:: bash
		     
      >>> a.genobs(filters='r,i,g,u', ovhtime=100, exptime=45, mode=2)   
      >>> a.obs
      <Table length=2772>
      n      ra      dec    fovra   fovdec filter   time
      int64 float64  float64 float64 float64  str1  float64
      ----- -------- ------- ------- ------- ------ --------
      0  2.13491    19.0     1.0     1.0      g    100.0
      0  3.19253    19.0     1.0     1.0      g    145.0
      0  4.25015    19.0     1.0     1.0      g    190.0
      0  4.25671    20.0     1.0     1.0      g    235.0
      0  3.19253    20.0     1.0     1.0      g    280.0
      0  3.19253    21.0     1.0     1.0      g    325.0
      ...      ...     ...     ...     ...    ...      ...
      76 39.33957    38.0     1.0     1.0      u 141455.0
      76 38.07055    38.0     1.0     1.0      u 141500.0
      76 36.80153    38.0     1.0     1.0      u 141545.0
      76 36.78379    39.0     1.0     1.0      u 141590.0
      76 38.07055    39.0     1.0     1.0      u 141635.0
      76 39.35731    39.0     1.0     1.0      u 141680.0
   
Assessment and optimization
-------------------------------------------------

:ref:`Previous <kbtrigger>`
