.. _kbtel:
   
:mod:`Telescope and Observatory`
===========================================

| :ref:`Previous <kbpoint>`
| :ref:`Next <kbcand>`
| :ref:`1. Telescope  <kbtel1>`
|      :ref:`1.1 set pointings <kbtel11>`
|      :ref:`1.2 set telescope name <kbtel12>`
|      :ref:`1.3 check telescope <kbtel13>`
| :ref:`2. Observatory  <kbtel2>`
|      :ref:`2.1 set telescopes <kbtel21>`
|      :ref:`2.2 set observatory location <kbtel22>`
|      :ref:`2.3 set observatory name <kbtel23>`
|      :ref:`2.4 check observatory <kbtel24>`
| :ref:`3. Observatory Frame <kbtel3>`
|      :ref:`3.1 parse a frame <kbtel31>`
|      :ref:`3.2 visibility calculations for observatory <kbtel32>`
|      :ref:`3.3 visibility calculations for a list of coordinates <kbtel33>`
|      :ref:`3.4 visibility calculations for pointings and candidates <kbtel34>`
|      :ref:`3.5 visibility plots <kbtel35>`
|      :ref:`3.6 rank pointings based on visibility <kbtel36>`
| :ref:`4. Observatories  <kbtel4>`
|      :ref:`4.1 set a list of observatories <kbtel41>`
|      :ref:`4.2 plot the sky <kbtel42>`

On the basis of telescope pointings setup (which is described in the :ref:`previous chapter <kbpoint>`), here in this chapter we will describe the construction and usage of `KOBE telescope` and `observatory`.
     
.. _kbtel1:

Telescope
----------------------------------------

`KOBE telescope` object contains 2 attributes:

* `pointings` - telescope pointings
* `telname` - telescope name

.. code-block:: bash
		
   >>> from kobe import telescope
   >>> a=telescope()
   >>> a.telname
   >>> a.pointings

.. _kbtel11:
   
1. set pointings

   The pointing settings are described in :ref:`previous chapter<kbpoint>`,
   here we show again in general:

   .. code-block:: bash

      # initialize a `telescope` object
      >>> from kobe import telescope
      >>> a=telescope()
      
      # now pointings is Nonetype
      >>> a.pointings

      # set tiling schedule
      >>> a.set_pointings(strategy='T')
      >>> a.pointings
      <kobe.KBpointings.tilings object at 0x113557fd0>

      # generate tilings
      >>> a.pointings.generatep(fovra=3,fovdec=3,limdec=[0,90])      
      >>> a.pointings.data
      <Table length=2337>
      n       ra      dec    fovra   fovdec
      int64  float64  float64 float64 float64
      ----- --------- ------- ------- -------
      0       3.0     0.0     3.0     3.0
      1       6.0     0.0     3.0     3.0
      2       9.0     0.0     3.0     3.0
      3      12.0     0.0     3.0     3.0
      4      15.0     0.0     3.0     3.0
      5      18.0     0.0     3.0     3.0
      ...       ...     ...     ...     ...
      2331  57.32197    87.0     3.0     3.0
      2332 114.64394    87.0     3.0     3.0
      2333  171.9659    87.0     3.0     3.0
      2334 229.28787    87.0     3.0     3.0
      2335 286.60984    87.0     3.0     3.0
      2336 343.93181    87.0     3.0     3.0
 
      # visualize
      >>> a.locshow_tel()
      >>> a.savefig('tel1')
		   
   .. image:: ../static/tel1.png
      :width: 800
      :align: center
      :alt: tel1.png

.. _kbtel12:

2. set telescope name

   Considering we have to distinguish different telescopes, to define a name for each telescope is needed, thus `KOBE` provide `set_telname` method:
   
   If one give the `telname` directly, it would be applied as the telescope name:

   .. code-block:: bash

      >>> a.set_telname(telname='VST')
      >>> a.telname
      'VST'

   Otherwise, `KOBE` will check all cashed telescope names, select the largest one from those integers, and plus one as the `telname`.
   If no integer name found, use 0 instead:
   
   .. code-block:: bash
		   
      >>> a.set_telname()
      >>> a.telname
      0

.. _kbtel13:
      
3. check telescope

   If `telname` and `pointings` are parsed correctly, `check_tel` function will work, otherwise it will return an `AssertionError`:
   
   .. code-block:: bash
		
      >>> a.check_tel()


.. _kbtel2:
      
Observatory
----------------------------------------

After having `telescope`, the next step is to set up an `observatory` object.
It includes 3 attributes:

* telescopes - a series of telescopes located at the observaroty
* location - location of observaroty 
* obsname - name of observaroty

.. code-block:: bash
		
   >>> from kobe import observatory
   >>> a=observatory() 
   >>> a.telescopes
   {}
   >>> a.location
   >>> a.obsname

.. _kbtel21:

1. set telescopes

   The telescope settings are described :ref:`above <kbtel1>`.
   As an example, we try to define 2 telescopes for one observatory:
   
   .. code-block:: bash

      # we set one telescope with tiling strategy
      >>> a.set_pointings(strategy='T')   
      >>> a.pointings.generatep(limdec=[-20,90],fovra=3, fovdec=3)
      >>> a.set_telname()

      # add the telescope to observatory
      >>> a.add_telescope()
      INFO:kobe.KBobservatory:telescope 0 added

      # we add one more telescope with galaxy strategy
      >>> a.set_pointings(strategy='G')      
      >>> a.pointings.generatep(limdist=[0,40])
      >>> a.set_telname()

      # add the telescope to observatory
      >>> a.add_telescope()
      INFO:kobe.KBobservatory:telescope 1 added

      # now we have 2 telescopes at the observatory
      >>> a.telescopes
      {0: <kobe.KBpointings.tilings object at 0x11b86b390>, 1: <kobe.KBpointings.galaxies object at 0x11b84db38>}
      # and both of them have already generated a list of pointings
      >>> a.telescopes[0].data
      <Table length=3110>
      n       ra      dec    fovra   fovdec
      int64  float64  float64 float64 float64
      ----- --------- ------- ------- -------
      0   3.19253   -20.0     3.0     3.0
      1   6.38507   -20.0     3.0     3.0
      2    9.5776   -20.0     3.0     3.0
      3  12.77013   -20.0     3.0     3.0
      4  15.96267   -20.0     3.0     3.0
      5   19.1552   -20.0     3.0     3.0
      ...       ...     ...     ...     ...
      3104 309.79026    85.0     3.0     3.0
      3105  344.2114    85.0     3.0     3.0
      3106  85.96113    88.0     3.0     3.0
      3107 171.92225    88.0     3.0     3.0
      3108 257.88338    88.0     3.0     3.0
      3109  343.8445    88.0     3.0     3.0

      # we can visualize them
      >>> a.locshow_obs(marker='.')
      >>> a.savefig(filepath='obs',wdir='./')
   
   .. image:: ../static/obs.png
      :width: 800
      :align: center
      :alt: obs.png

   Then, one can remove any of them:
   
   .. code-block:: bash
		   
      >>> a.del_telescope(1)
      >>> a.telescopes
      {0: <kobe.KBpointings.tilings object at 0x11a41cfd0>}

   Or, copy one of them:
   
   .. code-block:: bash

      # to_telescope will adopt telecsope 0 settings to the current telname and pointings
      >>> a.to_telescope(0)      
      
      >>> a.set_telname(3)
      >>> a.add_telescope()
      INFO:kobe.KBobservatory:telescope 3 added
      >>> a.telescopes
      {0: <kobe.KBpointings.tilings object at 0x11d541fd0>, 3: <kobe.KBpointings.galaxies object at 0x11d57a278>}


.. _kbtel22:

2. set observatory location

   Despite a list of telescopes, another important ingrediant for `observatory` is its location.
   For `KOBE`, there're 3 approaches to parse the location:

   2.1 via an observatory name

   When user ask `KOBE` to parse a location via an obsname,
   `KOBE` will firstly check if it's included by `observatory.EarthLocation_from_kobe`.
   If not, `KOBE` will then check with `astropy.coordinates.EarthLocation.of_site`.   
   Since astropy query will sometimes take time (depends on the network), we recommend user to add your observatory in `observatory.EarthLocation_from_kobe` so that to save time.
   
   For instance, if `keck` is called frequently, one could add the following line to `observatory.EarthLocation_from_kobe`:
   
   .. code-block:: python
		      
      keck = EarthLocation.from_geodetic(-155.47833333333332*u.deg, 19.828333333333326*u.deg, 4160*u.m)

   and then parse via obsname will save time:
   
   .. code-block:: bash
		
      >>> a.parse_location(obsname='keck')
      >>> a.location
      <EarthLocation (-5464487.81759887, -2492806.59108569, 2151240.19451846) m>   
      
   2.2 via coordinates

   user could parsed location via coordinates inetead:
   
   .. code-block:: bash
		
      >>> a.parse_location(longitude=-155.48, latitude=19.83, elevation=4160)
      >>> a.location
      <EarthLocation (-5464503.34830818, -2492621.64331625, 2151413.87241013) m>
      
   2.3 via astropy.coordinates.EarthLocation

   one could also provide an `astropy.coordinates.EarthLocation` object directly to `KOBE`:

   .. code-block:: bash

      # suppose b is an `astropy.coordinates.EarthLocation` object
      a.parse_location(location=b)

.. _kbtel23:

3. set observatory name

   Similar to `telescope`, to define a name for `observatory` is essential for a recognition purpose.
   One can set the `obsname` via `set_obsname` as followed:
   
   .. code-block:: bash
		   
      >>> a.set_obsname(obsname='keck')
      >>> a.obsname
      'keck'

   There's no need to set `obsname`, if the `observatory` location is parsed via a name or coordinates, e.g.
            
   .. code-block:: bash
		
      >>> a.parse_location(obsname='keck')      
      >>> a.obsname
      'keck'
      
      >>> a.parse_location(longitude=-155.48, latitude=19.83, elevation=4160)
      >>> a.obsname
      'T-155.48-19.83'

.. _kbtel24:

4. check observatory

   If `obsname` is parsed correctly,
   `check_obs` method will work normally, otherwise it will return an `AssertionError`.
   Meanwhile, if `telescopes` attribute is empty, it will alert a warning.
   
   .. code-block:: bash
		
      >>> a.check_obs()

.. _kbtel3:

Observatory Frame
----------------------------------------

After the `location` of an `observatory` set correctly (check :ref:`above<kbtel2>`), 
once observation time parsed, `KOBE` could generate an `astropy` frame via `KOBE frame` class.

The `frame` object contains 2 attributes:

* frame - an astropy frame
* obstime - observation time or times 

.. _kbtel31:
  
1. parse a frame
   
   .. code-block:: bash

      # initialize a frame and set location as keck for example
      >>> from kobe import frame
      >>> a=frame()
      >>> a.parse_location(obsname='keck')

      # then, let's try to parse observing time(s), the options include:
      #     1. an astropy.time.Time
      #     2. `str` now - present currently
      #     3. a time `str`, e.g. 1999-01-01T00:00:00.123456789
      #     4. a `float` or `int` number, means how much time before or later than currently
      #        e.g. time=-3600 (timeu=sec) stands for 1 hour before
      #     5. sequence of above options
      >>> a.parse_time(time='now')
      >>> a.obstime
      <Time object: scale='utc' format='datetime' value=2020-03-16 07:57:08.115352>
      >>> a.parse_frame()
      >>> a.frame
      <AltAz Frame (obstime=2020-03-16 07:57:08.115352, location=(-5464487.81759887, -2492806.59108569, 2151240.19451846) m, pressure=0.0 hPa, temperature=100.0 deg_C, relative_humidity=0.0, obswl=1.0 micron)>

      >>> a.parse_time(time=['now',1000,'1999-01-01T00:00:00.123456789'])
      >>> a.obstime
      array([<Time object: scale='utc' format='datetime' value=2020-03-24 08:13:09.122735>,
             <Time object: scale='utc' format='datetime' value=2020-03-24 08:29:49.123365>,
	     <Time object: scale='utc' format='isot' value=1999-01-01T00:00:00.123>],
	     dtype=object)	  
      >>> a.parse_frame()
      >>> a.frame
      <AltAz Frame (obstime=[datetime.datetime(2020, 3, 24, 8, 13, 9, 122735)
		   datetime.datetime(2020, 3, 24, 8, 29, 49, 123365)
		   datetime.datetime(1999, 1, 1, 0, 0, 0, 123457)],
		   location=(-5464487.81759887, -2492806.59108569, 2151240.19451846) m,
		   pressure=0.0 hPa, temperature=100.0 deg_C, relative_humidity=0.0, obswl=1.0 micron)>

.. _kbtel32:

2. visibility calculations for observatory
   
   Once a frame parsed, we could then calculate several visibilities for the observaroty,
   e.g. sun constrain, moon illumination and so on:

   .. code-block:: bash

      # calculate the height of sun
      >>> from kobe import frame
      >>> a=frame()
      >>> a.parse_location(obsname='keck')
      >>> a.calc_sun_height(time='2020-03-16 07:57:08.115352')
      <Latitude -48.30077726 deg>

      # calculate moon illumination (ranged from 0 to 1)   
      >>> a.calc_moon_illumination(time='2020-01-01')
      0.2995534608971628

      # estimate the start and end of the following observing night from an input time
      >>> a.at_night(time='2020-01-01')
      (<Time object: scale='utc' format='iso' value=2020-01-01 05:20:00.000>, <Time object: scale='utc' format='iso' value=2020-01-01 15:40:00.000>)


.. _kbtel33:

3. visibility calculations for a list of coordinates
   
   One could also input a series of coordinates to ask for their visibilities, e.g. airmass, seperation to solar object, and so on.

   .. code-block:: bash

      # obtain alt, az angle
      >>> a.altaz(ra=[0,45,90],dec=20,time=1,timeu='jd')
      <SkyCoord (AltAz: obstime=2020-03-26 04:04:04.836131, location=(-5464487.81759887, -2492806.59108569, 2151240.19451846) m, pressure=0.0 hPa, temperature=100.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt) in deg
      [(288.81014475,  7.26921235), (278.23873761, 48.46060971),
      ( 75.07976199, 89.33053841)]>

      # calculate the airmass
      >>> a.calc_airmass(ra=[0,45,90],dec=[20,-20,0],time='2020-03-25 03:44:03.391660')
      <Quantity [4.56849008, 1.730989  , 1.07024843]>

      # calculate specific fields seperation to the sun
      >>> a.calc_solar(solarobj='sun', time='2020-03-25 00:00:00', ra=30, dec=60)
      <Angle [74.76922248] deg>

      # calculate specific fields seperation to the galactic plane
      >>> a.calc_galactic(time='2020-03-25 00:00:00', ra=30, dec=60)
      <Angle [1.73839012] deg>

.. _kbtel34:

4. visibility calculations for `pointings` and `candidates`
   
   Instead of input coordinates manually, above functions are also applicable for `pointings` and `candidates` attributes if available:

   .. code-block:: bash

      # since for this case, we need both `frame` and `pointings` (or `candidates`)
      # we should import their son class, i.e. schedule instead
      >>> from kobe import schedule
      >>> a=schedule()
   
      # set pointings
      >>> a.set_pointings(strategy='T')
      >>> a.pointings.generatep(fovra=3,fovdec=3,limdec=[0,90])

      # parse observatory location
      >>> a.parse_location(obsname='keck')

      # get alt, az angle for all pointings
      >>> a.altaz_pointings(time='2020-01-01')
      <SkyCoord (AltAz: obstime=2020-01-01 00:00:00.000, location=(-5464487.81759887, -2492806.59108569, 2151240.19451846) m, pressure=0.0 hPa, temperature=100.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt) in deg
      [(101.57698197, 29.37889098), (100.27515512, 26.60745015),
      ( 99.03508821, 23.82472534), ..., (356.8031899 , 20.50632869),
      (358.90534896, 22.61682362), (  1.95809236, 22.06590241)]>

      # for `candidates` case
      # It will be described in the next chapter
      >>> a.readc(ra=[1,10,30], dec=[-30,20,57])
      >>> a.altaz_candidates(time='2020-01-01')
      <SkyCoord (AltAz: obstime=2020-01-01 00:00:00.000, location=(-5464487.81759887, -2492806.59108569, 2151240.19451846) m, pressure=0.0 hPa, temperature=100.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt) in deg
      [(131.06454937, 16.24572341), ( 77.35122346, 28.76716843),
      ( 34.91524674, 18.84961134)]>   


.. _kbtel35:

5. visibility plots

   To make visibility result shown directly, we provide methods to visualize them:

   .. code-block:: bash

      # import frame class
      >>> from kobe import frame
      >>> a=frame()
      >>> a.parse_location(obsname='keck')

      # sun height plot
      >>> import numpy as np
      # show the height of sun at keck, every half an hour from now to 10 hours later
      >>> a.calc_sun_height(time=np.arange(20)*1800,show=True,marker='o',color='r')
      >>> a.savefig('sunh')
      
   .. image:: ../static/sunh.png
      :width: 800
      :align: center
      :alt: sunh.png

   .. code-block:: bash
		   
      # moon illumination
      # show the moon illumination at keck, every day from now till 20 days later
      >>> a.calc_moon_illumination(time=np.arange(20),timeu='jd',show=True,marker='o',color='r')
      >>> a.savefig('mooni')

   .. image:: ../static/mooni.png
      :width: 800
      :align: center
      :alt: mooi.png      

   .. code-block:: bash

      # for asiago, calculate airmass for 3 coordinates, every half an hour from now to 10 hours later
      >>> a.parse_location(obsname='asiago')
      >>> a.calc_airmass(ra=[0,20,40],dec=[20,20,20],time=np.arange(20)*1800,show=True, marker='.')
      >>> import pylab as pl
      >>> pl.legend()
      >>> a.savefig('airmass')

   .. image:: ../static/airmass.png
      :width: 800
      :align: center
      :alt: airmass.png

   Meanwhile, visualization functions is also applicable to `pointings` and `candidates`:

   .. code-block:: bash

      >>> from kobe import schedule
      >>> a=schedule()

      >>> a.parse_location('asiago')
      >>> a.set_pointings(strategy='T')
      >>> a.pointings.generatep(fovra=3,fovdec=3,limdec=[30,40], limra=[40,50])
      >>> a.calc_airmass_pointings(time=np.arange(20)*1800,show=True)
      >>> import pylab as pl
      >>> pl.legend()
      >>> a.savefig('airmassp')

   .. image:: ../static/airmassp.png
      :width: 800
      :align: center
      :alt: airmassp.png

	    
   .. code-block:: bash

      >>> from kobe import schedule
      >>> a=schedule()
      >>> a.parse_location('asiago')
      >>> a.readc(ra=[1,10,30], dec=[30,20,57])
      >>> a.calc_airmass_candidates(time=np.arange(20)*1800,show=True)
      >>> import pylab as pl
      >>> pl.legend()
      >>> a.savefig('airmassc')

   .. image:: ../static/airmassc.png
      :width: 800
      :align: center
      :alt: airmassc.png


.. _kbtel36:
   
6. rank pointings based on visibility
   

   
.. _kbtel4:

Observatories
----------------------------------------

.. _kbtel41:

1. set a list of observatories
   
Then one can define a list of observatories, as shown below:

.. code-block:: bash

   # import `KOBE chedule`, a
   >>> from kobe import schedule
   >>> a=schedule()

   # add one `observatory`, lapalma
   >>> a.add_observatory('lapalma')

   # at lapalma, define a telescope, VST
   # whose field of view is 1 by 1 deg
   # here, in order to show tilings clearly (1 by 1 deg is too small)
   # we show OBs which is composed of 3 by 3 tilings
   >>> a.set_pointings(strategy='T')
   >>> a.pointings.generatep(fovra=3,fovdec=3,limdec=[-90,0])  
   >>> a.add_telescope('VST')

   # then we could define another observatory, asiago
   >>> a.add_observatory('asiago')

   # at asiago, define a telescope, schmidt   
   >>> a.set_pointings(strategy='T')
   >>> a.pointings.generatep(fovra=3,fovdec=3,limdec=[0,90])
   >>> a.add_telescope('schmidt')

   # then we have the observatories
   >>> a.observatories
   {'lapalma': {'VST': <kobe.KBpointings.tilings object at 0x11f3481d0>}, 'asiago': {'schmidt': <kobe.KBpointings.tilings object at 0x11f348240>}}
   >>> a.observatories['lapalma']['VST'].data
   <Table length=27666>
   n              ra             dec    fovra   fovdec
   int64         float64         float64 float64 float64
   ----- ----------------------- ------- ------- -------
   0 -1.6331239353195196e+16   -90.0     1.0     1.0
   0               171.89607   -90.0     1.0     1.0
   0   1.633123935319554e+16   -90.0     1.0     1.0
   0               114.59738   -89.0     1.0     1.0
   0               171.89607   -89.0     1.0     1.0
   0               229.19476   -89.0     1.0     1.0
   ...                     ...     ...     ...     ...
   3073               357.47579    19.0     1.0     1.0
   3073               358.53341    19.0     1.0     1.0
   3073               359.59103    19.0     1.0     1.0
   3073               357.46923    20.0     1.0     1.0
   3073               358.53341    20.0     1.0     1.0
   3073               359.59759    20.0     1.0     1.0

   # and we can visualize them
   >>> a.locshow_all()
   >>> a.savefig('observatories')

.. image:: ../static/observatories.png
   :width: 800
   :align: center
   :alt: observatories.png

.. _kbtel42:

2. plot the sky

.. code-block:: bash

   # import `KOBE chedule`, a
   >>> from kobe import schedule
   >>> a=schedule()

   # plot the sky
   >>> a.plot_sky()

   # paser a time
   >>> a.parse_time('20200101')
   >>> a.plot_sky(clear=True)

   # paser observatories
   
   
Check more detials in the corresponding API page.

| :ref:`Previous <kbpoint>`
| :ref:`Next <kbcand>`
| :ref:`Top <kbtel>`
