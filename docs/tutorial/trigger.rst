.. _kbtrigger:
   
:mod:`Trigger`
===========================================
| :ref:`Previous <kbcand>`
| :ref:`Next <kbprio>`
| :ref:`1. Parse trigger <kbtrigger1>`
| :ref:`2. Visualization <kbtrigger2>`
| :ref:`3. Trigger evaluation <kbtrigger3>`
| :ref:`4. Candidates evaluation based on trigger <kbtrigger4>`
| :ref:`5. Pointings evaluation based on trigger <kbtrigger5>`
|      :ref:`5.1 rank pointings based on probability <kbtrigger51>`
| :ref:`6. Trigger report and Circulate <kbtrigger6>`


In this chapter, we describe `KOBE trigger` usage, which would be used for followup search.

.. _kbtrigger1:

Parse healpix map
----------------------------------------

In `KOBE`, we have defined 5 possible trigger sources, and adopt `Healpy/HEALPix <https://healpy.readthedocs.io/en/latest/index.html>`_ to help parse trigger information from those alert filesm, which includes:

* healpix fits file
* healpix fits url
* VOevent XML file
  
* the root object (convinient for `pygcn <https://github.com/lpsinger/pygcn>`_ users)
* the coordinates.
  
The examples are shown below:

.. code-block:: bash

   # import trigger module
   >>> from kobe import trigger
   >>> a = trigger()

   # parsed healpix map via healpix fits url
   >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
   NSIDE = 512
   ORDERING = NESTED in fits file
   INDXSCHM = IMPLICIT
   Ordering converted to RING
   Ordering converted to RING
   Ordering converted to RING
   Ordering converted to RING
   
   >>> a.hpmap
   array([1.25419392e-07, 1.62144031e-07, 1.79526856e-07, ...,
   2.46590458e-10, 9.75865543e-10, 6.87730424e-10])
   
   >>> a.hpd1
   array([144.84488626, 154.92359484, 140.96139722, ..., 272.89617891,
   221.54107632, 282.52651738])
   
   >>> a.hphdr
   {'XTENSION': 'BINTABLE', 'BITPIX': 8, 'NAXIS': 2, 'NAXIS1': 32, 'NAXIS2': 3145728, 'PCOUNT': 0, 'GCOUNT': 1, 'TFIELDS': 4, 'TTYPE1': 'PROB', 'TFORM1': 'D', 'TUNIT1': 'pix-1', 'TTYPE2': 'DISTMU', 'TFORM2': 'D', 'TUNIT2': 'Mpc', 'TTYPE3': 'DISTSIGMA', 'TFORM3': 'D', 'TUNIT3': 'Mpc', 'TTYPE4': 'DISTNORM', 'TFORM4': 'D', 'TUNIT4': 'Mpc-2', 'MOC': True, 'PIXTYPE': 'HEALPIX', 'ORDERING': 'NESTED', 'COORDSYS': 'C', 'NSIDE': 512, 'INDXSCHM': 'IMPLICIT', 'OBJECT': 'G331903', 'REFERENC': 'https://gracedb.ligo.org/events/G331903', 'INSTRUME': 'H1,L1,V1', 'DATE-OBS': '2019-05-10T02:59:39.292500', 'MJD-OBS': 58613.12476032978, 'DATE': '2019-05-10T03:00:47.000000', 'CREATOR': 'BAYESTAR', 'ORIGIN': 'LIGO/Virgo', 'RUNTIME': 18.0, 'DISTMEAN': 268.8566049372629, 'DISTSTD': 108.0709050006497, 'LOGBCI': 0.6949211109947058, 'LOGBSN': 7.032293281836687, 'VCSVERS': 'ligo.skymap 0.1.6', 'VCSREV': '79504ec9fb1890fa91665bd69d7aa66cdaf11184', 'DATE-BLD': '2019-03-26T18:11:21', 'HISTORY': 'gwcelery worker -l info -n gwcelery-openmp-worker -Q openmp -c 1'}

   # an equivalent approach is to parse url via parse_trigger method   
   >>> a = trigger(url='https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
   >>> a.parse_trigger()

   # If one had a voevent file, LVC#S190510g-2-Initial.xml, downloaded in current folder
   >>> from kobe import trigger
   >>> a=trigger()
   >>> a.xml('LVC#S190510g-2-Initial.xml',wdir='./')

   # or try parse_trigger
   >>> a = trigger()
   >>> a.parse_trigger(vofile='https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz', wdir='./')   


.. _kbtrigger2:

Visualization
----------------------------------------

Once trigger healpix map parsed successfully, we could afterwards visualize it:

.. code-block:: bash

   >>> from kobe import trigger
   >>> a = trigger()
   >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')

   # show healpix map, together with 50 percent contours, via mollview scheme
   >>> a.locshow(cls=[.5], showhp=True, showtype='m', showgrid=True)
   >>> a.savefig('trigger1')

.. image:: ../static/trigger1.png
   :width: 800
   :align: center
   :alt: trigger1.png
	    
.. code-block:: bash

   # show 50 and 90 percent contours, rotating 90 deg along theta direction
   # set clear True to remove previous plot
   >>> a.locshow(cls=[.5, .9], rot_theta=90, clear=True)
   >>> a.savefig('trigger2')

.. image:: ../static/trigger2.png
   :width: 800
   :align: center
   :alt: trigger2.png
	 
.. code-block:: bash

   # use healpix.cartview scheme instead of mollview
   # check more in healpy page, or corresponding KOBE API page
   >>> a.locshow(cls=[.5, .9], showtype='c', clear=True)
   >>> a.savefig('trigger3')

.. image:: ../static/trigger3.png
   :width: 800
   :align: center
   :alt: trigger3.png

.. code-block:: bash

   # show contours only, witout healpix map
   >>> a.locshow(cls=[.5, .9], showhp=False, clear=True)
   >>> a.savefig('trigger4')

.. image:: ../static/trigger4.png
   :width: 800
   :align: center
   :alt: trigger4.png
	 
.. code-block:: bash

   # set visualization threshold for healpix map
   >>> a.locshow(clear=True,min=1e-8,max=1e-5)
   >>> a.savefig('trigger5')

.. image:: ../static/trigger5.png
   :width: 800
   :align: center
   :alt: trigger5.png

	 
.. _kbtrigger3:
	 
Trigger evaluation
----------------------------------------

Since not every trigger is that interesting, `KOBE` would try gather all possible informations from sources, so that user could assess them and decide if it's worth follow or not.

1. trigger information from healpix fits header or VOevent file

   Most informations are already recorded in the fits header or VOevent:
   
   .. code-block:: bash

      >>> from kobe import trigger
      >>> a = trigger()

      # if healpix map parsed from a VOevent xml file
      # both hphdr and voinf should be available
      >>> a.parse_trigger(vofile='https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz', wdir='./')      
      >>> a.hphdr
      {'XTENSION': 'BINTABLE', 'BITPIX': 8, 'NAXIS': 2, 'NAXIS1': 32, 'NAXIS2': 3145728, 'PCOUNT': 0, 'GCOUNT': 1, 'TFIELDS': 4, 'TTYPE1': 'PROB', 'TFORM1': 'D', 'TUNIT1': 'pix-1', 'TTYPE2': 'DISTMU', 'TFORM2': 'D', 'TUNIT2': 'Mpc', 'TTYPE3': 'DISTSIGMA', 'TFORM3': 'D', 'TUNIT3': 'Mpc', 'TTYPE4': 'DISTNORM', 'TFORM4': 'D', 'TUNIT4': 'Mpc-2', 'MOC': True, 'PIXTYPE': 'HEALPIX', 'ORDERING': 'NESTED', 'COORDSYS': 'C', 'NSIDE': 512, 'INDXSCHM': 'IMPLICIT', 'OBJECT': 'G331903', 'REFERENC': 'https://gracedb.ligo.org/events/G331903', 'INSTRUME': 'H1,L1,V1', 'DATE-OBS': '2019-05-10T02:59:39.292500', 'MJD-OBS': 58613.12476032978, 'DATE': '2019-05-10T03:00:47.000000', 'CREATOR': 'BAYESTAR', 'ORIGIN': 'LIGO/Virgo', 'RUNTIME': 18.0, 'DISTMEAN': 268.8566049372629, 'DISTSTD': 108.0709050006497, 'LOGBCI': 0.6949211109947058, 'LOGBSN': 7.032293281836687, 'VCSVERS': 'ligo.skymap 0.1.6', 'VCSREV': '79504ec9fb1890fa91665bd69d7aa66cdaf11184', 'DATE-BLD': '2019-03-26T18:11:21', 'HISTORY': 'gwcelery worker -l info -n gwcelery-openmp-worker -Q openmp -c 1'}
      >>> a.voinf
      {'Packet_Type': '151', 'internal': '0', 'Pkt_Ser_Num': '2', 'GraceID': 'S190510g', 'AlertType': 'Initial', 'HardwareInj': '0', 'OpenAlert': '1', 'EventPage': 'https://gracedb.ligo.org/superevents/S190510g/view/', 'Instruments': 'H1,L1,V1', 'FAR': '8.42945108717e-10', 'Group': 'CBC', 'Pipeline': 'gstlal', 'Search': 'AllSky', 'skymap_fits': 'https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz', 'BNS': '0.979698536622', 'NSBH': '0.0', 'BBH': '0.0', 'MassGap': '0.0', 'Terrestrial': '0.0203014633781', 'HasNS': '1.0', 'HasRemnant': '1.0'}

      # obviously, if healpix map is parsed from a healpix fits file
      # the voinf attribute would be an Nonetype
      >>> from kobe import trigger
      >>> a = trigger()
      >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
      >>> a.voinf

2. trigger localizarion area and distance statistics

   Since currently the localization area is been not recorded in healpix fits header or VOevent xml file yet for most cases, which is however important for the trigger evaluation.
   `KOBE` provide methods for the area estimations (**note** that for CBC type GW triggers, LVC would provide also distance assessment: the average estimation is recorded in the header, while informations at each directions are stored in the healpix fits):

   For 2d area calculations:
   
   .. code-block:: bash

      # get healpix indices for specific confidence levels
      >>> a.calc_loc_contours(cls=[.1,.5,.9,.99])
      {0.1: array([2414072, 2418168, 2416119, ..., 1570953, 1573001, 1573000]), 0.5: array([2414072, 2418168, 2416119, ..., 2018356, 2022450, 2024499]), 0.9: array([2414072, 2418168, 2416119, ...,  783552,  734398,  771264]), 0.99: array([2414072, 2418168, 2416119, ..., 1309038, 1309052, 1309051])}

      # estimate the erea of error boxes at specific confidence levels
      # unit in sq deg
      >>> a.calc_loc_area(cls=[.1,.5,.9,.99])
      {0.1: 36.351906008208665, 0.5: 575.5456172035578, 0.9: 3463.7648737864856, 0.99: 11508.39446313552}

.. _kbtrigger4:
      
Candidates evaluation based on trigger
------------------------------------------------------------

Then, one could assess and select candidates depending on trigger map:

.. _kbtrigger5:
      
Pointings evaluation based on trigger
------------------------------------------------------------

.. _kbtrigger51:

1. rank pointings based on probability
   
   Then, let try greedy appracoh for which a trigger is needed
   
   .. code-block:: bash

      # parse a ligo trigger and visualize it
      >>> a=schedule()
      >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
      >>> a.locshow(max=1e-5, cls=[.5, .9, .99])

      # set telescope pointings
      >>> a.set_pointings(strategy='T')
      >>> a.pointings.generatep(fovra=3,fovdec=3,limra=[20, 100], limdec=[0,30])

      # rank with GW probability
      >>> a.rankp(approach=2, mode=2, threshold=5,sort=1)

      # show pointings and the routes
      >>> a.locshow_tel(color='r')
      >>> a.routeshow_tel(color='g')     
      >>> a.zoomin([-1.2,0],[-.5,.6])
      >>> a.savefig('routeshow1')

   .. image:: ../static/routeshow1.png
      :width: 800
      :align: center
      :alt: routeshow1.png

   Sometimes, one might prefer to cover trigger probability as faster as possible,
   and for such case, we could adopte strict mode, namely `mode` = 1:
   
   .. code-block:: bash
      
      >>> a.rankp(approach=2, mode=1, threshold=5,sort=1)      
      >>> a.savefig('routeshow2')
      
   .. image:: ../static/routeshow2.png
      :width: 800
      :align: center
      :alt: routeshow2.png





	    
One could assess and select pointings depending on trigger:

   .. code-block:: bash

      # ask 2d trigger probablity at specific positions
      >>> a.calc_loc_prob(ra=[1,10], dec=[1,-20])
      array([7.36631063e-07, 8.50340927e-06])
      >>> a.calc_loc_prob(ra=[1,10], dec=[1,-20],fovra=1,fovdec=1)
      [5.8577539241927334e-05, 0.0006069629043450725]
      >>> a.calc_loc_prob(ra=[1,10], dec=[1,-20],fovra=[1,10],fovdec=[1,10])
      [5.8577539241927334e-05, 0.042806591674792185]

      # ask confidence level at specific positions, e.g. for such trigger:
      # position at ra=1, dec=1 is located at the top 87.18% region
      # while ra=10, deg=-20 is located at the top 34.51% region, which is more likely to have the contourpart in principle
      >>> a.calc_loc_cl(ra=[1,10], dec=[1,-20])
      [0.8718188163741416, 0.34519511496241106]
   
      # ask distance mean, variance and normalization factor at specific positions
      >>> a.calc_dis_inf(ra=[1,10], dec=[1,-20])
      (array([233.86346031, 273.32868364]), array([121.25377326, 112.45704833]), array([1.44320589e-05, 1.14500239e-05]))

      # for the trigger probability at the distance direction, despite coordinates, distance of each dots or tiles is also needed. We provide 2 approaches for the distance probablity estimation.

      # 1. user provide distance, and we provide how many sigma at specific coordinates:
      >>> a.calc_dis_sigma(ra=[1,10], dec=[1,-20], dist=[100,200])
      array([1.10399418, 0.65205947])
      >>> a.calc_dis_sigma(ra=[1,10], dec=[1,-20], dist=100)
      array([1.10399418, 1.54128786])
      
      # 2. or, one provide how much deeper a telescope could reach, i.e. the limiting magnitude at one filter, and parsed a targeting lightcurve via `KOBE candidates`. Thus `KOBE` will try to find the highest probablity that one telescope could reach at specific positions:
      >>> a.readlc('tmplc')
      >>> a.calc_dis_prob(ra=[1,10], dec=[1,-20], limmag=20)
      [2.8391830707199967e-06, 1.6860153679224023e-06]


      # similar to `frame`, here user could adopt above calculations to `candidates` or `pointings`, by adding '_candidates', '_pointings' correspondingly. But note that we should import their son classes instead:
      >>> from kobe import schedule
      >>> a = schedule()   

      # generate a list of candidates first
      >>> a.readc(ra=[1,10,30], dec=[-30,20,57])
      >>> a.calc_loc_prob_candidates(add=True)
      >>> a.candidates
      <Table length=3>
      n     ra   dec  ...   mag            prob
      int64 int64 int64 ... float64        float64
      ----- ----- ----- ... ------- ----------------------
      0     1   -30 ...     0.0 1.9284269010978466e-08
      1    10    20 ...     0.0 2.0231011501525045e-09
      2    30    57 ...     0.0 1.8087733560399342e-10
   
      # or, we simulate a list of candidates depending on trigger informations
      >>> a.sim_candidates(100)
      >>> a.candidates
      <Table length=100>
      n           ra         ...        dist
      int64      float64       ...      float64
      ----- ------------------ ... ------------------
      0       23.115234375 ...  375.9515624386379
      1       24.169921875 ... 329.81103124515795
      2 202.93945312499997 ...  298.0457228059061
      3 20.390624999999996 ...  340.7921527021853
      4 29.179687500000004 ... 348.23222429754685
      5 219.46289062499997 ...  209.4316853009612
      ...                ... ...                ...
      94  34.89257812499999 ...  312.6570100894105
      95         234.140625 ... 194.50830434378472
      96       189.66796875 ...  314.2610304124869
      97        4.306640625 ...   336.507412048297
      98  90.08789062499999 ...  136.7218608009585
      99 223.76953124999997 ... 182.44403569924526

      # then, we calculate the probability for those candidates
      # with add True, will append probablities to `candidates` attribitue
      >>> a.calc_loc_prob_candidates(add=True)
      >>> a.candidates
      <Table length=100>
      n           ra         ...          prob
      int64      float64       ...        float64
      ----- ------------------ ... ----------------------
      0      198.720703125 ... 1.0871144359444234e-06
      1        34.27734375 ... 3.6993686880535465e-06
      2       354.90234375 ... 1.5092874203885056e-06
      3  87.36328124999999 ... 3.1338349257123985e-05
      4       30.673828125 ...  1.952836299467096e-05
      5 204.78515624999997 ... 2.6864287804351886e-05
      ...                ... ...                    ...
      94        222.5390625 ...  7.778463626082537e-06
      95       37.177734375 ... 1.3224485371845156e-06
      96        88.06640625 ...  5.666885303976507e-05
      97        220.4296875 ...  5.495307908097957e-06
      98       87.978515625 ... 4.8930073608119987e-05
      99 345.49804687499994 ...  3.423411778829105e-07

      # otherwise, will return a probability sequence
      >>> a.calc_loc_prob_candidates(add=False)
      array([1.08711444e-06, 3.69936869e-06, 1.50928742e-06, 3.13383493e-05,
      1.95283630e-05, 2.68642878e-05, 1.10427040e-05, 6.71512096e-06, ...])

      # working the same for `calc_dis_sigma`:
      >>> a.calc_dis_sigma_candidates()
      <Column name='dist' dtype='float64' length=100>
      0.14880546273025794
      0.5665491638766624
      0.8015084423389162
      0.7534019306893466
      0.5606431781802166
      0.4155206143818324
      1.9909015138436261
      0.1306552537277863
      ...
      2.358471744440794
      0.5990096215373739
      0.9449567534390637
      1.154412315546838
      0.002524681517209831
      0.8895511354394909
      0.47304138903496756
      0.5663346170205388

      # working the same for `pointings`:
      >>> from kobe import schedule
      >>> a = schedule()       
      >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
      >>> a.set_pointings(strategy='T')
      >>> a.pointings.generatep(limdec=[-20,90])
      >>> a.calc_loc_prob_pointings()

      # One could also cut pointings depending on trigger, say, remaining only those pointings at high probability region
      >>> from kobe import schedule
      >>> a = schedule()

      # parse trigger
      >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')

      # parse tilings with dec above -20 deg, fov as 5 by 5 deg
      >>> a.set_pointings(strategy='T')      
      >>> a.pointings.generatep(limdec=[-20,90], fovra=5, fovdec=5)

      # select those tilings at 50% and 90% confidence level region
      >>> b=a.cut_pointings_hpmap(cls=[.5,.9])
      >>> b[.5]
      <Table length=331>
      n      ra      dec    fovra   fovdec
      int64 float64  float64 float64 float64
      ----- -------- ------- ------- -------
      5  6.38507   -20.0     1.0     1.0
      6  7.44924   -20.0     1.0     1.0
      7  8.51342   -20.0     1.0     1.0
      8   9.5776   -20.0     1.0     1.0
      9 10.64178   -20.0     1.0     1.0
      10 11.70596   -20.0     1.0     1.0
      ...      ...     ...     ...     ...
      14462 36.41893    21.0     1.0     1.0
      14463 37.49007    21.0     1.0     1.0
      14464 38.56122    21.0     1.0     1.0
      14465 39.63236    21.0     1.0     1.0
      14798 36.67018    22.0     1.0     1.0
      14799 37.74872    22.0     1.0     1.0

      # show the trigger and tilings at 90% C.L.
      >>> a.locshow(cls=[.9])
      >>> a.plot_box_data(b[.9],color='k')
      >>> a.savefig('cutp')
   
   .. image:: ../static/cutp.png
      :width: 800
      :align: center
      :alt: cutp.png

.. _kbtrigger6:
	    
Trigger report and Circulate
----------------------------------------

.. code-block:: bash
		
   >>> a.trigger_report(cls=[.5, .9], style='sms')
   >>> a.texts
   'DISTMEAN:268.86 DISTSTD:108.07 DATE-OBS:2019-05-10T02:59:39.292500 MJD-OBS:58613.12 OBJECT:G331903 INSTRUME:H1,L1,V1 CREATOR:BAYESTAR 0.5:575.55 0.9:3463.76 '

:ref:`Previous <kbcand>`
:ref:`Next <kbprio>`
