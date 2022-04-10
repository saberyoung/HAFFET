.. _kbpoint:
   
:mod:`Pointings - Tilings, Galaxies`
===========================================

| :ref:`Previous <kbscheme>`
| :ref:`Next <kbtel>`    
| :ref:`1. tiling search  <tiles>`
|      :ref:`1.1 generate pointings <tiles1>`
|      :ref:`1.2 show pointings <tiles2>`
|      :ref:`1.3 skip portion of sky <tiles3>`
|      :ref:`1.4 group and divide pointings <tiles4>`
|      :ref:`1.5 rank pointings <tiles5>`
|      :ref:`1.6 report and circulate <tiles6>`     
| :ref:`2. galaxy search  <galaxies>`
|      :ref:`2.1 generate pointings <galaxies1>`
|      :ref:`2.2 show pointings <galaxies2>`
|      :ref:`2.3 skip portion of sky <galaxies3>`
|      :ref:`2.4 group and divide pointings <galaxies4>`
|      :ref:`2.5 rank pointings <galaxies5>`
|      :ref:`2.6 report and circulate <galaxies6>`
|      :ref:`2.7 make galaxy healpix map <galaxies7>`
| :ref:`3. pointing  <point>`
     
In this chapter, we describe how we construct pointings for telescopes. 

Currently there're 2 strategies adopted by astronomers to survey the sky, namely, tiling and galaxy search.
Tiling search is suitable for the telescopes with relatively large field of view.
Astronomers tile the trigger localization (or portion of the sky) into a series of telescope shots, observe them in order.
While as for the galaxy strategy which is more suitable for small field of view telescopes, the galaxies inside the trigger region (or portion of the sky) are monitored with a certain frequency.
The crucial point of such search is a good knowledge of the local galaxies,
while currently the most used public galaxy catalog, i.e. GLADE, would have a good completeness up to several tens (or hundreds, depends on how much completeness one would like to cover) of Mpc.
The detection of the kilonova has proved the success of galaxy search, however, we should know that tiling search would be the only choice if one trigger locate in further volume, where we had few knowledges about the mass distributions.

.. _tiles:

Tiling search
----------------------------------------

.. _tiles1:

1. generate pointings
      
   There're several approaches to generate tile pointings after a `tilings` object was initialized.
   Once generation running successfully, there will be one `data` attribute available, that contains a list of tilings, namely ra, dec at the center, field of view of tiling box, and an index.

   1.1 One approach for tiling generation is to input pointings directly as shown below.
   
   .. code-block:: bash

      # import module and initialize a tiling object, a
      >>> from kobe import tilings
      >>> a=tilings()

      # fovra/fovdec should be a number when ra/dec is a number
      >>> a.readp_coo(ra=20,dec=30,fovra=5,fovdec=5)
      >>> a.data
      <Table length=1>
      n     ra   dec  fovra fovdec
      int64 int64 int64 int64 int64
      ----- ----- ----- ----- ------
      0    20    30     5      5

      # when ra/dec is sequence
      # fovra/fovdec could be either a number, or a same length sequence
      >>> a.readp_coo(ra=[0,10,20,30],dec=[0,10,20,30],fovra=5,fovdec=5)            
      >>> a.data
      <Table length=4>
      n     ra   dec  fovra fovdec
      int64 int64 int64 int64 int64
      ----- ----- ----- ----- ------
      0     0     0     5      5
      1    10    10     5      5
      2    20    20     5      5
      3    30    30     5      5
      
      >>> a.readp_coo(ra=[0,10,20,30],dec=[0,10,20,30],fovra=[1,2,3,4],fovdec=[1,2,3,4])      
      >>> a.data
      <Table length=4>
      n     ra   dec  fovra fovdec
      int64 int64 int64 int64 int64
      ----- ----- ----- ----- ------
      0     0     0     1      1
      1    10    10     2      2
      2    20    20     3      3
      3    30    30     4      4

      # if ra/dec is hh:mm:ss format
      >>> a.readp_coo(ra=['12:20:30','0:20:30'],dec=['10:30:00','-5:00:20'],fovra=1, fovdec=1,hms=True,clobber=True)
      >>> a.data
      <Table length=2>
      n           ra                 dec         fovra fovdec
      int64      float64             float64       int64 int64
      ----- ------------------ ------------------- ----- ------
      0 185.12499999999997                10.5     1      1
      1  5.124999999999999 -5.0055555555555555     1      1
      
   .. _rg:
      
   1.2 Then, tilings could be also generated via a region query, e.g. an example is shown below to create a tiling list varied dec from -20 deg to north pole, with 3 by 3 sq deg field of view.
   
   .. code-block:: bash
		
      >>> from kobe import tilings
      >>> a=tilings()
      >>> a.generatep(limdec=[-20, 90], fovra=3, fovdec=3)      
      >>> a.data
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
      6  22.34773   -20.0     3.0     3.0
      ...       ...     ...     ...     ...
      3103 275.36912    85.0     3.0     3.0
      3104 309.79026    85.0     3.0     3.0
      3105  344.2114    85.0     3.0     3.0
      3106  85.96113    88.0     3.0     3.0
      3107 171.92225    88.0     3.0     3.0
      3108 257.88338    88.0     3.0     3.0
      3109  343.8445    88.0     3.0     3.0

   1.3 Since the reference images are important for transient search, telescope pointings are always fixed. Thus, to save and read tilings could be useful for a long term consideration.

     .. code-block:: bash
	
	# Save tilings into a npz file	
	>>> a.savep('tmp_tiles', filetype='npz')
	
	# Read tilings from a npz file
	>>> a=tilings()
	>>> a.readp('tmp_tiles', filetype='npz')
	>>> a.data
	<Table length=3110>
	n        ra      dec    fovra   fovdec
	float64  float64  float64 float64 float64
	------- --------- ------- ------- -------
	0.0   3.19253   -20.0     3.0     3.0
	1.0   6.38507   -20.0     3.0     3.0
	2.0    9.5776   -20.0     3.0     3.0
	3.0  12.77013   -20.0     3.0     3.0
	4.0  15.96267   -20.0     3.0     3.0
	5.0   19.1552   -20.0     3.0     3.0
	6.0  22.34773   -20.0     3.0     3.0
	...       ...     ...     ...     ...
	3103.0 275.36912    85.0     3.0     3.0
	3104.0 309.79026    85.0     3.0     3.0
	3105.0  344.2114    85.0     3.0     3.0
	3106.0  85.96113    88.0     3.0     3.0
	3107.0 171.92225    88.0     3.0     3.0
	3108.0 257.88338    88.0     3.0     3.0
	3109.0  343.8445    88.0     3.0     3.0

     If one had already pre-defined pointing list, check KOBE products and adopt the same scheme.

   1.4. As supplementary of :ref:`1.2 region generation <rg>`, `KOBE` could generate tilings depends on trigger informations if available.

      Let's start with an example of tiling generation based on trigger localization.
      Since for such case, we need both `trigger` and `pointings` utilities,
      thus we should import their son class, i.e. `schedule` instead.
      For more detials of `trigger` utilities, please check the :ref:`trigger chapter <kbtrigger>`.
      
      .. code-block:: bash
	 
	 >>> from kobe import schedule
	 >>> a=schedule()

	 # parse a trigger first		      
	 >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')

	 # initialize a tiling network
	 >>> a.set_pointings(strategy='T')

	 # `KOBE` will run a monte carlo on `shiftra` and `shiftdec`,
	 # in order to maximize the probability coverage of top 100 ranked tilings
	 >>> a.generatep_trigger(100,limra=[10,40], limdec=[20, 40], fovra=3, fovdec=3)
	 <Table length=60>
	 n      ra      dec     fovra   fovdec
	 int64 float64  float64  float64 float64
	 ----- -------- -------- ------- -------
	 0 16.02716 22.96553     3.0     3.0
	 1 19.28541 22.96553     3.0     3.0
	 2 22.54366 22.96553     3.0     3.0
	 3 25.80191 22.96553     3.0     3.0
	 4 29.06016 22.96553     3.0     3.0
	 5 32.31841 22.96553     3.0     3.0
	 ...      ...      ...     ...     ...
	 54 22.85897 40.96553     3.0     3.0
	 55 26.83193 40.96553     3.0     3.0
	 56  30.8049 40.96553     3.0     3.0
	 57 34.77786 40.96553     3.0     3.0
	 58 38.75082 40.96553     3.0     3.0
	 59 42.72378 40.96553     3.0     3.0

      More methods concerning both `trigger` and `pointings`, e.g. check how much trigger probabilities inside a tiling, would be described in the :ref:`trigger chapter <kbtrigger>`.
      
.. _tiles2:
		   
2. show pointings

   After tilings parsed successfully, the `data` attribute is not `nonetype` anymore.
   Therefore, the tiles could be visualized via `KOBE visualization` utils:
   
   .. code-block:: bash

      # do plot
      >>> a.locshow()
      
      # save figure to the current directory (make sure one had the read-write permission), with filename as tiling1.png      
      >>> a.savefig(filepath='tiling1',wdir='./')

   .. image:: ../static/tiling1.png
      :width: 800
      :align: center
      :alt: tiling1.png

   In case one would like to zoom in the sight to somewhere specific, e.g. the north pole,
   you could firstly set a rotation scheme in `locshow` so that the north pole is in the center of plot, and then zoomin to the region:
   
   .. code-block:: bash

      # set clear to clean previous plot
      # set rot_phi to rotate map along dec direction
      # set network color to red
      >>> a.locshow(clear=True,rot_phi=90,color='r')

      # zoom in the region
      >>> a.zoomin([-.4,.4], [-.4,.4])

      # save plot
      >>> a.savefig(filepath='tiling2',wdir='./')

   .. image:: ../static/tiling2.png
      :width: 800
      :align: center
      :alt: tiling2.png

.. _tiles3:
	    
3. skip portion of sky

   Due to various of reasons, e.g. apart from bright sources, or avoid duplicated observation for same field, a portion of sky is sometimes needed to be skipped.
   Here, we show an example of the tiling removement process:

   First, let's generate a series of tiling with 1 by 1 sq deg, ranging ra from 10 to 40, dec from 20 to 40 deg, and zoom in:
   
   .. code-block:: bash
		   
      >>> from kobe import tilings
      >>> a=tilings()
      >>> a.generatep(limra=[10,40], limdec=[20, 40], fovra=1, fovdec=1)
      >>> a.locshow()
      >>> a.zoomin([-.6, 0], [.1,.8])
      >>> a.savefig('skip1')
   
   .. image:: ../static/skip1.png
      :width: 400
      :align: center
      :alt: skip1.png

   Then, assuming a bright source is located at ra=20 deg, dec=30 deg.
   We could generate a 5 by 5 deg rectangle at the position, as the forbiden region.   

   .. code-block:: bash

      >>> b=tilings()
      >>> b.readp_coo(ra=20,dec=30,fovra=5,fovdec=5)
      >>> b.locshow(color='r')
      >>> b.savefig('skip11')
      
   .. image:: ../static/skip11.png
      :width: 400
      :align: center
      :alt: skip11.png

   We would like to delete all fields that have overlaps with the forbiden box.
   Here, `skipfrac` is defined as the threshold on the fraction of one field overlapped by the forbiden region.
   All the fields whose overlapped fraction larger than `skipfrac` would be removed.
      
   .. code-block:: bash

      # set `skipfrac` = 0, to remove tilings that have any overlaps with the forbiden region
      >>> a.removep_coo([20],[30],[5],[5],skipfrac=0,nest=False)

      # show tilings
      >>> a.locshow(clear=True)

      # show forbiden box
      >>> b.locshow(color='r')

      # zoomin and save plot
      >>> a.zoomin([-.6, 0], [.1,.8])
      >>> a.savefig('skip2')

   .. image:: ../static/skip2.png
      :width: 400
      :align: center
      :alt: skip2.png

   .. code-block:: bash

      # set `skipfrac` = 0.5, to remove tilings which was 50% overlapped
      >>> a.removep_coo([20],[30],[5],[5],skipfrac=.5,nest=False)
      >>> a.locshow(clear=True)
      >>> b.locshow(color='r')
      >>> a.zoomin([-.6, 0], [.1,.8])
      >>> a.savefig('skip22')

   .. image:: ../static/skip22.png
      :width: 400
      :align: center
      :alt: skip22.png

   If there're lots of fields are needed to be skipped:
   
   * either, input sequences with `removep_coo` as above
   * or one could put all forbiden regions in a file (with demanding scheme) and adopt `removep_file` method.     

.. _tiles4:
     
4. group and divide pointings

   In `data`, there's one column named `n`, which is the index.
   `KOBE` would treat tilings (e.g. assign observing time to each pointing) with the same index as the same observing block, that means they're assumed to be excuted together, thus, the overhat time is consumed only once.
   The crucial point is to make sure the pointings are spatially close to each other, so that the autofocus and autoguide is reliable.   
   Therefore, we provide functions to group and divide `KOBE` pointings.
   
   For instance, one would like to use a 1 by 1 sq deg FoV telescope to tile the sky.
   He could build up a series of 3 by 3 deg pointings first, and then divide each of them into 9 sub-cells.
   The sub-pointings are then considered as in the same OBs by `KOBE`.   
   
   .. code-block:: bash

      # generate 3 by 3 deg boxes in specific portion of sky
      >>> from kobe import tilings
      >>> a=tilings()
      >>> a.generatep(limra=[0,40], limdec=[20, 40], fovra=3, fovdec=3)
      >>> a.data
      n     ra    dec  fovra fovdec
      --- -------- ---- ----- ------
      0  3.19253 20.0   3.0    3.0
      1  6.38507 20.0   3.0    3.0
      2   9.5776 20.0   3.0    3.0
      3 12.77013 20.0   3.0    3.0
      ...      ...  ...   ...    ...
      73 26.64938 38.0   3.0    3.0
      74 30.45644 38.0   3.0    3.0
      75 34.26349 38.0   3.0    3.0
      76 38.07055 38.0   3.0    3.0
      Length = 77 rows

      # show the tiling network
      >>> a.locshow()      
      >>> a.zoomin([-.6, 0], [.1,.8])
      >>> a.savefig('divide')
      
   .. image:: ../static/divide.png
      :width: 400
      :align: center
      :alt: divide.png

   .. code-block:: bash

      # then, divide each of them into a 3 by 3 deg sub-tilings      
      >>> a.dividep(3,3)

      # as shown, pointings within same OB would have the same index
      >>> a.data
      n     ra    dec  fovra fovdec
      --- -------- ---- ----- ------
      0  2.13491 19.0   1.0    1.0
      0  3.19253 19.0   1.0    1.0
      0  4.25015 19.0   1.0    1.0
      0  2.12835 20.0   1.0    1.0
      ...      ...  ...   ...    ...
      76 39.33957 38.0   1.0    1.0
      76 36.78379 39.0   1.0    1.0
      76 38.07055 39.0   1.0    1.0
      76 39.35731 39.0   1.0    1.0

      # show tiling network again
      Length = 693 rows
      >>> a.locshow(clear=True)      
      >>> a.savefig('divide1')
      
   .. image:: ../static/divide1.png
      :width: 400
      :align: center
      :alt: divide1.png
   
   .. code-block:: bash

      # here, there's a `shown` option provided, to visualize pointings in different OBs with different colors.
      # check more options, e.g. how to set colors and so on, in the corresponding API page.
      >>> a.locshow(shown=True)      
      >>> a.savefig('divide2')
      
   .. image:: ../static/divide2.png
      :width: 400
      :align: center
      :alt: divide2.png

   Meanwhile as an equivalent approach, one could build up sub-cells firstly and then ask `KOBE` to group them:
   
   .. code-block:: bash

      # generate 1 by 1 deg boxes in specific portion of sky
      >>> from kobe import tilings
      >>> a=tilings()
      >>> a.generatep(limra=[0,40], limdec=[20, 40], fovra=1, fovdec=1)

      # group them
      # in such case, `KOBE` would split tiling network into a series of 3 by 3 deg boxes, and consider all pointings inside the same box as in the same OB
      >>> a.groupp(3,3)
      >>> a.data
      n     ra    dec  fovra fovdec
      --- -------- ---- ----- ------
      0  2.13491 19.0   1.0    1.0
      0  3.19253 19.0   1.0    1.0
      0  4.25015 19.0   1.0    1.0
      0  2.12835 20.0   1.0    1.0
      ...      ...  ...   ...    ...
      88 39.33957 38.0   1.0    1.0
      87 36.78379 39.0   1.0    1.0
      87 38.07055 39.0   1.0    1.0
      88 39.35731 39.0   1.0    1.0
      Length = 693 rows

      # show them in different colors
      >>> a.locshow(shown=True)      
      >>> a.savefig('divide3')

   .. image:: ../static/divide3.png
      :width: 400
      :align: center
      :alt: divide3.png


.. _tiles5:
	    
5. rank pointings

   After we have a list of pointings, another important issue that need to be considered
   is to arange them so that telescope could observe them in order.
   Since now, we have just generated pointings (without :ref:`trigger <kbtrigger>` and
   :ref:`observaroty <kbtel2>` objects),
   we would discuss only ranking pointings spatially here
   (i.e. rank from west to east since the west pointings would descend earlier),
   and provide more ranking approaches, e.g. ranking based on trigger probability,
   later when they were parsed.

   In `KOBE`, we use index `n` to represent their priority, so the nature of
   pointing priorization is to change the index according to specific tracers.      
   There're 2 modes to rank pointings from west to east:

   * strict mode: from west to east strictly
   * adjacent mode: start from the pointing which is the westest, and the next pointing is set either adjacent or closest to the previous one, and so on. The aim is to take full advantage of every movements of telescope, in case telescope is not that easy to rotate.

   Let's start with an example of adjacent mode.

   .. code-block:: bash

      # initialize a tilings, a
      >>> from kobe import tilings     
      >>> a=tilings()      
      
      # generate pointings
      >>> a.generatep(fovra=3,fovdec=3,limra=[20, 100], limdec=[0,30])

      # rank pointings
      # approach 1 will start from the westest pointings
      # mode 2 will make pointings adjacent to each other      
      >>> a.rankp(mode=2)

      # show the routes
      >>> a.locshow(color='r')
      >>> a.routeshow(color='k')
      >>> a.zoomin([-1.2,0],[-.5,.6])
      >>> a.savefig('routeshow')

   .. image:: ../static/routeshow.png
      :width: 800
      :align: center
      :alt: routeshow.png
	    

.. _tiles6:
	    
6. report and circulate

   In `KOBE`, there're 3 types of informations defined, i.e. `texts`, `images` and `attachments`.
   One could use various of functions (e.g. `reportp` shown below) to append them, and call `circulate` functions to circulate to user/API/etc at any stage.
   
   .. code-block:: bash
		   
      >>> from kobe import tilings
      >>> a=tilings()
      >>> a.generatep(limra=[10,40], limdec=[20, 40], fovra=3, fovdec=3)

      # append `KOBE texts`
      >>> a.reportp(split=' ')
      >>> a.texts      
      'n ra dec fovra fovdec \n0 10.64178 20.0 1.0 1.0 \n1 11.70596 20.0 1.0 1.0 \n2 12.77013 20.0 1.0 1.0 \n3 13.83431 20.0 1.0 1.0 \n4 14.89849 20.0 1.0 1.0 \n5 15.96267 20.0 1.0 1.0 \n6 17.02684 20.0 1.0 1.0 \n7 18.09102 20.0 1.0 1.0 \n8 19.1552 20.0 1.0 1.0 \n
      ...'

      # call `send_email` to circulate texts to list of customers
      >>> a.send_email(emailpass='123', emailsmtp='smtp.qq.com',
          fromaddr='1@xx.com', toaddrs=['2@xx.com', '3@xx.com'])      
      
.. _galaxies:

Galaxy search
----------------------------------------

For galaxy search, the main difference compared to `tilings` is that the pointings are requested via Vizier instead.

.. _galaxies1:

1. generate pointings

   For pointing generation, `galaxies` is quite similar to `tilings` object:
   
   .. code-block:: bash

      # ask from GLADE for galaxies from 0 to 40 Mpc, with dec above -20 deg
      >>> from kobe import galaxies
      >>> a=galaxies()
      >>> a.generatep(catalog='GLADE', limdec=[-20, 90], limdist=[0,40])      
      >>> a.data
      <Table length=2492>
      n                    name                  ...      dist
      int64                 str62                  ...    float64
      ----- -------------------------------------- ... -------------
      0 28655:NGC3034:NGC3034:09555243+6940469 ... 4.70228466231
      1 42407:NGC4594:NGC4594:12395949-1137230 ... 3.65995814653
      2            --:---:---:12564369+2140575 ... 3.87868462045
      3 41220:NGC4472:NGC4472:12294679+0800014 ... 14.6955181559
      ...                                    ... ...           ...
      2488     --:---:SDSSJ123626.88+255738.2:--- ... 20.7004974885
      2489     --:---:SDSSJ122250.38+155056.9:--- ... 24.8768244613
      2490     --:---:SDSSJ124211.24+074016.0:--- ... 29.7521819633
      2491                     --:---:NGC5496:--- ... 23.4249708242

      # compared to `tiling data`, here each row includes a galaxy name, magnitude and distance, insetead of field of view
      >>> a.data[0]
      <Row index=0>
      n                    name                      ra       dec      mag         dist
      int64                 str62                   float64   float64  float64     float64
      ----- -------------------------------------- --------- --------- -------- -------------
      0 28655:NGC3034:NGC3034:09555243+6940469 148.96846 69.679703 -19.2115 4.70228466231

   Considering it's weird that one could manually provide galaxy informations, thus we don't provide `readp_coo` methods for `galaxies`.
   Meanwhile in case one would like to go for further galaxies, the query of larger catalog will take time. Also, working offline is sometimes unavoidable, thus `savep` and `readp` functions are essential and recommended.
   
   .. code-block:: bash
		
      >>> a.savep('tmp_galaxies', filetype='txt')     
      >>> a.readp('tmp_tiles', filetype='txt')
      
   `KOBE` could generate galaxies depends on trigger informations also:
      
   .. code-block:: bash
		
      >>> from kobe import schedule
      >>> a=schedule()

      # parse a trigger first		      
      >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')

      # initialize a tiling network
      >>> a.set_pointings(strategy='G')

      # `KOBE` will run a monte carlo on `shiftra` and `shiftdec`,
      # in order to maximize the probability coverage of top 100 ranked tilings
      >>> a.generatep_trigger(100,limra=[10,40], limdec=[20, 40], fovra=3, fovdec=3)
	 

.. _galaxies2:

2. show pointings

   The visualization call is similar to `tilings`.
         
   .. code-block:: bash
		
      >>> a.locshow(marker='o',ms=1,clear=True)
      >>> a.savefig(filepath='galaxy1',wdir='./')

   .. image:: ../static/galaxy1.png
      :width: 800
      :align: center
      :alt: galaxy1.png
	    
   .. note::

      for `galaxy.locshow`, marker should be setted anyhow, otherwise, the method will show nothing
     
   In case one would like to konw how galaxy properties distributed, `KOBE` provide methods to show some statistic plots:
   
   .. code-block:: bash

      # galaxy distance distribution
      >>> a.distshow(nbin=100)
      >>> a.savefig('dist.png')

   .. image:: ../static/dist.png
      :width: 800
      :align: center
      :alt: dist.png

   .. code-block:: bash

      # galaxy luminosity distribution
      >>> a.lumshow(nbin=1)
      >>> a.savefig('lum.png')

   .. image:: ../static/lum.png
      :width: 800
      :align: center
      :alt: lum.png

.. _galaxies3:

3. skip portion of sky

   It's quite similar to `tilings` case.
   
   First, we query galaxies from GLADE as pointings and show them.
   
   .. code-block:: bash

      >>> from kobe import galaxies
      >>> a=galaxies()
      >>> a.generatep(catalog='GLADE', limdec=[-20, 90], limdist=[0,40])
      >>> a.locshow(marker='x')
      >>> a.savefig('skipg1')
   
   .. image:: ../static/skipg1.png
      :width: 800
      :align: center
      :alt: skipg1.png

   Then, we generate a 10 by 10 deg polygon at ra=0 deg, dec=60 deg, as the forbiden region.   

   .. code-block:: bash

      >>> b=tilings()
      >>> b.readp_coo(ra=0,dec=60,fovra=10,fovdec=10)
      >>> b.locshow(color='r')
      >>> b.savefig('skipg11')
      
   .. image:: ../static/skipg11.png
      :width: 800
      :align: center
      :alt: skipg11.png

   We would like to delete all fields that inside the forbiden box.
   Considering galaxies have no area, `skipfrac` is not an option for galaxy case.
      
   .. code-block:: bash

      >>> a.removep_coo([0],[60],[10],[10],nest=False)
      >>> a.locshow(clear=True, marker='x')
      >>> b.locshow(color='r')
      >>> a.savefig('skipg2')

   .. image:: ../static/skipg2.png
      :width: 800
      :align: center
      :alt: skipg2.png

.. _galaxies4:

4. group and divide pointings

   Since `galaxy` is considered as a dot by `KOBE`, there will be no `divide` function.
   Here we show how to group the galaxies, which is again more or less the same as `tiling` case.
   
   .. code-block:: bash

      # query galaxies
      >>> from kobe import galaxies
      >>> a=galaxies()
      >>> a.generatep(catalog='GLADE', limdec=[-20, 90], limdist=[0,40])

      # group with every 10 by 10 sq deg region
      >>> a.groupp(10,10)
      >>> a.data
      <Table length=2492>
      n   ...      dist
      int64 ...    float64
      ----- ... -------------
      265 ... 4.70228466231
      18 ... 3.65995814653
      147 ... 3.87868462045
      82 ... 14.6955181559
      206 ... 2.49371678269
      35 ... 18.4142043388
      ... ...           ...
      35 ... 21.2227303453
      146 ... 25.2571985761
      146 ... 20.7004974885
      114 ... 24.8768244613
      82 ... 29.7521819633
      53 ... 23.4249708242

      # show them in different colors
      >>> a.locshow(marker='.', shown=True)      
      >>> a.savefig('divideg')

   .. image:: ../static/divideg.png
      :width: 800
      :align: center
      :alt: divideg.png

.. _galaxies5:

5. rank pointings

   For galaxy search:

   .. code-block:: bash

      # initialize a schedule pointings
      >>> from kobe import pointings     
      >>> a=pointings('G')

      # generate pointings
      >>> a.generatep(limdist=[0,40], limra=[20, 100], limdec=[0,30])

      # rank pointings
      >>> a.rankp(mode=2,threshold=10)

      # show
      >>> a.locshow(color='r')
      >>> a.routeshow(color='k')
      >>> a.zoomin([-1.2,0],[-.5,.6])

   .. image:: ../static/routegshow.png
      :width: 800
      :align: center
      :alt: routegshow.png
	    

.. _galaxies6:

6. report and circulate

   The circulate usage of `galaxies` is the same as `tilings`:
   
   .. code-block:: bash
		
      >>> from kobe import galaxies
      >>> a=galaxies()
      >>> a.generatep(catalog='GLADE', limdec=[-20, 90], limdist=[0,40])
      >>> a.reportp(split=' ')
      >>> a.texts

.. _galaxies7:

7. make healpix map with galaxies

   For the followup search, besides trigger probability, mass distribution is also an important tracer, since it's believed that the violent explosion should be accompanied with galaxies.
   Here, we provide one function that could generate a healpix map for mass distribution, illustrating how local galaxies distributes, which can be then used to weight the trigger map/serve as trigger to telescopes.
   
   .. code-block:: bash

      # construct a list of galaxies
      >>> from kobe import galaxies
      >>> a=galaxies()
      >>> a.generatep(catalog='GLADE', limdec=[-20, 90], limdist=[0,40])

      # construct a healpix map
      # here, tracer l, converts galaxy mass to relative galaxy luminosities
      >>> a.galaxies2hpmap(tracer='l', nside=512, nest=False)
      >>> a.hpmapm

      # show 50% contours of mass healpix map
      array([0., 0., 0., ..., 0., 0., 0.])
      >>> a.hplocshow(cls=[.5])
      >>> a.savefig('galhp')

   .. image:: ../static/galhp.png
      :width: 800
      :align: center
      :alt: galhp.png

	    
.. _point:
	  
Pointings
----------------------------------------

`pointings` is a combination of `tilings` and `galaxies`.
One should initialize a `pointings` object with either option [T]iling, or [G]alaxy.
If option T adopted, then `pointings` object became `tilings` object.

We show a tile pointing case below for example:
   
   .. code-block:: bash

      # initialize a tiling pointings
      >>> from kobe import pointings
      >>> a=pointings('T')
  
      # generate pointings
      >>> a.generatep(fovra=3,fovdec=3,limdec=[-20,90])
      
      # show
      >>> a.locshow()   
      >>> a.savefig(filepath='tiling1',wdir='./')

   .. image:: ../static/tiling1.png
      :width: 800
      :align: center
      :alt: tiling1.png

More examples, parameters and usages are described in the dedicated API page.
   
| :ref:`Previous <kbscheme>`
| :ref:`Next <kbtel>`
| :ref:`Top <kbpoint>`
