.. _kbcand:
   
:mod:`Candidates`
===========================================

| :ref:`Previous <kbtel>`
| :ref:`Next <kbtrigger>`
| :ref:`1. Candidates <kbcand1>`
|      :ref:`1.1 from list <kbcand11>`
|      :ref:`1.2 from text/npz/fits file <kbcand12>`
|      :ref:`1.3 from avro file <kbcand13>`
|      :ref:`1.4 simulate candidates based on trigger <kbcand13>`
|      :ref:`1.5 rank pointings based on candidates <kbcand14>`
| :ref:`2. Targeting ligutCurve <kbcand2>`

Besides observatories, the other important elements of observations include: candidates and triggers.
In this chapter, we show how `KOBE` would deal with `candidates` objects.

.. note::
   
   Candidates and triggers are all external alerts.
   Currently, main GW/GRB/neutrino/etc detectors, e.g. LIGO, Fermi and so on, distribute triggers via `VoEvent <https://en.wikipedia.org/wiki/VOEvent>`_, while major optical facilities, e.g. ZTF/LSST, would brodcast their data stream via `Avro <https://en.wikipedia.org/wiki/Apache_Avro>`_ format.
   Thus, `KOBE` provide methods to parse triggers from VOevent file, and candidares from Avro file, which might be slightly supplemented later.   

.. _kbcand1:

Candidate list
----------------------------------------

.. _kbcand11:

1. from list

   `KOBE` could generate a list of candidates with lists via `readc`.

   .. code-block:: bash

      # initialize candidates object
      >>> from kobe import candidates
      >>> a=candidates()

      # input ra, dec
      >>> a.readc(ra=[1,10,30], dec=[-30,20,57])

      # visualize them
      >>> a.candshow(marker='x',color='r')
      >>> a.savefig('cand')

   .. image:: ../static/cand.png
      :width: 800
      :align: center
      :alt: cand.png

.. _kbcand12:

2. from text file
   
   As same as `pointings`, one could generate candidates via a file.
   
   .. code-block:: bash
		   
      >>> a.readc_file(filename)   

.. _kbcand13:
      
3. from avro file

   `KOBE` could parse candidates from avro files, and here we show an example of ZTF alerts.

   .. code-block:: bash

      # initialize candidates object
      >>> from kobe import candidates
      >>> a=candidates()

      # download a zipped ZTF alert file and parse it
      >>> a.readc_avro('ztf_public_20190519.tar.gz',rbt=.8)
      >>> a.candidates
      <Table length=1130>
      n        ra        dec            mag               time
      int64   float64    float64        float64           float64
      ----- ----------- ---------- ------------------ ---------------
      0  162.876382 15.3103797 17.297101974487305 2458622.6578356
      1 167.8154547 34.7418207 17.316316604614258 2458622.6592014
      2 170.8729163 35.9226149 17.473085403442383 2458622.6587384
      3 154.7517604 37.8779819 18.452634811401367 2458622.6601273
      4 159.1361093 45.2796218 17.956443786621094 2458622.6610648
      5  158.419986 40.9514874 16.885969161987305 2458622.6601273
      ...         ...        ...                ...             ...
      1123 199.0804272  73.472593  16.56187629699707 2458622.6830787
      1124 206.5294535 73.5271186  18.14190101623535 2458622.6830787
      1125 210.1924989 73.5482038  16.03114891052246 2458622.6830787
      1126 215.3502885 72.8003903 16.690505981445312 2458622.6830787
      1127  213.335147 73.8961953 16.415159225463867 2458622.6830787
      1128 191.5811877 77.1084053 15.541926383972168 2458622.6830787
      1129 192.2536266 77.7498254 16.618961334228516 2458622.6830787
      >>> a.candshow(marker='.',ms=1)
      >>> a.savefig('candavro')

   .. image:: ../static/candavro.png
      :width: 800
      :align: center
      :alt: candavro.png

.. _kbcand14:
	    
4. simulate candidates based on trigger
   
   `KOBE` could simulate a series of candidates depending on trigger informations (which would be described in detail in the :ref:`next chapter <kbtrigger>`):
   
   .. code-block:: bash

      >>> from kobe import schedule
      >>> a = schedule()       
      >>> a.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
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
      >>> a.locshow(cls=[.9])
      >>> a.candshow(marker='x',color='k')
      >>> a.savefig('candl')
      
   .. image:: ../static/candl.png
      :width: 800
      :align: center
      :alt: candl.png

.. _kbcand15:
	    
5. rank pointings based on candidates
   
.. _kbcand2:

Targeting lightcurve
----------------------------------------

`KOBE` generate targeting lightcurves via OAC API:
   
   .. code-block:: bash
		
      >>> from kobe import candidates     
      >>> a=candidates()

      # suppose we're aiming to detect Kilonova
      # query via open supernova catalog
      >>> a.generatelc(tname='AT2017gfo', timeout=60)
      >>> a.lc
      <Table length=825>
      time        magnitude ... e_upper_magnitude
      float64        float64  ...        str4
      ------------------- --------- ... -----------------
      -239.64600000000064     20.44 ...
      -192.62200000000303     21.39 ...
      -191.65700000000652     21.34 ...
      -190.65400000000227     21.26 ...
      -189.64500000000407      21.1 ...
      -188.65600000000268     20.58 ...
      ...       ... ...               ...
      20.98899999999412     21.46 ...
      21.008999999998196     21.48 ...
      24.98899999999412     22.06 ...
      25.008999999998196     20.21 ...
      27.98899999999412     19.96 ...
      28.98899999999412      20.6 ...

      # one can save the lightcurves, for the use of next time
      >>> a.savelc('tmplc')
      
   .. code-block:: bash
		
      >>> from kobe import candidates      
      >>> a=candidates()

      # read lightcurve locally
      >>> a.readlc('tmplc')

      # show lightcurves in r,g,i and u band without upperlimit
      >>> a.lcshow(clear=True,showlim=False,filts='rgiu')
      <Figure size 1280x960 with 1 Axes>      
      >>> a.savefig('lc')

   .. image:: ../static/lc.png
      :width: 800
      :align: center
      :alt: lc.png
	    
:ref:`Previous <kbtel>`
:ref:`Next <kbtrigger>`
