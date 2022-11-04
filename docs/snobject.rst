.. _snobject:

:mod:`snobject` -- class deal with one single object
=================================================================================

.. currentmodule:: sdapy.snerun
.. autosummary::
   :toctree: generated/

   snobject   
   
Below provide various of functions for `snobject`:

- General

.. currentmodule:: sdapy.snerun
.. autosummary::
   :toctree: generated/   

   snobject.__init__
   snobject.version
   snobject.urllist
   snobject.keys2query_lc
   snobject.keys2query_spec   
   snobject.read_kwargs   
   snobject.run
   
- Lightcurves:

.. currentmodule:: sdapy.snerun
.. autosummary::
   :toctree: generated/

   snobject.parse_coo
   snobject.mjd_now
   
   snobject.query_fp_atlas
   snobject.query_fp_ztf
   snobject.query_alert_ztf
   snobject.query_oac
   
   snobject.get_alert_ztf
   snobject.get_fp_ztf
   snobject.get_fp_atlas
   snobject.get_oac
   snobject.get_external_phot
   
   snobject.bin_lc
   snobject.clip_lc
   snobject.cut_lc
   
   snobject.correct_baseline
   snobject.calibrate_baseline
   snobject.rapid
   snobject.rapid_plot

   snobject._nepochs
   snobject._ncolors
   snobject._peak_accuracy
   snobject._earlypoints
   
- Spectra:

.. currentmodule:: sdapy.snerun
.. autosummary::
   :toctree: generated/

   snobject.query_spectra   
   snobject.get_local_spectra
   snobject.get_external_spectra
   
- SED ann bolometric LCs:

.. currentmodule:: sdapy.snerun
.. autosummary::
   :toctree: generated/
   
   snobject.dm_error
   snobject.read_c10
   snobject.est_hostebv_with_colours
   snobject.calc_colors
   snobject.lyman_bol
   snobject.bb_colors
   snobject.bb_bol   
   snobject._to_luminosity
   snobject._cal_spectrum

- Fittings:

.. currentmodule:: sdapy.snerun
.. autosummary::
   :toctree: generated/

   snobject.run_gp
   snobject.run_fit
   
   snobject.set_peak_gp
   snobject.set_peak_bol_main
   snobject.set_peak_multiband_main

   snobject.set_vexp
   snobject.set_texp_pl
   snobject.set_texp_bol_main
   snobject.set_texp_midway   

   snobject._flux_at
   snobject._flux_at_list
   snobject._mag_at
   snobject._mag_at_list
   snobject._absmag_at
   snobject._absmag_at_list
   snobject._color_at
   snobject._rate_at

   snobject.all_fittings
   snobject.get_model
   snobject.get_par   
   
- Plots:

.. currentmodule:: sdapy.snerun
.. autosummary::
   :toctree: generated/
   
   snobject.summarize_plot
   snobject._ax
   snobject._ax1
   snobject._ax2
   snobject._ax3
   snobject._ax4
   snobject._ax5
   snobject._ax6
   
   snobject.show_corner
   snobject.savefig
   snobject.showfig
   
- Others

.. currentmodule:: sdapy.snerun
.. autosummary::
   :toctree: generated/

   snobject.add_lc
   snobject.add_flux
   snobject._add_flux
   snobject.add_mag
   snobject._add_mag
   snobject.mag_to_flux
   snobject.flux_to_mag
   snobject.bin_fp_atlas
   snobject._get_fp_atlas
   snobject.oac_phot_url
   snobject._get_par
   snobject.sym_mag
   snobject.match_colors
   snobject.combine_multi_obs
   snobject.bin_df
   snobject.merge_df_cols
   read_kwargs
   print_logo
   check_dir

There're some other functions used by the `snobject`:

.. currentmodule:: sdapy
.. autosummary::
   :toctree: generated/
   
   read_default.read_default
   read_default.get_keypairs
   read_default.get_parameters   
   gaussian_process.fit_gp
   model_fitters.get_engine
   model_fitters.get_model
   model_fitters.get_pars
   model_fitters.fit_model
   specline_handler.handle_spectra
   specline_handler.handle_spectrum
   image_tool.handle_image
