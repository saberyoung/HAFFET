# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
"""
from .__version__ import version, description
__version__ = version
__description__ = description

__all__ = [
    'snerun',           # main programmer
    'sneGui',           # a GUI for the main programmer
    'functions',        # different functions, e.g. blackbody, etc      
    'gaussian_process', # gaussian process
    'corner_hack',      # for corner plots
    'read_default',     # read default parameters or authorization
    'plot_atlas_fp',    # from ATLAS page, codes to bin atlas forced photometry
    'models',           # different models        
    'image_tool',       # class to handle images  
    'trigger_tool',     # class to handle triggers
    'model_fitters',    # class for fitting models        
    'specline_handler', # class for spectra
    'pbar',             # a progress bar class
    'filters',          # filter informations
    'constants',        # define constants and boundries
]

from . import functions, image_tool, gaussian_process, corner_hack, read_default, plot_atlas_fp, model_fitters, sneGui, snerun, specline_handler, pbar, filters, constants
