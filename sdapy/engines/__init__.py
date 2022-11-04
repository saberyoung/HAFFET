"""Initilization modules."""
import os
import numpy as np

path = os.path.dirname(os.path.abspath(__file__))
__all__ = []

for model in [
        f for f in os.listdir(path)
        if f.endswith('.py') and not f.startswith('__')
]:   
    __all__.append(model[:-3])

from .multiband_early import *
from .multiband_main import *
from .bol_early import *
from .bol_main import *
from .bol_tail import *
from .bol_full import *
from .specline import *
from .specv_evolution import *
from .sed import *
