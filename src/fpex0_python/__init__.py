"""
This file is executed when the module is imported (i.e. "import fpex0"). Nothing else.
That also means that no files are loaded contained in the folder.
"""

from .fpex0 import *
from .setup import *
from .InitialDistribution import *
from .fokker_planck import FokkerPlanck

from . import process
from . import example