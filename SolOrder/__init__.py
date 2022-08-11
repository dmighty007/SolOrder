# This is so that you can import ppack or import average from ppack
# in stead of from ppack.functions import average

import numpy as np_
import pandas as pd_
import MDAnalysis as mda_
from MDAnalysis.lib.nsgrid import FastNS as FastNS_
from MDAnalysis.lib.distances import calc_angles, calc_bonds
import tqdm.auto as tqdm_
from numba import njit as _njit

from .dihedral import OrderParameter
from .entropy import eFingerprint

# Version of the package
__version__ = "1.0.3"
