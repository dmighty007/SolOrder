# This is so that you can import ppack or import average from ppack
# in stead of from ppack.functions import average

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.lib.distances import calc_angles, calc_bonds
import tqdm.auto as tqdm
from numba import njit

from .dihedral import OrderParameter
from .entropy import eFingerprint
