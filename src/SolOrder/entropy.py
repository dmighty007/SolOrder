"""Parallel local pair-entropy fingerprints, expressed in units of k_B."""

import numpy as np
from MDAnalysis.lib.distances import minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors

from . import _kernels as kernels
from .neighbors import build_neighbors, positive
from .structure import StructuralOrder


def dist(a, b, box):
    """Minimum-image distance; accepts lengths alone or six box components."""
    box = np.asarray(box, dtype=np.float64)
    if box.shape == (3,):
        box = np.r_[box, [90., 90., 90.]]
    return np.linalg.norm(minimize_vectors(np.asarray([np.asarray(b)-a], dtype=float), box)[0])


def neighbours(i, arr):
    values = np.unique(arr[np.any(arr == i, axis=1)])
    return values[values != i]


class eFingerprint(StructuralOrder):
    """Entropy of oxygen centers using oxygen-only or all-atom radial density.

    ``local`` selects local versus global number-density normalization. Results
    are in units of k_B. An isolated center has undefined local density and gives
    NaN in local mode; global mode has the finite g(r)=0 limit. ``localEntropy``
    averages only oxygen-center fingerprints, even when hydrogen contributes to
    the radial density. Uniform neighbor averaging includes the central oxygen.
    """

    def __init__(self, universe, frame=0, include_hydrogen=False, cutoff=5.,
                 sigma=.15, local=True, *, selection="name OW", n_jobs=1,
                 backend="threading"):
        self.sigma = positive(sigma, "sigma")
        self.local = bool(local)
        self.include_hydrogen = bool(include_hydrogen)
        self.oxygen_selection = selection
        super().__init__(universe, cutoff, frame,
                         selection="all" if include_hydrogen else selection,
                         n_jobs=n_jobs, backend=backend)
        if self.sigma > self.cutoff:
            raise ValueError("sigma must not exceed cutoff")

    def refresh(self, frame=None):
        super().refresh(frame)
        oxygen_indices = self.u.select_atoms(self.oxygen_selection).indices
        self.ox_ids = np.flatnonzero(np.isin(self.atom_indices, oxygen_indices))
        if not len(self.ox_ids):
            raise ValueError("Oxygen selection contains no atoms")
        cell = triclinic_vectors(self.dimensions, dtype=np.float64)
        heights = self.volume / np.linalg.norm(np.cross(cell[[1, 2, 0]], cell[[2, 0, 1]]), axis=1)
        if self.cutoff > .5*np.min(heights):
            raise ValueError("Entropy cutoff must fit within half the smallest periodic cell height")
        self.global_rho = self.density()
        self.r = np.linspace(0., self.cutoff, int(self.cutoff/self.sigma)+1)
        self.rsq = self.r**2
        self.e = self.le = None
        return self

    def Entropy(self):
        g = self.neighbor_list()
        self.e = self._run(kernels.entropy_chunk, len(self.ox_ids), self.ox_ids,
                           g.offsets, g.distances, self.r, self.sigma,
                           self.global_rho, self.local, self.cutoff)
        return self.e

    def localEntropy(self):
        e = self.Entropy()
        g = build_neighbors(self.pos[self.ox_ids], self.dimensions, self.cutoff)
        self.le = self._run(kernels.average_chunk, len(e), g.offsets, g.indices, e[:, None])[:, 0]
        return self.le
