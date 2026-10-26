"""Water-specific analysis and compatibility with the original SolOrder API."""

import numpy as np
from MDAnalysis.lib.distances import minimize_vectors

from . import _kernels as kernels
from .structure import StructuralOrder


class OrderParameter(StructuralOrder):
    """StructuralOrder with the original calculation names and result attributes.

    Calls snapshot one frame. Use ``refresh(frame)`` for subsequent frames.
    Default ``n_jobs=1``; use -1 for all available joblib workers.
    """

    def __init__(self, universe, cutoff=4.5, frame=0, **kwargs):
        super().__init__(universe, cutoff, frame, **kwargs)

    @property
    def arr(self):
        """Unique neighbor pairs in selected-oxygen row indices."""
        g = self.neighbor_list()
        source = np.repeat(np.arange(self.n), g.counts)
        mask = source < g.indices
        return np.column_stack((source[mask], g.indices[mask]))

    def dist(self, a, b):
        return np.linalg.norm(minimize_vectors(np.asarray(b, dtype=float)[None]-a, self.dimensions)[0])

    def calculate_dihedral(self, a, b, c, d):
        points = np.asarray([a, b, c, d], dtype=np.float64)
        bonds = minimize_vectors(np.diff(points, axis=0), self.dimensions)
        points = np.vstack((np.zeros(3), np.cumsum(bonds, axis=0)))
        return kernels.dihedral_value(*points)

    def neighbours(self, i, arr):
        values = np.unique(arr[np.any(arr == i, axis=1)])
        return values[values != i]

    def _filter_(self, item, arr, pos):
        neighbors = self.neighbours(item, arr)
        distances = np.array([self.dist(pos[item], pos[j]) for j in neighbors])
        order = np.argsort(distances)[:4]
        return neighbors[order], distances[order]

    def maping(self, i):
        return self.u.atoms[self.atom_indices[i:i+1]].resids

    def get_atoms(self, i, j):
        """Legacy residue-ID lookup; duplicated IDs are explicitly ambiguous."""
        self._water_geometry()
        resids = self.u.atoms[self.atom_indices].resids
        rows = []
        for resid in (i, j):
            candidates = np.flatnonzero(resids == resid)
            if len(candidates) != 1:
                raise ValueError("Residue ID lookup must identify exactly one selected oxygen")
            rows.append(candidates[0])
        a, b = rows
        return (self.pos[a]+self._oh[a, 0], self.pos[a]+self._oh[a, 1], self.pos[a],
                self.pos[b]+self._oh[b, 0], self.pos[b]+self._oh[b, 1], self.pos[b])

    def _sort_dihedral(self, h1a, h2a, oa, h1b, h2b, ob):
        oo = minimize_vectors(np.asarray([ob-oa], dtype=float), self.dimensions)[0]
        ha = minimize_vectors(np.asarray([h1a-oa, h2a-oa]), self.dimensions)
        hb = minimize_vectors(np.asarray([h1b-ob, h2b-ob]), self.dimensions)
        ia = np.argmax(np.linalg.norm(ha-oo, axis=1))
        ib = np.argmax(np.linalg.norm(hb+oo, axis=1))
        return np.array([ha[ia], np.zeros(3), oo, oo+hb[ib]])

    def _directed_dihedrals(self, cutoff=None):
        oh = self._water_geometry()
        g = self.neighbor_list(cutoff)
        source = np.repeat(np.arange(self.n), g.counts)
        values = self._run(kernels.dihedrals_chunk, len(g.indices), source, g.indices, g.vectors, oh)
        return g, source, values

    def get_dihedrals(self):
        g, source, values = self._directed_dihedrals()
        self.dihedrals = values[source < g.indices]
        return self.dihedrals

    def F4(self, cutoff=None):
        g, _, angles = self._directed_dihedrals(cutoff)
        self.f4 = self._run(kernels.row_mean_chunk, self.n, g.offsets, np.cos(3.*angles))
        return self.f4

    def OTO(self):
        self.tetra_orient = self.tetrahedral_order()
        return self.tetra_orient

    def TTO(self):
        self.tetra_trans = self.translational_order()
        return self.tetra_trans

    def LSI(self):
        self.lsi = self.local_structure_index()
        return self.lsi

    def MinimumAngle(self):
        g = self.neighbor_list()
        self.minAngles = self._run(kernels.minimum_angle_chunk, self.n, g.offsets,
                                  g.indices, g.vectors, g.distances, self._water_geometry())
        return self.minAngles

    # Retain the original per-particle helpers for callers using them directly.
    def singleF4(self, i, arr):
        values = []
        oh = self._water_geometry()
        for j in self.neighbours(i, arr):
            oo = minimize_vectors(self.pos[j:j+1]-self.pos[i], self.dimensions)[0]
            values.append(np.cos(3.*kernels.water_dihedral(oo, oh[i], oh[j])))
        return np.mean(values) if values else np.nan

    def singleOTO(self, item, arr):
        neighbors, distances = self._filter_(item, arr, self.pos)
        if len(neighbors) < 4 or np.any(distances == 0):
            return np.nan
        vectors = minimize_vectors(self.pos[neighbors]-self.pos[item], self.dimensions)
        units = vectors / distances[:, None]
        cosines = np.clip(units @ units.T, -1., 1.)
        return 1. - 3./8.*np.sum((cosines[np.triu_indices(4, 1)]+1./3.)**2)

    def singleTTO(self, item):
        return self.translational_order()[item]

    def singleLSI(self, item):
        return self.local_structure_index()[item]

    def singleMinimumAngle(self, item):
        g = self.neighbor_list()
        if g.counts[item] == 0:
            return np.nan, -1
        return self.MinimumAngle()[item], g.indices[g.offsets[item]]
