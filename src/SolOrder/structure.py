"""Local structural order on an immutable snapshot of selected oxygen atoms."""

from dataclasses import dataclass
from functools import lru_cache
from math import exp, fsum, lgamma
from numbers import Integral

import numpy as np
from MDAnalysis.lib.distances import minimize_vectors
from scipy.special import sph_harm_y

from . import _kernels as kernels
from .neighbors import build_neighbors, periodic_box, positive
from .parallel import map_chunks, validate_parallel
from .extended import ExtendedOrder
from .descriptors import DescriptorOrder


@dataclass(frozen=True)
class HydrogenBonds:
    """Per-oxygen counts and (donor row, hydrogen topology index, acceptor row).

    A donor may contribute one bond per hydrogen. ``total`` is donated plus
    accepted, so summing it over all waters counts every bond twice.
    """

    donated: np.ndarray
    accepted: np.ndarray
    total: np.ndarray
    bonds: np.ndarray


def _harmonics_chunk(start, stop, offsets, vectors, distances, l):
    lo, hi = offsets[start], offsets[stop]
    v, d = vectors[lo:hi], distances[lo:hi]
    theta = np.arccos(np.clip(v[:, 2] / d, -1., 1.))
    phi = np.mod(np.arctan2(v[:, 1], v[:, 0]), 2.*np.pi)
    y = sph_harm_y(l, np.arange(-l, l+1)[None, :], theta[:, None], phi[:, None])
    return kernels.harmonic_mean(offsets[start:stop+1]-lo, y)


@lru_cache(maxsize=32)
def _wigner_terms(l):
    """Integer (l,l,l) Wigner-3j coefficients using the Racah factorial sum."""
    terms, coefficients = [], []
    log_delta = 3*lgamma(l+1) - lgamma(3*l+2)
    for m1 in range(-l, l+1):
        for m2 in range(-l, l+1):
            m3 = -m1-m2
            if abs(m3) > l:
                continue
            log_prefactor = .5 * (log_delta + sum(
                lgamma(l+m+1)+lgamma(l-m+1) for m in (m1, m2, m3)
            ))
            summands = []
            for z in range(max(0, -m1, m2), min(l, l-m1, l+m2)+1):
                denominator = sum(lgamma(v+1) for v in
                                  (z, l-z, l-m1-z, l+m2-z, m1+z, -m2+z))
                summands.append((-1.)**z * exp(log_prefactor-denominator))
            coefficients.append((-1.)**m3 * fsum(summands))
            terms.append((m1+l, m2+l, m3+l))
    return np.array(terms, dtype=np.int64), np.array(coefficients)


class StructuralOrder(ExtendedOrder, DescriptorOrder):
    """Snapshot analysis with periodic neighbors and configurable joblib workers.

    Distances are in Angstrom, angles in radians unless stated otherwise.
    Every returned per-atom array follows ``atom_indices`` (selection order).
    ``refresh(frame)`` explicitly replaces the snapshot after trajectory changes.
    Threads share numerical arrays; workers never read a live Universe.
    """

    def __init__(self, universe, cutoff=5., frame=0, *, selection="name OW",
                 hydrogen_selection="name HW1 HW2", n_jobs=1, backend="threading"):
        validate_parallel(n_jobs, backend)
        self.u = universe
        self.cutoff = positive(cutoff, "cutoff")
        self.selection = selection
        self.hydrogen_selection = hydrogen_selection
        self.n_jobs, self.backend = n_jobs, backend
        self.refresh(frame)

    def refresh(self, frame=None):
        """Copy one trajectory frame and invalidate cached geometry/results."""
        if frame is not None:
            self.u.trajectory[frame]
        atoms = self.u.select_atoms(self.selection)
        if not len(atoms):
            raise ValueError(f"Selection {self.selection!r} contains no atoms")
        self.dimensions, self.volume = periodic_box(self.u.dimensions)
        self.box = self.dimensions
        self.pos = np.array(atoms.positions, dtype=np.float64, copy=True)
        if not np.isfinite(self.pos).all():
            raise ValueError("Coordinates must be finite")
        self.atom_indices = atoms.indices.copy()
        self._resindices = atoms.resindices.copy()
        self.n = self.n_atoms = len(atoms)
        self._neighbors = {}
        self._voronoi = None
        self._wedge_cache = {}
        self._harmonics = {}
        self._oh = self._hydrogen_indices = None
        # Capture hydrogens on the caller too: later Universe changes cannot
        # mix oxygen and hydrogen coordinates from different frames.
        hydrogens = self.u.select_atoms(self.hydrogen_selection)
        self._h_positions = np.array(hydrogens.positions, dtype=np.float64, copy=True)
        self._h_resindices = hydrogens.resindices.copy()
        self._h_indices = hydrogens.indices.copy()
        for array in (self.pos, self.dimensions, self.atom_indices, self._resindices,
                      self._h_positions, self._h_resindices, self._h_indices):
            array.flags.writeable = False
        for name in ("f4", "dihedrals", "tetra_orient", "tetra_trans", "lsi", "minAngles"):
            self.__dict__.pop(name, None)
        return self

    def _run(self, kernel, size, *args):
        return map_chunks(kernel, size, *args, n_jobs=self.n_jobs, backend=self.backend)

    def neighbor_list(self, cutoff=None):
        """Return a cached, directed, distance-sorted neighbor list."""
        radius = self.cutoff if cutoff is None else positive(cutoff, "cutoff")
        if radius not in self._neighbors:
            graph = build_neighbors(self.pos, self.dimensions, radius)
            for array in (graph.offsets, graph.indices, graph.distances, graph.vectors):
                array.flags.writeable = False
            self._neighbors[radius] = graph
        return self._neighbors[radius]

    def density(self):
        """Global selected-particle number density in Angstrom^-3 (not mass density)."""
        return self.n_atoms / self.volume

    def coordination_number(self, cutoff=None):
        """Number of other selected atoms at distance <= cutoff; self excluded."""
        radius = self.cutoff if cutoff is None else positive(cutoff, "cutoff")
        g = self.neighbor_list(radius)
        return self._run(kernels.coordination_chunk, self.n, g.offsets, g.distances,
                         radius, False, 12, 28).astype(np.int64)

    def smooth_coordination(self, r0=3.5, *, n=12, m=28, cutoff=None):
        """Sum (1-(r/r0)^n)/(1-(r/r0)^m) over the search cutoff.

        At r=r0 the analytic limit n/m is used. The switching function has an
        infinite tail: ``cutoff`` truncates it and is independent of ``r0``.
        """
        r0 = positive(r0, "r0")
        if any(isinstance(v, bool) or not isinstance(v, Integral) for v in (n, m)) or not 0 < n < m:
            raise ValueError("Switch exponents must be integers with 0 < n < m")
        g = self.neighbor_list(cutoff)
        return self._run(kernels.coordination_chunk, self.n, g.offsets, g.distances,
                         r0, True, n, m)

    def kth_neighbor_distance(self, k=5, *, cutoff=None):
        """Distance to kth neighbor; NaN if fewer than k lie inside the cutoff."""
        if isinstance(k, bool) or not isinstance(k, Integral) or k < 1:
            raise ValueError("k must be a positive integer")
        g = self.neighbor_list(cutoff)
        return self._run(kernels.kth_chunk, self.n, g.offsets, g.distances, k)

    def d5(self, *, cutoff=None):
        """Distance to the fifth neighbor in Angstrom; NaN if fewer than five lie inside the cutoff."""
        return self.kth_neighbor_distance(5, cutoff=cutoff)

    def local_average(self, values, *, cutoff=None):
        """Average each value with its neighbors, including the central atom."""
        values = np.asarray(values)
        if values.ndim not in (1, 2) or values.shape[0] != self.n:
            raise ValueError("values must have one row per selected atom")
        values = np.asarray(values, dtype=np.complex128 if np.iscomplexobj(values) else np.float64)
        g = self.neighbor_list(cutoff)
        out = self._run(kernels.average_chunk, self.n, g.offsets, g.indices,
                        values[:, None] if values.ndim == 1 else values)
        return out[:, 0] if values.ndim == 1 else out

    def tetrahedral_order(self, *, averaged=False, cutoff=None):
        """Orientational tetrahedral order Q; NaN for fewer than four neighbors."""
        g = self.neighbor_list(cutoff)
        q = self._run(kernels.tetrahedral_chunk, self.n, g.offsets, g.distances, g.vectors, True)
        return self.local_average(q, cutoff=cutoff) if averaged else q

    def translational_order(self, *, cutoff=None):
        """1 - sum((r_k-mean(r))^2)/(12*mean(r)^2) for four nearest neighbors."""
        g = self.neighbor_list(cutoff)
        return self._run(kernels.tetrahedral_chunk, self.n, g.offsets, g.distances, g.vectors, False)

    def local_structure_index(self, shell=3.7, *, cutoff=None):
        """Variance of radial gaps through the first neighbor outside ``shell``.

        Search cutoff must exceed shell. Missing inner/outer neighbors give NaN.
        """
        shell = positive(shell, "shell")
        radius = self.cutoff if cutoff is None else positive(cutoff, "cutoff")
        if radius <= shell:
            raise ValueError("LSI search cutoff must exceed the shell radius")
        g = self.neighbor_list(radius)
        return self._run(kernels.lsi_chunk, self.n, g.offsets, g.distances, shell)

    def _qlm(self, l, cutoff):
        if isinstance(l, bool) or not isinstance(l, Integral) or not 0 <= l <= 32:
            raise ValueError("l must be an integer between 0 and 32")
        radius = self.cutoff if cutoff is None else positive(cutoff, "cutoff")
        key = (l, radius)
        if key not in self._harmonics:
            g = self.neighbor_list(radius)
            if np.any(g.distances == 0):
                raise ValueError("Bond order is undefined for coincident atoms")
            self._harmonics[key] = self._run(_harmonics_chunk, self.n, g.offsets,
                                              g.vectors, g.distances, l)
        return self._harmonics[key]

    def steinhardt(self, l=6, *, averaged=False, invariant="q", normalized=False, cutoff=None):
        """Steinhardt q_l or w_l; averaging acts on q_lm before the invariant.

        ``normalized=True`` computes w_l / (sum |q_lm|^2)^(3/2). No-neighbor
        environments give NaN. For numerical stability w_l is limited to l<=12.
        """
        if invariant not in ("q", "w"):
            raise ValueError("invariant must be 'q' or 'w'")
        if normalized and invariant != "w":
            raise ValueError("normalized is only applicable to w_l")
        qlm = self._qlm(l, cutoff)
        if averaged:
            qlm = self.local_average(qlm, cutoff=cutoff)
        if invariant == "q":
            return np.sqrt(4.*np.pi/(2*l+1) * np.sum(np.abs(qlm)**2, axis=1))
        if l > 12:
            raise ValueError("w_l is supported for l <= 12")
        terms, coefficients = _wigner_terms(l)
        return self._run(kernels.third_order_chunk, self.n, qlm, terms, coefficients, normalized)

    def crystalline_connections(self, l=6, threshold=.5, *, cutoff=None):
        """Count neighbors with Re(sum(q_i*conj(q_j)))/(|q_i||q_j|) > threshold.

        Zero/undefined harmonic vectors contribute no connections. The threshold
        is an analysis choice, not a universal ice-phase classification rule.
        """
        if not np.isfinite(threshold) or not -1 <= threshold <= 1:
            raise ValueError("threshold must lie in [-1, 1]")
        qlm = self._qlm(l, cutoff)
        g = self.neighbor_list(cutoff)
        return self._run(kernels.connections_chunk, self.n, g.offsets, g.indices, qlm, threshold)

    def _water_geometry(self):
        if self._oh is None:
            if len(np.unique(self._resindices)) != self.n:
                raise ValueError("Water analyses require exactly one selected oxygen per residue")
            hpos = np.empty((self.n, 2, 3))
            hindices = np.empty((self.n, 2), dtype=np.int64)
            by_residue = {}
            for index, residue in enumerate(self._h_resindices):
                by_residue.setdefault(residue, []).append(index)
            for i, residue in enumerate(self._resindices):
                rows = by_residue.get(residue, [])
                if len(rows) != 2:
                    raise ValueError(f"Oxygen row {i} requires exactly two selected hydrogens in its residue")
                hpos[i] = self._h_positions[rows]
                hindices[i] = self._h_indices[rows]
            oh = minimize_vectors((hpos-self.pos[:, None, :]).reshape(-1, 3), self.dimensions).reshape(-1, 2, 3)
            if not np.isfinite(oh).all() or np.any(np.linalg.norm(oh, axis=2) == 0):
                raise ValueError("O-H bond vectors must be finite and nonzero")
            self._oh, self._hydrogen_indices = oh, hindices
        return self._oh

    def hydrogen_bonds(self, cutoff=3.5, angle=30.):
        """Luzar-Chandler geometry: O--O <= cutoff, H--O--O < angle degrees.

        The angle is measured at the donor oxygen between O->H and O->acceptor.
        Hydrogens are assigned by residue membership, not guessed topology bonds.
        """
        cutoff = positive(cutoff, "cutoff")
        angle = positive(angle, "angle")
        if angle > 180:
            raise ValueError("angle must be <= 180 degrees")
        oh = self._water_geometry()
        g = self.neighbor_list(cutoff)
        source = np.repeat(np.arange(self.n), g.counts)
        mask = self._run(kernels.hbond_chunk, len(g.indices), source, g.indices,
                         g.vectors, g.distances, oh, cutoff, np.cos(np.deg2rad(angle)))
        edges, h = np.nonzero(mask)
        donors, acceptors = source[edges], g.indices[edges]
        bonds = np.column_stack((donors, self._hydrogen_indices[donors, h], acceptors))
        donated = np.bincount(donors, minlength=self.n)
        accepted = np.bincount(acceptors, minlength=self.n)
        return HydrogenBonds(donated, accepted, donated+accepted, bonds)
