"""Periodic neighbor geometry shared by the order-parameter kernels."""

from dataclasses import dataclass

import numpy as np
from MDAnalysis.lib.distances import minimize_vectors, self_capped_distance
from MDAnalysis.lib.mdamath import triclinic_vectors


@dataclass(frozen=True)
class NeighborList:
    """Directed CSR neighbors sorted by distance, then atom index.

    ``vectors`` point from the central atom to its neighbor. Rows refer to the
    selected oxygen order, rather than topology IDs or residue numbers.
    """

    offsets: np.ndarray
    indices: np.ndarray
    distances: np.ndarray
    vectors: np.ndarray

    @property
    def counts(self):
        return np.diff(self.offsets)


def positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def periodic_box(dimensions):
    if dimensions is None:
        raise ValueError("A periodic box with lengths and angles is required")
    box = np.array(dimensions, dtype=np.float64, copy=True)
    if box.shape != (6,) or not np.isfinite(box).all():
        raise ValueError("Box must contain six finite lengths and angles")
    if np.any(box[:3] <= 0) or np.any(box[3:] <= 0) or np.any(box[3:] >= 180):
        raise ValueError("Box lengths and angles are invalid")
    volume = abs(np.linalg.det(triclinic_vectors(box, dtype=np.float64)))
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError("Box must have positive volume")
    return box, volume


def build_neighbors(positions, box, cutoff):
    cutoff = positive(cutoff, "cutoff")
    pairs, _ = self_capped_distance(positions, max_cutoff=cutoff, box=box)
    if len(pairs):
        source = np.concatenate([pairs[:, 0], pairs[:, 1]])
        target = np.concatenate([pairs[:, 1], pairs[:, 0]])
        vectors = minimize_vectors(positions[target] - positions[source], box)
        distances = np.linalg.norm(vectors, axis=1)
        order = np.lexsort((target, distances, source))
        source, target = source[order], target[order]
        vectors, distances = vectors[order], distances[order]
    else:
        source = target = np.empty(0, dtype=np.int64)
        vectors = np.empty((0, 3), dtype=np.float64)
        distances = np.empty(0, dtype=np.float64)
    offsets = np.zeros(len(positions) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(np.bincount(source, minlength=len(positions)))
    return NeighborList(offsets, np.asarray(target, dtype=np.int64), distances, vectors)
