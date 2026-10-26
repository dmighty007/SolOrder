"""Optional non-neural descriptors with explicit species and periodic cells."""

from dataclasses import dataclass
import importlib
from copy import deepcopy
from numbers import Integral

import numpy as np
from MDAnalysis.lib.mdamath import triclinic_vectors


def _local_descriptor_chunk(start, stop, descriptor, atoms, centers):
    descriptor = deepcopy(descriptor)
    if stop == start:
        return np.empty((0, descriptor.get_number_of_features()))
    return np.asarray(descriptor.create(atoms, centers=centers[start:stop].tolist(), n_jobs=1))


@dataclass(frozen=True)
class GraphDescriptor:
    """Directed periodic graph with invariant distance-based features.

    Node columns: degree, mean bond distance, bond-distance standard deviation.
    Edge column: minimum-image distance. Isolated distance features are NaN.
    """
    atom_indices: np.ndarray
    edge_index: np.ndarray
    node_features: np.ndarray
    edge_features: np.ndarray


@dataclass(frozen=True)
class TemplateMatch:
    labels: np.ndarray
    rmsd: np.ndarray
    type_names: dict


class DescriptorOrder:
    def atomic_descriptor(self, kind, *, symbols, centers=None, **settings):
        """DScribe descriptor of selected snapshot atoms, with optional joblib.

        kinds: soap, mbtr, lmbtr, acsf, coulomb, sine, ewald. ``symbols`` is one
        explicit chemical symbol for all atoms or a sequence in selection order.
        Local outputs follow requested center row indices (all rows by default).
        Global matrix/MBTR descriptors return a single flattened vector.
        Settings are passed to DScribe, including required basis/grid choices.
        """
        try:
            from ase import Atoms
            module = importlib.import_module("dscribe.descriptors")
        except ImportError as exc:
            raise ImportError('Install SolOrder[descriptors] for DScribe/ASE methods') from exc
        names = {"soap": "SOAP", "mbtr": "MBTR", "lmbtr": "LMBTR", "acsf": "ACSF",
                 "coulomb": "CoulombMatrix", "sine": "SineMatrix", "ewald": "EwaldSumMatrix"}
        if kind not in names:
            raise ValueError(f"Unknown descriptor {kind!r}; choose from {tuple(names)}")
        if isinstance(symbols, str):
            symbols = [symbols]*self.n
        if len(symbols) != self.n:
            raise ValueError("Provide one chemical symbol per selected atom")
        atoms = Atoms(symbols=symbols, positions=self.pos.copy(),
                      cell=triclinic_vectors(self.dimensions, dtype=float), pbc=True)
        local = kind in ("soap", "lmbtr", "acsf")
        options = dict(settings)
        if kind in ("soap", "mbtr", "lmbtr", "acsf"):
            options.setdefault("species", sorted(set(symbols)))
            options.setdefault("periodic", True)
            options.setdefault("sparse", False)
            if local and options.get("average", "off") != "off":
                raise ValueError("Local descriptor output requires average='off'")
        else:
            options.setdefault("n_atoms_max", self.n)
            if options["n_atoms_max"] < self.n:
                raise ValueError("n_atoms_max must cover all selected atoms")
        descriptor = getattr(module, names[kind])(**options)
        if local:
            center_rows = np.arange(self.n) if centers is None else np.asarray(centers)
            if center_rows.ndim == 1 and center_rows.size == 0:
                center_rows = center_rows.astype(np.int64)
            if center_rows.ndim != 1 or not np.issubdtype(center_rows.dtype, np.integer) or np.any(center_rows < 0) or np.any(center_rows >= self.n):
                raise ValueError("centers must contain valid selected-atom row indices")
            # Sparse objects do not satisfy the ordered numeric chunk contract.
            if options.get("sparse", False):
                raise ValueError("Local parallel descriptors require sparse=False")
            return self._run(_local_descriptor_chunk, len(center_rows), descriptor, atoms, center_rows)
        if centers is not None:
            raise ValueError("centers only applies to local descriptors")
        return descriptor.create(atoms, n_jobs=1)

    def soap(self, *, symbols="O", r_cut=5., n_max=8, l_max=6, sigma=.5, centers=None, **settings):
        return self.atomic_descriptor("soap", symbols=symbols, centers=centers,
                                      r_cut=r_cut, n_max=n_max, l_max=l_max, sigma=sigma, **settings)

    def graph_descriptor(self, cutoff=None):
        """Return a periodic graph; no learned model or node ordering assumption."""
        graph = self.neighbor_list(cutoff)
        source = np.repeat(np.arange(self.n), graph.counts)
        mean = np.full(self.n, np.nan)
        std = mean.copy()
        for i in range(self.n):
            values = graph.distances[graph.offsets[i]:graph.offsets[i+1]]
            if len(values):
                mean[i], std[i] = values.mean(), values.std()
        return GraphDescriptor(self.atom_indices.copy(), np.vstack((source, graph.indices)),
                               np.column_stack((graph.counts, mean, std)), graph.distances[:, None].copy())

    def template_matching(self, *, templates=("CUBIC_DIAMOND", "HEX_DIAMOND"), rmsd_cutoff=.1):
        """Optional OVITO PTM on the selected periodic snapshot.

        Returns OVITO structure codes and scale-normalized RMSD. Diamond
        templates classify oxygen networks; these are not proton-order labels.
        OVITO controls its own compiled worker pool, independently of n_jobs.
        """
        try:
            from ovito.data import DataCollection, Particles, SimulationCell
            from ovito.pipeline import Pipeline, StaticSource
            from ovito.modifiers import PolyhedralTemplateMatchingModifier
        except ImportError as exc:
            raise ImportError('Install SolOrder[ptm] for polyhedral template matching') from exc
        if not np.isfinite(rmsd_cutoff) or rmsd_cutoff <= 0:
            raise ValueError("rmsd_cutoff must be finite and positive")
        modifier = PolyhedralTemplateMatchingModifier(output_rmsd=True, rmsd_cutoff=rmsd_cutoff)
        ids = []
        for name in templates:
            if not isinstance(name, str) or not hasattr(modifier.Type, name) or name == "OTHER":
                raise ValueError(f"Unknown PTM template {name!r}")
            ids.append(int(getattr(modifier.Type, name)))
        if not ids:
            raise ValueError("Select at least one PTM template")
        for item in modifier.structures:
            item.enabled = item.id in ids
        data = DataCollection()
        particles = Particles(count=self.n)
        particles.create_property("Position", data=self.pos.copy())
        data.objects.append(particles)
        cell = SimulationCell(pbc=(True, True, True))
        cell[...] = np.column_stack((triclinic_vectors(self.dimensions, dtype=float).T, np.zeros(3)))
        data.objects.append(cell)
        pipeline = Pipeline(source=StaticSource(data=data))
        pipeline.modifiers.append(modifier)
        result = pipeline.compute()
        return TemplateMatch(np.array(result.particles["Structure Type"]),
                             np.array(result.particles["RMSD"]),
                             {item.id: item.name for item in modifier.structures})
