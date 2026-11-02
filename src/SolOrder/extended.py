"""Additional manuscript parameters and periodic ice classifiers."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from MDAnalysis.lib.distances import minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors
from scipy.special import sph_harm_y

from . import _extended_kernels as kernels
from . import _kernels as base_kernels


@dataclass(frozen=True)
class AngularDistribution:
    angles: np.ndarray
    pairs: np.ndarray
    counts: np.ndarray
    bin_edges: np.ndarray
    probability_density: np.ndarray


@dataclass(frozen=True)
class VoronoiGeometry:
    volume: np.ndarray
    area: np.ndarray
    density: np.ndarray
    asphericity: np.ndarray


@dataclass(frozen=True)
class IceClassification:
    """Codes: 0 other, 1 Ih, 2 Ic, 3 ice interface, 4 hydrate, 5 hydrate interface."""
    labels: np.ndarray
    staggered: np.ndarray
    eclipsed: np.ndarray
    correlations: np.ndarray
    offsets: np.ndarray
    neighbors: np.ndarray

    @property
    def names(self):
        return np.array(["other", "hexagonal_ice", "cubic_ice", "interfacial_ice",
                         "hydrate", "interfacial_hydrate"])[self.labels]


@dataclass(frozen=True)
class RingStatistics:
    """Smallest ring through every bond angle j-i-k at each center.

    ``counts[i, s-3]`` counts angles whose smallest ring has ``sizes[s-3] = s``
    members; ``unclosed`` counts angles with no ring up to the maximum size.
    """
    sizes: np.ndarray
    counts: np.ndarray
    unclosed: np.ndarray

    @property
    def fractions(self):
        """Fraction of each center's angles per ring size; NaN without angles."""
        total = self.counts.sum(axis=1)+self.unclosed
        with np.errstate(invalid="ignore", divide="ignore"):
            return self.counts/total[:, None]

    @property
    def mean_size(self):
        """Mean smallest-ring size over closed angles; NaN if none closed."""
        closed = self.counts.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return self.counts @ self.sizes / closed


def radial_values(function, radius):
    """Evaluate a caller's radial function or integer power outside JIT code."""
    if callable(function):
        values = np.asarray(function(radius), dtype=float)
        values = np.broadcast_to(values, radius.shape)
    elif isinstance(function, Integral) and not isinstance(function, bool):
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            values = radius**int(function)
    else:
        raise ValueError("Radial functions must be callables or integer powers")
    if not np.isfinite(values).all():
        raise ValueError("Radial function returned non-finite values")
    return values


def _cell_geometry(dimensions):
    """Box rows, their inverse and half the smallest cell height (for minimum images)."""
    cell = triclinic_vectors(dimensions, dtype=np.float64)
    volume = abs(np.linalg.det(cell))
    heights = [volume/np.linalg.norm(np.cross(cell[(i+1) % 3], cell[(i+2) % 3])) for i in range(3)]
    return cell, np.linalg.inv(cell), .5*min(heights)


class ExtendedOrder:
    """Mixin whose workers operate only on snapshot arrays."""

    def water_angular_distribution(self, cutoff=3.5, *, bins=90):
        """Minimum of four H--O--O angles per unique O--O pair, in radians.

        All pairs within cutoff are sampled, irrespective of hydrogen bonding.
        Histogram spans [0,pi]; density integrates to one for nonempty samples.
        An empty sample returns zero counts and NaN probability densities.
        """
        graph = self.neighbor_list(cutoff)
        source = np.repeat(np.arange(self.n), graph.counts)
        selected = source < graph.indices
        pairs = np.column_stack((source[selected], graph.indices[selected]))
        angles = self._run(kernels.pair_angle_chunk, int(selected.sum()),
                           source[selected], graph.indices[selected], graph.vectors[selected],
                           graph.distances[selected], self._water_geometry())
        if isinstance(bins, Integral) and not isinstance(bins, bool):
            if bins <= 0:
                raise ValueError("bins must be positive")
            edges = np.linspace(0, np.pi, bins+1)
        else:
            edges = np.asarray(bins, dtype=float)
            if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
                raise ValueError("Histogram edges must be finite and strictly increasing")
            if edges[0] != 0 or not np.isclose(edges[-1], np.pi, rtol=0, atol=1e-12):
                raise ValueError("Histogram edges must cover [0, pi]")
        finite = angles[np.isfinite(angles)]
        counts, _ = np.histogram(finite, bins=edges)
        density = counts/(len(finite)*np.diff(edges)) if len(finite) else np.full(len(edges)-1, np.nan)
        return AngularDistribution(angles, pairs, counts, edges, density)

    def centrosymmetry(self, n_neighbors=12, *, pairing="conventional", averaged=False, cutoff=None):
        """CSP in Angstrom^2. Missing neighbors give NaN.

        ``conventional`` sums the N/2 smallest pair weights, as in LAMMPS.
        ``rank`` implements the manuscript's first-half/second-half rank pairs.
        """
        if isinstance(n_neighbors, bool) or not isinstance(n_neighbors, Integral) or n_neighbors < 2 or n_neighbors % 2:
            raise ValueError("n_neighbors must be a positive even integer >= 2")
        if pairing not in ("conventional", "rank"):
            raise ValueError("pairing must be 'conventional' or 'rank'")
        graph = self.neighbor_list(cutoff)
        out = self._run(kernels.centrosymmetry_chunk, self.n, graph.offsets,
                        graph.vectors, n_neighbors, int(pairing == "rank"))
        return self.local_average(out, cutoff=cutoff) if averaged else out

    def _wedges(self, cutoff=None):
        graph = self.neighbor_list(cutoff)
        key = self.cutoff if cutoff is None else float(cutoff)
        if key in self._wedge_cache:
            return self._wedge_cache[key]
        if np.any(graph.distances == 0):
            raise ValueError("Neighbor angles are undefined for coincident atoms")
        left, right = [], []
        offsets = np.zeros(self.n+1, dtype=np.int64)
        for i, n in enumerate(graph.counts):
            a, b = np.triu_indices(n, 1)
            left.append(a+graph.offsets[i])
            right.append(b+graph.offsets[i])
            offsets[i+1] = offsets[i]+len(a)
        left, right = np.concatenate(left), np.concatenate(right)
        data = self._run(kernels.wedge_chunk, len(left), left, right, graph.vectors, graph.distances)
        opposite = minimize_vectors(graph.vectors[left]-graph.vectors[right], self.dimensions)
        result = (offsets, data[:, 0], data[:, 1], data[:, 2], np.linalg.norm(opposite, axis=1))
        for array in result:
            array.flags.writeable = False
        self._wedge_cache[key] = result
        return result

    def _wedge_average(self, offsets, values, averaged, cutoff):
        out = self._run(base_kernels.row_mean_chunk, self.n, offsets, np.asarray(values, dtype=float))
        return self.local_average(out, cutoff=cutoff) if averaged else out

    def neighbor_distance_parameter(self, f_alpha=1, f_beta=1, f_gamma=1, *,
                                    symmetrize=True, averaged=False, cutoff=None):
        """Mean f_alpha(r_ij) f_beta(r_ik) f_gamma(r_jk) over neighbor pairs.

        Integer functions specify radial powers; callables receive arrays in A.
        By default alpha/beta products are symmetrized under neighbor exchange.
        Set symmetrize=False for the manuscript's distance-ranked pair ordering.
        """
        offsets, ra, rb, _, rc = self._wedges(cutoff)
        values = radial_values(f_alpha, ra)*radial_values(f_beta, rb)
        if symmetrize:
            values = .5*(values+radial_values(f_alpha, rb)*radial_values(f_beta, ra))
        values *= radial_values(f_gamma, rc)
        return self._wedge_average(offsets, values, averaged, cutoff)

    def angular_fourier(self, a=1., radial=0, *, radial_beta=None, weighting="min",
                        symmetrize=True, averaged=False, cutoff=None):
        """Angular Fourier descriptor with explicit radial-weight convention.

        weighting='min': f(min(r_ij,r_ik))*cos(a*theta), as in the manuscript.
        weighting='product': f(r_ij)*g(r_ik)*cos(a*theta), as in MALIO; the
        product is symmetrized by default. Integers specify radial powers.
        """
        if not np.isfinite(a):
            raise ValueError("Angular frequency must be finite")
        offsets, ra, rb, angles, _ = self._wedges(cutoff)
        if weighting == "min":
            if radial_beta is not None:
                raise ValueError("radial_beta requires weighting='product'")
            weights = radial_values(radial, np.minimum(ra, rb))
        elif weighting == "product":
            beta = radial if radial_beta is None else radial_beta
            weights = radial_values(radial, ra)*radial_values(beta, rb)
            if symmetrize:
                weights = .5*(weights+radial_values(radial, rb)*radial_values(beta, ra))
        else:
            raise ValueError("weighting must be 'min' or 'product'")
        return self._wedge_average(offsets, weights*np.cos(a*angles), averaged, cutoff)

    def bond_angle_parameter(self, n1=1, n2=1, phi=0., *, averaged=False, cutoff=None):
        """Mean cos(n2*theta+phi)^n1; radians, following Takahashi (2023)."""
        if any(isinstance(n, bool) or not isinstance(n, Integral) or n < 1 for n in (n1, n2)):
            raise ValueError("n1 and n2 must be positive integers")
        if not np.isfinite(phi):
            raise ValueError("phi must be finite")
        offsets, _, _, angles, _ = self._wedges(cutoff)
        return self._wedge_average(offsets, np.cos(n2*angles+phi)**n1, averaged, cutoff)

    def angle_histogram_parameter(self, frequencies=(1, 2, 3, 4), *, bins=180,
                                  power=False, averaged=False, cutoff=None):
        """Normalized DFT amplitudes of per-atom neighbor-angle histograms.

        Equal bins cover [0,pi]. Divide histogram counts by the pair count before
        rFFT, so the DC channel equals one. Frequencies are integer DFT indices.
        No-pair environments give NaN. This explicitly fixes the manuscript's
        unspecified histogram/FFT conventions. power=True returns squared
        amplitudes, matching the author's MALIO implementation.
        """
        if isinstance(bins, bool) or not isinstance(bins, Integral) or bins < 2:
            raise ValueError("bins must be an integer >= 2")
        freq = np.asarray(frequencies)
        if freq.ndim != 1 or not np.issubdtype(freq.dtype, np.integer) or np.any(freq < 0) or np.any(freq > bins//2):
            raise ValueError("Frequencies must be integer rFFT indices in [0, bins//2]")
        offsets, _, _, angles, _ = self._wedges(cutoff)
        out = self._run(kernels.angle_histogram_chunk, self.n, offsets, angles, bins, freq)
        if power:
            out = out**2
        return self.local_average(out, cutoff=cutoff) if averaged else out

    def neighborhood_parameters(self, m=1, *, averaged=False, cutoff=None):
        """A1/A2/A3 (A/P/N) in A^2, following the author's MALIO convention.

        For each bond i--j choose up to m other neighbors k of i ranked by
        r_ik+r_jk. A1 averages |sum_k(r_ik+r_jk)|^2 over j. A2 averages
        |sum_k(r_ij+r_kj)|^2. A3 instead squares the combined sum over j before
        dividing by N. This resolves the manuscript's duplicated A2/A3 formula.
        The neighborhood subset follows MALIO's radial ranking, rather than
        imposing an additional common-cutoff-neighbor intersection.
        """
        if isinstance(m, bool) or not isinstance(m, Integral) or m < 1:
            raise ValueError("m must be a positive integer")
        graph = self.neighbor_list(cutoff)
        out = self._run(kernels.neighborhood_chunk, self.n, graph.offsets, graph.vectors,
                        graph.distances, *_cell_geometry(self.dimensions), m)
        return self.local_average(out, cutoff=cutoff) if averaged else out

    def quadrupole_parameter(self, dipole, *, negative_site_distance=0.):
        """Three-charge water Q_T in supplied dipole units times Angstrom.

        Q_T = 3*(r_OH*sin(theta))^2*mu / (4*(r_OH*cos(theta)-d_OM)).
        theta is half the H-O-H angle; mean O-H length approximates flexible
        geometries. Supply dipoles and the negative-site offset for the chosen
        water model. This molecular moment is not a lattice-order classifier.
        """
        dipole = np.broadcast_to(np.asarray(dipole, dtype=float), (self.n,))
        offset = np.broadcast_to(np.asarray(negative_site_distance, dtype=float), (self.n,))
        if not np.isfinite(dipole).all() or not np.isfinite(offset).all() or np.any(offset < 0):
            raise ValueError("Dipoles and nonnegative site offsets must be finite")
        oh = self._water_geometry()
        lengths = np.linalg.norm(oh, axis=2)
        cosine = np.sum(oh[:, 0]*oh[:, 1], axis=1)/(lengths[:, 0]*lengths[:, 1])
        theta = .5*np.arccos(np.clip(cosine, -1, 1))
        r = lengths.mean(axis=1)
        denominator = r*np.cos(theta)-offset
        if np.any(denominator <= 1e-12):
            raise ValueError("Water geometry/site offset gives a nonpositive dipole lever arm")
        return .75*(r*np.sin(theta))**2*dipole/denominator

    def voronoi_geometry(self):
        """Periodic Voro++ cells via freud, using all selected positions.

        Returns volume (A^3), area (A^2), number density (A^-3), and A^3/(36pi V^2).
        The total cell volume must equal the periodic box volume.
        """
        if self._voronoi is None:
            import freud
            cell = triclinic_vectors(self.dimensions, dtype=np.float64)
            fractional = np.mod(self.pos @ np.linalg.inv(cell), 1.)
            points = (fractional-.5) @ cell
            if len(np.unique(points.astype(np.float32), axis=0)) != self.n:
                raise ValueError("Periodic Voronoi cells require distinct selected positions")
            box = freud.box.Box.from_matrix(cell.T)
            tess = freud.locality.Voronoi().compute((box, points))
            volume = np.array(tess.volumes, dtype=float)
            area = np.bincount(tess.nlist.query_point_indices, weights=tess.nlist.weights, minlength=self.n)
            if not np.isfinite(volume).all() or np.any(volume <= 0) or not np.isclose(volume.sum(), self.volume, rtol=1e-5):
                raise RuntimeError("Periodic Voronoi tessellation failed volume conservation")
            arrays = (volume, area, 1./volume, area**3/(36*np.pi*volume**2))
            for array in arrays:
                array.flags.writeable = False
            self._voronoi = VoronoiGeometry(*arrays)
        return self._voronoi

    def voronoi_density(self):
        return self.voronoi_geometry().density

    def voronoi_asphericity(self):
        return self.voronoi_geometry().asphericity

    def chill(self, *, plus=False, cutoff=3.5):
        """CHILL / CHILL+ ice classification with explicit neighbor conventions.

        CHILL uses the four nearest oxygens inside cutoff and eclipsed band
        (-.2,-.05). CHILL+ follows OVITO's cutoff implementation: all oxygens
        inside cutoff contribute to q_3m, but only four-coordinated centers
        receive a structure label. Its eclipsed band is (-.35,.25). Both use
        staggered correlations < -.8. Undercoordinated neighbors still have
        harmonic vectors and may contribute to another center's correlations.
        No hydrogen coordinates are needed. Labels are not confidence scores.
        """
        graph = self.neighbor_list(cutoff)
        counts = graph.counts if plus else np.minimum(graph.counts, 4)
        offsets = np.r_[0, np.cumsum(counts)].astype(np.int64)
        edges = np.concatenate([np.arange(graph.offsets[i], graph.offsets[i]+n) for i, n in enumerate(counts)])
        indices, vectors, distances = graph.indices[edges], graph.vectors[edges], graph.distances[edges]
        if np.any(distances == 0):
            raise ValueError("CHILL is undefined for coincident oxygens")
        theta = np.arccos(np.clip(vectors[:, 2]/distances, -1, 1))
        phi = np.mod(np.arctan2(vectors[:, 1], vectors[:, 0]), 2*np.pi)
        harmonics = sph_harm_y(3, np.arange(-3, 4)[None, :], theta[:, None], phi[:, None])
        qlm = base_kernels.harmonic_mean(offsets, harmonics)
        # One output row per directed edge; keep map_chunks' row contract.
        source = np.repeat(np.arange(self.n), counts)
        correlations = self._run(kernels.edge_correlation_chunk, len(indices), source, indices, qlm)
        staggered = np.bincount(source[correlations < -.8], minlength=self.n)
        if plus:
            mask = (correlations > -.35) & (correlations < .25)
        else:
            mask = (correlations > -.2) & (correlations < -.05)
        eclipsed = np.bincount(source[mask], minlength=self.n)
        labels = self._run(kernels.classify_chunk, self.n, offsets, indices, staggered, eclipsed, plus)
        return IceClassification(labels, staggered, eclipsed, correlations, offsets, indices)

    def chill_plus(self, cutoff=3.5):
        return self.chill(plus=True, cutoff=cutoff)

    def _hbond_edges(self, graph, cutoff, angle):
        """Mask of directed edges in ``graph`` that are hydrogen bonds either way."""
        bonds = self.hydrogen_bonds(cutoff, angle).bonds
        keys = np.r_[bonds[:, 0]*self.n+bonds[:, 2], bonds[:, 2]*self.n+bonds[:, 0]]
        source = np.repeat(np.arange(self.n), graph.counts)
        return np.isin(source*self.n+graph.indices, keys)

    def zeta(self, cutoff=3.5, angle=30., *, search=None):
        """Russo-Tanaka zeta: first non-H-bonded minus last H-bonded O--O distance.

        Hydrogen bonds use the Luzar-Chandler geometry of ``hydrogen_bonds``,
        counted in either direction. Neighbors are searched within ``search``
        (default: the snapshot cutoff), which must exceed ``cutoff``. Centers
        without a hydrogen bond, or without a non-bonded oxygen inside
        ``search``, give NaN. Units: Angstrom.
        """
        radius = self.cutoff if search is None else float(search)
        if not np.isfinite(radius) or radius <= cutoff:
            raise ValueError("zeta search radius must exceed the hydrogen-bond cutoff")
        graph = self.neighbor_list(radius)
        bonded = self._hbond_edges(graph, cutoff, angle)
        source = np.repeat(np.arange(self.n), graph.counts)
        last_bonded = np.full(self.n, -np.inf)
        first_free = np.full(self.n, np.inf)
        np.maximum.at(last_bonded, source[bonded], graph.distances[bonded])
        np.minimum.at(first_free, source[~bonded], graph.distances[~bonded])
        out = first_free - last_bonded
        out[~np.isfinite(out)] = np.nan
        return out

    def ring_statistics(self, max_size=10, *, network="hbond", cutoff=3.5, angle=30.):
        """Smallest-ring (King/Guttman) statistics of the oxygen network.

        ``network='hbond'`` links waters joined by a hydrogen bond in either
        direction; ``network='distance'`` links all oxygens within ``cutoff``
        and needs no hydrogens. Paths that wrap around the periodic box are
        not rings. Search cost grows quickly with ``max_size``.
        """
        if isinstance(max_size, bool) or not isinstance(max_size, Integral) or not 3 <= max_size <= 16:
            raise ValueError("max_size must be an integer in [3, 16]")
        graph = self.neighbor_list(cutoff)
        if network == "hbond":
            mask = self._hbond_edges(graph, cutoff, angle)
        elif network == "distance":
            mask = np.ones(len(graph.indices), dtype=bool)
        else:
            raise ValueError("network must be 'hbond' or 'distance'")
        source = np.repeat(np.arange(self.n), graph.counts)
        offsets = np.r_[0, np.cumsum(np.bincount(source[mask], minlength=self.n))].astype(np.int64)
        counts = self._run(kernels.smallest_ring_chunk, self.n, offsets, graph.indices[mask],
                           np.ascontiguousarray(graph.vectors[mask]), max_size)
        return RingStatistics(np.arange(3, max_size+1), counts[:, 3:], counts[:, 0])
