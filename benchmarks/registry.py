"""Order-parameter registry for the review benchmark.

Settings are fixed here before any separability is computed (design principle
1 in docs/benchmark_design.md). Each entry maps one SolOrder call to one or
more per-oxygen columns. ``op`` is an OrderParameter snapshot built with a
5 A neighbor cutoff; ``ent`` is the matching eFingerprint snapshot.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np

FIRST_SHELL = 3.5  # A; O-O first-shell cutoff used by every first-shell parameter


@dataclass(frozen=True)
class Entry:
    name: str
    group: str
    compute: Callable  # (op, ent) -> array, or dict of column -> array
    needs_h: bool = False
    smooth: bool = False  # smooth functional form (switching function, Gaussian); not a validated CV
    kind: str = "scalar"  # "scalar" or "label"
    source: str = ""

    def columns(self, op, ent):
        out = self.compute(op, ent)
        return out if isinstance(out, dict) else {self.name: out}


def _q_tet_averaged(op, ent):
    return op.local_average(op.tetrahedral_order(), cutoff=FIRST_SHELL)


def _angle_histogram(op, ent):
    values = op.angle_histogram_parameter((1, 2, 3), cutoff=FIRST_SHELL)
    return {f"ah{k}": values[:, i] for i, k in enumerate((1, 2, 3))}


def _neighborhood(op, ent):
    values = op.neighborhood_parameters(m=1, cutoff=FIRST_SHELL)
    return {f"A{k}": values[:, k-1] for k in (1, 2, 3)}


def _steinhardt(op, ent):
    out = {}
    for l in (3, 4, 6):
        out[f"q{l}"] = op.steinhardt(l, cutoff=FIRST_SHELL)
        out[f"q{l}_avg"] = op.steinhardt(l, averaged=True, cutoff=FIRST_SHELL)
    for l in (4, 6):
        out[f"w{l}_hat"] = op.steinhardt(l, invariant="w", normalized=True, cutoff=FIRST_SHELL)
    return out


def _voronoi(op, ent):
    return {"vor_density": op.voronoi_density(), "vor_asphericity": op.voronoi_asphericity()}


def _entropy(op, ent):
    return {"entropy": ent.Entropy(), "entropy_avg": ent.localEntropy()}


RING_SIZES = (4, 5, 6, 7, 8)


def _rings(op, ent):
    rings = op.ring_statistics(10, network="hbond", cutoff=FIRST_SHELL)
    fractions = rings.fractions
    out = {"ring_mean": rings.mean_size}
    for size in RING_SIZES:
        out[f"ring_f{size}"] = fractions[:, size-3]
    return out


def _ptm(op, ent):
    match = op.template_matching()
    code = {v: k for k, v in match.type_names.items()}
    labels = np.zeros(op.n, dtype=np.int8)  # common codes: 0 other, 1 Ih, 2 Ic
    labels[match.labels == code["Hexagonal diamond"]] = 1
    labels[match.labels == code["Cubic diamond"]] = 2
    return labels


SCALARS = [
    Entry("cn", "radial", lambda op, ent: op.coordination_number(FIRST_SHELL),
          source="cutoff count"),
    Entry("cn_smooth", "radial", lambda op, ent: op.smooth_coordination(r0=FIRST_SHELL, n=12, m=28),
          smooth=True, source="rational switching function"),
    Entry("d5", "radial", lambda op, ent: op.d5(), source="Cuthbertson & Poole 2011"),
    Entry("lsi", "radial", lambda op, ent: op.local_structure_index(3.7),
          source="Shiratani & Sasai 1996"),
    Entry("voronoi", "radial", _voronoi, source="periodic Voronoi (freud)"),
    Entry("entropy", "radial", _entropy, smooth=True, source="Piaggi & Parrinello 2017"),
    Entry("q_tet", "angular", lambda op, ent: op.tetrahedral_order(), source="Errington & Debenedetti 2001"),
    Entry("q_tet_avg", "angular", _q_tet_averaged),
    Entry("s_k", "angular", lambda op, ent: op.translational_order(), source="Chau & Hardwick 1998"),
    Entry("bond_angle", "angular", lambda op, ent: op.bond_angle_parameter(1, 1, cutoff=FIRST_SHELL),
          source="Takahashi 2023"),
    Entry("afs3", "angular", lambda op, ent: op.angular_fourier(a=3, radial=0, cutoff=FIRST_SHELL)),
    Entry("angle_hist", "angular", _angle_histogram),
    Entry("neighborhood", "angular", _neighborhood, source="MALIO A/P/N"),
    Entry("f4", "hbond", lambda op, ent: op.F4(cutoff=FIRST_SHELL), needs_h=True,
          source="Rodger et al. 1996"),
    Entry("n_hb", "hbond", lambda op, ent: op.hydrogen_bonds(FIRST_SHELL, 30.).total.astype(float),
          needs_h=True, source="Luzar & Chandler 1996"),
    Entry("zeta", "hbond", lambda op, ent: op.zeta(FIRST_SHELL, 30.), needs_h=True,
          source="Russo & Tanaka 2014"),
    Entry("rings", "hbond", _rings, needs_h=True, source="King 1967; Guttman 1990"),
    Entry("steinhardt", "bond_orientational", _steinhardt, source="Steinhardt 1983; Lechner & Dellago 2008"),
    Entry("xconn6", "bond_orientational",
          lambda op, ent: op.crystalline_connections(6, .5, cutoff=FIRST_SHELL).astype(float),
          source="ten Wolde et al. 1996"),
    Entry("csp", "template", lambda op, ent: op.centrosymmetry(12), source="Kelchner et al. 1998"),
    Entry("ndp", "template", lambda op, ent: op.neighbor_distance_parameter(cutoff=FIRST_SHELL)),
]

LABELS = [
    Entry("chill", "label", lambda op, ent: op.chill(cutoff=FIRST_SHELL).labels, kind="label",
          source="Moore et al. 2010"),
    Entry("chill_plus", "label", lambda op, ent: op.chill_plus(FIRST_SHELL).labels, kind="label",
          source="Nguyen & Molinero 2015"),
    Entry("ptm", "label", _ptm, kind="label", source="Larsen et al. 2016 (OVITO)"),
]

REGISTRY = SCALARS + LABELS
