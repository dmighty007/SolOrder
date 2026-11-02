"""Regression tests on raw, seeded GenIce2 ice and empty-clathrate structures."""

import hashlib
from itertools import product
import json
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals, distance_array, minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors
import numpy as np
import pytest
from scipy.special import eval_legendre

from SolOrder import OrderParameter, eFingerprint


DATA = Path(__file__).parent/"data/genice"
MANIFEST = json.loads((DATA/"manifest.json").read_text())
PHASES = list(MANIFEST["structures"])


@pytest.fixture(params=PHASES)
def ice(request):
    phase = request.param
    u = mda.Universe(str(DATA/f"{phase}.gro"))
    return phase, u


def analyze(u, **kwargs):
    return OrderParameter(u, cutoff=5, selection="name O", hydrogen_selection="name H", **kwargs)


def reference_geometry(u):
    oxygen = u.select_atoms("name O").positions.astype(float)
    distances = distance_array(oxygen, oxygen, box=u.dimensions)
    np.fill_diagonal(distances, np.inf)
    return oxygen, distances


def reference_f4(u, cutoff=3.5):
    """Choose each outer hydrogen independently; use MDAnalysis torsions."""
    oxygen, distances = reference_geometry(u)
    donor, acceptor = np.where(distances <= cutoff)
    oh = np.array([minimize_vectors(res.atoms.select_atoms("name H").positions-o,
                                   u.dimensions)
                   for res, o in zip(u.residues, oxygen)])
    oo = minimize_vectors(oxygen[acceptor]-oxygen[donor], u.dimensions)
    ha, hb = oh[donor], oh[acceptor]
    # The outer H of each water is the one farthest from the opposite oxygen.
    ia = np.argmax(np.linalg.norm(ha-oo[:, None, :], axis=2), axis=1)
    ib = np.argmax(np.linalg.norm(hb+oo[:, None, :], axis=2), axis=1)
    a = ha[np.arange(len(donor)), ia]
    b = oo + hb[np.arange(len(donor)), ib]
    angles = calc_dihedrals(a, np.zeros_like(oo), oo, b)
    sums = np.bincount(donor, weights=np.cos(3*angles), minlength=len(oxygen))
    return sums/np.bincount(donor, minlength=len(oxygen))


def test_genice_provenance_and_ice_rules(ice):
    phase, u = ice
    entry = MANIFEST["structures"][phase]
    assert hashlib.sha256((DATA/entry["file"]).read_bytes()).hexdigest() == entry["sha256"]
    assert len(u.residues) == entry["molecules"]
    assert set(u.atoms.names) == {"O", "H"}
    hb = analyze(u, n_jobs=2).hydrogen_bonds()
    np.testing.assert_array_equal(hb.donated, 2)
    np.testing.assert_array_equal(hb.accepted, 2)
    assert len(hb.bonds) == 2*entry["molecules"]
    assert len(set(map(tuple, hb.bonds))) == len(hb.bonds)


def test_genice_distances_and_tetrahedral_reference(ice):
    phase, u = ice
    op = analyze(u, n_jobs=2)
    oxygen, distances = reference_geometry(u)
    np.testing.assert_array_equal(op.coordination_number(3.5), np.sum(distances <= 3.5, axis=1))
    np.testing.assert_allclose(op.d5(), np.sort(distances, axis=1)[:, 4], atol=3e-6)
    cell = triclinic_vectors(u.dimensions, dtype=float)
    assert op.density() == pytest.approx(op.n/abs(np.linalg.det(cell)), rel=1e-7)

    # Independently enumerate periodic images for distributed central atoms.
    # Wrapping fractional coordinates permits multiple-cell input translations.
    wrapped = np.mod(oxygen @ np.linalg.inv(cell), 1) @ cell
    rows = np.linspace(0, op.n-1, 12, dtype=int)
    delta = wrapped[None, :, :]-wrapped[rows, None, :]
    brute = np.full(delta.shape[:2], np.inf)
    for shift in product(range(-2, 3), repeat=3):
        brute = np.minimum(brute, np.linalg.norm(delta + np.array(shift) @ cell, axis=2))
    brute[np.arange(len(rows)), rows] = np.inf
    np.testing.assert_allclose(distances[rows], brute, atol=6e-6)

    # Check four-neighbor formulas directly, outside the Numba kernels.
    nearest = np.argsort(distances, axis=1)[:, :4]
    radii = np.take_along_axis(distances, nearest, axis=1)
    mean = radii.mean(axis=1)
    expected_tto = 1-((radii-mean[:, None])**2).sum(axis=1)/(12*mean**2)
    np.testing.assert_allclose(op.TTO(), expected_tto, atol=2e-7)
    expected_q = []
    for i in range(op.n):
        vectors = minimize_vectors(oxygen[nearest[i]]-oxygen[i], u.dimensions)
        unit = vectors/np.linalg.norm(vectors, axis=1)[:, None]
        cosine = unit @ unit.T
        expected_q.append(1-3/8*((cosine[np.triu_indices(4, 1)]+1/3)**2).sum())
    # Degenerate kth shells can have multiple equally valid four-neighbor sets;
    # test exact identities only when the fourth neighbor is unambiguous.
    sorted_r = np.sort(distances, axis=1)
    unambiguous = sorted_r[:, 4]-sorted_r[:, 3] > 2e-5
    assert np.any(unambiguous) or phase in ("VII", "VIII")
    np.testing.assert_allclose(op.OTO()[unambiguous], np.array(expected_q)[unambiguous], atol=3e-6)
    if phase in ("Ih", "Ic"):
        assert np.min(op.OTO()) > .99
        np.testing.assert_array_equal(op.coordination_number(3.5), 4)


def test_genice_f4_independent_outer_hydrogens(ice):
    _, u = ice
    np.testing.assert_allclose(analyze(u, n_jobs=2).F4(cutoff=3.5), reference_f4(u), atol=2e-6)


def test_genice_bond_order_and_lsi_reference(ice):
    _, u = ice
    op = analyze(u, n_jobs=2)
    oxygen, distances = reference_geometry(u)
    rows = np.linspace(0, op.n-1, 12, dtype=int)
    for l in (3, 4, 6):
        expected = []
        for i in rows:
            neighbors = np.flatnonzero(distances[i] <= 3.5)
            vectors = minimize_vectors(oxygen[neighbors]-oxygen[i], u.dimensions)
            unit = vectors/np.linalg.norm(vectors, axis=1)[:, None]
            expected.append(np.sqrt(max(0, eval_legendre(l, np.clip(unit @ unit.T, -1, 1)).mean())))
        np.testing.assert_allclose(op.steinhardt(l, cutoff=3.5)[rows], expected, atol=3e-7)
    expected = []
    for i in rows:
        radii = np.sort(distances[i][distances[i] <= 5])
        n = np.sum(radii < 3.7)
        expected.append(np.var(np.diff(radii[:n+1])))
    np.testing.assert_allclose(op.LSI()[rows], expected, atol=2e-6)


def test_genice_serial_parallel_and_periodic_translation(ice):
    _, u = ice
    a, b = analyze(u), analyze(u, n_jobs=2)
    for method in ("coordination_number", "smooth_coordination", "d5", "OTO", "TTO", "LSI", "MinimumAngle"):
        serial, parallel = getattr(a, method)(), getattr(b, method)()
        assert np.isfinite(serial).all()
        np.testing.assert_allclose(serial, parallel, atol=1e-14)
    np.testing.assert_allclose(a.F4(cutoff=3.5), b.F4(cutoff=3.5), atol=1e-14)
    for kwargs in ({"l": 3}, {"l": 4, "averaged": True},
                   {"l": 6, "invariant": "w", "normalized": True}):
        values = a.steinhardt(**kwargs)
        assert np.isfinite(values).all()
        np.testing.assert_allclose(values, b.steinhardt(**kwargs), atol=1e-14)
    for local in (False, True):
        e1 = eFingerprint(u, selection="name O", local=local)
        e2 = eFingerprint(u, selection="name O", local=local, n_jobs=2)
        assert np.isfinite(e1.Entropy()).all()
        np.testing.assert_allclose(e1.e, e2.Entropy(), atol=1e-14)
        np.testing.assert_allclose(e1.localEntropy(), e2.localEntropy(), atol=1e-14)
    shifted = u.copy()
    cell = triclinic_vectors(u.dimensions, dtype=float)
    shift = 2*cell[0]-cell[1]+cell[2]
    # Translate whole waters across several cells, keeping the atom-name input.
    shifted.atoms.positions += shift
    other = analyze(shifted, n_jobs=2)
    np.testing.assert_array_equal(a.coordination_number(3.5), other.coordination_number(3.5))
    np.testing.assert_allclose(a.d5(), other.d5(), atol=2e-5)
    np.testing.assert_allclose(a.steinhardt(6), other.steinhardt(6), atol=2e-5)
    np.testing.assert_array_equal(a.hydrogen_bonds().total, other.hydrogen_bonds().total)


@pytest.mark.parametrize("phase", ["Ih", "V"])
def test_genice_process_backend(phase):
    u = mda.Universe(str(DATA/f"{phase}.gro"))
    serial = analyze(u)
    process = analyze(u, n_jobs=2, backend="loky")
    for calculation in (lambda x: x.steinhardt(6), lambda x: x.F4(cutoff=3.5),
                        lambda x: x.hydrogen_bonds().total):
        np.testing.assert_allclose(calculation(serial), calculation(process), atol=1e-14)
    a = eFingerprint(u, selection="name O")
    b = eFingerprint(u, selection="name O", n_jobs=2, backend="loky")
    np.testing.assert_allclose(a.Entropy(), b.Entropy(), atol=1e-14)
