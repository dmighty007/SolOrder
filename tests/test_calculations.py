from importlib.metadata import version

import numpy as np
import pytest
from joblib import parallel_config
from MDAnalysis.lib.distances import calc_dihedrals

from SolOrder import OrderParameter, __version__, eFingerprint
from SolOrder.entropy import dist


def test_public_api_and_version():
    assert __version__ == version("SolOrder") == "1.1.0"
    assert callable(OrderParameter)
    assert callable(eFingerprint)


def test_tetrahedral_order(water_universe):
    op = OrderParameter(water_universe, cutoff=5)
    with parallel_config(backend="threading"):
        op.OTO()
    op.TTO()
    assert op.tetra_orient.shape == (5,)
    assert op.tetra_orient[0] == pytest.approx(1, abs=1e-6)
    assert op.tetra_trans[0] == pytest.approx(1, abs=1e-6)


def test_dihedrals_and_parallel_f4(water_universe):
    op = OrderParameter(water_universe, cutoff=5)
    op.get_dihedrals()
    assert len(op.dihedrals) == 10  # Unique pairs, rather than directed neighbors.
    assert np.isfinite(op.dihedrals).all()
    expected = []
    for i in range(op.n):
        values = []
        for j in op.neighbours(i, op.arr):
            atoms = op.get_atoms(op.maping(i)[0], op.maping(j)[0])
            points = op._sort_dihedral(*atoms)
            phi = calc_dihedrals(*points)
            values.append(np.cos(3 * phi))
        expected.append(np.mean(values))
    with parallel_config(backend="threading"):
        op.F4(cutoff=5)
    np.testing.assert_allclose(op.f4, expected, atol=1e-6)


def test_entropy_and_neighbor_average(water_universe):
    fingerprint = eFingerprint(water_universe, cutoff=5)
    fingerprint.Entropy()
    assert fingerprint.e.shape == (5,)
    assert np.isfinite(fingerprint.e).all()
    # All five oxygens lie inside each other's cutoff in this fixture.
    expected_average = fingerprint.e.mean()
    fingerprint.localEntropy()
    np.testing.assert_allclose(fingerprint.le, expected_average)


def test_periodic_distance():
    a = np.array([0.5, 1., 1.])
    b = np.array([9.5, 1., 1.])
    assert dist(a, b, np.array([10., 10., 10.])) == pytest.approx(1.)
