"""Optional-backend integration checks, skipped when extras are not installed."""

import numpy as np
import pytest
from MDAnalysis.lib.mdamath import triclinic_vectors

from SolOrder import StructuralOrder
from test_structure import oxygen_universe
from test_extended import ice_analysis, PHASES


def test_local_dscribe_descriptors_and_parallel_centers():
    pytest.importorskip("dscribe")
    op = StructuralOrder(oxygen_universe([[5, 5, 5], [7, 5, 5], [5, 7, 5]]), n_jobs=2)
    serial = StructuralOrder(op.u)
    for kind, settings in (
        ("soap", dict(r_cut=4., n_max=2, l_max=2, sigma=.4)),
        ("lmbtr", dict(geometry={"function": "distance"}, grid={"min": 0, "max": 5, "sigma": .1, "n": 30}, weighting={"function": "exp", "scale": 1., "threshold": 1e-3})),
        ("acsf", dict(r_cut=4., g2_params=[[1., 0.], [.5, 1.]])),
    ):
        a = op.atomic_descriptor(kind, symbols="O", centers=[2, 0], **settings)
        b = serial.atomic_descriptor(kind, symbols="O", centers=[2, 0], **settings)
        assert a.shape[0] == 2 and a.shape[1] > 0 and np.isfinite(a).all()
        np.testing.assert_allclose(a, b, atol=1e-12)
        all_centers = serial.atomic_descriptor(kind, symbols="O", **settings)
        np.testing.assert_allclose(a, all_centers[[2, 0]], atol=1e-12)
        empty = op.atomic_descriptor(kind, symbols="O", centers=[], **settings)
        assert empty.shape == (0, a.shape[1])


def test_global_descriptors_and_coulomb_reference():
    pytest.importorskip("dscribe")
    op = StructuralOrder(oxygen_universe([[5, 5, 5], [7, 5, 5]]))
    coulomb = op.atomic_descriptor("coulomb", symbols="O", permutation="none").reshape(2, 2)
    np.testing.assert_allclose(coulomb, [[.5*8**2.4, 32], [32, .5*8**2.4]])
    sine = op.atomic_descriptor("sine", symbols="O", permutation="none").reshape(2, 2)
    assert sine[0, 1] == pytest.approx(64/(20*np.sin(np.pi*2/20)**2))
    ewald = op.atomic_descriptor("ewald", symbols="O", permutation="none")
    assert ewald.shape == (4,) and np.isfinite(ewald).all()
    mbtr = op.atomic_descriptor("mbtr", symbols="O", geometry={"function": "distance"},
                               grid={"min": 0, "max": 5, "sigma": .1, "n": 30},
                               weighting={"function": "exp", "scale": 1., "threshold": 1e-3})
    assert mbtr.shape == (30,) and np.isfinite(mbtr).all()


def test_soap_rotation_and_process_backend():
    pytest.importorskip("dscribe")
    positions = np.array([[5, 5, 5], [7, 5, 5], [5, 7, 5]], dtype=float)
    a = StructuralOrder(oxygen_universe(positions))
    b = StructuralOrder(oxygen_universe(positions[:, [2, 0, 1]]), n_jobs=2, backend="loky")
    np.testing.assert_allclose(a.soap(r_cut=4, n_max=2, l_max=2), b.soap(r_cut=4, n_max=2, l_max=2), atol=1e-10)


def test_graph_descriptor_mapping_and_features():
    op = StructuralOrder(oxygen_universe([[5, 5, 5], [7, 5, 5], [5, 7, 5]]), cutoff=2.1)
    graph = op.graph_descriptor()
    np.testing.assert_array_equal(graph.node_features[:, 0], [2, 1, 1])
    np.testing.assert_allclose(graph.node_features[:, 1], 2)
    np.testing.assert_allclose(graph.node_features[:, 2], 0)
    assert graph.edge_index.shape == (2, 4)
    np.testing.assert_array_equal(graph.atom_indices, op.atom_indices)


def _ovito_pipeline(op, modifier):
    pytest.importorskip("ovito")
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    data = DataCollection()
    particles = Particles(count=op.n)
    particles.create_property("Position", data=op.pos.copy())
    data.objects.append(particles)
    cell = SimulationCell(pbc=(True, True, True))
    cell[...] = np.column_stack((triclinic_vectors(op.box, dtype=float).T, np.zeros(3)))
    data.objects.append(cell)
    pipeline = Pipeline(source=StaticSource(data=data))
    pipeline.modifiers.append(modifier)
    return pipeline


@pytest.mark.parametrize("phase", PHASES)
def test_chill_plus_matches_ovito(phase):
    pytest.importorskip("ovito")
    from ovito.modifiers import ChillPlusModifier
    op = ice_analysis(phase, n_jobs=2)
    expected = np.array(_ovito_pipeline(op, ChillPlusModifier(cutoff=3.5)).compute().particles["Structure Type"])
    np.testing.assert_array_equal(op.chill_plus().labels, expected)


@pytest.mark.parametrize("phase,expected", [("Ih", 7), ("Ic", 6)])
def test_ice_polyhedral_templates(phase, expected):
    pytest.importorskip("ovito")
    result = ice_analysis(phase).template_matching()
    np.testing.assert_array_equal(result.labels, expected)
    assert np.max(result.rmsd) < .1


@pytest.mark.parametrize("noise", [.05, .2, .5])
def test_chill_plus_matches_ovito_with_disorder_and_surface(noise):
    pytest.importorskip("ovito")
    from ovito.modifiers import ChillPlusModifier
    ideal = ice_analysis("Ih")
    positions = ideal.pos.copy()+np.random.default_rng(25).normal(0, noise, ideal.pos.shape)
    # Remove a slab to create undercoordinated surface environments.
    positions = positions[positions[:, 0] > 5.]
    op = StructuralOrder(oxygen_universe(positions, ideal.dimensions), n_jobs=2)
    expected = np.array(_ovito_pipeline(op, ChillPlusModifier(cutoff=3.5)).compute().particles["Structure Type"])
    np.testing.assert_array_equal(op.chill_plus().labels, expected)


@pytest.mark.parametrize("phase", ["Ih", "V"])
def test_soap_on_genice_periodic_cells(phase):
    pytest.importorskip("dscribe")
    op = ice_analysis(phase, n_jobs=2)
    reference = op.soap(r_cut=3.5, n_max=2, l_max=2, centers=[0, 1, 2])
    cell = triclinic_vectors(op.box, dtype=float)
    shifted = StructuralOrder(oxygen_universe(op.pos+cell[0]-2*cell[1], op.box), n_jobs=2)
    np.testing.assert_allclose(reference, shifted.soap(r_cut=3.5, n_max=2, l_max=2, centers=[0, 1, 2]), rtol=1e-4, atol=1e-6)
