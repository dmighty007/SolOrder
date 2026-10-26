import numpy as np
import MDAnalysis as mda
import pytest
from MDAnalysis.lib.distances import distance_array, minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors
from scipy.special import eval_legendre

from SolOrder import StructuralOrder, OrderParameter, eFingerprint
from SolOrder.parallel import map_chunks
from SolOrder import _kernels


def oxygen_universe(positions, dimensions=(20, 20, 20, 90, 90, 90)):
    u = mda.Universe.empty(len(positions), n_residues=len(positions),
                           atom_resindex=np.arange(len(positions)), trajectory=True)
    u.add_TopologyAttr("names", ["OW"]*len(positions))
    u.add_TopologyAttr("resids", np.arange(1, len(positions)+1))
    u.atoms.positions = positions
    u.dimensions = dimensions
    return u


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_coordination_and_kth_reference(n_jobs):
    rng = np.random.default_rng(13)
    u = oxygen_universe(rng.uniform(0, 20, (60, 3)))
    op = StructuralOrder(u, cutoff=8, n_jobs=n_jobs)
    distances = distance_array(u.atoms.positions, u.atoms.positions, box=u.dimensions)
    np.fill_diagonal(distances, np.inf)
    expected_cn = np.sum(distances <= 8, axis=1)
    np.testing.assert_array_equal(op.coordination_number(), expected_cn)
    expected_fifth = np.sort(distances, axis=1)[:, 4]
    expected_fifth[expected_fifth > 8] = np.nan
    np.testing.assert_allclose(op.d5(), expected_fifth, atol=2e-6)
    x = distances/3.5
    with np.errstate(invalid="ignore"):
        weights = (1-x**12)/(1-x**28)
    weights[distances > 8] = 0
    np.testing.assert_allclose(op.smooth_coordination(), weights.sum(axis=1), atol=2e-6)
    assert op.density() == pytest.approx(60/8000)


def test_switching_limit_and_no_self():
    u = oxygen_universe([[0, 0, 0], [3.5, 0, 0]])
    op = StructuralOrder(u, cutoff=4, n_jobs=2)
    np.testing.assert_allclose(op.smooth_coordination(), 12/28)
    np.testing.assert_array_equal(op.coordination_number(), [1, 1])
    assert np.isnan(op.d5()).all()


def test_triclinic_multiple_box_translations():
    box = np.array([15., 16., 17., 75., 80., 65.])
    cell = triclinic_vectors(box)
    fractional = np.random.default_rng(4).uniform(.1, .9, (30, 3))
    positions = fractional @ cell
    u = oxygen_universe(positions, box)
    op = StructuralOrder(u, cutoff=6, n_jobs=2)
    shifted = positions.copy()
    shifted[::2] += 2*cell[0] - 3*cell[1] + cell[2]
    other = StructuralOrder(oxygen_universe(shifted, box), cutoff=6)
    np.testing.assert_array_equal(op.coordination_number(), other.coordination_number())
    np.testing.assert_allclose(op.kth_neighbor_distance(2), other.kth_neighbor_distance(2), atol=1e-5)
    assert op.density() == pytest.approx(30/abs(np.linalg.det(cell)), rel=1e-6)
    expected = distance_array(u.atoms.positions, u.atoms.positions, box=box)
    for i in range(op.n):
        graph = op.neighbor_list()
        sl = slice(graph.offsets[i], graph.offsets[i+1])
        np.testing.assert_allclose(graph.distances[sl], expected[i, graph.indices[sl]], atol=2e-6)


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_fcc_steinhardt_known_values(n_jobs):
    neighbors = []
    for zero in range(3):
        for a in (-1, 1):
            for b in (-1, 1):
                point = [a, b]
                point.insert(zero, 0)
                neighbors.append(point)
    positions = np.vstack([np.zeros(3), np.array(neighbors)*2.]) + 10
    op = StructuralOrder(oxygen_universe(positions), cutoff=3.1, n_jobs=n_jobs)
    assert op.steinhardt(4)[0] == pytest.approx(np.sqrt(7/192), abs=1e-7)
    assert op.steinhardt(6)[0] == pytest.approx(np.sqrt(169/512), abs=1e-7)
    assert op.steinhardt(4, invariant="w", normalized=True)[0] == pytest.approx(-.1593173731, abs=1e-8)
    assert op.steinhardt(6, invariant="w", normalized=True)[0] == pytest.approx(-.01316060073, abs=1e-8)


def test_bond_order_independent_legendre_reference_and_rotation():
    rng = np.random.default_rng(7)
    positions = rng.normal(0, 1, (18, 3)) + 10
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    op = StructuralOrder(oxygen_universe(positions), cutoff=5, n_jobs=2)
    rotated = StructuralOrder(oxygen_universe((positions-10) @ rotation+10), cutoff=5)
    g = op.neighbor_list()
    for l in (3, 4, 6):
        expected = []
        for i in range(op.n):
            sl = slice(g.offsets[i], g.offsets[i+1])
            unit = g.vectors[sl]/g.distances[sl, None]
            # Addition theorem avoids using spherical harmonics as the oracle.
            q2 = np.mean(eval_legendre(l, np.clip(unit @ unit.T, -1, 1)))
            expected.append(np.sqrt(max(q2, 0)))
        np.testing.assert_allclose(op.steinhardt(l), expected, atol=1e-7)
        for averaged in (False, True):
            np.testing.assert_allclose(op.steinhardt(l, averaged=averaged),
                                       rotated.steinhardt(l, averaged=averaged), atol=2e-6)
            np.testing.assert_allclose(op.steinhardt(l, invariant="w", normalized=True, averaged=averaged),
                                       rotated.steinhardt(l, invariant="w", normalized=True, averaged=averaged), atol=2e-6)
    np.testing.assert_array_equal(op.crystalline_connections(0, threshold=.5), g.counts)


def test_average_harmonics_before_invariant(water_universe):
    op = StructuralOrder(water_universe, cutoff=5, n_jobs=2)
    qlm = op._qlm(3, None)
    averaged = np.broadcast_to(qlm.mean(axis=0), qlm.shape)
    expected = np.sqrt(4*np.pi/7 * (np.abs(averaged)**2).sum(axis=1))
    np.testing.assert_allclose(op.steinhardt(3, averaged=True), expected, atol=1e-14)
    np.testing.assert_allclose(op.tetrahedral_order(averaged=True), op.tetrahedral_order().mean())


def test_lsi_includes_first_outside_shell():
    positions = [[10, 10, 10], [11, 10, 10], [12.5, 10, 10], [13, 10, 10], [14, 10, 10]]
    op = StructuralOrder(oxygen_universe(positions), cutoff=5, n_jobs=2)
    assert op.local_structure_index()[0] == pytest.approx(np.var([1.5, .5, 1.]))
    with pytest.raises(ValueError, match="exceed"):
        op.local_structure_index(cutoff=3)


def test_undercoordinated_and_isolated():
    op = StructuralOrder(oxygen_universe([[0, 0, 0], [5, 5, 5]]), cutoff=2, n_jobs=2)
    np.testing.assert_array_equal(op.coordination_number(), 0)
    for result in (op.tetrahedral_order(), op.translational_order(), op.d5(),
                   op.steinhardt(), op.local_structure_index(shell=1)):
        assert np.isnan(result).all()
    np.testing.assert_array_equal(op.crystalline_connections(), 0)


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_hydrogen_bonds_periodic_and_duplicate_resids(water_universe, n_jobs):
    u = water_universe
    # Place one acceptor across the boundary along the donor's first O-H bond.
    u.atoms.positions = np.array([
        [19, 10, 10], [19.95, 10, 10], [18.76, 10.93, 10],
        [1.8, 10, 10], [2.75, 10, 10], [1.56, 10.93, 10],
        [5, 5, 5], [5.95, 5, 5], [4.76, 5.93, 5],
        [10, 5, 5], [10.95, 5, 5], [9.76, 5.93, 5],
        [15, 5, 5], [15.95, 5, 5], [14.76, 5.93, 5],
    ])
    u.residues.resids = [1, 1, 1, 1, 1]
    result = StructuralOrder(u, n_jobs=n_jobs).hydrogen_bonds()
    np.testing.assert_array_equal(result.donated, [1, 0, 0, 0, 0])
    np.testing.assert_array_equal(result.accepted, [0, 1, 0, 0, 0])
    np.testing.assert_array_equal(result.total, [1, 1, 0, 0, 0])
    np.testing.assert_array_equal(result.bonds, [[0, 1, 1]])
    np.testing.assert_array_equal(StructuralOrder(u).hydrogen_bonds(cutoff=2).total, 0)


def test_hydrogen_bonds_against_direct_geometry(water_universe):
    u = water_universe
    op = StructuralOrder(u, n_jobs=2)
    result = op.hydrogen_bonds(cutoff=5, angle=75)
    bonds = []
    for i in range(5):
        for j in range(5):
            if i == j:
                continue
            oo = minimize_vectors(op.pos[j:j+1]-op.pos[i], u.dimensions)[0]
            if np.linalg.norm(oo) > 5:
                continue
            for h in (1, 2):
                oh = minimize_vectors(u.atoms.positions[3*i+h:3*i+h+1]-op.pos[i], u.dimensions)[0]
                angle = np.rad2deg(np.arccos(np.clip(np.dot(oh, oo)/(np.linalg.norm(oh)*np.linalg.norm(oo)), -1, 1)))
                if angle < 75:
                    bonds.append((i, 3*i+h, j))
    assert sorted(map(tuple, result.bonds)) == sorted(bonds)


def test_snapshot_and_refresh(water_universe):
    op = OrderParameter(water_universe, n_jobs=2)
    old = op.F4().copy()
    water_universe.atoms.positions += [7, 3, 1]
    np.testing.assert_allclose(op.F4(), old)
    op.refresh()
    assert "f4" not in op.__dict__
    np.testing.assert_allclose(op.F4(), old, atol=1e-5)


@pytest.mark.parametrize("kwargs", [{"n_jobs": 0}, {"backend": "unknown"}, {"cutoff": -1}])
def test_invalid_configuration(water_universe, kwargs):
    with pytest.raises(ValueError):
        StructuralOrder(water_universe, **kwargs)


def test_invalid_calculation_inputs(water_universe):
    op = StructuralOrder(water_universe)
    for call in (lambda: op.d5(cutoff=0), lambda: op.kth_neighbor_distance(0),
                 lambda: op.smooth_coordination(n=28, m=12),
                 lambda: op.steinhardt(-1), lambda: op.steinhardt(4, normalized=True),
                 lambda: op.crystalline_connections(threshold=2),
                 lambda: op.hydrogen_bonds(angle=181)):
        with pytest.raises(ValueError):
            call()


def test_global_entropy_reference_and_all_atom_average(water_universe):
    obj = eFingerprint(water_universe, local=False, n_jobs=2)
    result = obj.Entropy()
    graph = obj.neighbor_list()
    expected = []
    for i in range(5):
        radii = graph.distances[graph.offsets[i]:graph.offsets[i+1]]
        r = obj.r[1:]
        g = np.exp(-(r[:, None]-radii)**2/(2*obj.sigma**2)).sum(axis=1)
        g /= 4*np.pi*obj.global_rho*r**2*np.sqrt(2*np.pi*obj.sigma**2)
        with np.errstate(divide="ignore", invalid="ignore"):
            integrand = np.where(g > 1e-10, g*np.log(g)-g+1, 1)*r**2
        expected.append(-2*np.pi*obj.global_rho*np.trapz(np.r_[0, integrand], obj.r))
    np.testing.assert_allclose(result, expected, atol=1e-10)
    assert not obj.local
    all_atom = eFingerprint(water_universe, include_hydrogen=True, n_jobs=2)
    values = all_atom.Entropy()
    assert values.shape == (5,)
    np.testing.assert_allclose(all_atom.localEntropy(), values.mean())


def test_entropy_honors_frame(water_universe):
    positions = water_universe.atoms.positions.copy()
    trajectory = np.stack([positions, positions*1.1])
    water_universe.load_new(trajectory, order="fac", dimensions=np.tile([20, 20, 20, 90, 90, 90], (2, 1)))
    obj = eFingerprint(water_universe, frame=1)
    np.testing.assert_allclose(obj.pos, trajectory[1, ::3])
    obj.refresh(0)
    np.testing.assert_allclose(obj.pos, trajectory[0, ::3])


@pytest.mark.parametrize("backend", ["threading", "loky"])
def test_parallel_chunk_backend(backend):
    offsets = np.arange(17, dtype=np.int64)
    distances = np.arange(1, 17, dtype=float)
    actual = map_chunks(_kernels.kth_chunk, 16, offsets, distances, 1, n_jobs=2, backend=backend)
    np.testing.assert_array_equal(actual, distances)


def test_serial_parallel_all_calculations(water_universe):
    serial = OrderParameter(water_universe, cutoff=5)
    parallel = OrderParameter(water_universe, cutoff=5, n_jobs=2)
    for method in ("coordination_number", "smooth_coordination", "d5", "OTO", "TTO",
                   "LSI", "F4", "get_dihedrals", "MinimumAngle"):
        np.testing.assert_allclose(getattr(parallel, method)(), getattr(serial, method)(), equal_nan=True)
    for kwargs in ({"l": 3}, {"l": 4, "averaged": True},
                   {"l": 6, "invariant": "w", "normalized": True}):
        np.testing.assert_allclose(parallel.steinhardt(**kwargs), serial.steinhardt(**kwargs), atol=1e-14)
    for local in (False, True):
        a = eFingerprint(water_universe, local=local)
        b = eFingerprint(water_universe, local=local, n_jobs=2)
        np.testing.assert_allclose(a.Entropy(), b.Entropy())
        np.testing.assert_allclose(a.localEntropy(), b.localEntropy())


def test_process_backend_analysis(water_universe):
    a = OrderParameter(water_universe, cutoff=5)
    b = OrderParameter(water_universe, cutoff=5, n_jobs=2, backend="loky")
    for method in ("OTO", "F4", "steinhardt"):
        np.testing.assert_allclose(getattr(a, method)(), getattr(b, method)(), atol=1e-14)


def test_isolated_entropy_limits():
    u = oxygen_universe([[0, 0, 0], [5, 5, 5]])
    obj = eFingerprint(u, cutoff=2, local=False, sigma=.02)
    expected = -2*np.pi*obj.global_rho*np.trapz(obj.r**2, obj.r)
    np.testing.assert_allclose(obj.Entropy(), expected)
    assert np.isnan(eFingerprint(u, cutoff=2, local=True).Entropy()).all()


def test_missing_hydrogens_and_box():
    u = oxygen_universe([[1, 1, 1], [3, 1, 1]])
    with pytest.raises(ValueError, match="two selected hydrogens"):
        StructuralOrder(u).hydrogen_bonds()
    u.dimensions = None
    with pytest.raises(ValueError, match="periodic box"):
        StructuralOrder(u)
