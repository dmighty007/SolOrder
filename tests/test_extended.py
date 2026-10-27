"""Physical/reference checks for the remaining geometric manuscript methods."""

from itertools import combinations, product
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.lib.distances import calc_angles, minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors

from SolOrder import OrderParameter, StructuralOrder
from SolOrder import _extended_kernels as kernels
from test_structure import oxygen_universe


DATA = Path(__file__).parent/"data/genice"
PHASES = ("Ih", "Ic", "II", "III", "V", "VI", "VII", "VIII", "XII", "XVII", "sI")


def ice_analysis(phase, **kwargs):
    return OrderParameter(mda.Universe(str(DATA/f"{phase}.gro")), selection="name O",
                          hydrogen_selection="name H", cutoff=5, **kwargs)


def test_full_water_angle_distribution(water_universe):
    op = OrderParameter(water_universe, cutoff=5, n_jobs=2)
    result = op.water_angular_distribution(cutoff=5, bins=np.linspace(0, np.pi, 31))
    assert len(result.angles) == 10
    assert len(set(map(tuple, result.pairs))) == 10
    assert result.counts.sum() == 10
    assert np.sum(result.probability_density*np.diff(result.bin_edges)) == pytest.approx(1)
    for (i, j), angle in zip(result.pairs, result.angles):
        hi = water_universe.residues[i].atoms.select_atoms("name HW1 HW2").positions
        hj = water_universe.residues[j].atoms.select_atoms("name HW1 HW2").positions
        oi, oj = op.pos[i], op.pos[j]
        angles = [calc_angles(h, oi, oj, box=op.box) for h in hi]
        angles += [calc_angles(h, oj, oi, box=op.box) for h in hj]
        assert angle == pytest.approx(min(angles), abs=2e-6)
    empty = op.water_angular_distribution(cutoff=.1)
    assert len(empty.angles) == 0 and np.isnan(empty.probability_density).all()
    with pytest.raises(ValueError):
        op.water_angular_distribution(bins=[0, 1])


def test_centrosymmetry_reference_and_pairing():
    # Three opposite pairs give a perfect centrosymmetric central site.
    vectors = np.vstack((np.eye(3), -np.eye(3)))
    positions = np.vstack((np.zeros(3), vectors))+10
    op = StructuralOrder(oxygen_universe(positions), cutoff=2, n_jobs=2)
    assert op.centrosymmetry(6)[0] == pytest.approx(0)
    distorted = positions.copy()
    distorted[1] += [.1, .2, 0]
    op = StructuralOrder(oxygen_universe(distorted), cutoff=2)
    g = op.neighbor_list()
    v = g.vectors[g.offsets[0]:g.offsets[1]][:6]
    weights = sorted(np.dot(a+b, a+b) for a, b in combinations(v, 2))
    assert op.centrosymmetry(6)[0] == pytest.approx(sum(weights[:3]))
    assert op.centrosymmetry(6, pairing="rank")[0] == pytest.approx(np.sum((v[:3]+v[3:])**2))
    assert np.isnan(op.centrosymmetry(12)).all()
    with pytest.raises(ValueError):
        op.centrosymmetry(3)


def test_angular_and_distance_parameters_reference():
    positions = [[10, 10, 10], [11, 10, 10], [10, 12, 10], [10, 10, 13]]
    op = StructuralOrder(oxygen_universe(positions), cutoff=3.2, n_jobs=2)
    # Central atom has three mutually perpendicular bonds, lengths 1,2,3.
    assert op.bond_angle_parameter(n1=2, n2=2)[0] == pytest.approx(1)
    assert op.angular_fourier(a=2, radial=0)[0] == pytest.approx(-1)
    assert op.angular_fourier(a=2, radial=1)[0] == pytest.approx(-4/3)
    expected = np.mean([1*2*np.sqrt(5), 1*3*np.sqrt(10), 2*3*np.sqrt(13)])
    assert op.neighbor_distance_parameter()[0] == pytest.approx(expected)
    # Different alpha/beta functions are symmetric under exchange by default.
    expected = np.mean([.5*(a*b*b+b*a*a)*np.sqrt(a*a+b*b) for a, b in combinations((1, 2, 3), 2)])
    assert op.neighbor_distance_parameter(f_alpha=lambda x: x, f_beta=2)[0] == pytest.approx(expected)
    expected_ranked = np.mean([a*b*b*np.sqrt(a*a+b*b) for a, b in combinations((1, 2, 3), 2)])
    assert op.neighbor_distance_parameter(f_beta=2, symmetrize=False)[0] == pytest.approx(expected_ranked)
    np.testing.assert_allclose(op.bond_angle_parameter(averaged=True), op.local_average(op.bond_angle_parameter()))
    histogram = op.angle_histogram_parameter(frequencies=[0, 1, 2], bins=20)
    np.testing.assert_allclose(histogram[0], 1.)  # All central angles in the same bin.
    with pytest.raises(ValueError):
        op.neighbor_distance_parameter(f_alpha=lambda x: x*np.nan)


def test_neighborhood_distinct_variants_reference():
    op = StructuralOrder(oxygen_universe([[10, 10, 10], [11, 10, 10], [10, 12, 10], [10, 10, 13]]), cutoff=3.2, n_jobs=2)
    # Independent direct sum using source-verified M-neighbor ranking.
    v = np.diag([1., 2., 3.])
    a_terms, p_terms = [], []
    for j in range(3):
        k = min((k for k in range(3) if k != j), key=lambda k: np.linalg.norm(v[k])+np.linalg.norm(v[k]-v[j]))
        a_terms.append(2*v[k]-v[j])
        p_terms.append(2*v[j]-v[k])
    expected = [np.sum(np.array(a_terms)**2)/3, np.sum(np.array(p_terms)**2)/3,
                np.sum(np.sum(p_terms, axis=0)**2)/3]
    np.testing.assert_allclose(op.neighborhood_parameters()[0], expected)
    assert expected[1] != expected[2]
    np.testing.assert_allclose(op.neighborhood_parameters(averaged=True), op.local_average(op.neighborhood_parameters()))


def test_quadrupole_geometry(water_universe):
    op = OrderParameter(water_universe)
    oh = op._water_geometry()
    length = np.linalg.norm(oh, axis=2).mean(axis=1)
    theta = np.arccos(np.sum(oh[:, 0]*oh[:, 1], axis=1)/(np.linalg.norm(oh[:, 0], axis=1)*np.linalg.norm(oh[:, 1], axis=1)))/2
    expected = 3*(length*np.sin(theta))**2*2.35/(4*length*np.cos(theta))
    np.testing.assert_allclose(op.quadrupole_parameter(2.35), expected)
    assert np.all(op.quadrupole_parameter(2.35, negative_site_distance=.15) > expected)
    with pytest.raises(ValueError):
        op.quadrupole_parameter(2.35, negative_site_distance=2.)


def test_periodic_voronoi_exact_cubic_cells():
    positions = np.array(list(product((1., 5.), repeat=3)))
    op = StructuralOrder(oxygen_universe(positions, dimensions=[8, 8, 8, 90, 90, 90]))
    geometry = op.voronoi_geometry()
    np.testing.assert_allclose(geometry.volume, 64, rtol=1e-6)
    np.testing.assert_allclose(geometry.area, 96, rtol=1e-6)
    np.testing.assert_allclose(geometry.density, 1/64, rtol=1e-6)
    np.testing.assert_allclose(geometry.asphericity, 6/np.pi, rtol=1e-6)
    op.refresh()
    assert op._voronoi is None


@pytest.mark.parametrize("phase", PHASES)
def test_genice_voronoi_and_extended_serial_parallel(phase):
    a, b = ice_analysis(phase), ice_analysis(phase, n_jobs=2)
    voronoi = a.voronoi_geometry()
    assert voronoi.volume.sum() == pytest.approx(a.volume, rel=2e-6)
    assert np.isfinite(voronoi.asphericity).all() and np.all(voronoi.asphericity >= 1)
    for method, kwargs in (("centrosymmetry", {"n_neighbors": 4}),
                           ("neighbor_distance_parameter", {}), ("angular_fourier", {}),
                           ("bond_angle_parameter", {}), ("angle_histogram_parameter", {}),
                           ("neighborhood_parameters", {})):
        values = getattr(a, method)(**kwargs)
        assert np.isfinite(values).all()
        np.testing.assert_allclose(values, getattr(b, method)(**kwargs), atol=1e-12)
    angular = a.water_angular_distribution()
    np.testing.assert_allclose(angular.angles, b.water_angular_distribution().angles)
    assert angular.counts.sum() == len(angular.pairs)
    assert np.sum(angular.probability_density*np.diff(angular.bin_edges)) == pytest.approx(1)
    for plus in (False, True):
        np.testing.assert_array_equal(a.chill(plus=plus).labels, b.chill(plus=plus).labels)
    for phase_expected, code in (("Ih", 1), ("Ic", 2)):
        if phase == phase_expected:
            np.testing.assert_array_equal(a.chill().labels, code)
            np.testing.assert_array_equal(a.chill_plus().labels, code)
    if phase == "sI":
        np.testing.assert_array_equal(a.chill_plus().labels, 4)


def test_extended_process_backend(water_universe):
    a = OrderParameter(water_universe, cutoff=5)
    b = OrderParameter(water_universe, cutoff=5, n_jobs=2, backend="loky")
    for method in ("neighbor_distance_parameter", "angular_fourier", "bond_angle_parameter",
                   "angle_histogram_parameter", "neighborhood_parameters"):
        np.testing.assert_allclose(getattr(a, method)(), getattr(b, method)(), atol=1e-12)
    np.testing.assert_allclose(a.water_angular_distribution(5).angles, b.water_angular_distribution(5).angles)


def test_chill_interfacial_rules():
    offsets = np.arange(0, 25, 4)
    neighbors = np.array([1, 2, 3, 4]*6)
    neighbors[4] = 0  # The 3-staggered site has a 2-staggered neighbor.
    staggered = np.array([2, 3, 4, 3, 0, 0])
    eclipsed = np.array([0, 0, 0, 1, 4, 3])
    labels = kernels.classify_chunk(0, 6, offsets, neighbors, staggered, eclipsed, True)
    np.testing.assert_array_equal(labels, [3, 3, 2, 1, 4, 5])
    labels = kernels.classify_chunk(0, 6, offsets, neighbors, staggered, eclipsed, False)
    np.testing.assert_array_equal(labels, [3, 3, 2, 1, 0, 0])


def test_triclinic_voronoi_translation_and_duplicate_rejection():
    op = ice_analysis("V")
    before = op.voronoi_geometry()
    cell = triclinic_vectors(op.box, dtype=float)
    shifted = StructuralOrder(oxygen_universe(op.pos+2*cell[0]-cell[2], op.box))
    np.testing.assert_allclose(before.volume, shifted.voronoi_geometry().volume, atol=2e-4)
    duplicate = StructuralOrder(oxygen_universe([[0, 0, 0], [20, 0, 0]]))
    with pytest.raises(ValueError, match="distinct"):
        duplicate.voronoi_geometry()


def test_fourier_product_and_histogram_dft_reference():
    op = StructuralOrder(oxygen_universe([[10, 10, 10], [11, 10, 10], [10, 12, 10], [10, 10, 13]]), cutoff=3.2)
    expected = -np.mean([1*2, 1*3, 2*3])
    assert op.angular_fourier(a=2, radial=1, weighting="product")[0] == pytest.approx(expected)
    offsets, _, _, angles, _ = op._wedges()
    histogram, _ = np.histogram(angles[offsets[1]:offsets[2]], bins=25, range=(0, np.pi))
    reference = np.abs(np.fft.rfft(histogram/histogram.sum()))[:4]
    np.testing.assert_allclose(op.angle_histogram_parameter([0, 1, 2, 3], bins=25)[1], reference)
    np.testing.assert_allclose(op.angle_histogram_parameter([0, 1, 2, 3], bins=25, power=True)[1], reference**2)


def test_ring_statistics_square_lattice_and_periodic_wrapping():
    grid = np.array([(2.8*x, 2.8*y, 10.) for x in range(6) for y in range(6)])
    op = StructuralOrder(oxygen_universe(grid, (16.8, 16.8, 20, 90, 90, 90)), cutoff=3)
    rings = op.ring_statistics(8, network="distance", cutoff=3)
    # Four right angles close 4-rings; two straight angles need a 6-ring detour.
    np.testing.assert_array_equal(rings.counts[:, rings.sizes == 4], 4)
    np.testing.assert_array_equal(rings.counts[:, rings.sizes == 6], 2)
    assert rings.counts.sum() == 6*len(grid) and not rings.unclosed.any()
    np.testing.assert_allclose(rings.mean_size, 14/3)
    # A five-atom chain closes only by wrapping around the box: not a ring.
    chain = np.array([(2.8*x, 10., 10.) for x in range(5)])
    wrapped = StructuralOrder(oxygen_universe(chain, (14, 20, 20, 90, 90, 90)), cutoff=3)
    rings = wrapped.ring_statistics(10, network="distance", cutoff=3)
    assert not rings.counts.any()
    np.testing.assert_array_equal(rings.unclosed, 1)
    assert np.isnan(rings.mean_size).all()
    with pytest.raises(ValueError):
        op.ring_statistics(2)
    with pytest.raises(ValueError):
        op.ring_statistics(network="voronoi")


@pytest.mark.parametrize("phase", PHASES)
def test_genice_zeta_reference_and_rings(phase):
    a, b = ice_analysis(phase), ice_analysis(phase, n_jobs=2)
    oxygen = a.pos
    distances = np.linalg.norm(minimize_vectors(
        (oxygen[None]-oxygen[:, None]).reshape(-1, 3), a.dimensions), axis=1).reshape(len(oxygen), -1)
    np.fill_diagonal(distances, np.inf)
    bonds = a.hydrogen_bonds().bonds
    bonded = np.zeros_like(distances, dtype=bool)
    bonded[bonds[:, 0], bonds[:, 2]] = bonded[bonds[:, 2], bonds[:, 0]] = True
    free = ~bonded & (distances <= 5)
    expected = np.where(free, distances, np.inf).min(1) - np.where(bonded, distances, -np.inf).max(1)
    np.testing.assert_allclose(a.zeta(), expected, atol=1e-10)
    rings, parallel = a.ring_statistics(), b.ring_statistics()
    np.testing.assert_array_equal(rings.counts, parallel.counts)
    assert not rings.unclosed.any()
    present = set(rings.sizes[rings.counts.sum(0) > 0])
    known = {"Ih": {6}, "Ic": {6}, "VI": {4, 8}, "VII": {6}, "VIII": {6}, "XII": {7, 8}}
    if phase in known:
        assert present == known[phase]
    with pytest.raises(ValueError):
        a.zeta(search=3.5)


def test_minimum_image_matches_mdanalysis_in_skewed_cell():
    from SolOrder.extended import _cell_geometry
    box = np.array([20., 23., 26., 62., 75., 110.])
    cell, inverse, half_height = _cell_geometry(box)
    d = np.random.default_rng(3).uniform(-30, 30, (500, 3))
    ours = np.array([kernels.minimum_image(v, cell, inverse, half_height) for v in d])
    reference = np.linalg.norm(minimize_vectors(d, box), axis=1)  # MDAnalysis works in float32
    np.testing.assert_allclose(np.linalg.norm(ours, axis=1), reference, rtol=1e-6)
    assert np.all(np.linalg.norm(ours, axis=1) <= reference + 1e-5)
