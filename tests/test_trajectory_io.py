"""External-file readers are tested on tiny synthetic files, not user datasets."""

import zlib
import numpy as np
import pytest

from SolOrder.trajectory_io import sample_indices, scan_coordinate_chunks, read_coordinate_chunk, index_xyz, read_xyz_frame


def test_sampling_endpoints_and_extra_indices():
    assert sample_indices(10, 3, extra=(3, 99)) == [0, 3, 4, 9]
    assert sample_indices(2, 8) == [0, 1]
    with pytest.raises(ValueError):
        sample_indices(0, 1)


def test_bounded_complete_coordinate_recovery(tmp_path):
    a = np.arange(24, dtype="<f8").reshape(8, 3)
    b = a+.25
    first, second = zlib.compress(a.tobytes()), zlib.compress(b.tobytes())
    truncated = zlib.compress((a+1).tobytes())[:-3]
    path = tmp_path/"damaged.h5"
    path.write_bytes(b"unreadable HDF header\x78\x9c junk"+first+b"other metadata"+second+truncated)
    records = scan_coordinate_chunks(path, 8)
    assert len(records) == 2
    np.testing.assert_array_equal(read_coordinate_chunk(path, records[0], 8), a)
    np.testing.assert_array_equal(read_coordinate_chunk(path, records[1], 8), b)


def test_xyz_index_and_boxes(tmp_path):
    path = tmp_path/"movie.xyz"
    path.write_text("2\nframe 8 box 20.5\nO 1 2 3\nO 4 5 6\n2\nframe 10\nO 7 8 9\nO 1 3 4\n")
    records = index_xyz(path)
    assert len(records) == 2
    points, box, original = read_xyz_frame(path, records[0])
    np.testing.assert_array_equal(points, [[1, 2, 3], [4, 5, 6]])
    assert box == 20.5 and original == 8
    assert read_xyz_frame(path, records[1])[1:] == (None, 10)


def test_wedge_cache_invalidated(water_universe):
    from SolOrder import OrderParameter
    op = OrderParameter(water_universe)
    before = op._wedges()
    assert op._wedges() is before
    op.refresh()
    assert not op._wedge_cache
    after = op._wedges()
    assert after is not before


def test_pdb_snapshot_does_not_reload_topology_coordinates(tmp_path, water_universe):
    from SolOrder.trajectory_io import snapshot_template
    from SolOrder import OrderParameter, eFingerprint
    path = tmp_path/"water.pdb"
    with pytest.warns(UserWarning):
        water_universe.atoms.write(path)
        u = snapshot_template(path)
    positions = water_universe.atoms.positions+np.array([3, 2, 1], dtype=np.float32)
    u.atoms.positions = positions
    u.dimensions = [20, 20, 20, 90, 90, 90]
    op = OrderParameter(u)
    np.testing.assert_array_equal(op.pos, positions[::3])
    eFingerprint(u).Entropy()
    np.testing.assert_array_equal(u.atoms.positions, positions)
