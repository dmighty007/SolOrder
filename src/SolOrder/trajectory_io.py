"""Bounded, read-only recovery and sampling for external validation data."""

import mmap
import re
import zlib

import numpy as np


ZLIB_HEADER = re.compile(b"\x78[\x01\x5e\x9c\xda]")


def sample_indices(size, count, extra=()):
    """Evenly spaced indices including endpoints and valid requested extras."""
    if size < 1 or count < 1:
        raise ValueError("size and count must be positive")
    selected = np.linspace(0, size-1, min(size, count), dtype=int).tolist()
    selected += [int(i) for i in extra if 0 <= int(i) < size]
    return sorted(set(selected))


def scan_coordinate_chunks(path, n_atoms):
    """Index complete zlib streams encoding (n_atoms,3) little-endian float64.

    Used only when the checkpoint's HDF5 metadata is unreadable. Recovered order
    is physical storage order, not a guaranteed logical HDF5 trajectory order.
    Decompression is bounded and exact-size/eof/finite checks reject false hits.
    """
    expected = n_atoms*3*8
    records = []
    with open(path, "rb") as handle, mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
        cursor = 0
        while True:
            match = ZLIB_HEADER.search(mapped, cursor)
            if match is None:
                break
            offset = match.start()
            block = mapped[offset:min(len(mapped), offset+expected+4096)]
            decoder = zlib.decompressobj()
            try:
                decoded = decoder.decompress(block, expected+1)
                valid = decoder.eof and len(decoded) == expected
                if valid:
                    points = np.frombuffer(decoded, dtype="<f8")
                    valid = np.isfinite(points).all() and np.max(np.abs(points)) < 1e7
                if valid:
                    consumed = len(block)-len(decoder.unused_data)
                    records.append((offset, consumed))
                    cursor = offset+consumed
                    continue
            except zlib.error:
                pass
            cursor = offset+1
    if not records:
        raise ValueError(f"No valid coordinate chunks found in {path}")
    return records


def read_coordinate_chunk(path, record, n_atoms):
    offset, size = record
    with open(path, "rb") as handle:
        handle.seek(offset)
        data = handle.read(size)
    return np.frombuffer(zlib.decompress(data), dtype="<f8").reshape(n_atoms, 3).copy()


def index_xyz(path):
    """Index XYZ frame offsets and comments without retaining coordinates."""
    result = []
    with open(path, "rb") as handle:
        while True:
            offset = handle.tell()
            count_line = handle.readline()
            if not count_line:
                break
            count = int(count_line)
            comment = handle.readline().decode().strip()
            for _ in range(count):
                if not handle.readline():
                    raise ValueError("Truncated XYZ frame")
            result.append((offset, count, comment))
    return result


def read_xyz_frame(path, record):
    offset, count, comment = record
    with open(path, "rb") as handle:
        handle.seek(offset)
        handle.readline()
        handle.readline()
        rows = [handle.readline().split() for _ in range(count)]
    positions = np.array([[float(x) for x in row[1:4]] for row in rows])
    tokens = comment.split()
    box = float(tokens[tokens.index("box")+1]) if "box" in tokens else None
    original_frame = int(tokens[tokens.index("frame")+1]) if "frame" in tokens else None
    return positions, box, original_frame


def snapshot_template(topology):
    """Copy topology/coordinates into MemoryReader; no file reader can reload it."""
    import MDAnalysis as mda
    return mda.Merge(mda.Universe(topology).atoms)
