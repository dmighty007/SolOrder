import MDAnalysis as mda
import numpy as np
import pytest


@pytest.fixture
def water_universe():
    """Five waters with an ideal tetrahedral shell around the first oxygen."""
    directions = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
    oxygens = np.vstack([np.zeros(3), directions * (2.8 / np.sqrt(3))]) + 10
    offsets = np.array([[0, 0, 0], [0.9572, 0, 0], [-0.239, 0.927, 0]])
    universe = mda.Universe.empty(
        15, n_residues=5, atom_resindex=np.repeat(np.arange(5), 3), trajectory=True
    )
    universe.add_TopologyAttr("names", ["OW", "HW1", "HW2"] * 5)
    universe.add_TopologyAttr("ids", np.arange(1, 16))
    universe.add_TopologyAttr("resids", np.arange(1, 6))
    universe.atoms.positions = (oxygens[:, None, :] + offsets).reshape(-1, 3)
    universe.dimensions = [20, 20, 20, 90, 90, 90]
    return universe
