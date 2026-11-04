"""Write the small example structures used by the beginner's guide (docs/guide/data).

Takes thermal TIP4P/Ice frames at 220 K from the benchmark cache: ice Ih, ice Ic
and supercooled liquid (O and H only; the massless M site is dropped), a short
ice Ih trajectory, and one frame of a 1600-molecule ice Ih seed in liquid
(oxygens only, gzipped) with its constructed-seed mask.

  python scripts/make_guide_data.py
"""

import gzip
import shutil
from pathlib import Path

import MDAnalysis as mda
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPLICAS = ROOT/"benchmarks/cache/replicas"
INTERFACE = ROOT/"benchmarks/cache/interface"
OUT = ROOT/"docs/guide/data"
STRIDE_PS = 5.  # saved-frame spacing of the replica runs


def water(gro, frames, box, out, trajectory=None):
    u = mda.Universe(str(gro))
    u.atoms.positions = frames[-1]
    u.dimensions = [*box, 90., 90., 90.]
    keep = u.select_atoms("name OW HW1 HW2")
    keep.write(str(out))
    if trajectory:
        with mda.Writer(str(trajectory), keep.n_atoms) as w:
            for i, x in enumerate(frames):
                u.atoms.positions = x
                u.trajectory.ts.time = STRIDE_PS*i
                w.write(keep)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for phase, name in (("Ih", "ice_Ih"), ("Ic", "ice_Ic"), ("liquid", "liquid")):
        d = np.load(REPLICAS/f"{phase}_p1_220.npz")
        traj = OUT/"ice_Ih.xtc" if phase == "Ih" else None
        frames = d["x"][-10:] if traj else d["x"]
        water(REPLICAS/f"{phase}_p1.gro", frames, d["box"], OUT/f"{name}.gro", traj)

    seed = "F_seed_Ih_1600_p1"
    d = np.load(INTERFACE/f"{seed}.npz")
    u = mda.Universe(str(INTERFACE/f"{seed}.gro"))
    u.atoms.positions = d["x"][-1]
    u.dimensions = [*d["box"], 90., 90., 90.]
    u.select_atoms("name OW").write(str(OUT/"seed.gro"))
    with open(OUT/"seed.gro", "rb") as src, gzip.open(OUT/"seed.gro.gz", "wb", compresslevel=9) as dst:
        shutil.copyfileobj(src, dst)
    (OUT/"seed.gro").unlink()
    truth = np.load(INTERFACE/f"{seed}.truth.npz")
    np.savez_compressed(OUT/"seed_truth.npz", seed=truth["seed"], centre=truth["centre"], R=truth["R"])


if __name__ == "__main__":
    main()
