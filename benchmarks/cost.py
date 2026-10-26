"""Wall time of every registry entry against system size.

Systems are periodic tilings of one thermal frame of run C1 at 220 K (bulk Ih,
1920 molecules; liquid, 4096 molecules), from about 2,000 to 65,000
molecules. Each measurement starts from a fresh OrderParameter/eFingerprint
snapshot, as compute.py does, so it includes neighbor-list construction:

  entry       one registry entry alone on a fresh snapshot
  snapshot    building the snapshot (coordinates and box copy)
  all         every entry on one shared snapshot (shared neighbor caches)

Numba kernels are compiled first on the smallest system; reported times are
the median of three warm repeats. ``--jobs 1`` limits freud (Voronoi) to one
thread as well; OVITO (PTM) always uses its own thread pool.

  OMP_NUM_THREADS=1 python cost.py --jobs 1   # results/cost_jobs1.csv
  python cost.py --jobs 8                     # results/cost_jobs8.csv
"""

import argparse
from pathlib import Path
import platform
import time
import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

import MDAnalysis as mda
import numpy as np
import pandas as pd

from SolOrder import OrderParameter, eFingerprint

from compute import PREPS
from registry import REGISTRY

RESULTS = Path(__file__).resolve().parent/"results"
TILINGS = {"Ih": ((1, 1, 1), (2, 1, 1), (2, 2, 1), (2, 2, 2), (4, 2, 2), (4, 4, 2)),
           "liquid": ((1, 1, 1), (2, 1, 1), (2, 2, 1), (2, 2, 2), (4, 2, 2))}
REPEATS = 3


def tiled(phase, reps):
    base = mda.Universe(str(PREPS/f"{phase}_p1.gro"))
    frames = np.load(PREPS/f"{phase}_p1_220.npz")
    box = frames["box"]
    shifts = np.array([(i, j, k) for i in range(reps[0]) for j in range(reps[1]) for k in range(reps[2])])*box
    u = mda.Merge(*[base.atoms]*len(shifts))
    u.atoms.positions = np.concatenate([frames["x"][0] + s for s in shifts])
    u.dimensions = [*(box*reps), 90., 90., 90.]
    return u


def snapshot(u, jobs):
    op = OrderParameter(u, cutoff=5., selection="name OW", hydrogen_selection="name HW1 HW2",
                        n_jobs=jobs)
    ent = eFingerprint(u, cutoff=5., selection="name OW", n_jobs=jobs)
    return op, ent


def timed(function):
    times = []
    for _ in range(REPEATS):
        start = time.perf_counter()
        function()
        times.append(time.perf_counter() - start)
    return float(np.median(times))


def measure(jobs):
    if jobs == 1:
        import freud
        freud.parallel.set_num_threads(1)
    small = tiled("Ih", (1, 1, 1))
    op, ent = snapshot(small, jobs)
    for entry in REGISTRY:  # compile every kernel before timing
        entry.columns(op, ent)
    rows = []
    for phase, tilings in TILINGS.items():
        for reps in tilings:
            u = tiled(phase, reps)
            n = len(u.select_atoms("name OW"))
            base = dict(phase=phase, molecules=n, jobs=jobs)
            rows.append(dict(base, entry="snapshot", seconds=timed(lambda: snapshot(u, jobs))))

            def run_all():
                op, ent = snapshot(u, jobs)
                for entry in REGISTRY:
                    entry.columns(op, ent)
            rows.append(dict(base, entry="all", seconds=timed(run_all)))
            for entry in REGISTRY:
                rows.append(dict(base, entry=entry.name,
                                 seconds=timed(lambda: entry.columns(*snapshot(u, jobs)))))
            done = rows[-len(REGISTRY)-2:]
            print(phase, n, " ".join(f"{r['entry']}={r['seconds']:.2f}" for r in done), flush=True)
    out = pd.DataFrame(rows)
    out["us_per_molecule"] = 1e6*out.seconds/out.molecules
    out["cpu"] = platform.processor() or platform.machine()
    RESULTS.mkdir(exist_ok=True)
    out.to_csv(RESULTS/f"cost_jobs{jobs}.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    measure(args.jobs)
