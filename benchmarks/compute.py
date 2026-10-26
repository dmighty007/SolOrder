"""Compute every registry column for the pilot datasets, one row per water.

Pilot datasets (docs/benchmark_design.md, milestone 1):
  genice    19 ideal GenIce phases x 5 GenIce seeds (cache/genice/seed<N>,
            from scripts/generate_genice_fixtures.py), each jittered once
            (Gaussian 0.02 A on every atom, about twice the .gro
            coordinate precision) so that "separable" is judged at a stated
            resolution rather than at rounding noise.
  seedtest  TIP4P/Ice 220 K bulk Ih/Ic/stacking-disordered crystals and
            spherical seeds in liquid, thermal and quenched frames, with
            per-molecule ground truth (ice_nucleation_pg/seedtest). Run A.
  replica   Independent TIP4P/Ice NVT runs from the same starting structures
            (md_train.py): bulk Ih/Ic/Isd and bulk liquid at 220 K, Ih/Ic at
            250 K, 60 thermal and 60 quenched frames 5 ps apart. Run B; used
            as held-out data for thresholds fitted on run A, and vice versa.
  prep      Independent preparations (md_replicas.py, run C1-C5): new proton
            configurations (and stacking for Isd), re-melted liquid, same
            protocol and frame spacing as run B.
  model     Second and third water models (run D): TIP4P/2005 on run C
            structures at 220 K and 200 K and a TIP4P/Ice control with the
            same NPT-then-NVT protocol (md_conditions.py); oxygen-only mW at
            220 K (md_mw.py; hydrogen-dependent columns are NaN).
  pressure  Thermal high-pressure phases (md_conditions.py, run E): TIP4P/Ice
            ices II, III, V, VI, VII, VIII and liquid at 0.3 and 1.0 GPa, 220 K.

Large outputs go to benchmarks/cache/ (not committed).
"""

import argparse
import hashlib
import json
from pathlib import Path
import time
import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

from joblib import Parallel, delayed
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_box
import numpy as np
import pandas as pd

import SolOrder
from SolOrder import OrderParameter, eFingerprint

from registry import REGISTRY

ROOT = Path(__file__).resolve().parents[1]
GENICE = ROOT/"benchmarks/cache/genice"
SEEDTEST = Path("/home/dm/Dibyendu/Projects/ice_nucleation_pg/data/seedtest")
CACHE = ROOT/"benchmarks/cache"
JITTER = .02  # A
SEEDS = (1, 2, 3, 4, 5)  # GenIce seeds; proton-ordered phases repeat
# Ice XVI is the emptied sII framework; GenIce gives identical oxygen positions.
DUPLICATES = {"XVI": "sII"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# Columns of the hydrogen-dependent entries; NaN for oxygen-only models (mW).
H_COLUMNS = ("f4", "n_hb", "zeta", "ring_mean", *(f"ring_f{k}" for k in (4, 5, 6, 7, 8)))


def evaluate(u, oxygen, hydrogen):
    oxygen_only = hydrogen is None
    op = OrderParameter(u, cutoff=5., selection=oxygen, hydrogen_selection=hydrogen or oxygen)
    ent = eFingerprint(u, cutoff=5., selection=oxygen)
    out = {}
    for entry in REGISTRY:
        if oxygen_only and entry.needs_h:
            continue
        for name, values in entry.columns(op, ent).items():
            out[name] = np.asarray(values, dtype=np.float32)
    if oxygen_only:
        n = len(next(iter(out.values())))
        out.update({name: np.full(n, np.nan, dtype=np.float32) for name in H_COLUMNS})
    return out


def genice_task(phase, seed):
    u = mda.Universe(str(GENICE/f"seed{seed}/{phase}.gro"))
    rng = np.random.default_rng([12345, seed, sum(map(ord, phase))])
    u.atoms.positions = u.atoms.positions + rng.normal(0., JITTER, u.atoms.positions.shape)
    columns = evaluate(u, "name O", "name H")
    n = len(next(iter(columns.values())))
    meta = dict(dataset="genice", system=phase, kind="ideal", frame=seed,
                mol=np.arange(n), phase=phase)
    return pd.DataFrame({**meta, **columns})


def seed_centroid_distance(oxygen, seed, box):
    ref = oxygen[np.flatnonzero(seed)[0]]
    offsets = minimize_vectors(oxygen[seed]-ref, box)
    centre = ref + offsets.mean(axis=0)
    return np.linalg.norm(minimize_vectors(oxygen-centre, box), axis=1)


def seedtest_task(system, kind, frame):
    u = mda.Universe(str(SEEDTEST/f"{system}.gro"))
    frames = np.load(SEEDTEST/f"{'frames' if kind == 'thermal' else 'framesq'}_{system}.npz")
    u.atoms.positions = frames["x"][frame]
    u.dimensions = [*frames["box"], 90., 90., 90.]
    columns = evaluate(u, "name OW", "name HW1 HW2")
    truth = np.load(SEEDTEST/"truth.npz")
    seed = truth[f"{system}__seed"].astype(bool)
    oxygen = u.select_atoms("name OW").positions
    meta = dict(dataset="seedtest", system=system, kind=kind, frame=frame,
                mol=np.arange(len(seed)), phase=system.split("_")[1],
                seed=seed, interior=truth[f"{system}__interior"].astype(bool),
                stack=truth[f"{system}__stack"].astype(np.int8),
                r_seed=seed_centroid_distance(oxygen, seed, u.dimensions).astype(np.float32),
                R=float(truth[f"{system}__R"]) if f"{system}__R" in truth.files else np.inf)
    return pd.DataFrame({**meta, **columns})


def replica_task(system, kind, frame):
    structure = system.removesuffix("_250K")
    u = mda.Universe(str(SEEDTEST/f"{structure}.gro"))
    frames = np.load(SEEDTEST/f"train_{system}.npz")
    u.atoms.positions = frames["x" if kind == "thermal" else "xq"][frame]
    u.dimensions = [*frames["box"], 90., 90., 90.]
    columns = evaluate(u, "name OW", "name HW1 HW2")
    n = len(next(iter(columns.values())))
    truth = np.load(SEEDTEST/"truth.npz")
    stack = truth[f"{structure}__stack"].astype(np.int8) if f"{structure}__stack" in truth.files \
        else np.zeros(n, dtype=np.int8)
    meta = dict(dataset="replica", system=structure, kind=kind, frame=frame, mol=np.arange(n),
                phase=structure.removeprefix("bulk_"), stack=stack,
                temperature=250. if system.endswith("_250K") else 220.)
    return pd.DataFrame({**meta, **columns})


PREPS = Path(__file__).resolve().parent/"cache/replicas"


def prep_task(run, kind, frame):
    phase, prep, temperature = run.split("_")
    u = mda.Universe(str(PREPS/f"{phase}_{prep}.gro"))
    frames = np.load(PREPS/f"{run}.npz")
    u.atoms.positions = frames["x" if kind == "thermal" else "xq"][frame]
    u.dimensions = [*frames["box"], 90., 90., 90.]
    columns = evaluate(u, "name OW", "name HW1 HW2")
    n = len(next(iter(columns.values())))
    stack_file = PREPS/f"{phase}_{prep}.stack.npy"
    stack = np.load(stack_file) if stack_file.exists() else np.zeros(n, dtype=np.int8)
    system = "liquid" if phase == "liquid" else f"bulk_{phase}"
    meta = dict(dataset="prep", system=system, kind=kind, frame=frame, mol=np.arange(n), phase=phase,
                stack=stack, temperature=float(temperature), run=f"C{prep[1:]}")
    return pd.DataFrame({**meta, **columns})


def prep_runs():
    return sorted(p.stem for p in PREPS.glob("*_p*_*.npz"))


CONDITIONS = Path(__file__).resolve().parent/"cache/conditions"
MW = CONDITIONS/"mw_inputs"  # md_mw.py inputs and outputs (run on solaris)


def conditions_task(run, kind, frame):
    outputs = MW if "_mW_" in run else CONDITIONS
    info = json.loads((outputs/f"{run}.json").read_text())
    structure = info["structure"]
    oxygen_only = info["model"] == "mW"
    folder = MW if oxygen_only else CONDITIONS if structure.startswith("E_") else PREPS
    u = mda.Universe(str(folder/structure))
    frames = np.load(outputs/f"{run}.npz")
    u.atoms.positions = frames["x" if kind == "thermal" else "xq"][frame]
    u.dimensions = triclinic_box(*frames["cell"])
    columns = evaluate(u, "name OW", None if oxygen_only else "name HW1 HW2")
    n = len(next(iter(columns.values())))
    stack_file = folder/structure.replace(".gro", ".stack.npy")
    stack = np.load(stack_file) if stack_file.exists() else np.zeros(n, dtype=np.int8)
    phase = info["phase"]
    system = "liquid" if phase == "liquid" else f"bulk_{phase}"
    meta = dict(dataset=info["dataset"], system=system, kind=kind, frame=frame, mol=np.arange(n), phase=phase,
                stack=stack, temperature=info["temperature"], pressure=info["pressure_gpa"], model=info["model"],
                run=f"p{info['prep']}")
    return pd.DataFrame({**meta, **columns})


def conditions_runs(prefix):
    runs = [p.stem for p in CONDITIONS.glob(f"{prefix}_*.npz")]
    if prefix == "D":
        runs += [p.stem for p in MW.glob("D_mW_*.npz")]
    return sorted(runs)


REPLICA_SYSTEMS = ("bulk_Ih", "bulk_Ic", "bulk_Isd", "liquid", "bulk_Ih_250K", "bulk_Ic_250K")
REPLICA_FRAMES = 60


def seedtest_systems():
    systems = [f"bulk_{p}" for p in ("Ih", "Ic", "Isd")]
    systems += [f"seed_{p}_{n}" for p in ("Ih", "Ic", "Isd") for n in (50, 100, 200)]
    return systems


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", default="genice,seedtest,replica,prep,model,pressure")
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    manifest = dict(solorder=SolOrder.__version__,
                    registry_sha256=sha256(Path(__file__).with_name("registry.py")),
                    jitter_angstrom=JITTER, genice_seeds=SEEDS, inputs={})
    previous = CACHE/"manifest.json"
    if previous.exists():  # keep hashes of datasets not recomputed in this run
        manifest["inputs"] = json.loads(previous.read_text()).get("inputs", {})
    for dataset in args.datasets.split(","):
        start = time.time()
        if dataset == "genice":
            phases = [p for p in json.loads((GENICE/f"seed{SEEDS[0]}/manifest.json").read_text())["structures"]
                      if p not in DUPLICATES]
            tasks = [delayed(genice_task)(p, s) for p in phases for s in SEEDS]
            files = [GENICE/f"seed{s}/{p}.gro" for s in SEEDS for p in phases]
        elif dataset == "seedtest":
            tasks = [delayed(seedtest_task)(s, k, f) for s in seedtest_systems()
                     for k in ("thermal", "quenched") for f in range(10)]
            files = [SEEDTEST/"truth.npz"]
            for s in seedtest_systems():
                files += [SEEDTEST/f"{s}.gro", SEEDTEST/f"frames_{s}.npz", SEEDTEST/f"framesq_{s}.npz"]
        elif dataset == "prep":
            runs = prep_runs()
            tasks = [delayed(prep_task)(r, k, f) for r in runs for k in ("thermal", "quenched")
                     for f in range(len(np.load(PREPS/f"{r}.npz")["energy"]))]
            files = [PREPS/f"{r}.npz" for r in runs] + sorted(PREPS.glob("*_p?.gro"))
        elif dataset in ("model", "pressure"):
            runs = conditions_runs("D" if dataset == "model" else "E")
            tasks = [delayed(conditions_task)(r, k, f) for r in runs for k in ("thermal", "quenched")
                     for f in range(len(np.load((MW if "_mW_" in r else CONDITIONS)/f"{r}.npz")["energy"]))]
            files = [(MW if "_mW_" in r else CONDITIONS)/f"{r}.npz" for r in runs]
            files += sorted(CONDITIONS.glob("E_*_p?.gro"))
        elif dataset == "replica":
            tasks = [delayed(replica_task)(s, k, f) for s in REPLICA_SYSTEMS
                     for k in ("thermal", "quenched") for f in range(REPLICA_FRAMES)]
            files = [SEEDTEST/"truth.npz"]
            for s in REPLICA_SYSTEMS:
                files += [SEEDTEST/f"{s.removesuffix('_250K')}.gro", SEEDTEST/f"train_{s}.npz"]
        else:
            raise ValueError(f"Unknown dataset {dataset!r}")
        tables = Parallel(n_jobs=args.workers)(tasks)
        table = pd.concat(tables, ignore_index=True)
        table.to_pickle(CACHE/f"{dataset}.pkl")
        manifest["inputs"][dataset] = {str(f.relative_to(f.parents[1])): sha256(f) for f in files}
        print(f"{dataset}: {len(table)} rows, {len(tasks)} frames, {time.time()-start:.0f} s")
    (CACHE/"manifest.json").write_text(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
