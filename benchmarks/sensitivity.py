"""Sensitivity of separability to the neighbor convention (run C subset).

Registry defaults (registry.py) stay the canonical benchmark. Here each factor
is varied on its own, everything else at its default:

  shell      first-shell O-O cutoff for every cutoff-based column (CN, smooth
             CN r0, angular, H-bond, ring and bond-orientational columns,
             averaging of Q, and the CHILL/CHILL+ labels): 3.2-3.8 A
  lsi        LSI shell radius: 3.2-3.8 A (default 3.7)
  entropy    pair-entropy integration radius r_m, which is also its averaging
             radius in eFingerprint: 4.0-6.0 A (default 5.0)
  average    averaging radius of q3/q4/q6 and Q with q_lm and Q computed at
             3.5 A: 3.2-5.0 A (default 3.5)

Rank-based columns (d5, Q, S_k, centrosymmetry) and Voronoi columns have no
cutoff and are not scanned.

Data: preparations C1-C3 at 220 K, bulk Ih, bulk Ic and liquid, twelve thermal
frames 25 ps apart and the matching quenched frames. Scores per variant:
oracle (pooled C1-C3); held-out, threshold fitted on C1 and applied to C2 and
C3 (worse of the two); and cross-convention, threshold fitted on C1 at the
default setting and applied to C2 and C3 computed with the variant.

  python sensitivity.py compute   # cache/sensitivity.pkl
  python sensitivity.py metrics   # results/sensitivity.csv
"""

import argparse
from pathlib import Path
import sys
import time
import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

from joblib import Parallel, delayed
import MDAnalysis as mda
import numpy as np
import pandas as pd

from SolOrder import OrderParameter, eFingerprint

from compute import PREPS
from metrics import fit_threshold, hits, ks
from registry import FIRST_SHELL, RING_SIZES

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
SHELLS = (3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8)
ENTROPY_RADII = (4.0, 4.5, 5.0, 5.5, 6.0)
AVERAGE_RADII = (3.2, 3.5, 4.0, 4.5, 5.0)
DEFAULTS = dict(shell=FIRST_SHELL, lsi=3.7, entropy=5.0, average=FIRST_SHELL)
PREPARATIONS = ("p1", "p2", "p3")
PHASES = ("Ih", "Ic", "liquid")
FRAMES = tuple(range(0, 60, 5))
TASKS = (("Ih", "liquid"), ("Ic", "liquid"), ("Ih", "Ic"))


def shell_columns(op, r):
    """Every cutoff-based registry column at first-shell cutoff r."""
    out = {"cn": op.coordination_number(r),
           "cn_smooth": op.smooth_coordination(r0=r, n=12, m=28),
           "q_tet_avg": op.local_average(op.tetrahedral_order(), cutoff=r),
           "bond_angle": op.bond_angle_parameter(1, 1, cutoff=r),
           "afs3": op.angular_fourier(a=3, radial=0, cutoff=r),
           "f4": op.F4(cutoff=r),
           "n_hb": op.hydrogen_bonds(r, 30.).total.astype(float),
           "zeta": op.zeta(r, 30.),
           "xconn6": op.crystalline_connections(6, .5, cutoff=r).astype(float),
           "ndp": op.neighbor_distance_parameter(cutoff=r)}
    ah = op.angle_histogram_parameter((1, 2, 3), cutoff=r)
    out.update({f"ah{k}": ah[:, i] for i, k in enumerate((1, 2, 3))})
    a = op.neighborhood_parameters(m=1, cutoff=r)
    out.update({f"A{k}": a[:, k-1] for k in (1, 2, 3)})
    for l in (3, 4, 6):
        out[f"q{l}"] = op.steinhardt(l, cutoff=r)
        out[f"q{l}_avg"] = op.steinhardt(l, averaged=True, cutoff=r)
    for l in (4, 6):
        out[f"w{l}_hat"] = op.steinhardt(l, invariant="w", normalized=True, cutoff=r)
    rings = op.ring_statistics(10, network="hbond", cutoff=r)
    out["ring_mean"] = rings.mean_size
    out.update({f"ring_f{k}": rings.fractions[:, k-3] for k in RING_SIZES})
    out["chill"] = op.chill(cutoff=r).labels.astype(float)
    out["chill_plus"] = op.chill_plus(r).labels.astype(float)
    return out


def average_columns(op, radius):
    out = {f"q{l}_avg": op.local_average(op._qlm(l, FIRST_SHELL), cutoff=radius) for l in (3, 4, 6)}
    for l in (3, 4, 6):
        out[f"q{l}_avg"] = np.sqrt(4.*np.pi/(2*l+1)*np.sum(np.abs(out[f"q{l}_avg"])**2, axis=1))
    out["q_tet_avg"] = op.local_average(op.tetrahedral_order(), cutoff=radius)
    return out


def frame_task(phase, prep, kind, frame):
    u = mda.Universe(str(PREPS/f"{phase}_{prep}.gro"))
    frames = np.load(PREPS/f"{phase}_{prep}_220.npz")
    u.atoms.positions = frames["x" if kind == "thermal" else "xq"][frame]
    u.dimensions = [*frames["box"], 90., 90., 90.]
    op = OrderParameter(u, cutoff=6., selection="name OW", hydrogen_selection="name HW1 HW2")
    out = {}
    for r in SHELLS:
        out.update({f"{c}|shell|{r}": v for c, v in shell_columns(op, r).items()})
        out[f"lsi|lsi|{r}"] = op.local_structure_index(r)
    for r in AVERAGE_RADII:
        out.update({f"{c}|average|{r}": v for c, v in average_columns(op, r).items()})
    for r in ENTROPY_RADII:
        ent = eFingerprint(u, cutoff=r, selection="name OW")
        out[f"entropy|entropy|{r}"] = ent.Entropy()
        out[f"entropy_avg|entropy|{r}"] = ent.localEntropy()
    n = len(next(iter(out.values())))
    meta = dict(phase=phase, run=prep, kind=kind, frame=frame, mol=np.arange(n))
    return pd.DataFrame({**meta, **{k: np.asarray(v, dtype=np.float32) for k, v in out.items()}})


def compute(workers):
    start = time.time()
    tasks = [delayed(frame_task)(phase, prep, kind, frame) for phase in PHASES for prep in PREPARATIONS
             for kind in ("thermal", "quenched") for frame in FRAMES]
    table = pd.concat(Parallel(n_jobs=workers)(tasks), ignore_index=True)
    table.to_pickle(CACHE/"sensitivity.pkl")
    print(f"sensitivity: {len(table)} rows, {len(tasks)} frames, {time.time()-start:.0f} s")


def ice_fraction(labels, phase):
    """Labeler score as a balanced-accuracy analogue: correct polytype (ice) or not ice (liquid)."""
    code = {"Ih": 1, "Ic": 2}.get(phase)
    return np.mean(labels == code) if code else np.mean(~np.isin(labels, (1, 2)))


def metrics():
    table = pd.read_pickle(CACHE/"sensitivity.pkl")
    variants = [c for c in table.columns if "|" in c]
    rows = []
    for kind in ("thermal", "quenched"):
        part = table[table.kind == kind]
        masks = {(p, r): ((part.phase == p) & (part.run == r)).to_numpy() for p in PHASES for r in PREPARATIONS}
        pooled = {p: (part.phase == p).to_numpy() for p in PHASES}
        for variant in variants:
            column, factor, value = variant.split("|")
            v = part[variant].to_numpy(float)
            default = part[f"{column}|{factor}|{DEFAULTS[factor]}"].to_numpy(float)
            if column in ("chill", "chill_plus"):
                for phase in PHASES:
                    rows.append(dict(kind=kind, column=column, factor=factor, value=float(value),
                                     a=phase, b="", label_score=ice_fraction(v[pooled[phase]], phase)))
                continue
            for a, b in TASKS:
                row = dict(kind=kind, column=column, factor=factor, value=float(value), a=a, b=b,
                           oracle=(1 + ks(v[pooled[a]], v[pooled[b]]))/2,
                           nan_fraction=float(np.mean(~np.isfinite(v[pooled[a] | pooled[b]]))))
                for name, train in (("heldout", v), ("cross", default)):
                    t, below = fit_threshold(train[masks[a, "p1"]], train[masks[b, "p1"]])
                    scores = [.5*(hits(v[masks[a, r]], t, below).mean() + hits(v[masks[b, r]], t, not below).mean())
                              for r in PREPARATIONS[1:]] if np.isfinite(t) else [np.nan]
                    row[name] = min(scores)
                rows.append(row)
    out = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    out.to_csv(RESULTS/"sensitivity.csv", index=False)
    print(f"sensitivity.csv: {len(out)} rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=("compute", "metrics"))
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    compute(args.workers) if args.step == "compute" else metrics()
