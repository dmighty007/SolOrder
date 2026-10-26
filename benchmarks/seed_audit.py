"""Seed-label audit (W2 in docs/research_workflow.md).

Construction labels say which molecules were built as seed and with which
stacking. This script checks what is still true in each analyzed frame and
where the labelers fail, without using any benchmarked order parameter:

  survival   per seed molecule: displacement from its constructed site (after
             removing the seed's mean drift) and how many of its constructed
             lattice bonds are still shorter than 3.5 A. "Intact" means
             displacement < 1 A and every constructed bond kept.
  depth      R - r: distance inside the constructed seed radius, measured from
             the current seed centroid.
  labels     CHILL, CHILL+ and PTM by depth, split into correct polytype, wrong
             polytype, CHILL/CHILL+ interfacial ice, and other.
  PTM RMSD   best diamond-template RMSD without a cutoff, by depth (NaN where
             no diamond template fits at all).
  CHILL bonds in bulk ice, where each 4-coordinated molecule has its bond
             correlations: does the least negative one (the eclipsed bond in
             a hexagonal layer) fall in CHILL's window (-0.2, -0.05) or only
             in CHILL+'s (-0.35, 0.25)?

Reads cache/seedtest.pkl (compute.py) and the seed-test frames; writes
results/seed_audit_*.csv.
"""

from pathlib import Path
import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

from joblib import Parallel, delayed
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from SolOrder import OrderParameter

from compute import SEEDTEST, seedtest_systems

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
LATTICE_BOND = 2.76*1.1  # A; constructed O-O bonds (build.py)
KEPT_BOND = 3.5  # A
INTACT_DISPLACEMENT = 1.  # A, about three times the thermal mean in bulk ice at 220 K
DEPTH_BINS = np.arange(-2., 13., 1.)
CHILL_WINDOW, CHILL_PLUS_WINDOW = (-.2, -.05), (-.35, .25)


def oxygen_sites(path):
    lines = Path(path).read_text().splitlines()
    xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in lines[2:-1]])*10
    return xyz[::4]


def survival(system, kind, frame):
    truth = np.load(SEEDTEST/"truth.npz")
    seed = truth[f"{system}__seed"].astype(bool)
    built = oxygen_sites(SEEDTEST/f"{system}.gro")[seed]
    frames = np.load(SEEDTEST/f"{'frames' if kind == 'thermal' else 'framesq'}_{system}.npz")
    box = np.r_[frames["box"], 90., 90., 90.]
    now = frames["x"][frame][::4][seed]
    displacement = minimize_vectors(now-built, box)
    displacement -= displacement.mean(axis=0)
    pairs = cKDTree(built).query_pairs(LATTICE_BOND, output_type="ndarray")  # seed is whole and centred
    kept = np.linalg.norm(minimize_vectors(now[pairs[:, 0]]-now[pairs[:, 1]], box), axis=1) < KEPT_BOND
    bonds = np.bincount(pairs.ravel(), minlength=len(built))
    lost = np.bincount(pairs[~kept].ravel(), minlength=len(built))
    out = pd.DataFrame(dict(system=system, kind=kind, frame=frame, mol=np.flatnonzero(seed),
                            displacement=np.linalg.norm(displacement, axis=1), bonds=bonds, bonds_lost=lost))
    out["intact"] = (out.displacement < INTACT_DISPLACEMENT) & (out.bonds_lost == 0)
    return out


def universe(system, kind, frame):
    u = mda.Universe(str(SEEDTEST/f"{system}.gro"))
    frames = np.load(SEEDTEST/f"{'frames' if kind == 'thermal' else 'framesq'}_{system}.npz")
    u.atoms.positions = frames["x"][frame]
    u.dimensions = [*frames["box"], 90., 90., 90.]
    return u


def ptm_rmsd(system, kind, frame):
    op = OrderParameter(universe(system, kind, frame), cutoff=5., selection="name OW")
    match = op.template_matching(rmsd_cutoff=10.)
    rmsd = np.where(match.labels == 0, np.nan, match.rmsd)  # OVITO reports 0 when no template fits
    return pd.DataFrame(dict(system=system, kind=kind, frame=frame, mol=np.arange(op.n), ptm_rmsd=rmsd))


def chill_bonds(system, kind, frame):
    """Least negative bond correlation of each 4-coordinated molecule, CHILL neighbor rule."""
    op = OrderParameter(universe(system, kind, frame), cutoff=5., selection="name OW")
    result = op.chill(cutoff=3.5)
    rows = []
    for i in range(op.n):
        c = result.correlations[result.offsets[i]:result.offsets[i+1]]
        if len(c) == 4:
            c = np.sort(c)
            rows.append((i, c[-1], int(np.sum(c < -.8))))
    rows = np.array(rows)
    return pd.DataFrame(dict(system=system, kind=kind, frame=frame, mol=rows[:, 0].astype(int),
                             least_negative=rows[:, 1], staggered=rows[:, 2].astype(int)))


def label_classes(labels, stack, method):
    correct = np.isin(labels, (1, 2)) & (labels == stack)
    wrong = np.isin(labels, (1, 2)) & (labels != stack)
    interfacial = labels == 3 if method != "ptm" else np.zeros(len(labels), bool)
    return dict(correct=correct, wrong_polytype=wrong, interfacial=interfacial,
                other=~(correct | wrong | interfacial))


def depth_table(seeds):
    rows = []
    seeds = seeds.assign(depth_bin=pd.cut(seeds.depth, DEPTH_BINS))
    for (kind, phase, system, depth_bin), part in seeds.groupby(["kind", "phase", "system", "depth_bin"],
                                                              observed=True):
        base = dict(kind=kind, phase=phase, system=system, depth=depth_bin.mid, n=len(part),
                    intact=part.intact.mean(), ptm_matched=part.ptm_rmsd.notna().mean(),
                    ptm_rmsd_median=part.ptm_rmsd.median(), ptm_below_cutoff=(part.ptm_rmsd < .1).mean())
        for method in ("chill", "chill_plus", "ptm"):
            classes = label_classes(part[method].to_numpy(), part["stack"].to_numpy(), method)
            for name, hit in classes.items():
                rows.append(dict(base, method=method, outcome=name, fraction=hit.mean(),
                                 fraction_intact=hit[part.intact.to_numpy()].mean() if part.intact.any() else np.nan))
    return pd.DataFrame(rows)


def system_table(seeds):
    rows = []
    for (kind, system), part in seeds.groupby(["kind", "system"], sort=False):
        row = dict(kind=kind, system=system, n_seed=part.groupby("frame").size().iloc[0],
                   R=part.R.iloc[0], intact=part.intact.mean(), intact_interior=part.intact[part.interior].mean(),
                   displacement_median=part.displacement.median(), bonds_lost_mean=part.bonds_lost.mean(),
                   deep_fraction=(part.depth > 5.).mean())
        for method in ("chill", "chill_plus", "ptm"):
            classes = label_classes(part[method].to_numpy(), part["stack"].to_numpy(), method)
            intact = part.intact.to_numpy()
            deep = (part.depth > 5.).to_numpy()
            row[f"{method}_correct"] = classes["correct"].mean()
            row[f"{method}_correct_intact"] = classes["correct"][intact].mean()
            row[f"{method}_correct_deep"] = classes["correct"][deep].mean() if deep.any() else np.nan
            row[f"{method}_interfacial"] = classes["interfacial"].mean()
            row[f"{method}_wrong_polytype"] = classes["wrong_polytype"].mean()
        rows.append(row)
    return pd.DataFrame(rows)


def bond_table(bonds, table):
    rows = []
    merged = bonds.merge(table[["system", "kind", "frame", "mol", "stack"]], on=["system", "kind", "frame", "mol"])
    for (kind, system, stack), part in merged.groupby(["kind", "system", "stack"]):
        c = part.least_negative.to_numpy()
        rows.append(dict(kind=kind, system=system, layer={1: "hexagonal", 2: "cubic"}[stack], n=len(part),
                         three_staggered=(part.staggered == 3).mean(),
                         four_staggered=(part.staggered == 4).mean(),
                         in_chill_window=np.mean((c > CHILL_WINDOW[0]) & (c < CHILL_WINDOW[1])),
                         in_chill_plus_window=np.mean((c > CHILL_PLUS_WINDOW[0]) & (c < CHILL_PLUS_WINDOW[1])),
                         least_negative_median=np.median(c),
                         least_negative_p05=np.percentile(c, 5), least_negative_p95=np.percentile(c, 95)))
    return pd.DataFrame(rows)


def main():
    table = pd.read_pickle(CACHE/"seedtest.pkl")
    seeded = [s for s in seedtest_systems() if s.startswith("seed")]
    bulk = [s for s in seedtest_systems() if s.startswith("bulk")]
    tasks = [(s, k, f) for s in seeded for k in ("thermal", "quenched") for f in range(10)]
    with Parallel(n_jobs=10) as parallel:
        alive = pd.concat(parallel(delayed(survival)(*t) for t in tasks))
        rmsd = pd.concat(parallel(delayed(ptm_rmsd)(*t) for t in tasks))
        bonds = pd.concat(parallel(delayed(chill_bonds)(s, k, f) for s in bulk
                                   for k in ("thermal", "quenched") for f in range(10)))
    keys = ["system", "kind", "frame", "mol"]
    seeds = (table[table.seed][keys + ["phase", "stack", "interior", "r_seed", "R", "chill", "chill_plus", "ptm"]]
             .merge(alive, on=keys).merge(rmsd, on=keys))
    seeds["depth"] = seeds.R - seeds.r_seed
    RESULTS.mkdir(exist_ok=True)
    depth_table(seeds).to_csv(RESULTS/"seed_audit_depth.csv", index=False)
    system_table(seeds).to_csv(RESULTS/"seed_audit_systems.csv", index=False)
    bond_table(bonds, table).to_csv(RESULTS/"seed_audit_chill_bonds.csv", index=False)
    print("seed audit done")


if __name__ == "__main__":
    main()
