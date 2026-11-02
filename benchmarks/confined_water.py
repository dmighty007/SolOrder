"""Application: detecting confined water in T4 lysozyme L99A (4W52) from water-only local order parameters.

System: T4L L99A with benzene (OPLS-AA, 15,735 SPC/E waters, 300 K), 89
PathGennie benzene-egress paths (benchmarks/cache/t4l, copied from the
campaign's Correct_Runs; see its README). Two ensembles: "dry" (replicas r1,
r4, 48 paths; no water in the hydrophobic cavity) and "wet" (r2, r3, 41
paths; the same three waters, residues 5250, 5264 and 6060, sit beside
benzene in the cavity from the start).

Labels (geometric, independent of every order parameter), per water and frame:
  cavity3   one of the three cavity waters, while it is buried (B >= BURIED)
  buried    any other water with B >= BURIED, B = protein heavy atoms within 8 A
  surface   a water within 4 A of a protein heavy atom with B < BURIED (sample)
  bulk      a water with no protein heavy atom within 12 A (sample)
  near_lig  a water within 6 A of a benzene carbon
Order parameters: as in host_guest.py, from water oxygens and their hydrogens
only (protein, benzene and ions are not neighbors).

Analysis: held-out scores with thresholds fitted on one pair of replicas and
applied to the other (r1+r2 versus r3+r4, worse direction): cavity3 and
buried against surface and against bulk. Frozen thresholds (fitted on
buried versus surface water of one replica pair) are then applied near the
cavity in the other pair, to count waters a rule would call confined in the
dry and the wet ensemble, and cavity3 parameters
are followed against benzene's distance from the cavity.

Step export writes cache/applications/t4l_features.npz (archived; analyze
reads it when cache/t4l_ops.pkl is absent). The replica split is the column
'replica'.

  python confined_water.py compute | analyze | export
"""

import argparse
import json
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

from joblib import Parallel, delayed
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, distance_array
import numpy as np
import pandas as pd

from host_guest import OP_COLUMNS, export_features, heldout, load_features, oracle

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
DATA = CACHE/"t4l"
CAVITY_WATERS = (5250, 5264, 6060)
BURIED = 60
FRAMES, SAMPLE = 60, 200


def water_ops(u, frame):
    from SolOrder import OrderParameter, eFingerprint
    op = OrderParameter(u, frame=frame, cutoff=5., selection="resname SOL and name OW",
                        hydrogen_selection="resname SOL and (name HW1 or name HW2)")
    ent = eFingerprint(u, frame=frame, cutoff=5., selection="resname SOL and name OW")
    rings = op.ring_statistics(10, network="hbond", cutoff=3.5).fractions
    return dict(Q=op.tetrahedral_order(), Q_avg=op.local_average(op.tetrahedral_order(), cutoff=3.5),
                lsi=op.local_structure_index(3.7), d5=op.d5(), cn=op.coordination_number(3.5),
                cn_smooth=op.smooth_coordination(r0=3.5, n=12, m=28), zeta=op.zeta(3.5, 30.),
                n_hb=op.hydrogen_bonds(3.5, 30.).total.astype(float),
                q6_avg=op.steinhardt(6, averaged=True, cutoff=3.5), q3=op.steinhardt(3, cutoff=3.5),
                q3_avg=op.steinhardt(3, averaged=True, cutoff=3.5), f4=op.F4(cutoff=3.5),
                s=ent.Entropy(), s_avg=ent.localEntropy(), vor_density=op.voronoi_density(),
                vor_asphericity=op.voronoi_asphericity(), ring_f5=rings[:, 2], ring_f6=rings[:, 3],
                ring_f7=rings[:, 4])


def path_task(row):
    u = mda.Universe(str(DATA/"system/pbcmol.gro"), str(DATA/row.dcd))
    ow = u.select_atoms("resname SOL and name OW")
    prot = u.select_atoms("protein and not name H*")
    lig = u.select_atoms("resname BNZ and not name H*")
    resid = ow.resids
    cav = np.isin(resid, CAVITY_WATERS)
    n = len(u.trajectory)
    frames = np.unique(np.linspace(0, n-1, min(n, FRAMES)).round().astype(int))
    rng = np.random.default_rng([2028, sum(map(ord, row.case))])
    # cavity centre: benzene centre in the first frame (it starts in the L99A cavity)
    u.trajectory[0]
    centre = lig.center_of_geometry()
    rows = []
    for f in frames:
        u.trajectory[f]
        box = u.dimensions
        pairs = capped_distance(ow.positions, prot.positions, 8., box=box, return_distances=True)
        burial = np.bincount(pairs[0][:, 0], minlength=ow.n_atoms)
        contact = np.zeros(ow.n_atoms, bool)
        contact[pairs[0][pairs[1] < 4., 0]] = True
        near12 = np.zeros(ow.n_atoms, bool)
        near12[capped_distance(ow.positions, prot.positions, 12., box=box, return_distances=False)[:, 0]] = True
        d_lig = distance_array(ow.positions, lig.positions, box=box).min(axis=1)
        d_centre = distance_array(ow.positions, centre[None], box=box)[:, 0]
        lig_shift = float(np.linalg.norm(lig.center_of_geometry() - centre))
        region = np.full(ow.n_atoms, "", dtype=object)
        surf = np.flatnonzero(contact & (burial < BURIED))
        region[rng.choice(surf, size=min(SAMPLE, len(surf)), replace=False)] = "surface"
        bulk = np.flatnonzero(~near12)
        region[rng.choice(bulk, size=min(SAMPLE, len(bulk)), replace=False)] = "bulk"
        region[d_lig < 6.] = "near_lig"
        region[burial >= BURIED] = "buried"
        region[cav & (burial >= BURIED)] = "cavity3"
        keep = np.flatnonzero(region != "")
        ops = water_ops(u, f)
        rows.append(pd.DataFrame(dict(case=row.case, state=row.state, replica=row.replica, frame=f,
                                      lig_shift=lig_shift, water=resid[keep], region=region[keep],
                                      burial=burial[keep], d_lig=d_lig[keep].astype(np.float32),
                                      d_centre=d_centre[keep].astype(np.float32),
                                      **{c: np.asarray(ops[c])[keep].astype(np.float32) for c in OP_COLUMNS})))
    return pd.concat(rows, ignore_index=True)


def compute(workers):
    index = pd.read_csv(DATA/"mapping/path_index.csv")
    tables = Parallel(n_jobs=workers)(delayed(path_task)(r) for r in index.itertuples())
    table = pd.concat(tables, ignore_index=True)
    table.to_pickle(CACHE/"t4l_ops.pkl")
    print(f"t4l_ops: {len(table)} rows")
    print(table.groupby(["region", "state"]).size().unstack().to_string())


def analyze(workers):
    from metrics import fit_threshold, hits
    t = load_features("t4l")
    t["fold"] = t.replica.isin(("r3", "r4")).astype(int)  # r1+r2 (one dry, one wet) versus r3+r4
    rows = []
    for target in ("cavity3", "buried", "near_lig"):
        for reference in ("surface", "bulk"):
            a_part, b_part = t[t.region == target], t[t.region == reference]
            if target == "near_lig":
                a_part = a_part[a_part.state == "wet"]
            for c in OP_COLUMNS:
                a, b = a_part[c].to_numpy(float), b_part[c].to_numpy(float)
                rows.append(dict(target=target, reference=reference, column=c, n=len(a_part),
                                 oracle=oracle(a, b), heldout=heldout(a, b, a_part.fold.to_numpy(),
                                                                       b_part.fold.to_numpy()),
                                 median=float(np.nanmedian(a)), reference_median=float(np.nanmedian(b))))
    scores = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    scores.to_csv(RESULTS/"t4l_scores.csv", index=False)
    # frozen rules: fitted on buried versus surface water of one replica pair, applied near the
    # cavity in the other pair (each pair has one dry and one wet replica); mean of both directions
    frozen = []
    for c in OP_COLUMNS:
        directions = []
        for train in (0, 1):
            fit, test = t[t.fold == train], t[t.fold != train]
            th, below = fit_threshold(fit.loc[fit.region == "buried", c].to_numpy(float),
                                      fit.loc[fit.region == "surface", c].to_numpy(float))
            if not np.isfinite(th):
                continue
            near = test[test.d_centre < 8.]
            flag = hits(near[c].to_numpy(float), th, below)
            for state in ("dry", "wet"):
                m = (near.state == state).to_numpy()
                frames = near[m].groupby(["case", "frame"]).ngroups
                directions.append(dict(state=state, flagged=float(flag[m].sum()/max(frames, 1)),
                                       recall=float(hits(test.loc[test.region == "cavity3", c].to_numpy(float),
                                                         th, below).mean())))
        if not directions:
            continue
        d = pd.DataFrame(directions)
        for state, g in d.groupby("state"):
            frozen.append(dict(column=c, state=state, flagged_per_frame=float(g.flagged.mean()),
                               flagged_min=float(g.flagged.min()), flagged_max=float(g.flagged.max()),
                               cavity3_recall=float(g.recall.min())))
    pd.DataFrame(frozen).to_csv(RESULTS/"t4l_frozen.csv", index=False)
    # cavity waters as benzene leaves
    cav = t[t.region == "cavity3"].copy()
    cav["stage"] = pd.cut(cav.lig_shift, [0, 3, 6, 10, 100], labels=["<3", "3-6", "6-10", ">10"])
    track = cav.groupby("stage", observed=True)[["cn_smooth", "vor_density", "n_hb", "q6_avg", "Q", "burial"]].median()
    track["n"] = cav.groupby("stage", observed=True).size()
    track.to_csv(RESULTS/"t4l_cavity_track.csv")
    pd.set_option("display.width", 220)
    print(scores.pivot_table(index="column", columns=["target", "reference"], values="heldout").round(3).to_string())
    print(pd.DataFrame(frozen).pivot(index="column", columns="state", values="flagged_per_frame").round(2)
          .join(pd.DataFrame(frozen).groupby("column").cavity3_recall.first().round(3)).to_string())
    print(track.round(3).to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=("compute", "analyze", "export"))
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    {"compute": compute, "analyze": analyze,
     "export": lambda workers: export_features("t4l")}[args.step](args.workers)


if __name__ == "__main__":
    main()
