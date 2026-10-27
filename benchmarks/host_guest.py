"""Application: key water molecules in OAMe-G2 unbinding paths, from water-only local order parameters.

System: octa-acid (OAMe, residue OCB) with guest G2 (residue MOL), 9 Na+ and
2100 TIP3P waters, 300 K; 208 PathGennie escape paths (GROMACS, accepted
states saved in order, no time axis). Data: benchmarks/cache/oame/ (box.gro,
config/conf_template.pdb, trajectory_<k>.xtc).

Frame of reference. Each frame is rotated onto the host template that the
path-sampling bias used (PLUMED FIT_TO_TEMPLATE on 40 host atoms); the cavity
axis is then x = y = 20.136 A, and z runs from the bottom of the bowl (about
14 A) through the portal (host atoms end at about 28 A); the bound guest has
its centre of mass at z = 25.7 A.

Definitions (geometric, independent of every order parameter):
  cavity water    oxygen with 16 < z < 26 A and within 4.5 A of the axis
  bridging water  oxygen within 3.5 A of a guest polar atom (O1, O2, N1) and of
                  a host oxygen, with the guest not yet out of the portal
  portal frame    first frame with the guest centre of mass at z > PORTAL_Z
  dry / wet path  mean number of cavity waters over the PRE frames before the
                  portal frame: 0 to DRY_MAX dry, at least WET_MIN wet
  completed path  guest centre of mass beyond z = 35 A at the end of the path
  entering water  a water that becomes a cavity water while the guest is
                  still in the cavity or crossing the portal

Steps:
  scan      geometry of every frame of every path -> cache/oame_scan.pkl
  compute   water-only order parameters (SolOrder, water oxygens and their
            hydrogens; host, guest and ions are not neighbors) on every
            frame within WINDOW frames of the portal frame and on up to 40
            further frames per path -> cache/oame_ops.pkl
  analyze   oracle and held-out scores of every parameter for the three kinds
            of key water against bulk water, dry and wet paths separately, and
            cavity water against other water near the portal, and the
            prediction of entering waters from their parameters before
            entry -> results/oame_*.csv
  export    the per-water table as cache/applications/oame_features.npz
            (archived; analyze reads it when cache/oame_ops.pkl is absent)

  python host_guest.py scan | compute | analyze | export
"""

import argparse
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

from joblib import Parallel, delayed
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import capped_distance, distance_array
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
DATA = CACHE/"oame"
AXIS = np.array([20.136, 20.136])
CAVITY_Z, CAVITY_R = (16., 26.), 4.5
PORTAL_Z, OUT_Z = 28., 35.
BRIDGE = 3.5
PRE, DRY_MAX, WET_MIN = 10, .5, 2.
WINDOW, EXTRA = 30, 40
GUEST_POLAR = ("O1", "O2", "N1")


def paths():
    return sorted((p.stem for p in DATA.glob("trajectory_*.xtc")), key=lambda s: int(s.split("_")[1]))


class Geometry:
    """Rotation of each frame onto the host template, and the groups used by the definitions."""

    def __init__(self, path):
        self.u = mda.Universe(str(DATA/"box.gro"), str(DATA/f"{path}.xtc"))
        template = mda.Universe(str(DATA/"config/conf_template.pdb"))
        self.fit = self.u.atoms[template.atoms.ids - 1]
        assert (self.fit.names == template.atoms.names).all()
        self.ref = template.atoms.positions.astype(float)
        self.ref_centre = self.ref.mean(axis=0)
        self.guest = self.u.select_atoms("resname MOL")
        self.polar = self.guest.select_atoms("name " + " ".join(GUEST_POLAR))
        self.host_o = self.u.select_atoms("resname OCB and name O*")
        self.water_o = self.u.select_atoms("resname WAT and name O")

    def frame(self):
        box = self.u.dimensions
        x = self.fit.positions.astype(float)
        x = x[0] + (x - x[0]) - box[:3]*np.round((x - x[0])/box[:3])  # whole around its first atom
        centre = x.mean(axis=0)
        R, _ = rotation_matrix(x - centre, self.ref - self.ref_centre)

        def to_template(p):
            d = np.asarray(p, float) - centre
            d -= box[:3]*np.round(d/box[:3])
            return d @ R.T + self.ref_centre
        w = to_template(self.water_o.positions)
        g = to_template(self.guest.center_of_mass()[None])[0]
        r = np.linalg.norm(w[:, :2] - AXIS, axis=1)
        cavity = (w[:, 2] > CAVITY_Z[0]) & (w[:, 2] < CAVITY_Z[1]) & (r < CAVITY_R)
        near_guest = distance_array(self.water_o.positions, self.polar.positions, box=box).min(axis=1) < BRIDGE
        near_host = distance_array(self.water_o.positions, self.host_o.positions, box=box).min(axis=1) < BRIDGE
        return dict(guest_z=g[2], guest_r=float(np.linalg.norm(g[:2] - AXIS)), cavity=np.flatnonzero(cavity),
                    bridge=np.flatnonzero(near_guest & near_host & (g[2] < PORTAL_Z + 2.)),
                    water_z=w[:, 2].astype(np.float32), water_r=r.astype(np.float32))


def scan_path(path):
    geo = Geometry(path)
    rows = []
    for ts in geo.u.trajectory:
        f = geo.frame()
        rows.append(dict(path=path, frame=ts.frame, guest_z=f["guest_z"], guest_r=f["guest_r"],
                         n_cavity=len(f["cavity"]), cavity=f["cavity"], bridge=f["bridge"]))
    return pd.DataFrame(rows)


def classify(scan):
    out = []
    for path, p in scan.groupby("path", sort=False):
        z = p.guest_z.to_numpy()
        crossed = np.flatnonzero(z > PORTAL_Z)
        portal = int(crossed[0]) if len(crossed) else -1
        before = p.n_cavity.to_numpy()[max(portal - PRE, 0):portal] if portal > 0 else np.array([np.nan])
        occupancy = float(np.mean(before)) if len(before) else np.nan
        label = "dry" if occupancy <= DRY_MAX else "wet" if occupancy >= WET_MIN else "intermediate"
        out.append(dict(path=path, frames=len(p), portal_frame=portal, completed=bool(z[-1] > OUT_Z),
                        cavity_before_portal=occupancy, max_cavity=int(p.n_cavity.max()),
                        label=label if portal > 0 else "no_exit"))
    return pd.DataFrame(out)


def scan(workers):
    tables = Parallel(n_jobs=workers)(delayed(scan_path)(p) for p in paths())
    table = pd.concat(tables, ignore_index=True)
    table.to_pickle(CACHE/"oame_scan.pkl")
    labels = classify(table)
    RESULTS.mkdir(exist_ok=True)
    labels.to_csv(RESULTS/"oame_paths.csv", index=False)
    print(labels.groupby(["label", "completed"]).size().to_string())
    print(labels.describe().round(2).to_string())


MAX_FRAMES, NEAR_PORTAL, BULK_FROM_HOST, BULK_SAMPLE = 300, 10., 15., 60
PORTAL = np.array([*AXIS, PORTAL_Z])
OP_COLUMNS = ("Q", "Q_avg", "lsi", "d5", "cn", "cn_smooth", "zeta", "n_hb", "q6_avg", "q3", "q3_avg", "f4",
              "s", "s_avg", "vor_density", "vor_asphericity", "ring_f5", "ring_f6", "ring_f7")


def water_ops(u, frame):
    from SolOrder import OrderParameter, eFingerprint
    op = OrderParameter(u, frame=frame, cutoff=5., selection="resname WAT and name O",
                        hydrogen_selection="resname WAT and (name H1 or name H2)")
    ent = eFingerprint(u, frame=frame, cutoff=5., selection="resname WAT and name O")
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


def next_entry(scan_path):
    """For each scanned frame and water: frames until it next becomes a cavity water (inf if never)."""
    cav = [set(c.tolist()) for c in scan_path.cavity]
    entries = {}
    prev = set()
    for f, cur in zip(scan_path.frame, cav):
        for w in cur - prev:
            entries.setdefault(w, []).append(f)
        prev = cur
    return {w: np.array(sorted(fs)) for w, fs in entries.items()}


def compute_path(path, scan_path, label):
    geo = Geometry(path)
    n = len(geo.u.trajectory)
    frames = np.unique(np.linspace(0, n-1, min(n, MAX_FRAMES)).round().astype(int))
    entries = next_entry(scan_path)
    rng = np.random.default_rng([2027, int(path.split("_")[1])])
    host = geo.u.select_atoms("resname OCB")
    rows = []
    for f in frames:
        geo.u.trajectory[f]
        g = geo.frame()
        ops = water_ops(geo.u, f)
        box = geo.u.dimensions
        d_portal = np.hypot(g["water_r"], g["water_z"] - PORTAL_Z)
        d_host = distance_array(geo.water_o.positions, host.center_of_mass()[None], box=box)[:, 0]
        region = np.full(len(d_portal), "", dtype=object)
        region[d_portal < NEAR_PORTAL] = "near_portal"
        region[g["bridge"]] = "bridge"
        region[g["cavity"]] = "cavity"
        bulk = np.flatnonzero((d_host > BULK_FROM_HOST) & (region == ""))
        region[rng.choice(bulk, size=min(BULK_SAMPLE, len(bulk)), replace=False)] = "bulk"
        keep = np.flatnonzero(region != "")
        until = np.full(len(keep), np.inf)
        for i, w in enumerate(keep):
            later = entries.get(w)
            if later is not None:
                after = later[later > f]
                if len(after):
                    until[i] = after[0] - f
        rows.append(pd.DataFrame(dict(path=path, label=label, frame=f, guest_z=g["guest_z"], water=keep,
                                      region=region[keep], z=g["water_z"][keep], r=g["water_r"][keep],
                                      d_portal=d_portal[keep].astype(np.float32), frames_to_entry=until,
                                      **{c: np.asarray(ops[c])[keep].astype(np.float32) for c in OP_COLUMNS})))
    return pd.concat(rows, ignore_index=True)


def compute(workers):
    scan_table = pd.read_pickle(CACHE/"oame_scan.pkl")
    labels = pd.read_csv(RESULTS/"oame_paths.csv").set_index("path")
    groups = dict(tuple(scan_table.groupby("path")))
    tables = Parallel(n_jobs=workers)(delayed(compute_path)(p, groups[p], labels.loc[p, "label"]) for p in paths())
    table = pd.concat(tables, ignore_index=True)
    table.to_pickle(CACHE/"oame_ops.pkl")
    print(f"oame_ops: {len(table)} rows, {table.frame.groupby(table.path).nunique().sum()} frames")
    print(table.groupby(["region", "label"]).size().unstack().to_string())


HORIZON = 10  # original path frames
FEATURES = CACHE/"applications"


def export_features(name):
    """Write the per-water table as a compressed, pickle-free archive."""
    t = pd.read_pickle(CACHE/f"{name}_ops.pkl")
    FEATURES.mkdir(parents=True, exist_ok=True)
    arrays = {c: (t[c].to_numpy(str) if t[c].dtype.kind in "OUT" or str(t[c].dtype) == "str" else t[c].to_numpy())
              for c in t.columns}
    np.savez_compressed(FEATURES/f"{name}_features.npz", **arrays)
    print(f"{FEATURES/name}_features.npz: {len(t)} rows, {len(t.columns)} columns")


def load_features(name):
    """The per-water table: cache/<name>_ops.pkl if computed here, else the archived npz."""
    pkl = CACHE/f"{name}_ops.pkl"
    if pkl.exists():
        return pd.read_pickle(pkl)
    with np.load(FEATURES/f"{name}_features.npz") as z:
        return pd.DataFrame({c: z[c] for c in z.files})


def oracle(a, b):
    from scipy.stats import ks_2samp
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    return (1 + ks_2samp(a, b).statistic)/2 if len(a) > 20 and len(b) > 20 else np.nan


def heldout(a, b, fold_a, fold_b):
    """Worse of the two directions: threshold fitted on one half of the paths, applied to the other."""
    from metrics import fit_threshold, hits
    scores = []
    for train in (0, 1):
        ta, tb = a[fold_a == train], b[fold_b == train]
        sa, sb = a[fold_a != train], b[fold_b != train]
        if min(len(ta), len(tb), len(sa), len(sb)) < 20:
            return np.nan
        t, below = fit_threshold(ta, tb)
        scores.append(.5*(hits(sa, t, below).mean() + hits(sb, t, not below).mean()))
    return min(scores)


def stage(z):
    return np.where(z < PORTAL_Z, "inside", np.where(z < OUT_Z, "crossing", "out"))


def analyze(workers):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    t = load_features("oame")
    t["fold"] = t.path.str.split("_").str[1].astype(int) % 2
    t["stage"] = stage(t.guest_z.to_numpy())
    t["path_type"] = np.where(t.label == "dry", "dry", np.where(t.label.isin(("wet", "intermediate")), "wet", "none"))
    rows = []
    for (region, reference) in (("cavity", "bulk"), ("bridge", "bulk"), ("cavity", "near_portal")):
        bulk = t[t.region == reference]
        for path_type in ("dry", "wet", "all"):
            for st in ("inside", "crossing", "out", "all"):
                if reference == "near_portal" and (path_type, st) != ("all", "all"):
                    continue
                part = t[t.region == region]
                ref = bulk
                if path_type != "all":
                    part, ref = part[part.path_type == path_type], ref[ref.path_type == path_type]
                if st != "all":
                    part, ref = part[part.stage == st], ref[ref.stage == st]
                if len(part) < 50:
                    continue
                for c in OP_COLUMNS:
                    a, b = part[c].to_numpy(float), ref[c].to_numpy(float)
                    rows.append(dict(region=region, reference=reference, path_type=path_type, stage=st,
                                     column=c, n=len(part),
                                     paths=part.path.nunique(), oracle=oracle(a, b),
                                     heldout=heldout(a, b, part.fold.to_numpy(), ref.fold.to_numpy()),
                                     median=float(np.nanmedian(a)), bulk_median=float(np.nanmedian(b))))
    regions = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    regions.to_csv(RESULTS/"oame_regions.csv", index=False)

    # dry versus wet: cavity waters of the two kinds of path
    cav = t[t.region == "cavity"]
    drywet = [dict(column=c, oracle=oracle(cav.loc[cav.path_type == "dry", c].to_numpy(float),
                                            cav.loc[cav.path_type == "wet", c].to_numpy(float)),
                   median_dry=float(cav.loc[cav.path_type == "dry", c].median()),
                   median_wet=float(cav.loc[cav.path_type == "wet", c].median())) for c in OP_COLUMNS]
    pd.DataFrame(drywet).to_csv(RESULTS/"oame_drywet.csv", index=False)

    # entry prediction among waters near the portal, guest not yet out
    cand = t[(t.region.isin(("near_portal", "bridge"))) & (t.guest_z < OUT_Z)].copy()
    cand["y"] = (cand.frames_to_entry <= HORIZON).astype(int)
    position = ["d_portal", "z", "r"]
    out = []
    for c in ["d_portal", "z", "r", *OP_COLUMNS]:
        x = cand[c].to_numpy(float)
        ok = np.isfinite(x)
        auc = roc_auc_score(cand.y[ok], x[ok])
        out.append(dict(model=c, features=1, auc=max(auc, 1-auc), direction="+" if auc >= .5 else "-"))

    def grouped(features):
        aucs = []
        for train in (0, 1):
            tr, te = cand[cand.fold == train], cand[cand.fold != train]
            X_tr, X_te = tr[features].to_numpy(float), te[features].to_numpy(float)
            mu, sd = np.nanmean(X_tr, axis=0), np.nanstd(X_tr, axis=0)
            X_tr = np.nan_to_num((X_tr - mu)/sd)
            X_te = np.nan_to_num((X_te - mu)/sd)
            model = LogisticRegression(class_weight="balanced", max_iter=500).fit(X_tr, tr.y)
            aucs.append(roc_auc_score(te.y, model.predict_proba(X_te)[:, 1]))
        return min(aucs)
    out.append(dict(model="position (logistic)", features=3, auc=grouped(position), direction="held-out"))
    for c in OP_COLUMNS:
        out.append(dict(model=f"position + {c}", features=4, auc=grouped(position + [c]), direction="held-out"))
    out.append(dict(model="position + all", features=3 + len(OP_COLUMNS),
                    auc=grouped(position + list(OP_COLUMNS)), direction="held-out"))
    entry = pd.DataFrame(out)
    entry.to_csv(RESULTS/"oame_entry.csv", index=False)
    pd.set_option("display.width", 220)
    pd.set_option("display.max_rows", 200)
    print(f"candidates {len(cand)}, entering within {HORIZON} frames: {cand.y.sum()} ({cand.y.mean():.3%}), "
          f"paths {cand.path.nunique()}")
    print(entry.round(3).to_string(index=False))
    for region, reference in (("cavity", "bulk"), ("bridge", "bulk"), ("cavity", "near_portal")):
        part = regions[(regions.region == region) & (regions.reference == reference) & (regions.stage == "all")]
        print(region, "vs", reference)
        print(part.pivot(index="column", columns="path_type", values="heldout").round(3).to_string())
    print(pd.DataFrame(drywet).round(3).sort_values("oracle", ascending=False).to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=("scan", "compute", "analyze", "export"))
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    {"scan": scan, "compute": compute, "analyze": analyze,
     "export": lambda workers: export_features("oame")}[args.step](args.workers)


if __name__ == "__main__":
    main()
