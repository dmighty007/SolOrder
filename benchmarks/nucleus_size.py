"""Planar-interface calibration of nucleus-size undercounting: forward and inverse (run F).

Convention. A crystal of N molecules has the equimolar (dividing-surface)
radius R_s = (3N / 4 pi n)^(1/3), with n the number density of the
constructed lattice. A labeler with bulk recall f misses a surface layer of
effective thickness delta_s, measured from the equimolar surface:

  calibration (planar slabs only):
      f        = recall of slab molecules deeper than 8 A
      delta_s  = integral over depth of (1 - p(z)/f) dz,
                 p(z) the recall at depth z below the outermost oxygen plane
  The planar dead layer is measured from the outermost molecular plane, which
  on these rough, cut spheres coincides with the equimolar radius (within
  0.4 A). Measuring it instead from the slab's own equimolar surface, which
  lies a = 0.92 A outside a flat basal facet, is reported as the alternative
  convention ("equimolar"); the two differ by about 20% in N at R = 15 A.
  forward:  N_lab = f n (4 pi/3) (R_s - delta_s)^3
  inverse:  R_s   = (3 N_lab / (4 pi f n))^(1/3) + delta_s,  N = (4 pi/3) n R_s^3

Every seed (50-1600 molecules, Ih/Ic/Isd, three preparations) is held out:
nothing in the calibration comes from seeds. Observables, each calibrated
with its own per-molecule recall on the slabs:
  correct          molecules of the constructed crystal labeled with their
                   correct polytype (CHILL, CHILL+, PTM)
  cluster          largest connected cluster of molecules labeled Ih or Ic,
                   which can include liquid molecules (CHILL+, PTM)
  cluster_interf   CHILL+ largest cluster including the interfacial class
Errors are reported per system (one seed, one preparation, mean over its 10
frames) and pooled (mean over preparations and polytypes per size).

  python nucleus_size.py   # results/nucleus_size.csv, nucleus_size_summary.csv, nucleus_size_calibration.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

from interface_depth import outcomes

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
INTERFACE = CACHE/"interface"
B = 2.76  # A, O-O distance of the constructed lattice (md_replicas.py)
A = B*np.sqrt(8/3)
N_DENSITY = 3/(A*A*np.sqrt(3)*B)  # molecules per A^3 of the constructed lattice (0.0309)
BIN, DEEP = .25, 8.
OBSERVABLES = {"correct": ("chill", "chill_plus", "ptm"), "cluster": ("chill_plus", "ptm"),
               "cluster_interf": ("chill_plus",)}
CLUSTER_COLUMN = {"cluster": "{}_largest", "cluster_interf": "{}_largest_with_interface"}


def per_molecule(table, method, observable):
    o = outcomes(table, method)
    if observable == "correct":
        return o["correct"]
    if observable == "cluster":
        return o["correct"] | o["wrong_polytype"]
    return o["correct"] | o["wrong_polytype"] | o["interfacial"]


def slab_offset(name):
    """Distance from the outermost oxygen plane to the slab's equimolar surface (A)."""
    truth = np.load(INTERFACE/f"{name}.truth.npz")
    meta = pd.read_json(INTERFACE/f"{name}.prep.json", typ="series")
    n = int(truth["seed"].sum())
    area = meta["box"][0]*meta["box"][1]
    span = float(np.nanmax(truth["depth"]))*2  # outermost plane to outermost plane
    return (n/(N_DENSITY*area) - span)/2


def calibrate(slabs, method, observable):
    p = per_molecule(slabs, method, observable).astype(float)
    depth = slabs.depth.to_numpy()
    f = float(p[depth > DEEP].mean())
    bins = np.floor(depth/BIN).astype(int)
    recall = pd.Series(p).groupby(bins).mean()
    edges = np.arange(recall.index.max()+1)
    curve = recall.reindex(edges).ffill().to_numpy()
    offset = float(np.mean([slab_offset(n) for n in slabs.name.unique()]))
    delta = float(np.sum(1 - curve/f)*BIN) + offset
    return f, delta, offset


def main():
    table = pd.read_pickle(CACHE/"interface.pkl")
    clusters = pd.read_pickle(CACHE/"interface_clusters.pkl")
    seeds, slabs = table[table.geometry == "seed"], table[table.geometry == "slab"]
    calibration, rows = [], []
    for observable, methods in OBSERVABLES.items():
        for method in methods:
            f, delta_equimolar, offset = calibrate(slabs, method, observable)
            calibration.append(dict(observable=observable, method=method, bulk_recall=f,
                                    dead_layer=delta_equimolar - offset, dead_layer_equimolar=delta_equimolar,
                                    equimolar_offset=offset))
            for (convention, delta), (name, sub) in ((c, g) for c in (("molecular", delta_equimolar - offset),
                                                                       ("equimolar", delta_equimolar))
                                                     for g in seeds.groupby("name")):
                n = int(sub["size"].iloc[0])
                frames = sub.frame.nunique()
                if observable == "correct":
                    observed = per_molecule(sub, method, observable).sum()/frames
                else:
                    c = clusters[clusters.name == name]
                    observed = float(c[CLUSTER_COLUMN[observable].format(method)].mean())
                Rs = (3*n/(4*np.pi*N_DENSITY))**(1/3)
                forward = f*N_DENSITY*4/3*np.pi*max(Rs - delta, 0)**3
                Rs_inv = (3*observed/(4*np.pi*f*N_DENSITY))**(1/3) + delta
                inverse = 4/3*np.pi*N_DENSITY*Rs_inv**3
                rows.append(dict(convention=convention, observable=observable, method=method, name=name,
                                 kind=sub.kind.iloc[0],
                                 size=n, prep=int(sub.prep.iloc[0]), R_equimolar=Rs, observed=observed,
                                 forward=forward, forward_error=(forward-observed)/observed if observed else np.nan,
                                 inverse=inverse, inverse_error=(inverse-n)/n))
    out = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    out.to_csv(RESULTS/"nucleus_size.csv", index=False)
    pd.DataFrame(calibration).to_csv(RESULTS/"nucleus_size_calibration.csv", index=False)

    def summary(group):
        pooled = group.groupby("size")[["observed", "forward", "inverse"]].mean()
        return pd.Series(dict(
            systems=len(group),
            pooled_forward=float(((pooled.forward - pooled.observed)/pooled.observed).abs().max()),
            pooled_inverse=float(((pooled.inverse - pooled.index)/pooled.index).abs().max()),
            individual_forward_median=float(group.forward_error.abs().median()),
            individual_forward_max=float(group.forward_error.abs().max()),
            individual_inverse_median=float(group.inverse_error.abs().median()),
            individual_inverse_max=float(group.inverse_error.abs().max())))
    tested = out[out["size"] >= 200]
    table_ = tested.groupby(["convention", "observable", "method"]).apply(summary, include_groups=False).reset_index()
    table_.to_csv(RESULTS/"nucleus_size_summary.csv", index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(calibration).round(3).to_string(index=False))
    print("\nseeds of 200-1600 molecules (pooled: max over sizes of |mean error|; individual: per system)")
    print(table_.round(3).to_string(index=False))
    print("\ninverse error by size (mean, min, max over systems)")
    primary = out[out.convention == "molecular"]
    print(primary.groupby(["observable", "method", "size"])[["forward_error", "inverse_error"]]
          .agg(["mean", "min", "max"]).round(3).to_string())


if __name__ == "__main__":
    main()
