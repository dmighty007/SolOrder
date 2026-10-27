"""Planar-interface calibration of nucleus-size undercounting: forward and inverse (run F).

Definitions. A crystal of N molecules has the equimolar radius
R = (3N / 4 pi n)^(1/3), with n the number density of the constructed lattice.
A labeler has bulk recall f (recall of slab molecules deeper than DEEP) and
misses, per unit area of a planar interface, Gamma = sum_i (1 - p_i/f) / area
molecules, p_i the labeled outcome of crystal molecule i (one face, mean over
frames). Its dead layer is the count

    delta = Gamma / n,

which needs no reference plane: to leading order in 1/R a sphere misses
4 pi R^2 Gamma molecules whatever plane its depth is measured from.

  forward:  N_lab = f n (4 pi/3) (R - delta)^3
  inverse:  R     = (3 N_lab / (4 pi f n))^(1/3) + delta,  N = (4 pi/3) n R^3

Sensitivity (reported, never used to choose): the integral definition
delta = int (1 - p(z)/f) dz below the outermost oxygen plane, on BIN-wide
bins. A basal slab has molecules only at discrete depths (0, 2.76, 3.68,
6.44 A, ...), so the integral depends on how empty bins are filled: forward
fill (the definition of an earlier draft), linear interpolation or backward
fill; "integral_equimolar" adds the 0.92 A from the outermost plane to the
slab's equimolar surface.

Every seed (50-1600 molecules, Ih/Ic/Isd, three preparations) is held out:
nothing in the calibration comes from seeds. Observables, each calibrated
with its own per-molecule recall on the slabs:
  correct          molecules of the constructed crystal labeled with their
                   correct polytype (CHILL, CHILL+, PTM)
  cluster          largest connected cluster of molecules labeled Ih or Ic,
                   which can include liquid molecules (CHILL+, PTM)
  cluster_interf   CHILL+ largest cluster including the interfacial class
Each inverse estimate is compared with the constructed size N and with the
intact size (constructed molecules still within 1 A of their site with all
lattice bonds, interface_depth.py), mean over the system's frames. Errors are
reported per system (one seed, one preparation, mean over its 10 frames) and
pooled (mean over preparations and polytypes per size).

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
TESTED = 200  # smallest seed size in the validated range, fixed before the confirmation set (run G)
OBSERVABLES = {"correct": ("chill", "chill_plus", "ptm"), "cluster": ("chill_plus", "ptm"),
               "cluster_interf": ("chill_plus",)}
CLUSTER_COLUMN = {"cluster": "{}_largest", "cluster_interf": "{}_largest_with_interface"}
PRIMARY = "count"
DEFINITIONS = ("count", "integral_ffill", "integral_linear", "integral_bfill", "integral_equimolar")


def per_molecule(table, method, observable):
    o = outcomes(table, method)
    if observable == "correct":
        return o["correct"]
    if observable == "cluster":
        return o["correct"] | o["wrong_polytype"]
    return o["correct"] | o["wrong_polytype"] | o["interfacial"]


def slab_area(name):
    box = pd.read_json(INTERFACE/f"{name}.prep.json", typ="series")["box"]
    return box[0]*box[1]


def slab_offset(name):
    """Distance from the outermost oxygen plane to the slab's equimolar surface (A)."""
    truth = np.load(INTERFACE/f"{name}.truth.npz")
    n = int(truth["seed"].sum())
    span = float(np.nanmax(truth["depth"]))*2  # outermost plane to outermost plane
    return (n/(N_DENSITY*slab_area(name)) - span)/2


def calibrate(slabs, method, observable):
    """Bulk recall and the dead layer under every definition, from planar slabs only."""
    p = per_molecule(slabs, method, observable).astype(float)
    depth = slabs.depth.to_numpy()
    f = float(p[depth > DEEP].mean())
    per_slab = []
    for name, g in slabs.groupby("name"):
        missing = np.sum(1 - per_molecule(g, method, observable).astype(float)/f)/g.frame.nunique()/2
        per_slab.append(missing/(N_DENSITY*slab_area(name)))
    bins = np.floor(depth/BIN).astype(int)
    recall = pd.Series(p).groupby(bins).mean().reindex(np.arange(bins.max()+1))
    integral = {fill: float(np.sum(1 - getattr(recall, fill)().to_numpy()/f)*BIN)
                for fill in ("ffill", "interpolate", "bfill")}
    offset = float(np.mean([slab_offset(n) for n in slabs.name.unique()]))
    return dict(bulk_recall=f, count=float(np.mean(per_slab)), count_min=float(np.min(per_slab)),
                count_max=float(np.max(per_slab)), integral_ffill=integral["ffill"],
                integral_linear=integral["interpolate"], integral_bfill=integral["bfill"],
                integral_equimolar=integral["ffill"] + offset, equimolar_offset=offset)


def forward(n, f, delta):
    R = (3*n/(4*np.pi*N_DENSITY))**(1/3)
    return f*N_DENSITY*4/3*np.pi*max(R - delta, 0)**3


def inverse(observed, f, delta):
    return 4/3*np.pi*N_DENSITY*((3*observed/(4*np.pi*f*N_DENSITY))**(1/3) + delta)**3


def observed_size(sub, clusters, method, observable):
    if observable == "correct":
        return per_molecule(sub, method, observable).sum()/sub.frame.nunique()
    c = clusters[clusters.name == sub.name.iloc[0]]
    return float(c[CLUSTER_COLUMN[observable].format(method)].mean())


def apply(seeds, clusters, calibration, definitions=DEFINITIONS):
    """Forward and inverse estimates for every seed system, observable, labeler and definition."""
    rows = []
    for c in calibration.itertuples():
        for definition in definitions:
            delta = getattr(c, definition)
            for name, sub in seeds.groupby("name"):
                n = int(sub["size"].iloc[0])
                intact = float(sub.intact.sum()/sub.frame.nunique())
                observed = observed_size(sub, clusters, c.method, c.observable)
                fwd, inv = forward(n, c.bulk_recall, delta), inverse(observed, c.bulk_recall, delta)
                rows.append(dict(definition=definition, observable=c.observable, method=c.method, name=name,
                                 kind=sub.kind.iloc[0], size=n, prep=int(sub.prep.iloc[0]),
                                 R_equimolar=(3*n/(4*np.pi*N_DENSITY))**(1/3), intact=intact, observed=observed,
                                 forward=fwd, forward_error=(fwd-observed)/observed if observed else np.nan,
                                 inverse=inv, inverse_error=(inv-n)/n, inverse_error_intact=(inv-intact)/intact))
    return pd.DataFrame(rows)


def summary(group):
    pooled = group.groupby("size")[["observed", "forward", "inverse", "intact"]].mean()
    by_prep = group.groupby("prep").inverse_error.mean()
    return pd.Series(dict(
        systems=len(group),
        pooled_forward=float(((pooled.forward - pooled.observed)/pooled.observed).abs().max()),
        pooled_inverse=float(((pooled.inverse - pooled.index)/pooled.index).abs().max()),
        individual_forward_median=float(group.forward_error.abs().median()),
        individual_forward_max=float(group.forward_error.abs().max()),
        individual_inverse_median=float(group.inverse_error.abs().median()),
        individual_inverse_max=float(group.inverse_error.abs().max()),
        signed_inverse_mean=float(group.inverse_error.mean()),
        signed_inverse_prep_min=float(by_prep.min()), signed_inverse_prep_max=float(by_prep.max()),
        pooled_inverse_intact=float(((pooled.inverse - pooled.intact)/pooled.intact).abs().max()),
        individual_inverse_intact_median=float(group.inverse_error_intact.abs().median()),
        individual_inverse_intact_max=float(group.inverse_error_intact.abs().max()),
        signed_inverse_intact_mean=float(group.inverse_error_intact.mean())))


def calibration_table(slabs):
    return pd.DataFrame([dict(observable=o, method=m, **calibrate(slabs, m, o))
                         for o, methods in OBSERVABLES.items() for m in methods])


def main():
    table = pd.read_pickle(CACHE/"interface.pkl")
    clusters = pd.read_pickle(CACHE/"interface_clusters.pkl")
    seeds, slabs = table[table.geometry == "seed"], table[table.geometry == "slab"]
    calibration = calibration_table(slabs)
    out = apply(seeds, clusters, calibration)
    RESULTS.mkdir(exist_ok=True)
    out.to_csv(RESULTS/"nucleus_size.csv", index=False)
    calibration.to_csv(RESULTS/"nucleus_size_calibration.csv", index=False)
    tested = out[out["size"] >= TESTED]
    table_ = tested.groupby(["definition", "observable", "method"]).apply(summary, include_groups=False).reset_index()
    table_.to_csv(RESULTS/"nucleus_size_summary.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(calibration.round(3).to_string(index=False))
    print(f"\nseeds of {TESTED}-1600 molecules (pooled: max over sizes of |mean error|; individual: per system)")
    print(table_.round(3).to_string(index=False))
    print("\ninverse error by size, primary definition (mean, min, max over systems)")
    primary = out[out.definition == PRIMARY]
    print(primary.groupby(["observable", "method", "size"])[["forward_error", "inverse_error", "inverse_error_intact"]]
          .agg(["mean", "min", "max"]).round(3).to_string())


if __name__ == "__main__":
    main()
