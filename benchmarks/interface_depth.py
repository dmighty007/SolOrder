"""Polytype recall against depth, and the nucleus-size bias it predicts (run F).

Step "compute": label every molecule of every run F frame with CHILL, CHILL+
and PTM, and audit whether each constructed crystal molecule is still intact
(within 1 A of its constructed site after removing the crystal's drift, all
constructed bonds shorter than 3.5 A). Per frame we also record the largest
connected cluster of molecules labeled ice (Ih or Ic; for CHILL and CHILL+ also
with the interfacial class), with 3.5 A O-O edges, as nucleation studies
measure it.

Step "metrics":
  recall      P(outcome | depth bin) per labeler and polytype, per preparation
              and pooled, for seeds and slabs
  prediction  the expected number of correctly labeled molecules of every seed,
              summed over its constructed molecules, from four models fitted
              without the seed sizes they predict (50, 400, 800, 1600):
                depth_seeds          P(correct | depth) from seeds of 100 and 200
                depth_slab           P(correct | depth) from the planar slabs
                support_seeds        P(correct | support) from seeds of 100 and
                                     200, where support is the fraction of a
                                     probe sphere (radius and margin fitted)
                                     around the molecule that lies in the crystal
                support_seeds_slab   the same, fitted on those seeds and the slabs
  dead layer  integral over depth of (1 - P(correct | depth))

  python interface_depth.py compute   # cache/interface.pkl, cache/interface_clusters.pkl
  python interface_depth.py metrics   # results/interface_*.csv
  python interface_depth.py compute --confirmation   # run G -> cache/confirmation.pkl, confirmation_clusters.pkl
"""

import argparse
from pathlib import Path
import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

from joblib import Parallel, delayed
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from SolOrder import OrderParameter

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
INTERFACE = CACHE/"interface"
LABELERS = ("chill", "chill_plus", "ptm")
TRAIN_SIZES = (100, 200)
BIN = .5  # A
INTACT_DISTANCE, BOND = 1., 3.5


def parse(name):
    """F_seed_<kind>_<N>_p<k>, F_slab_<kind>_p<k>; confirmation set G_seed_..., G_prism_..., G_long_..."""
    parts = name.split("_")
    if parts[1] == "seed":
        return dict(geometry="seed", kind=parts[2], size=int(parts[3]), prep=int(parts[4][1:]), facet="")
    return dict(geometry="slab", kind=parts[2], size=0, prep=int(parts[3][1:]),
                facet="prism" if parts[1] == "prism" else "basal")


def labels(op):
    out = {"chill": op.chill(cutoff=3.5).labels, "chill_plus": op.chill_plus(3.5).labels}
    match = op.template_matching()
    code = {v: k for k, v in match.type_names.items()}
    ptm = np.zeros(op.n, dtype=np.int8)
    ptm[match.labels == code["Hexagonal diamond"]] = 1
    ptm[match.labels == code["Cubic diamond"]] = 2
    out["ptm"] = ptm
    return out


def largest_cluster(mask, pairs):
    """Size of the largest connected set of masked oxygens (pairs: O-O edges within 3.5 A)."""
    if not mask.any():
        return 0
    keep = mask[pairs[:, 0]] & mask[pairs[:, 1]]
    p = pairs[keep]
    n = len(mask)
    graph = coo_matrix((np.ones(len(p)), (p[:, 0], p[:, 1])), shape=(n, n))
    _, component = connected_components(graph, directed=False)
    return int(np.bincount(component[mask]).max())


def system_task(name):
    meta = parse(name)
    truth = np.load(INTERFACE/f"{name}.truth.npz")
    frames = np.load(INTERFACE/f"{name}.npz")
    u = mda.Universe(str(INTERFACE/f"{name}.gro"))
    oxygens = u.select_atoms("name OW")
    crystal = truth["seed"]
    start = oxygens.positions[crystal].astype(float)
    box = np.r_[frames["box"], 90., 90., 90.]
    # constructed lattice bonds inside the crystal
    within = mda.lib.distances.self_capped_distance(start, 3.0, box=box)[0]
    rows, clusters = [], []
    times = frames["time_ps"] if "time_ps" in frames.files else 2.*(np.arange(len(frames["x"])) + 1)
    for f, x in enumerate(frames["x"]):
        u.atoms.positions = x
        u.dimensions = box
        op = OrderParameter(u, cutoff=5., selection="name OW", hydrogen_selection="name HW1 HW2")
        lab = labels(op)
        o = oxygens.positions.astype(float)
        moved = minimize_vectors(o[crystal] - start, box)
        moved -= moved.mean(axis=0)
        bond = np.linalg.norm(minimize_vectors(o[crystal][within[:, 0]] - o[crystal][within[:, 1]], box), axis=1)
        broken = np.bincount(within[bond > BOND].ravel(), minlength=crystal.sum())
        intact = (np.linalg.norm(moved, axis=1) < INTACT_DISTANCE) & (broken == 0)
        table = dict(frame=f, time_ps=float(times[f]), depth=truth["depth"][crystal], stack=truth["stack"][crystal], intact=intact,
                     **{k: v[crystal].astype(np.int8) for k, v in lab.items()})
        rows.append(pd.DataFrame(table))
        pairs = op.neighbor_list(3.5)
        edges = np.c_[np.repeat(np.arange(op.n), np.diff(pairs.offsets)), pairs.indices]
        cluster = dict(name=name, **meta, frame=f, time_ps=float(times[f]), constructed=int(crystal.sum()))
        for k, v in lab.items():
            ice = np.isin(v, (1, 2))
            cluster[f"{k}_largest"] = largest_cluster(ice, edges)
            cluster[f"{k}_in_crystal_correct"] = int(np.sum(ice[crystal] & (v[crystal] == truth["stack"][crystal])))
            cluster[f"{k}_outside"] = int(np.sum(ice[~crystal]))
            if k != "ptm":
                with_interface = np.isin(v, (1, 2, 3))
                cluster[f"{k}_largest_with_interface"] = largest_cluster(with_interface, edges)
        clusters.append(cluster)
    out = pd.concat(rows, ignore_index=True)
    for key, value in meta.items():
        out[key] = value
    out["name"] = name
    return out, pd.DataFrame(clusters)


def compute(names, workers, confirmation=False):
    """Run F -> cache/interface*.pkl; the confirmation set (run G) -> cache/confirmation*.pkl."""
    prefix, out = ("G_", "confirmation") if confirmation else ("F_s", "interface")
    names = names or sorted(p.name.removesuffix(".npz") for p in INTERFACE.glob(f"{prefix}*.npz")
                            if ".truth" not in p.name and "liquid" not in p.name)
    results = Parallel(n_jobs=workers)(delayed(system_task)(n) for n in names)
    pd.concat([r[0] for r in results], ignore_index=True).to_pickle(CACHE/f"{out}.pkl")
    pd.concat([r[1] for r in results], ignore_index=True).to_pickle(CACHE/f"{out}_clusters.pkl")
    print(f"{out}: {len(names)} systems")


def outcomes(table, method):
    v = table[method].to_numpy()
    ice = np.isin(v, (1, 2))
    return {"correct": ice & (v == table["stack"].to_numpy()), "wrong_polytype": ice & (v != table["stack"].to_numpy()),
            "interfacial": v == 3, "not_ice": ~ice & (v != 3)}


def curve(table, method):
    """P(correct | depth bin) as an array indexed by bin."""
    bins = np.floor(table.depth.to_numpy()/BIN).astype(int)
    correct = outcomes(table, method)["correct"]
    n = np.bincount(bins, minlength=bins.max()+1)
    k = np.bincount(bins, weights=correct, minlength=bins.max()+1)
    with np.errstate(invalid="ignore"):
        return k/np.where(n > 0, n, np.nan), n


def predict(p, depth):
    bins = np.minimum(np.floor(depth/BIN).astype(int), len(p)-1)
    filled = pd.Series(p).ffill().fillna(0).to_numpy()  # deeper than any training bin: use the deepest value
    return float(filled[bins].sum())


def support_fraction(depth, R, rho, margin):
    """Fraction of a probe sphere (radius rho) around each molecule that lies inside the crystal.

    Seeds: the crystal is a sphere of radius R + margin (R: outermost constructed oxygen) and the
    molecule sits at r = R - depth from its centre. Slabs (R infinite): a half-space whose face is
    margin beyond the outermost oxygen plane.
    """
    full = 4/3*np.pi*rho**3
    depth = np.asarray(depth, float)
    R = np.broadcast_to(np.asarray(R, float), depth.shape)
    planar = np.isinf(R)
    h = np.clip(rho - (depth + margin), 0, 2*rho)  # height of the probe cap outside a planar crystal
    slab = 1 - np.pi*h**2*(3*rho - h)/3/full
    Rc = np.where(planar, 1., R + margin)
    r = np.maximum(np.where(planar, 1., R - depth), 1e-9)
    lens = np.pi*(Rc + rho - r)**2*(r**2 + 2*r*rho - 3*rho**2 + 2*r*Rc + 6*rho*Rc - 3*Rc**2)/(12*r)
    v = np.where(r + rho <= Rc, full, np.where(r >= Rc + rho, 0., lens))
    v = np.where(r + Rc <= rho, 4/3*np.pi*Rc**3, v)  # probe swallows the crystal
    return np.clip(np.where(planar, slab, v/full), 0, 1)


RHOS = np.arange(2., 14.01, .5)
MARGINS = (0., 1.38)  # A; 1.38 is half the O-O bond length


def fit_support(train, method, radius):
    """Probe radius and margin maximizing the likelihood of a monotone recall(support) curve."""
    from sklearn.isotonic import IsotonicRegression
    y = outcomes(train, method)["correct"].astype(float)
    best = None
    for margin in MARGINS:
        for rho in RHOS:
            phi = support_fraction(train.depth.to_numpy(), radius, rho, margin)
            model = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip").fit(phi, y)
            p = np.clip(model.predict(phi), 1e-6, 1-1e-6)
            loglik = float(np.sum(y*np.log(p) + (1-y)*np.log(1-p)))
            if best is None or loglik > best[0]:
                best = (loglik, rho, margin, model)
    return best


def metrics():
    table = pd.read_pickle(CACHE/"interface.pkl")
    radius = {n: float(np.load(INTERFACE/f"{n}.truth.npz")["R"]) for n in table.name.unique()}
    table["R"] = table.name.map(radius)
    clusters = pd.read_pickle(CACHE/"interface_clusters.pkl")
    seeds = table[table.geometry == "seed"]
    rows = []
    pooled = table.assign(kind="all")
    for (geometry, kind), part in pd.concat([table, pooled]).groupby(["geometry", "kind"]):
        for method in LABELERS:
            sets = {"all": part}
            if geometry == "seed":
                sets["train"] = part[part["size"].isin(TRAIN_SIZES)]
                for size, sub in part.groupby("size"):
                    sets[f"n{size}"] = sub
            for prep, sub in part.groupby("prep"):
                sets[f"p{prep}"] = sub
            for label, sub in sets.items():
                bins = np.floor(sub.depth.to_numpy()/BIN).astype(int)
                result = outcomes(sub, method)
                for b in np.unique(bins):
                    m = bins == b
                    row = dict(geometry=geometry, kind=kind, method=method, subset=label, depth=(b+.5)*BIN,
                               n=int(m.sum()), intact=float(sub.intact.to_numpy()[m].mean()))
                    row.update({o: float(v[m].mean()) for o, v in result.items()})
                    rows.append(row)
    recall = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    recall.to_csv(RESULTS/"interface_recall.csv", index=False)
    predictions, layers, slab_rows = [], [], []
    slabs = table[table.geometry == "slab"]
    for kind, part in seeds.groupby("kind"):
        train = part[part["size"].isin(TRAIN_SIZES)]
        slab = slabs[slabs.kind == kind]
        both = pd.concat([train, slab])
        for method in LABELERS:
            models = {"depth_seeds": ("depth", curve(train, method)[0])}
            if len(slab):
                models["depth_slab"] = ("depth", curve(slab, method)[0])
            fits = {"support_seeds": fit_support(train, method, train.R.to_numpy())}
            if len(slab):
                fits["support_seeds_slab"] = fit_support(both, method, both.R.to_numpy())
            for name, (_, rho, margin, model) in fits.items():
                models[name] = ("support", (rho, margin, model))
            p = models["depth_seeds"][1]
            row = dict(kind=kind, method=method, bulk_recall=np.nanmean(p[-4:]),
                       dead_layer_seeds=float(np.nansum(1 - pd.Series(p).ffill().to_numpy())*BIN))
            if len(slab):
                ps = models["depth_slab"][1]
                row["dead_layer_slab"] = float(np.nansum(1 - pd.Series(ps).ffill().to_numpy())*BIN)
            for name in fits:
                row[f"{name}_probe_radius"], row[f"{name}_margin"] = fits[name][1], fits[name][2]
            layers.append(row)
            if len(slab):  # planar control: measured curve against curves fitted on seeds only
                correct = outcomes(slab, method)["correct"]
                bins = np.floor(slab.depth.to_numpy()/BIN).astype(int)
                rho, margin, model = models["support_seeds"][1]
                phi = support_fraction(slab.depth.to_numpy(), np.inf, rho, margin)
                filled = pd.Series(p).ffill().fillna(0).to_numpy()
                for b_ in np.unique(bins):
                    m = bins == b_
                    slab_rows.append(dict(kind=kind, method=method, depth=(b_+.5)*BIN, n=int(m.sum()),
                                          measured=float(correct[m].mean()),
                                          predicted_support_seeds=float(model.predict(phi[m]).mean()),
                                          predicted_depth_seeds=float(filled[min(b_, len(filled)-1)])))
            for name, sub in part.groupby("name"):
                frames = sub.frame.nunique()
                first = sub[sub.frame == sub.frame.min()]
                depth, R = first.depth.to_numpy(), first.R.iloc[0]
                c = clusters[clusters.name == name]
                out = dict(kind=kind, method=method, name=name, size=int(sub["size"].iloc[0]),
                           prep=int(sub.prep.iloc[0]), training=int(sub["size"].iloc[0]) in TRAIN_SIZES,
                           R=R, measured=outcomes(sub, method)["correct"].sum()/frames,
                           largest_cluster=float(c[f"{method}_largest"].mean()),
                           largest_with_interface=float(c[f"{method}_largest_with_interface"].mean())
                           if f"{method}_largest_with_interface" in c else np.nan,
                           intact=float(sub.intact.mean()))
                for mname, (form, params) in models.items():
                    if form == "depth":
                        out[mname] = predict(params, depth)
                    else:
                        rho, margin, model = params
                        out[mname] = float(model.predict(support_fraction(depth, R, rho, margin)).sum())
                predictions.append(out)
    predictions = pd.DataFrame(predictions)
    predictions.to_csv(RESULTS/"interface_prediction.csv", index=False)
    pd.DataFrame(layers).to_csv(RESULTS/"interface_dead_layer.csv", index=False)
    pd.DataFrame(slab_rows).to_csv(RESULTS/"interface_slab_prediction.csv", index=False)
    names = [m for m in ("depth_seeds", "depth_slab", "support_seeds", "support_seeds_slab") if m in predictions]
    test = predictions[~predictions.training]
    errors = test.assign(**{m: (test[m] - test.measured)/test.measured for m in names})
    errors.groupby(["method", "size"])[names].mean().to_csv(RESULTS/"interface_prediction_errors.csv")
    pd.set_option("display.width", 220)
    print(predictions.groupby(["method", "kind", "size"])[names + ["measured", "largest_with_interface"]]
          .mean().round(1).to_string())
    print("\nmean relative error on unseen sizes:")
    print(errors.groupby(["method", "size"])[names].mean().round(3).to_string())
    print(pd.DataFrame(layers).round(2).to_string(index=False))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=("compute", "metrics"))
    parser.add_argument("names", nargs="*")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--confirmation", action="store_true", help="the run G confirmation set")
    args = parser.parse_args()
    compute(args.names, args.workers, args.confirmation) if args.step == "compute" else metrics()
