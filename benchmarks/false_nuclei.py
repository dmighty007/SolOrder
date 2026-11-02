"""False crystalline clusters in liquid references, with frozen thresholds.

A balanced per-molecule accuracy does not say whether a labeling rule reports a
nucleus where there is none. Here every rule is frozen before it sees the test
data: scalar parameters get the threshold and direction that separate thermal
Ih from liquid in preparations C1-C3 (TIP4P/Ice, 220 K); CHILL, CHILL+ and PTM
keep their published settings. The rules are then applied to liquids not used
for fitting: preparations C4-C5, TIP4P/2005 at 220 and 200 K, mW at 220 K (the
TIP4P/Ice threshold, unchanged), and liquid at 0.3-3.0 GPa. Per frame we record
the fraction of molecules labeled ice and the largest connected cluster of
them (O-O edges below 3.5 A). Persistence: of the molecules in the largest
cluster of one frame, the fraction that are again in an ice-labeled cluster of
at least 5 molecules 5 ps later. Spanning: the largest cluster connects to its
own periodic image (an edge closes a loop with a nonzero lattice translation).

Two controls. (1) Thresholds chosen for balanced accuracy accept many liquid
molecules by construction, so every scalar rule is also applied with a
threshold that labels 1% of the fitting liquid as ice (rule "<column>@fp1").
(2) System size: the largest false cluster grows with the number of molecules
inspected, so the rules are also applied to the three independent
32,768-molecule TIP4P/Ice liquids of run F (220 K, one configuration each,
results/false_nuclei_size.csv).

The price of the stricter threshold, recall of thermal Ih, Ic and Isd in the
held-out preparations C4-C5, is in results/false_nuclei_thresholds.csv.

  python false_nuclei.py   # results/false_nuclei*.csv
"""

import json
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_box
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components

from metrics import fit_threshold, hits

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
SCALARS = ("q6_avg", "entropy_avg", "q_tet_avg", "q3_avg", "q4_avg", "ring_f6", "xconn6", "cn_smooth", "d5")
LABELERS = {"chill": (1, 2), "chill_plus": (1, 2), "chill_plus_with_interface": (1, 2, 3), "ptm": (1, 2)}
EDGE, MIN_CLUSTER = 3.5, 5


def components(x, box, mask):
    if not mask.any():
        return np.full(len(mask), -1)
    idx = np.flatnonzero(mask)
    pairs = capped_distance(x[idx], x[idx], EDGE, box=box, return_distances=False)
    pairs = pairs[pairs[:, 0] < pairs[:, 1]]
    n = len(idx)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, comp = connected_components(graph, directed=False)
    out = np.full(len(mask), -1)
    out[idx] = comp
    return out


def spans(x, box, members):
    """True if the cluster of oxygens `members` connects to its own periodic image."""
    if len(members) < 3:
        return False
    y = x[members]
    pairs = capped_distance(y, y, EDGE, box=box, return_distances=False)
    pairs = pairs[pairs[:, 0] < pairs[:, 1]]
    n = len(members)
    vectors = minimize_vectors(y[pairs[:, 1]] - y[pairs[:, 0]], box)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)).tocsr()
    order, parent = breadth_first_order(graph, 0, directed=False)
    unwrapped = np.zeros_like(y)
    lookup = {(a, b): v for (a, b), v in zip(map(tuple, pairs), vectors)}
    for node in order[1:]:
        a = parent[node]
        step = lookup[(a, node)] if (a, node) in lookup else -lookup[(node, a)]
        unwrapped[node] = unwrapped[a] + step
    closing = unwrapped[pairs[:, 1]] - unwrapped[pairs[:, 0]] - vectors
    return bool(np.any(np.linalg.norm(closing, axis=1) > 1.))


def references():
    out = []
    for prep in ("C4", "C5"):
        out.append(dict(name=f"TIP4P/Ice 220 K ({prep})", table="prep", npz=CACHE/f"replicas/liquid_p{prep[1]}_220.npz",
                        key=dict(run=prep, phase="liquid", temperature=220.)))
    for path in sorted((CACHE/"conditions").glob("*liquid*.npz")) + sorted((CACHE/"conditions/mw_inputs").glob("*liquid*.npz")):
        info = json.loads(path.with_suffix(".json").read_text())
        if info["model"] == "tip4pice" and info["dataset"] == "D":
            continue  # protocol control, same model and state as run C
        if info["dataset"] == "E" and info["prep"] == 3:
            continue  # held out for the ice V study
        model = {"tip4pice": "TIP4P/Ice", "tip4p2005": "TIP4P/2005", "mW": "mW"}[info["model"]]
        state = f"{info['pressure_gpa']} GPa" if info["dataset"] == "E" else f"{int(info['temperature'])} K"
        out.append(dict(name=f"{model} {state} (p{info['prep']})", table="pressure" if info["dataset"] == "E" else "model",
                        npz=path, key=dict(phase="liquid", model=info["model"], temperature=info["temperature"],
                                           pressure=info["pressure_gpa"], run=f"p{info['prep']}")))
    return out


def main():
    tables = {name: pd.read_pickle(CACHE/f"{name}.pkl") for name in ("prep", "model", "pressure")}
    prep = tables["prep"]
    fit = prep[(prep.kind == "thermal") & (prep.temperature == 220.) & prep.run.isin(("C1", "C2", "C3"))]
    rules = {}
    for c in SCALARS:
        liquid = fit.loc[fit.phase == "liquid", c].to_numpy(float)
        t, below = fit_threshold(fit.loc[fit.phase == "Ih", c].to_numpy(float), liquid)
        rules[c] = (t, below)
        rules[f"{c}@fp1"] = (float(np.nanquantile(liquid, .01 if below else .99)), below)
    threshold_table(rules, prep).to_csv(RESULTS/"false_nuclei_thresholds.csv", index=False)
    frames_out, summary = [], []
    for ref in references():
        table = tables[ref["table"]]
        mask = table.kind == "thermal"
        for column, value in ref["key"].items():
            if column in table:
                mask &= table[column] == value
        part = table[mask].sort_values(["frame", "mol"])
        data = np.load(ref["npz"])
        cell = data["cell"] if "cell" in data.files else np.diag(data["box"])
        box = triclinic_box(*cell)
        sites = 1 if "mw_inputs" in str(ref["npz"]) else 4
        labels = {}
        for c, (t, below) in rules.items():
            values = part[c.split("@")[0]].to_numpy(float)
            if np.isnan(values).all():
                continue
            labels[c] = hits(values, t, below)
        for name, codes in LABELERS.items():
            source = "chill_plus" if name.startswith("chill_plus") else name
            labels[name] = np.isin(part[source].to_numpy(), codes)
        frame_ids = part.frame.to_numpy()
        frames = np.unique(frame_ids)
        for rule, ice in labels.items():
            previous, largest, fractions, persistence, spanning = None, [], [], [], []
            for f in frames:
                m = ice[frame_ids == f]
                x = data["x"][f][::sites].astype(np.float32)
                comp = components(x, box, m)
                sizes = np.bincount(comp[comp >= 0]) if (comp >= 0).any() else np.array([0])
                biggest = int(sizes.max())
                in_cluster = np.isin(comp, np.flatnonzero(sizes >= MIN_CLUSTER))
                if previous is not None and previous.any():
                    persistence.append(float(in_cluster[previous].mean()))
                previous = (comp == sizes.argmax()) if biggest >= MIN_CLUSTER else None
                spanning.append(spans(x, box, np.flatnonzero(comp == sizes.argmax())) if biggest >= 20 else False)
                largest.append(biggest)
                fractions.append(float(m.mean()))
                frames_out.append(dict(reference=ref["name"], rule=rule, frame=int(f), fraction=float(m.mean()),
                                       largest=biggest, spans=spanning[-1]))
            summary.append(dict(reference=ref["name"], rule=rule, molecules=int((frame_ids == frames[0]).sum()),
                                false_positive=float(np.mean(fractions)), largest_mean=float(np.mean(largest)),
                                largest_max=int(np.max(largest)),
                                frames_with_cluster_20=float(np.mean(np.array(largest) >= 20)),
                                frames_spanning=float(np.mean(spanning)),
                                persistence_5ps=float(np.mean(persistence)) if persistence else np.nan))
        print(ref["name"], flush=True)
    RESULTS.mkdir(exist_ok=True)
    pd.DataFrame(frames_out).to_csv(RESULTS/"false_nuclei_frames.csv", index=False)
    summary = pd.DataFrame(summary)
    summary.to_csv(RESULTS/"false_nuclei.csv", index=False)
    size = size_check(rules)
    size.to_csv(RESULTS/"false_nuclei_size.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 400)
    print(summary.round(3).to_string(index=False))
    print(size.round(3).to_string(index=False))


def threshold_table(rules, prep):
    """Fraction of each held-out thermal phase (C4, C5, 220 K) that every scalar rule labels ice."""
    test = prep[(prep.kind == "thermal") & (prep.temperature == 220.) & prep.run.isin(("C4", "C5"))]
    rows = []
    for rule, (t, below) in rules.items():
        column = rule.split("@")[0]
        rows.append(dict(rule=rule, threshold=t, ice_below=below,
                         **{phase: float(hits(test.loc[test.phase == phase, column].to_numpy(float), t, below).mean())
                            for phase in ("Ih", "Ic", "Isd", "liquid")}))
    return pd.DataFrame(rows)


def size_check(rules):
    """Largest false cluster in the 32,768-molecule run F liquids, with the same frozen rules."""
    from compute import evaluate
    rows = []
    for path in sorted((CACHE/"interface").glob("F_liquid_p*.npz")):
        data = np.load(path)
        x, box3 = data["x"], np.asarray(data["box"], float)
        u = mda.Universe.empty(4*len(x), n_residues=len(x), atom_resindex=np.repeat(np.arange(len(x)), 4),
                               trajectory=True)
        u.add_TopologyAttr("names", np.tile(["OW", "HW1", "HW2", "MW"], len(x)))
        u.atoms.positions = x.reshape(-1, 3)
        u.dimensions = [*box3, 90., 90., 90.]
        box = u.dimensions
        columns = evaluate(u, "name OW", "name HW1 HW2")
        o = x[:, 0].astype(np.float32)
        labels = {c: hits(np.asarray(columns[c.split("@")[0]], float), t, below) for c, (t, below) in rules.items()}
        for name, codes in LABELERS.items():
            source = "chill_plus" if name.startswith("chill_plus") else name
            labels[name] = np.isin(columns[source], codes)
        for rule, ice in labels.items():
            comp = components(o, box, ice)
            sizes = np.bincount(comp[comp >= 0]) if (comp >= 0).any() else np.array([0])
            biggest = int(sizes.max())
            rows.append(dict(reference=f"TIP4P/Ice 220 K, 32,768 molecules ({path.stem[-2:]})", rule=rule,
                             molecules=len(x), false_positive=float(ice.mean()), largest=biggest,
                             spans=spans(o, box, np.flatnonzero(comp == sizes.argmax())) if biggest >= 20 else False))
        print(path.stem, flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
