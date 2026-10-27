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
at least 5 molecules 5 ps later.

  python false_nuclei.py   # results/false_nuclei.csv, false_nuclei_frames.csv
"""

import json
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.lib.mdamath import triclinic_box
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

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
        t, below = fit_threshold(fit.loc[fit.phase == "Ih", c].to_numpy(float), fit.loc[fit.phase == "liquid", c].to_numpy(float))
        rules[c] = (t, below)
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
            values = part[c].to_numpy(float)
            if np.isnan(values).all():
                continue
            labels[c] = hits(values, t, below)
        for name, codes in LABELERS.items():
            source = "chill_plus" if name.startswith("chill_plus") else name
            labels[name] = np.isin(part[source].to_numpy(), codes)
        frame_ids = part.frame.to_numpy()
        frames = np.unique(frame_ids)
        for rule, ice in labels.items():
            previous, largest, fractions, persistence = None, [], [], []
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
                largest.append(biggest)
                fractions.append(float(m.mean()))
                frames_out.append(dict(reference=ref["name"], rule=rule, frame=int(f), fraction=float(m.mean()),
                                       largest=biggest))
            summary.append(dict(reference=ref["name"], rule=rule, molecules=int((frame_ids == frames[0]).sum()),
                                false_positive=float(np.mean(fractions)), largest_mean=float(np.mean(largest)),
                                largest_max=int(np.max(largest)),
                                frames_with_cluster_20=float(np.mean(np.array(largest) >= 20)),
                                persistence_5ps=float(np.mean(persistence)) if persistence else np.nan))
        print(ref["name"], flush=True)
    RESULTS.mkdir(exist_ok=True)
    pd.DataFrame(frames_out).to_csv(RESULTS/"false_nuclei_frames.csv", index=False)
    summary = pd.DataFrame(summary)
    summary.to_csv(RESULTS/"false_nuclei.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 400)
    print(summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
