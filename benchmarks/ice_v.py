"""Why ice V is hard to tell from liquid at 0.5 GPa, and a two-parameter remedy.

Data: thermal frames of run E, ice V and liquid at 0.5 GPa, TIP4P/Ice 220 K,
preparations p1 and p2 (in cache/pressure.pkl) and a fresh preparation p3
(md_conditions.py E_V_p3_220_0.5GPa E_liquid_p3_220_0.5GPa) that no selection
step sees.

1. Structure of the errors. Each ice V molecule is characterized by its number
   of oxygens within 3.5 A in the time-averaged structure (thermal noise
   removed) and by whether it lies on a four-membered ring of the hydrogen-bond
   network. With the best single parameter (threshold fitted on p1), we report
   the misclassified fraction per class on p2.
2. Two-parameter remedy. For every pair of the scalar columns, a logistic
   regression on standardized values is fitted on frames 0-29 of p1 and scored
   (balanced accuracy) on frames 30-59 of p1; the best pair is refitted on all
   of p1 and tested once on p2 and on p3. The best single column is selected
   and tested the same way.

  python ice_v.py   # results/ice_v_classes.csv, ice_v_pairs.csv, ice_v_test.csv
"""

from itertools import combinations
from pathlib import Path
import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from compute import CONDITIONS, conditions_task
from metrics import fit_threshold, hits

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
EXCLUDE = {"chill", "chill_plus", "ptm", "n_hb", "cn"}  # labels and integer counts constant in these phases


def balanced(model, X, y):
    p = model.predict(X)
    return .5*(np.mean(p[y == 1] == 1) + np.mean(p[y == 0] == 0))


def threshold_score(train_a, train_b, test_a, test_b):
    t, below = fit_threshold(train_a, train_b)
    return .5*(hits(test_a, t, below).mean() + hits(test_b, t, not below).mean())


def load():
    pressure = pd.read_pickle(CACHE/"pressure.pkl")
    part = pressure[(pressure.kind == "thermal") & (pressure.run.isin(("p1", "p2")))
                    & (((pressure.phase == "V")) | ((pressure.phase == "liquid") & (pressure.pressure == .5)))]
    fresh = [delayed(conditions_task)(run, "thermal", f) for run in ("E_V_p3_220_0.5GPa", "E_liquid_p3_220_0.5GPa")
             for f in range(60)]
    p3 = pd.concat(Parallel(n_jobs=10)(fresh), ignore_index=True)
    return pd.concat([part, p3], ignore_index=True)


def structure_classes(prep):
    """Per ice V molecule: oxygens within 3.5 A in the time-averaged structure."""
    data = np.load(CONDITIONS/f"E_V_{prep}_220_0.5GPa.npz")
    cell, x = data["cell"], data["x"][:, ::4].astype(float)
    inv = np.linalg.inv(cell)
    d = x - x[0]
    frac = d @ inv
    mean = x[0] + ((frac - np.round(frac)) @ cell).mean(axis=0)
    diff = mean[:, None] - mean[None]
    frac = diff @ inv
    r = np.linalg.norm((frac - np.round(frac)) @ cell, axis=2)
    return (r < 3.5).sum(axis=1) - 1


def main():
    table = load()
    columns = [c for c in table.columns if c not in EXCLUDE and table[c].dtype.kind == "f"
               and c not in ("temperature", "pressure") and table[c].notna().all()]
    y = (table.phase == "V").to_numpy().astype(int)
    run, frame = table.run.to_numpy(), table.frame.to_numpy()
    train, valid = (run == "p1") & (frame < 30), (run == "p1") & (frame >= 30)
    p1, p2, p3 = run == "p1", run == "p2", run == "p3"

    # best single column, selected on p1 (first half -> second half), tested on p2 and p3
    singles = []
    for c in columns:
        v = table[c].to_numpy(float)
        singles.append(dict(column=c, valid=threshold_score(v[train & (y == 1)], v[train & (y == 0)],
                                                            v[valid & (y == 1)], v[valid & (y == 0)])))
    singles = pd.DataFrame(singles).sort_values("valid", ascending=False)
    best_single = singles.column.iloc[0]
    v = table[best_single].to_numpy(float)
    single_test = {name: threshold_score(v[p1 & (y == 1)], v[p1 & (y == 0)], v[m & (y == 1)], v[m & (y == 0)])
                   for name, m in (("p2", p2), ("p3", p3))}

    # every pair, logistic regression on standardized values
    X_all = table[columns].to_numpy(float)
    mu, sd = X_all[train].mean(axis=0), X_all[train].std(axis=0)
    Z = (X_all - mu)/np.where(sd > 0, sd, 1)
    idx = {c: i for i, c in enumerate(columns)}

    def pair_score(a, b):
        cols = [idx[a], idx[b]]
        model = LogisticRegression(class_weight="balanced", max_iter=200).fit(Z[train][:, cols], y[train])
        return dict(a=a, b=b, valid=balanced(model, Z[valid][:, cols], y[valid]))

    pairs = pd.DataFrame(Parallel(n_jobs=10)(delayed(pair_score)(a, b) for a, b in combinations(columns, 2)))
    pairs = pairs.sort_values("valid", ascending=False)
    a, b = pairs.a.iloc[0], pairs.b.iloc[0]
    cols = [idx[a], idx[b]]
    model = LogisticRegression(class_weight="balanced", max_iter=200).fit(Z[p1][:, cols], y[p1])
    pair_test = {name: balanced(model, Z[m][:, cols], y[m]) for name, m in (("p2", p2), ("p3", p3))}

    # where the errors are: misclassified fraction of ice V molecules per structural class (p2)
    rows = []
    for prep, mask in (("p2", p2), ("p3", p3)):
        cn = structure_classes(prep)
        part = table[mask & (y == 1)]
        cls = cn[part.mol.to_numpy()]
        four = part.ring_f4.to_numpy() > 0
        t, below = fit_threshold(v[p1 & (y == 1)], v[p1 & (y == 0)])
        single_miss = ~hits(part[best_single].to_numpy(float), t, below)
        pair_miss = model.predict(Z[mask & (y == 1)][:, cols]) == 0
        for label, groups in (("neighbors_3.5A", cls), ("on_4_ring", four)):
            for g in np.unique(groups):
                m = groups == g
                rows.append(dict(prep=prep, feature=label, value=str(g), fraction_of_ice_V=float(m.mean()),
                                 single_error=float(single_miss[m].mean()), pair_error=float(pair_miss[m].mean())))
    classes = pd.DataFrame(rows)
    liquid_cn = [dict(prep=p, mean_cn_liquid=float(table.loc[(run == p) & (y == 0), "cn_smooth"].mean()),
                      mean_cn_ice_V=float(table.loc[(run == p) & (y == 1), "cn_smooth"].mean())) for p in ("p1", "p2", "p3")]

    RESULTS.mkdir(exist_ok=True)
    classes.to_csv(RESULTS/"ice_v_classes.csv", index=False)
    pairs.to_csv(RESULTS/"ice_v_pairs.csv", index=False)
    test = pd.DataFrame([dict(model=f"single: {best_single}", **single_test, valid=singles.valid.iloc[0]),
                         dict(model=f"pair: {a} + {b}", **pair_test, valid=pairs.valid.iloc[0])])
    test.to_csv(RESULTS/"ice_v_test.csv", index=False)
    pd.set_option("display.width", 200)
    print(singles.head(6).round(3).to_string(index=False))
    print(pairs.head(10).round(3).to_string(index=False))
    print(test.round(3).to_string(index=False))
    print(classes.round(3).to_string(index=False))
    print(pd.DataFrame(liquid_cn).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
