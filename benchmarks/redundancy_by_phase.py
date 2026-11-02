"""Redundancy within single phases, compared with the pooled ice/liquid redundancy (F2).

F2 clusters parameters by Spearman |rho| pooled over thermal Ih, Ic and liquid,
so part of every correlation is the ice/liquid contrast itself. Here |rho| is
computed within one phase at a time, thermal frames at 220 K: Ih, Ic and
liquid from runs C1-C5, and ices V, VI and VII and liquid at 1.0 GPa from run
E. For every pair grouped together in the pooled analysis (|rho| > 0.8) we
report its within-phase |rho|; we also count the groups that average-linkage
clustering at |rho| > 0.8 gives in each phase.

  python redundancy_by_phase.py   # results/redundancy_by_phase.csv, redundancy_groups.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
PER_PHASE = 20000
CUT = .8
rng = np.random.default_rng(20261027)


def groups(rho):
    distance = np.clip(1. - np.abs(rho), 0, None)
    np.fill_diagonal(distance, 0.)
    return fcluster(linkage(squareform((distance+distance.T)/2), method="average"), 1-CUT, criterion="distance")


def spearman(table, columns):
    index = rng.choice(len(table), size=min(PER_PHASE, len(table)), replace=False)
    sample = table.iloc[index][columns]
    usable = [c for c in columns if sample[c].notna().all() and sample[c].nunique() > 1]
    rho = spearmanr(sample[usable].to_numpy(float)).statistic
    return pd.DataFrame(rho, index=usable, columns=usable)


def main():
    pooled = pd.read_csv(RESULTS/"redundancy_spearman.csv", index_col=0)
    pooled = pooled.loc[pooled.notna().all(), pooled.notna().all()]
    columns = pooled.columns.tolist()
    prep = pd.read_pickle(CACHE/"prep.pkl")
    prep = prep[(prep.kind == "thermal") & (prep.temperature == 220.)]
    pressure = pd.read_pickle(CACHE/"pressure.pkl")
    pressure = pressure[pressure.kind == "thermal"]
    phases = {name: prep[prep.phase == name] for name in ("Ih", "Ic", "liquid")}
    phases.update({name: pressure[pressure.phase == name] for name in ("V", "VI", "VII")})
    phases["liquid 1.0 GPa"] = pressure[(pressure.phase == "liquid") & (pressure.pressure == 1.)]
    within = {name: spearman(table, columns) for name, table in phases.items()}
    pooled_groups = dict(zip(columns, groups(pooled.to_numpy())))
    rows = []
    for i, a in enumerate(columns):
        for b in columns[i+1:]:
            row = dict(a=a, b=b, pooled=abs(pooled.loc[a, b]), same_pooled_group=pooled_groups[a] == pooled_groups[b])
            for name, rho in within.items():
                row[name] = abs(rho.loc[a, b]) if a in rho.index and b in rho.index else np.nan
            rows.append(row)
    pairs = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    pairs.to_csv(RESULTS/"redundancy_by_phase.csv", index=False)
    summary = [dict(phase="pooled (F2)", columns=len(columns), groups=len(set(pooled_groups.values())))]
    for name, rho in within.items():
        summary.append(dict(phase=name, columns=len(rho), groups=len(set(groups(rho.to_numpy())))))
    pd.DataFrame(summary).to_csv(RESULTS/"redundancy_groups.csv", index=False)
    grouped = pairs[pairs.same_pooled_group]
    print(pd.DataFrame(summary).to_string(index=False))
    print(f"\n{len(grouped)} pairs share a pooled group; within-phase |rho| > {CUT}:")
    print(grouped[list(within)].gt(CUT).mean().round(2).to_string())
    print("\nmedian within-phase |rho| of those pairs:")
    print(grouped[list(within)].median().round(2).to_string())


if __name__ == "__main__":
    main()
