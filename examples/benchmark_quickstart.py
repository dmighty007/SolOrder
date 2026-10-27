"""Reproduce one benchmark score end to end, from coordinates to the published number.

Computes averaged q6 for thermal ice Ih and supercooled liquid from preparation
C1 (TIP4P/Ice, 220 K, 60 frames each), the oracle score (1 + KS D)/2 for Ih
versus liquid, and compares it with the value in
benchmarks/results/separability_grouped.csv. Takes about a minute.

  python examples/benchmark_quickstart.py [--data benchmarks/cache/replicas]

With the archival release, point --data at coordinates/run_C_preparations.
"""

import argparse
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from SolOrder import OrderParameter

ROOT = Path(__file__).resolve().parents[1]


def averaged_q6(data, phase):
    u = mda.Universe(str(data/f"{phase}_p1.gro"))
    frames = np.load(data/f"{phase}_p1_220.npz")
    values = []
    for x in frames["x"]:
        u.atoms.positions = x
        u.dimensions = [*frames["box"], 90., 90., 90.]
        op = OrderParameter(u, cutoff=5., selection="name OW", hydrogen_selection="name HW1 HW2")
        values.append(op.steinhardt(6, averaged=True, cutoff=3.5))
    return np.concatenate(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, default=ROOT/"benchmarks/cache/replicas")
    args = parser.parse_args()
    ice, liquid = averaged_q6(args.data, "Ih"), averaged_q6(args.data, "liquid")
    oracle = (1 + ks_2samp(ice, liquid, method="asymp").statistic)/2
    print(f"averaged q6, Ih versus liquid, C1 thermal 220 K: oracle score {oracle:.4f}")
    candidates = (ROOT/"benchmarks/results", ROOT.parent/"results")  # repository, or code/ inside the release
    published = next((d/"separability_grouped.csv" for d in candidates if (d/"separability_grouped.csv").exists()), None)
    if published:
        table = pd.read_csv(published)
        row = table[(table.column == "q6_avg") & (table.a == "Ih") & (table.b == "liquid") & (table.kind == "thermal")]
        if len(row) and "oracle_C1" in row:
            print(f"published (separability_grouped.csv, C1): {float(row.oracle_C1.iloc[0]):.4f}")


if __name__ == "__main__":
    main()
