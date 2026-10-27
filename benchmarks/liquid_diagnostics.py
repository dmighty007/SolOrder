"""Relaxation and aging diagnostics for every liquid reference in the benchmark.

For each liquid run (C1-C5, run D models, run E pressures, mW) and its 60
thermal frames 5 ps apart (300 ps):
  msd_300ps      mean squared oxygen displacement over the run (A^2), drift removed
  fs_300ps       self-intermediate scattering function F_s(k, t) of oxygens at
                 k = 2.25 A^-1 and t = 295 ps; near 0 means structurally
                 relaxed within the run, near 1 means arrested
  tau_ps         time at which F_s first falls below 1/e (NaN if it does not)
  energy_slope   linear trend of the potential energy, kJ/mol per molecule per ns
  <op>_drift     trend of the frame mean of an order parameter over the run, in
                 units of the Ih - liquid separation of that parameter at 220 K
                 (difference of medians in run C); |drift| << 1 means no aging
                 toward ice on the scale that matters for classification
  ice_fraction   CHILL+ ice or interfacial-ice fraction, first and last 10 frames

  python liquid_diagnostics.py   # results/liquid_diagnostics.csv
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
CONDITIONS = CACHE/"conditions"
K = 2.25  # A^-1
STRIDE_PS = 5.
OPS = ("q6_avg", "entropy_avg", "q_tet_avg", "cn_smooth")


def wrap(d, cell):
    frac = d @ np.linalg.inv(cell)
    return (frac - np.round(frac)) @ cell


def dynamics(x, cell, sites):
    o = x[:, ::sites].astype(float)
    steps = [wrap(o[t] - o[t-1], cell) for t in range(1, len(o))]
    disp = np.cumsum([np.zeros_like(o[0])] + steps, axis=0)  # unwrapped from 5 ps steps
    disp -= disp.mean(axis=1, keepdims=True)  # remove drift of the whole box
    r = np.linalg.norm(disp, axis=2)
    kr = K*np.where(r > 0, r, 1e-12)
    fs = np.mean(np.sin(kr)/kr, axis=1)
    below = np.flatnonzero(fs < np.exp(-1))
    return float(np.mean(r[-1]**2)), float(fs[-1]), float(below[0]*STRIDE_PS) if len(below) else np.nan


def trend(values):
    t = np.arange(len(values))*STRIDE_PS/1000.  # ns
    return float(np.polyfit(t, values, 1)[0]*t[-1])  # change over the run


def runs():
    out = []
    for path in sorted((CACHE/"replicas").glob("liquid_p*_220.npz")):
        out.append(dict(run=f"C{path.stem.split('_')[1][1:]}", model="TIP4P/Ice", T=220., P=0., path=path,
                        sites=4, table="prep", key=dict(run=f"C{path.stem.split('_')[1][1:]}", phase="liquid",
                                                         temperature=220.)))
    for path in sorted(CONDITIONS.glob("*liquid*.npz")) + sorted((CONDITIONS/"mw_inputs").glob("*liquid*.npz")):
        info = json.loads(path.with_suffix(".json").read_text())
        model = {"tip4pice": "TIP4P/Ice", "tip4p2005": "TIP4P/2005", "mW": "mW"}[info["model"]]
        table = "pressure" if info["dataset"] == "E" else "model"
        key = dict(phase="liquid", model=info["model"], temperature=info["temperature"],
                   pressure=info["pressure_gpa"], run=f"p{info['prep']}")
        out.append(dict(run=path.stem, model=model, T=info["temperature"], P=info["pressure_gpa"], path=path,
                        sites=1 if info["model"] == "mW" else 4, table=table, key=key))
    return out


def main():
    tables = {name: pd.read_pickle(CACHE/f"{name}.pkl") for name in ("prep", "model", "pressure")}
    prep = tables["prep"]
    ref = prep[(prep.kind == "thermal") & (prep.temperature == 220.)]
    scale = {op: ref.loc[ref.phase == "Ih", op].median() - ref.loc[ref.phase == "liquid", op].median() for op in OPS}
    rows = []
    for spec in runs():
        data = np.load(spec["path"])
        cell = data["cell"] if "cell" in data.files else np.diag(data["box"])
        msd, fs, tau = dynamics(data["x"], cell, spec["sites"])
        n = data["x"].shape[1]//spec["sites"]
        table = tables[spec["table"]]
        mask = table.kind == "thermal"
        for column, value in spec["key"].items():
            if column in table:
                mask &= table[column] == value
        part = table[mask]
        means = part.groupby("frame")[list(OPS)].mean()
        ice = part.groupby("frame").chill_plus.apply(lambda v: np.isin(v, (1, 2, 3)).mean())
        row = dict(run=spec["run"], model=spec["model"], T=spec["T"], P=spec["P"], molecules=n,
                   msd_300ps=msd, fs_300ps=fs, tau_ps=tau,
                   energy_slope=trend(data["energy"]/n)/(len(data["energy"])*STRIDE_PS/1000.),
                   ice_fraction_first=float(ice.iloc[:10].mean()), ice_fraction_last=float(ice.iloc[-10:].mean()),
                   frames=len(means))
        for op in OPS:
            row[f"{op}_drift"] = trend(means[op].to_numpy())/scale[op] if means[op].notna().all() else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    out.to_csv(RESULTS/"liquid_diagnostics.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.drop(columns=["molecules"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
