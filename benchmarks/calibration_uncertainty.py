"""Uncertainty of the planar calibration itself, propagated into recovered nucleus sizes.

The calibration (f, delta) comes from nine run F basal slabs: three
polytypes times three independent preparations. Here f and delta are
re-estimated jointly on resampled sets of slabs and every re-estimate is
pushed through eq. 8 for every run F and run G seed:

  bootstrap   2000 resamples of the nine slabs with replacement
  jackknife   leave one preparation (three slabs) out

Reported separately from prediction errors: the spread of f and delta, and
the relative spread of the recovered size per seed (95% interval half-width
over its point estimate). The frozen calibration used for run G
(results/nucleus_size_calibration.csv) is not changed; the point estimate
here is checked against it.

  python calibration_uncertainty.py   # results/calibration_uncertainty.csv, calibration_uncertainty_sizes.csv
"""

import numpy as np
import pandas as pd

from nucleus_size import (CACHE, DEEP, N_DENSITY, OBSERVABLES, RESULTS, inverse, observed_size, per_molecule,
                          slab_area)

RESAMPLES, SEED = 2000, 2031


def slab_sums(slabs, method, observable):
    """Per slab: deep hits and counts (for f), crystal molecules and hits per frame and face area (for Gamma)."""
    rows = []
    for name, g in slabs.groupby("name"):
        p = per_molecule(g, method, observable).astype(float)
        deep = g.depth.to_numpy() > DEEP
        frames = g.frame.nunique()
        rows.append(dict(name=name, prep=int(g.prep.iloc[0]), deep_hits=p[deep].sum(), deep_n=deep.sum(),
                         crystal=len(p)/frames, hits=p.sum()/frames, area=slab_area(name)))
    return pd.DataFrame(rows)


def estimate(sums):
    """f pooled over molecules, delta = mean over slabs of Gamma/n (as nucleus_size.calibrate)."""
    f = sums.deep_hits.sum()/sums.deep_n.sum()
    gamma = (sums.crystal - sums.hits/f)/2/sums.area
    return f, float(np.mean(gamma/N_DENSITY))


def main():
    calibration = pd.read_csv(RESULTS/"nucleus_size_calibration.csv")
    run_f = pd.read_pickle(CACHE/"interface.pkl")
    slabs = run_f[run_f.geometry == "slab"]
    seeds = {"F": (run_f[run_f.geometry == "seed"], pd.read_pickle(CACHE/"interface_clusters.pkl"))}
    if (CACHE/"confirmation.pkl").exists():
        run_g = pd.read_pickle(CACHE/"confirmation.pkl")
        seeds["G"] = (run_g[run_g.geometry == "seed"], pd.read_pickle(CACHE/"confirmation_clusters.pkl"))
    rng = np.random.default_rng(SEED)
    summary, sizes = [], []
    for observable, methods in OBSERVABLES.items():
        for method in methods:
            sums = slab_sums(slabs, method, observable)
            f0, d0 = estimate(sums)
            frozen = calibration[(calibration.observable == observable) & (calibration.method == method)].iloc[0]
            assert abs(f0 - frozen.bulk_recall) < 1e-9 and abs(d0 - frozen["count"]) < 1e-9
            boot = np.array([estimate(sums.iloc[rng.integers(len(sums), size=len(sums))]) for _ in range(RESAMPLES)])
            jack = np.array([estimate(sums[sums.prep != k]) for k in sorted(sums.prep.unique())])
            summary.append(dict(observable=observable, method=method, f=f0, delta=d0,
                                f_lo=np.quantile(boot[:, 0], .025), f_hi=np.quantile(boot[:, 0], .975),
                                delta_lo=np.quantile(boot[:, 1], .025), delta_hi=np.quantile(boot[:, 1], .975),
                                delta_sd=boot[:, 1].std(), f_delta_corr=np.corrcoef(boot.T)[0, 1],
                                jack_delta_min=jack[:, 1].min(), jack_delta_max=jack[:, 1].max(),
                                jack_f_min=jack[:, 0].min(), jack_f_max=jack[:, 0].max()))
            for run, (table, clusters) in seeds.items():
                for name, sub in table.groupby("name"):
                    obs = observed_size(sub, clusters, method, observable)
                    point = inverse(obs, f0, d0)
                    draws = inverse(obs, boot[:, 0], boot[:, 1])
                    jd = inverse(obs, jack[:, 0], jack[:, 1])
                    lo, hi = np.quantile(draws, (.025, .975))
                    n = int(sub["size"].iloc[0])
                    sizes.append(dict(run=run, observable=observable, method=method, name=name, size=n,
                                      inverse=point, lo=lo, hi=hi, half_width=(hi - lo)/2/point,
                                      jack_spread=(jd.max() - jd.min())/2/point,
                                      error=(point - n)/n, covered=bool(lo <= n <= hi)))
    summary, sizes = pd.DataFrame(summary), pd.DataFrame(sizes)
    summary.to_csv(RESULTS/"calibration_uncertainty.csv", index=False)
    sizes.to_csv(RESULTS/"calibration_uncertainty_sizes.csv", index=False)
    pd.set_option("display.width", 220)
    print(summary.round(3).to_string(index=False))
    print(sizes.groupby(["observable", "method", "run", "size"])[["half_width", "jack_spread", "error"]]
          .agg(lambda x: np.abs(x).median()).round(3).to_string())
    print(sizes.groupby(["observable", "method", "run"]).covered.mean().round(2).to_string())


if __name__ == "__main__":
    main()
