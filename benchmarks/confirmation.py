"""Run G: test the frozen nucleus-size calibration on fresh seeds, and measure facet and time dependence.

Follows confirmation_protocol.md. Reads the calibration fixed before run G
(results/nucleus_size_calibration.csv, from the run F basal slabs) and never
refits it.

  seeds    eq. 8 (inverse, count dead layer) on every G_seed system: errors
           against the constructed and the intact size
  facets   f and the count dead layer of the prism slabs (G_prism) next to the
           frozen basal values and the fourth basal preparation (first ten
           frames of G_long)
  time     f and the count dead layer of the G_long basal slabs at every saved
           time, 2-500 ps after release

  python confirmation.py   # results/confirmation_seeds.csv, confirmation_slabs.csv, confirmation_time.csv
"""

import numpy as np
import pandas as pd

from nucleus_size import CACHE, PRIMARY, RESULTS, apply, calibrate

STANDARD_FRAMES = 10  # frames at 2-20 ps, the standard protocol


def main():
    table = pd.read_pickle(CACHE/"confirmation.pkl")
    clusters = pd.read_pickle(CACHE/"confirmation_clusters.pkl")
    calibration = pd.read_csv(RESULTS/"nucleus_size_calibration.csv")
    seeds = table[table.geometry == "seed"]
    out = apply(seeds, clusters, calibration, definitions=(PRIMARY,))
    out.to_csv(RESULTS/"confirmation_seeds.csv", index=False)

    slabs = table[table.geometry == "slab"]
    long_ = slabs[slabs.name.str.startswith("G_long")]
    standard = pd.concat([slabs[slabs.name.str.startswith("G_prism")], long_[long_.frame < STANDARD_FRAMES]])
    facet_rows = []
    for c in calibration.itertuples():
        facet_rows.append(dict(facet="basal (run F, frozen)", kind="all", observable=c.observable, method=c.method,
                               bulk_recall=c.bulk_recall, count=c.count, count_min=c.count_min, count_max=c.count_max))
        for facet, part in standard.groupby("facet"):
            label = "basal (p4)" if facet == "basal" else facet
            for kind, sub in [("all", part), *part.groupby("kind")]:
                cal = calibrate(sub, c.method, c.observable)
                facet_rows.append(dict(facet=label, kind=kind, observable=c.observable, method=c.method,
                                       bulk_recall=cal["bulk_recall"], count=cal["count"],
                                       count_min=cal["count_min"], count_max=cal["count_max"]))
    facets = pd.DataFrame(facet_rows)
    facets.to_csv(RESULTS/"confirmation_slabs.csv", index=False)

    time_rows = []
    for c in calibration.itertuples():
        for (kind, t), sub in long_.groupby(["kind", "time_ps"]):
            cal = calibrate(sub, c.method, c.observable)
            time_rows.append(dict(kind=kind, time_ps=t, observable=c.observable, method=c.method,
                                  bulk_recall=cal["bulk_recall"], count=cal["count"],
                                  intact=float(sub.intact.mean())))
    time = pd.DataFrame(time_rows)
    time.to_csv(RESULTS/"confirmation_time.csv", index=False)

    pd.set_option("display.width", 220)
    main_ = out[out.observable.isin(("correct", "cluster")) & out.method.isin(("chill_plus", "ptm"))]
    print(main_[["observable", "method", "name", "size", "observed", "inverse", "inverse_error",
                 "inverse_error_intact"]].round(3).to_string(index=False))
    print(main_.groupby(["observable", "method"])[["inverse_error", "inverse_error_intact"]]
          .agg(lambda x: np.abs(x).max()).round(3).to_string())
    print(facets[facets.observable == "correct"].round(2).to_string(index=False))
    t = time[(time.observable == "correct") & time.method.isin(("chill_plus", "ptm"))]
    print(t.pivot_table(index="time_ps", columns=["method", "kind"], values="count").round(2).to_string())


if __name__ == "__main__":
    main()
