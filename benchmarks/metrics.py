"""Separability, redundancy and label-accuracy metrics for the pilot.

Separability of two populations under one scalar parameter is the two-sample
Kolmogorov-Smirnov statistic D. It needs no binning, and (1 + D)/2 is the
balanced accuracy of the best single threshold chosen on the same data: an
*oracle* score, reported as ``oracle``, not a predictive accuracy. Confidence
intervals come from a bootstrap over frames (molecules within a frame are
correlated); frames of one short run are not independent, so they are too narrow.

Held-out scores fit the threshold and its direction on one run and apply it to
another (run A: seed test; run B: replica runs). NaN values count as errors
there, so missing values cannot improve a held-out score.
"""

import argparse

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, spearmanr

from registry import LABELS, SCALARS

ROOT = Path(__file__).resolve().parent
CACHE, RESULTS = ROOT/"cache", ROOT/"results"
BOOTSTRAP = 200
BLOCK = 5  # replica frames per bootstrap block (25 ps; see replica_autocorrelation.csv)
GRID = 4001  # pooled quantiles tried as thresholds
REPLICA_HALF = 30  # first 30 replica frames train, last 30 test (within-run temporal split)
LIQUID_MARGIN = 8.  # A beyond the seed radius: "far liquid" in seeded systems
rng = np.random.default_rng(20261026)


COLUMN_ENTRY = {
    "vor_density": "voronoi", "vor_asphericity": "voronoi",
    "entropy": "entropy", "entropy_avg": "entropy",
    **{f"ah{k}": "angle_hist" for k in (1, 2, 3)},
    **{f"A{k}": "neighborhood" for k in (1, 2, 3)},
    "ring_mean": "rings", **{f"ring_f{k}": "rings" for k in (4, 5, 6, 7, 8)},
    **{c: "steinhardt" for c in ("q3", "q3_avg", "q4", "q4_avg", "q6", "q6_avg", "w4_hat", "w6_hat")},
}


def column_info(column):
    name = COLUMN_ENTRY.get(column, column)
    entry = next(e for e in SCALARS if e.name == name)
    return dict(column=column, entry=entry.name, group=entry.group,
                needs_h=entry.needs_h, smooth=entry.smooth, source=entry.source)


def ks(a, b):
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if not len(a) or not len(b):
        return np.nan
    return ks_2samp(a, b, method="asymp").statistic


def by_frame(values, frames):
    return [values[frames == f] for f in np.unique(frames)]


def bootstrap_ks(a_frames, b_frames):
    stats = []
    for _ in range(BOOTSTRAP):
        a = np.concatenate([a_frames[i] for i in rng.integers(len(a_frames), size=len(a_frames))])
        b = np.concatenate([b_frames[i] for i in rng.integers(len(b_frames), size=len(b_frames))])
        stats.append(ks(a, b))
    return np.nanpercentile(stats, [2.5, 97.5])


def ideal_separability(table, columns):
    rows = []
    phases = list(dict.fromkeys(table.phase))
    for column in columns:
        values = {p: table.loc[table.phase == p, column].to_numpy(float) for p in phases}
        for a, b in combinations(phases, 2):
            d = ks(values[a], values[b])
            rows.append(dict(column=column, a=a, b=b, ks=d, oracle=(1+d)/2))
    return pd.DataFrame(rows)


def thermal_groups(table):
    """Reference populations from the seed test, per frame kind."""
    seeded = table.system.str.startswith("seed")
    far_liquid = seeded & ~table.seed & (table.r_seed > table.R + LIQUID_MARGIN)
    isd = table.system == "bulk_Isd"
    return {
        "Ih": table.system == "bulk_Ih",
        "Ic": table.system == "bulk_Ic",
        "liquid": far_liquid,
        "Isd_hex": isd & (table["stack"] == 1),
        "Isd_cubic": isd & (table["stack"] == 2),
    }


PAIRS = [("Ih", "liquid"), ("Ic", "liquid"), ("Ih", "Ic"), ("Isd_hex", "Isd_cubic")]


def thermal_separability(table, columns):
    rows = []
    for kind in ("thermal", "quenched"):
        part = table[table.kind == kind]
        groups = thermal_groups(part)
        # Frames from different systems are distinct sampling units.
        unit = (part.system + ":" + part.frame.astype(str)).to_numpy()
        for column in columns:
            values = part[column].to_numpy(float)
            for a, b in PAIRS:
                av, bv = values[groups[a]], values[groups[b]]
                d = ks(av, bv)
                lo, hi = bootstrap_ks(by_frame(av, unit[groups[a]]), by_frame(bv, unit[groups[b]]))
                rows.append(dict(kind=kind, column=column, a=a, b=b, ks=d, ks_lo=lo, ks_hi=hi,
                                 oracle=(1+d)/2, n_a=len(av), n_b=len(bv),
                                 nan_fraction=float(np.mean(~np.isfinite(np.r_[av, bv])))))
    return pd.DataFrame(rows)


def replica_groups(table):
    """Reference populations from the replica runs (run B), per frame kind."""
    cold, warm = table.temperature == 220., table.temperature == 250.
    isd = table.system == "bulk_Isd"
    return {
        "Ih": cold & (table.system == "bulk_Ih"),
        "Ic": cold & (table.system == "bulk_Ic"),
        "liquid": table.system == "liquid",
        "Isd_hex": isd & (table["stack"] == 1),
        "Isd_cubic": isd & (table["stack"] == 2),
        "Ih_250": warm & (table.system == "bulk_Ih"),
        "Ic_250": warm & (table.system == "bulk_Ic"),
    }


def fit_threshold(a, b):
    """Threshold t and direction (True: a is x <= t) maximizing balanced accuracy."""
    a, b = np.sort(a[np.isfinite(a)]), np.sort(b[np.isfinite(b)])
    if not len(a) or not len(b):
        return np.nan, True
    grid = np.unique(np.quantile(np.r_[a, b], np.linspace(0., 1., GRID)))
    diff = np.searchsorted(a, grid, side="right")/len(a) - np.searchsorted(b, grid, side="right")/len(b)
    best = np.argmax(np.abs(diff))
    return grid[best], bool(diff[best] > 0)


def hits(values, t, below):
    """Correct calls for one population; NaN is never correct."""
    return values <= t if below else values > t


def frame_hits(values, frames, t, below):
    """Correct calls and population size per frame, in frame order."""
    hit = hits(values, t, below)
    ids = np.unique(frames)
    return np.array([hit[frames == f].sum() for f in ids]), np.array([(frames == f).sum() for f in ids])


def block_interval(a_counts, b_counts):
    """Block bootstrap over contiguous test frames; both populations share the blocks."""
    (a_hit, a_n), (b_hit, b_n) = a_counts, b_counts
    blocks = np.arange(len(a_hit)) // BLOCK
    members = [np.flatnonzero(blocks == k) for k in range(blocks.max() + 1)]
    scores = []
    for _ in range(BOOTSTRAP):
        mask = np.concatenate([members[k] for k in rng.integers(len(members), size=len(members))])
        scores.append(.5*(a_hit[mask].sum()/a_n[mask].sum() + b_hit[mask].sum()/b_n[mask].sum()))
    return np.percentile(scores, [2.5, 97.5])


def transfer(train_a, train_b, test_a, test_b, frames_a=None, frames_b=None):
    """Balanced accuracy on (test_a, test_b) of the threshold fitted on (train_a, train_b)."""
    t, below = fit_threshold(train_a, train_b)
    if not np.isfinite(t):
        return dict(score=np.nan, lo=np.nan, hi=np.nan)
    score = .5*(hits(test_a, t, below).mean() + hits(test_b, t, not below).mean())
    lo = hi = np.nan
    if frames_a is not None:
        lo, hi = block_interval(frame_hits(test_a, frames_a, t, below),
                                frame_hits(test_b, frames_b, t, not below))
    return dict(score=score, lo=lo, hi=hi)


def heldout_separability(seedtest, replica, columns):
    """Oracle versus held-out single-threshold scores across independent runs."""
    rows = []
    for kind in ("thermal", "quenched"):
        run_a, run_b = seedtest[seedtest.kind == kind], replica[replica.kind == kind]
        groups_a, groups_b = thermal_groups(run_a), replica_groups(run_b)
        frames_b = run_b.frame.to_numpy()
        tasks = PAIRS + [("Ih_250", "Ic_250")]
        for column in columns:
            va, vb = run_a[column].to_numpy(float), run_b[column].to_numpy(float)
            for a, b in tasks:
                row = dict(kind=kind, column=column, a=a, b=b)
                ba, bb = vb[groups_b[a]], vb[groups_b[b]]
                fa, fb = frames_b[groups_b[a]], frames_b[groups_b[b]]
                row["oracle_B"] = (1 + ks(ba, bb))/2
                if a in groups_a:
                    aa, ab = va[groups_a[a]], va[groups_a[b]]
                    row["oracle_A"] = (1 + ks(aa, ab))/2
                    ab_ = transfer(aa, ab, ba, bb, fa, fb)
                    row.update(A_to_B=ab_["score"], A_to_B_lo=ab_["lo"], A_to_B_hi=ab_["hi"])
                    row["B_to_A"] = transfer(ba, bb, aa, ab)["score"]
                else:  # 250 K: thresholds fitted at 220 K in run B
                    ta, tb = vb[groups_b["Ih"]], vb[groups_b["Ic"]]
                    t = transfer(ta, tb, ba, bb, fa, fb)
                    row.update(T220_to_T250=t["score"], T220_to_T250_lo=t["lo"], T220_to_T250_hi=t["hi"])
                early_a, early_b = fa < REPLICA_HALF, fb < REPLICA_HALF
                row["B_early_to_late"] = transfer(ba[early_a], bb[early_b], ba[~early_a], bb[~early_b])["score"]
                row["nan_fraction_B"] = float(np.mean(~np.isfinite(np.r_[ba, bb])))
                rows.append(row)
    return pd.DataFrame(rows)


def grouped_separability(seedtest, prep, columns):
    """Scores across independent preparations C1..Ck (prep k pairs ice k with liquid k).

    Per column, pair and frame kind: oracle score in every preparation (mean,
    SD, min); leave-one-preparation-out held-out score (fit on the pooled
    other preparations, test on the left-out one); and thresholds fitted on
    run A (seed test) tested on each preparation. For Ih/Ic at 250 K the
    transfer fits on the pooled 220 K preparations instead of run A.
    """
    runs = sorted(prep.run.unique())
    rows = []
    for kind in ("thermal", "quenched"):
        run_a = seedtest[seedtest.kind == kind]
        groups_a = thermal_groups(run_a)
        parts = {r: prep[(prep.kind == kind) & (prep.run == r)] for r in runs}
        groups = {r: replica_groups(part) for r, part in parts.items()}
        for column in columns:
            va = run_a[column].to_numpy(float)
            values = {r: part[column].to_numpy(float) for r, part in parts.items()}
            for a, b in PAIRS + [("Ih_250", "Ic_250")]:
                pops = {r: (values[r][groups[r][a]], values[r][groups[r][b]]) for r in runs
                        if groups[r][a].any() and groups[r][b].any()}
                if len(pops) < 2:
                    continue
                oracle = {r: (1 + ks(*pops[r]))/2 for r in pops}
                loo = {}
                for r in pops:
                    train = [x for x in pops if x != r]
                    loo[r] = transfer(np.concatenate([pops[x][0] for x in train]),
                                      np.concatenate([pops[x][1] for x in train]), *pops[r])["score"]
                if a in groups_a:
                    fit = va[groups_a[a]], va[groups_a[b]]
                else:  # 250 K: fit on pooled 220 K Ih and Ic preparations
                    fit = (np.concatenate([values[r][groups[r]["Ih"]] for r in runs]),
                           np.concatenate([values[r][groups[r]["Ic"]] for r in runs]))
                moved = {r: transfer(*fit, *pops[r])["score"] for r in pops}
                o, l, m = (np.array(list(d.values())) for d in (oracle, loo, moved))
                rows.append(dict(kind=kind, column=column, a=a, b=b, n_preps=len(pops),
                                 oracle_mean=o.mean(), oracle_sd=o.std(ddof=1), oracle_min=o.min(),
                                 loo_mean=l.mean(), loo_min=l.min(), loo_gap_max=np.max(o-l),
                                 transfer_mean=m.mean(), transfer_min=m.min(), transfer_gap_max=np.max(o-m),
                                 transfer_from="A" if a in groups_a else "C 220 K",
                                 **{f"oracle_{r}": oracle[r] for r in pops}))
    return pd.DataFrame(rows)


MODEL_PREPS = ("C1", "C2", "C3")  # run C preparations that run D starts from


def model_separability(prep, model, columns):
    """Second water model (run D): thresholds fitted on TIP4P/Ice applied to TIP4P/2005.

    Conditions, each pooled over preparations p1-p3: C (run C, TIP4P/Ice, NVT in
    the constructed box), ice_npt (TIP4P/Ice, run D protocol), t2005_220 and
    t2005_200 (TIP4P/2005 at the same temperature and at the same undercooling),
    mW (oxygen-only, 220 K: about the TIP4P/Ice undercooling; hydrogen-dependent
    columns are NaN and left out).
    Transfers: C -> ice_npt is the protocol effect alone, ice_npt -> t2005_* the
    model effect alone, C -> t2005_* both. *_min is the worst preparation.
    """
    conditions = {
        "C": prep[prep.run.isin(MODEL_PREPS) & (prep.temperature == 220.)].assign(
            run=lambda t: "p" + t.run.str[1:]),
        "ice_npt": model[(model.model == "tip4pice") & (model.temperature == 220.)],
        "t2005_220": model[(model.model == "tip4p2005") & (model.temperature == 220.)],
        "t2005_200": model[(model.model == "tip4p2005") & (model.temperature == 200.)],
        "mW": model[model.model == "mW"],
    }
    moves = [("C", "ice_npt"), ("ice_npt", "t2005_220"), ("ice_npt", "t2005_200"),
             ("C", "t2005_220"), ("C", "t2005_200"), ("ice_npt", "mW"), ("C", "mW")]
    rows = []
    for kind in ("thermal", "quenched"):
        parts = {c: t[t.kind == kind] for c, t in conditions.items()}
        groups = {}
        for c, part in parts.items():
            g = replica_groups(part.assign(temperature=220.))  # one temperature per condition
            groups[c] = {k: v for k, v in g.items() if v.any()}
        for column in columns:
            values = {c: part[column].to_numpy(float) for c, part in parts.items()}
            for a, b in PAIRS:
                pops = {c: (values[c][groups[c][a]], values[c][groups[c][b]]) for c in parts
                        if a in groups[c] and b in groups[c]}
                row = dict(kind=kind, column=column, a=a, b=b)
                for c, (x, y) in pops.items():
                    row[f"oracle_{c}"] = (1 + ks(x, y))/2
                for src, dst in moves:
                    if src not in pops or dst not in pops or not np.isfinite(pops[dst][0]).any():
                        continue  # NaN throughout: hydrogen-dependent column on mW
                    row[f"{src}_to_{dst}"] = transfer(*pops[src], *pops[dst])["score"]
                    preps = parts[dst].run.to_numpy()
                    ga, gb = groups[dst][a].to_numpy(), groups[dst][b].to_numpy()
                    row[f"{src}_to_{dst}_min"] = min(
                        transfer(*pops[src], values[dst][ga & (preps == r)], values[dst][gb & (preps == r)])["score"]
                        for r in np.unique(preps))
                rows.append(row)
    return pd.DataFrame(rows)


HP_CRYSTALS = ("Ih", "Ic", "II", "III", "V", "VI", "VII", "VIII")
HP_LIQUIDS = {"liquid_0": 1e-4, "liquid_0.3": .3, "liquid_0.5": .5, "liquid_1.0": 1., "liquid_3.0": 3.}


def pressure_populations(prep, pressure, kind):
    """Thermal phases at 220 K, two preparations each (run C p1-p2 at ambient, run E)."""
    ambient = prep[(prep.kind == kind) & prep.run.isin(("C1", "C2")) & (prep.temperature == 220.)]
    ambient = ambient.assign(run="p" + ambient.run.str[1:], pressure=1e-4)
    table = pd.concat([ambient[ambient.phase.isin(("Ih", "Ic", "liquid"))], pressure[pressure.kind == kind]])
    names = np.where(table.phase == "liquid", "liquid_" + table.pressure.map(
        {v: k.split("_")[1] for k, v in HP_LIQUIDS.items()}).fillna(""), table.phase)
    return table.assign(population=names)


def pressure_separability(prep, pressure, ideal, columns):
    """Thermal high-pressure phases (run E): thermal oracle and held-out scores versus ideal.

    Pairs: every crystal pair, every crystal against every liquid, and the
    liquids against each other. Held-out: fit on one preparation, test on the
    other, both ways, mean (heldout) and worse (heldout_min). ideal: oracle
    score of the jittered GenIce structures (separability_ideal.csv).
    """
    ideal_score = {}
    for r in ideal.itertuples():
        ideal_score[r.column, r.a, r.b] = ideal_score[r.column, r.b, r.a] = r.oracle
    pairs = list(combinations(HP_CRYSTALS, 2)) + [(c, l) for c in HP_CRYSTALS for l in HP_LIQUIDS]
    pairs += list(combinations(HP_LIQUIDS, 2))
    rows = []
    for kind in ("thermal", "quenched"):
        table = pressure_populations(prep, pressure, kind)
        population, preps = table.population.to_numpy(), table.run.to_numpy()
        for column in columns:
            values = table[column].to_numpy(float)
            for a, b in pairs:
                ma, mb = population == a, population == b
                row = dict(kind=kind, column=column, a=a, b=b, oracle=(1 + ks(values[ma], values[mb]))/2,
                           ideal=ideal_score.get((column, a, b), np.nan))
                p1, p2 = (preps == "p1"), (preps == "p2")
                scores = [transfer(values[ma & tr], values[mb & tr], values[ma & te], values[mb & te])["score"]
                          for tr, te in ((p1, p2), (p2, p1))]
                row.update(heldout=np.mean(scores), heldout_min=np.min(scores),
                           nan_fraction=float(np.mean(~np.isfinite(np.r_[values[ma], values[mb]]))))
                rows.append(row)
    return pd.DataFrame(rows)


def replica_autocorrelation(replica, columns):
    """Lag correlation of per-frame means in run B (frames 5 ps apart), thermal frames."""
    rows = []
    part = replica[replica.kind == "thermal"]
    for (system, temperature), table in part.groupby(["system", "temperature"]):
        means = table.groupby("frame")[columns].mean()
        for column in columns:
            x = means[column].to_numpy() - means[column].mean()
            row = dict(system=system, temperature=temperature, column=column)
            for lag in (1, 2, 5):
                row[f"lag{lag}"] = float(np.dot(x[:-lag], x[lag:])/np.dot(x, x)) if np.dot(x, x) > 0 else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def redundancy(table, columns, per_group=20000):
    part = table[table.kind == "thermal"]
    groups = thermal_groups(part)
    picks = []
    for name in ("Ih", "Ic", "liquid"):
        index = np.flatnonzero(groups[name].to_numpy())
        picks.append(rng.choice(index, size=min(per_group, len(index)), replace=False))
    sample = part.iloc[np.concatenate(picks)][columns].to_numpy(float)
    rho = spearmanr(sample, nan_policy="omit").statistic
    return pd.DataFrame(rho, index=columns, columns=columns)


def label_scores(table):
    rows = []
    for kind in ("thermal", "quenched"):
        part = table[table.kind == kind]
        for system, sys_table in part.groupby("system", sort=False):
            frames = sys_table.frame.to_numpy()
            seed = sys_table.seed.to_numpy()
            interior = sys_table.interior.to_numpy()
            stack = sys_table["stack"].to_numpy()
            far = (~seed & (sys_table.r_seed > sys_table.R + LIQUID_MARGIN)).to_numpy()
            for entry in LABELS:
                labels = sys_table[entry.name].to_numpy()
                ice = np.isin(labels, (1, 2))
                correct = ice & (labels == stack)

                def score(mask, hit):
                    per_frame = [hit[mask & (frames == f)].mean() for f in np.unique(frames)
                                 if (mask & (frames == f)).any()]
                    if not per_frame:
                        return np.nan, np.nan, np.nan
                    boot = [np.mean(rng.choice(per_frame, len(per_frame))) for _ in range(BOOTSTRAP)]
                    return np.mean(per_frame), *np.percentile(boot, [2.5, 97.5])

                for metric, mask, hit in (("seed_correct", seed, correct),
                                          ("seed_ice_any", seed, ice),
                                          ("seed_wrong_polytype", seed, ice & ~correct),
                                          ("interior_correct", seed & interior, correct),
                                          ("liquid_false_positive", far, ice)):
                    mean, lo, hi = score(mask, hit)
                    rows.append(dict(kind=kind, system=system, method=entry.name, metric=metric,
                                     value=mean, lo=lo, hi=hi))
    return pd.DataFrame(rows)


STEPS = ("ideal", "thermal", "heldout", "grouped", "model", "pressure", "labels")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--steps", default=",".join(STEPS))
    steps = parser.parse_args().steps.split(",")
    RESULTS.mkdir(exist_ok=True)
    seedtest = pd.read_pickle(CACHE/"seedtest.pkl")
    columns = [c for c in seedtest.columns if c in COLUMN_ENTRY or c in {e.name for e in SCALARS}]
    pd.DataFrame([column_info(c) for c in columns]).to_csv(RESULTS/"columns.csv", index=False)
    if "ideal" in steps:
        genice = pd.read_pickle(CACHE/"genice.pkl")
        ideal_separability(genice, columns).to_csv(RESULTS/"separability_ideal.csv", index=False)
        print("ideal separability done")
    if "thermal" in steps:
        thermal_separability(seedtest, columns).to_csv(RESULTS/"separability_thermal.csv", index=False)
        print("thermal separability done")
    if "heldout" in steps:
        replica = pd.read_pickle(CACHE/"replica.pkl")
        replica_autocorrelation(replica, columns).to_csv(RESULTS/"replica_autocorrelation.csv", index=False)
        heldout_separability(seedtest, replica, columns).to_csv(RESULTS/"separability_heldout.csv", index=False)
        print("held-out separability done")
    if "grouped" in steps:
        prep = pd.read_pickle(CACHE/"prep.pkl")
        grouped_separability(seedtest, prep, columns).to_csv(RESULTS/"separability_grouped.csv", index=False)
        print("grouped separability done")
    if "model" in steps:
        prep, model = pd.read_pickle(CACHE/"prep.pkl"), pd.read_pickle(CACHE/"model.pkl")
        model_separability(prep, model, columns).to_csv(RESULTS/"separability_model.csv", index=False)
        print("model separability done")
    if "pressure" in steps:
        prep, pressure = pd.read_pickle(CACHE/"prep.pkl"), pd.read_pickle(CACHE/"pressure.pkl")
        ideal = pd.read_csv(RESULTS/"separability_ideal.csv")
        pressure_separability(prep, pressure, ideal, columns).to_csv(RESULTS/"separability_pressure.csv",
                                                                    index=False)
        print("pressure separability done")
    if "labels" in steps:
        redundancy(seedtest, columns).to_csv(RESULTS/"redundancy_spearman.csv")
        label_scores(seedtest).to_csv(RESULTS/"label_scores.csv", index=False)
        print("redundancy and label scores done")


if __name__ == "__main__":
    main()
