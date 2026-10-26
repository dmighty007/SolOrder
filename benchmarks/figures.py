"""Pilot figures F1 (separability), F2 (redundancy), F7 (seed-test labels) and
S2 (oracle versus held-out thresholds).

Reads only benchmarks/results/*.csv, so figures can be restyled without
recomputing anything. Writes PNG and PDF to benchmarks/results/figures/.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

RESULTS = Path(__file__).resolve().parent/"results"
FIGURES = RESULTS/"figures"

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SERIES = {"chill": "#2a78d6", "chill_plus": "#eb6834", "ptm": "#1baf7a"}
METHOD_NAMES = {"chill": "CHILL", "chill_plus": "CHILL+", "ptm": "PTM"}
# Single-hue sequential ramp: score 0.5 (chance) light, 1.0 dark.
ACCURACY = LinearSegmentedColormap.from_list(
    "accuracy", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
GROUP_ORDER = ["radial", "angular", "hbond", "bond_orientational", "template"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": GRID, "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
})


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES/f"{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def ordered_columns():
    info = pd.read_csv(RESULTS/"columns.csv")
    info["rank"] = info.group.map(GROUP_ORDER.index)
    return info.sort_values("rank", kind="stable")


def group_separators(ax, info, horizontal=True):
    edges = np.flatnonzero(info.group.to_numpy()[1:] != info.group.to_numpy()[:-1]) + .5
    for edge in edges:
        (ax.axhline if horizontal else ax.axvline)(edge, color=SURFACE, lw=2)


def heatmap(ax, matrix, info, xlabels, annotate):
    image = ax.imshow(matrix, cmap=ACCURACY, vmin=.5, vmax=1., aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(info)), info.column)
    ax.set_xticks(range(len(xlabels)), xlabels, rotation=90)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    group_separators(ax, info)
    if annotate:
        for (i, j), value in np.ndenumerate(matrix):
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=6,
                        color=SURFACE if value > .8 else INK)
    return image


def f1_ideal():
    info = ordered_columns()
    table = pd.read_csv(RESULTS/"separability_ideal.csv")
    table["pair"] = table.a + "–" + table.b
    matrix = table.pivot(index="column", columns="pair", values="oracle")
    pairs = list(dict.fromkeys(table.pair))
    matrix = matrix.loc[info.column, pairs].to_numpy()
    fig, ax = plt.subplots(figsize=(max(11, .11*len(pairs)), 7.5))
    image = heatmap(ax, matrix, info, pairs, annotate=False)
    ax.set_title("Ideal GenIce crystals: oracle score (1 + KS D)/2 per phase pair")
    fig.colorbar(image, ax=ax, fraction=.02, pad=.01, label="oracle score, (1 + KS D)/2")
    ax.tick_params(axis="x", labelsize=5)
    save(fig, "S1_separability_ideal_all_pairs")


def f1_ideal_counts():
    """Phase x phase: how many parameters separate each pair (accuracy >= 0.95)."""
    table = pd.read_csv(RESULTS/"separability_ideal.csv")
    phases = list(dict.fromkeys([*table.a, *table.b]))
    counts = (table.oracle >= .95).groupby([table.a, table.b]).sum()
    matrix = np.full((len(phases), len(phases)), np.nan)
    for (a, b), value in counts.items():
        i, j = phases.index(a), phases.index(b)
        matrix[i, j] = matrix[j, i] = value
    n_columns = table.column.nunique()
    ramp = LinearSegmentedColormap.from_list("count", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    image = ax.imshow(matrix, cmap=ramp, vmin=0, vmax=n_columns, interpolation="nearest")
    ax.set_xticks(range(len(phases)), phases, rotation=90)
    ax.set_yticks(range(len(phases)), phases)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for (i, j), value in np.ndenumerate(matrix):
        if np.isfinite(value):
            ax.text(j, i, f"{value:.0f}", ha="center", va="center", fontsize=6,
                    color=SURFACE if value > .55*n_columns else INK)
    ax.set_title(f"Ideal GenIce crystals: parameters (of {n_columns}) separating each pair at ≥ 0.95")
    fig.colorbar(image, ax=ax, fraction=.04, pad=.02, label="number of parameters")
    save(fig, "F1a_separability_ideal")


def f1_thermal():
    info = ordered_columns()
    table = pd.read_csv(RESULTS/"separability_thermal.csv")
    names = {"liquid": "liquid", "Ih": "Ih", "Ic": "Ic", "Isd_hex": "Isd hex", "Isd_cubic": "Isd cubic"}
    table["pair"] = table.a.map(names) + " vs " + table.b.map(names)
    pairs = list(dict.fromkeys(table.pair))
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 7.5), sharey=True)
    for ax, kind in zip(axes, ("thermal", "quenched")):
        part = table[table.kind == kind].pivot(index="column", columns="pair", values="oracle")
        image = heatmap(ax, part.loc[info.column, pairs].to_numpy(), info, pairs, annotate=True)
        ax.set_title(f"TIP4P/Ice 220 K, {kind}")
    fig.colorbar(image, ax=axes, fraction=.03, pad=.02, label="oracle score, (1 + KS D)/2")
    save(fig, "F1b_separability_thermal")


def f2_redundancy():
    rho = pd.read_csv(RESULTS/"redundancy_spearman.csv", index_col=0)
    keep = rho.notna().all().to_numpy()
    rho = rho.loc[keep, keep]
    distance = 1. - np.abs(rho.to_numpy())
    np.fill_diagonal(distance, 0.)
    tree = linkage(squareform(np.clip((distance+distance.T)/2, 0, None)), method="average")
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    dendrogram(tree, labels=rho.index.tolist(), ax=ax, color_threshold=0,
               above_threshold_color=INK2, leaf_font_size=7, leaf_rotation=90)
    ax.axhline(.2, color=GRID, lw=1, ls="--")
    ax.text(ax.get_xlim()[1], .21, "|ρ| = 0.8", ha="right", va="bottom", color=INK2, fontsize=7)
    ax.set_ylabel("1 − |Spearman ρ|")
    ax.set_title("Redundancy among parameters (TIP4P/Ice 220 K: Ih, Ic and liquid pooled)")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    save(fig, "F2_redundancy")


def f7_labels():
    table = pd.read_csv(RESULTS/"label_scores.csv")
    systems = [s for s in dict.fromkeys(table.system) if s.startswith("seed")]
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 5), sharex=True, sharey=True)
    width = .26
    for ax, kind in zip(axes, ("thermal", "quenched")):
        part = table[(table.kind == kind) & (table.metric == "interior_correct")]
        for k, method in enumerate(SERIES):
            rows = part[part.method == method].set_index("system").loc[systems]
            x = np.arange(len(systems)) + (k-1)*width
            ax.bar(x, rows.value, width-.03, color=SERIES[method], label=METHOD_NAMES[method])
            ax.errorbar(x, rows.value, yerr=[rows.value-rows.lo, rows.hi-rows.value],
                        fmt="none", ecolor=INK2, elinewidth=.8, capsize=1.5)
        ax.set_ylim(0, 1)
        ax.set_ylabel("interior seed molecules\nwith correct polytype")
        ax.set_title(f"{kind.capitalize()} frames", loc="left")
        ax.yaxis.grid(True, color=GRID, lw=.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].legend(frameon=False, ncol=3, loc="upper left")
    axes[1].set_xticks(range(len(systems)), [s.replace("seed_", "").replace("_", "\nN=") for s in systems])
    save(fig, "F7_seed_labels")


def s2_heldout():
    """Each point is one parameter and one population pair; on the diagonal, nothing is lost."""
    table = pd.read_csv(RESULTS/"separability_heldout.csv")
    tasks = {"ice vs liquid": table.b == "liquid", "Ih vs Ic stacking": table.b != "liquid"}
    colors = dict(zip(tasks, ("#2a78d6", "#eb6834")))
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6), sharex=True, sharey=True)
    panels = (("A_to_B", "(a) Fitted on run A, tested on run B", "oracle score, run B"),
              ("T220_to_T250", "(b) Ih vs Ic: fitted at 220 K, tested at 250 K", "oracle score at 250 K"))
    for ax, (column, title, xlabel) in zip(axes, panels):
        ax.plot([.5, 1.], [.5, 1.], color=GRID, lw=1, zorder=0)
        for task, mask in tasks.items():
            part = table[mask & table[column].notna()]
            if len(part):
                ax.scatter(part.oracle_B, part[column], s=14, color=colors[task], edgecolor=SURFACE,
                           linewidth=.6, label=task)
        ax.set_title(title, loc="left")
        ax.set_xlabel(xlabel)
        ax.set_xlim(.48, 1.01)
        ax.set_ylim(.48, 1.01)
        ax.set_aspect("equal")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("held-out balanced accuracy")
    axes[0].legend(frameon=False, loc="upper left")
    save(fig, "S2_heldout")


def s3_grouped():
    """Independent preparations C1-C5: worst held-out score against the mean oracle score."""
    table = pd.read_csv(RESULTS/"separability_grouped.csv")
    tasks = {"ice vs liquid": table.b == "liquid", "Ih vs Ic stacking": table.b != "liquid"}
    colors = dict(zip(tasks, ("#2a78d6", "#eb6834")))
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6), sharex=True, sharey=True)
    panels = (("loo_min", "(a) Leave one preparation out"),
              ("transfer_min", "(b) Fitted on run A (250 K: on C at 220 K)"))
    for ax, (column, title) in zip(axes, panels):
        ax.plot([.5, 1.], [.5, 1.], color=GRID, lw=1, zorder=0)
        for task, mask in tasks.items():
            part = table[mask]
            ax.scatter(part.oracle_mean, part[column], s=14, color=colors[task], edgecolor=SURFACE,
                       linewidth=.6, label=task)
        ax.set_title(title, loc="left")
        ax.set_xlabel("mean oracle score, C1–C5")
        ax.set_xlim(.43, 1.01)
        ax.set_ylim(.43, 1.01)
        ax.set_aspect("equal")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("lowest held-out balanced accuracy\nover the five preparations")
    axes[0].legend(frameon=False, loc="upper left")
    save(fig, "S3_grouped")


def f7b_depth():
    """Correct polytype by depth below the constructed seed surface, all seed sizes pooled."""
    table = pd.read_csv(RESULTS/"seed_audit_depth.csv")
    table = table[(table.kind == "thermal") & (table.outcome == "correct")]
    phases = {"Ih": "#2a78d6", "Ic": "#eb6834", "Isd": "#1baf7a"}
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.8), sharex=True, sharey=True)
    for ax, method in zip(axes, SERIES):
        part = table[table.method == method]
        intact = part.groupby("depth").apply(lambda g: np.average(g.intact, weights=g.n), include_groups=False)
        counts = part.groupby("depth").n.sum()
        keep = counts.index[counts >= 30]
        ax.plot(keep, intact.loc[keep], color=INK2, lw=1.2, ls="--", label="lattice intact")
        for phase, color in phases.items():
            rows = part[part.phase == phase]
            pooled = rows.groupby("depth").apply(lambda g: np.average(g.fraction, weights=g.n), include_groups=False)
            n = rows.groupby("depth").n.sum()
            pooled = pooled[n >= 30]
            ax.plot(pooled.index, pooled, color=color, lw=2, marker="o", ms=3.5, label=phase)
        ax.set_title(METHOD_NAMES[method], loc="left")
        ax.set_xlabel("depth below seed surface (Å)")
        ax.set_ylim(0, 1.02)
        ax.yaxis.grid(True, color=GRID, lw=.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("seed molecules with\ncorrect polytype")
    axes[0].legend(frameon=False, loc="lower right", fontsize=7)
    save(fig, "F7b_seed_depth")


def s4_model():
    """Second water model: thresholds moved between protocols and models (run D)."""
    table = pd.read_csv(RESULTS/"separability_model.csv")
    tasks = {"ice vs liquid": table.b == "liquid", "Ih vs Ic stacking": table.b != "liquid"}
    colors = dict(zip(tasks, ("#2a78d6", "#eb6834")))
    panels = (("C_to_ice_npt", "oracle_ice_npt", "(a) Protocol: NVT box → NPT", "TIP4P/Ice, NPT protocol"),
              ("ice_npt_to_t2005_220", "oracle_t2005_220", "(b) Model, same T (220 K)", "TIP4P/2005, 220 K"),
              ("ice_npt_to_t2005_200", "oracle_t2005_200", "(c) Model, same ΔT (200 K)", "TIP4P/2005, 200 K"),
              ("ice_npt_to_mW", "oracle_mW", "(d) Oxygen-only mW (220 K)", "mW, 220 K"))
    panels = [p for p in panels if p[0] in table]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.4), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, (column, oracle, title, xlabel) in zip(axes, panels):
        ax.plot([.5, 1.], [.5, 1.], color=GRID, lw=1, zorder=0)
        for task, mask in tasks.items():
            part = table[mask & table[column].notna()]
            ax.scatter(part[oracle], part[column], s=14, color=colors[task], edgecolor=SURFACE,
                       linewidth=.6, label=task)
        ax.set_title(title, loc="left")
        ax.set_xlabel(f"oracle score, {xlabel}")
        ax.set_xlim(.43, 1.01)
        ax.set_ylim(.43, 1.01)
        ax.set_aspect("equal")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in axes[::2]:
        ax.set_ylabel("held-out balanced accuracy\n(threshold from the source condition)")
    for ax in axes[len(panels):]:
        ax.set_visible(False)
    axes[0].legend(frameon=False, loc="upper left")
    save(fig, "S4_model")


PRESSURE_ORDER = ["Ih", "Ic", "II", "III", "V", "VI", "VII", "VIII", "liquid_0", "liquid_0.3", "liquid_0.5", "liquid_1.0",
                  "liquid_3.0"]
PRESSURE_NAMES = {"liquid_0": "liq 0 GPa", "liquid_0.3": "liq 0.3", "liquid_0.5": "liq 0.5", "liquid_1.0": "liq 1.0",
                  "liquid_3.0": "liq 3.0"}


def f8_pressure():
    """Thermal high-pressure phases (run E): best held-out score per pair, and ideal versus thermal."""
    table = pd.read_csv(RESULTS/"separability_pressure.csv")
    thermal = table[table.kind == "thermal"]
    best = thermal.groupby(["a", "b"]).heldout_min.max()
    n = len(PRESSURE_ORDER)
    matrix = np.full((n, n), np.nan)
    for (a, b), v in best.items():
        i, j = PRESSURE_ORDER.index(a), PRESSURE_ORDER.index(b)
        matrix[max(i, j), min(i, j)] = v
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), gridspec_kw=dict(width_ratios=(1.15, 1)))
    ax = axes[0]
    image = ax.imshow(matrix, cmap=ACCURACY, vmin=.5, vmax=1., interpolation="nearest")
    labels = [PRESSURE_NAMES.get(x, x) for x in PRESSURE_ORDER]
    ax.set_xticks(range(n), labels, rotation=90)
    ax.set_yticks(range(n), labels)
    ax.tick_params(length=0)
    for i in range(n):
        for j in range(i):
            if np.isfinite(matrix[i, j]):
                ax.text(j, i, f"{matrix[i, j]:.2f}".removeprefix("0"), ha="center", va="center", fontsize=6,
                        color=SURFACE if matrix[i, j] > .8 else INK)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title("(a) Best single parameter, held-out (thermal)", loc="left")
    fig.colorbar(image, ax=ax, shrink=.75, label="worse held-out balanced accuracy\nof the two preparations")
    ax = axes[1]
    crystals = thermal[thermal.ideal.notna()]
    ax.plot([.5, 1.], [.5, 1.], color=GRID, lw=1, zorder=0)
    ax.scatter(crystals.ideal, crystals.oracle, s=8, color="#2a78d6", edgecolor=SURFACE, linewidth=.4,
               alpha=.8)
    ax.set_xlim(.48, 1.01)
    ax.set_ylim(.48, 1.01)
    ax.set_aspect("equal")
    ax.set_xlabel("oracle score, ideal GenIce (jittered)")
    ax.set_ylabel("oracle score, thermal 220 K")
    ax.set_title("(b) Crystal pairs, every parameter", loc="left")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    save(fig, "F8_pressure")


# Fixed slot per parameter across all panels of F5 (validated categorical order).
SENSITIVITY_SERIES = {"cn": ("CN", "#1baf7a"), "cn_smooth": ("smooth CN", "#eda100"),
                      "q6_avg": ("q̄6", "#e87ba4"), "ring_f6": ("6-ring fraction", "#008300"),
                      "q3_avg": ("q̄3", "#6250d6"), "q4_avg": ("q̄4", "#e34948")}


def f5_sensitivity():
    """Neighbor-convention scan on run C1-C3 (sensitivity.py)."""
    table = pd.read_csv(RESULTS/"sensitivity.csv")
    thermal = table[table.kind == "thermal"]

    def series(factor, column, a, b, kind="thermal", score="oracle"):
        part = table[(table.kind == kind) & (table.factor == factor) & (table.column == column)
                     & (table.a == a) & (table.b == b)].sort_values("value")
        return part.value.to_numpy(), part[score].to_numpy()

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6))
    axes = axes.ravel()
    panels = ((axes[0], "(a) First-shell cutoff: Ih vs liquid", ("cn", "cn_smooth", "q6_avg", "ring_f6"),
               "Ih", "liquid", "oracle"),
              (axes[1], "(b) First-shell cutoff: Ih vs Ic", ("q3_avg", "q4_avg", "q6_avg"), "Ih", "Ic", "oracle"),
              (axes[2], "(c) Threshold fitted at 3.5 Å, applied\n     at another cutoff: Ih vs liquid",
               ("cn", "cn_smooth", "q6_avg", "ring_f6"), "Ih", "liquid", "cross"))
    for ax, title, columns, a, b, score in panels:
        for column in columns:
            name, color = SENSITIVITY_SERIES[column]
            x, y = series("shell", column, a, b, score=score)
            ax.plot(x, y, color=color, lw=2, marker="o", ms=4, label=name)
        ax.axvline(3.5, color=GRID, lw=1, zorder=0)
        ax.set_xlabel("first-shell cutoff (Å)")
        ax.set_title(title, loc="left")
        ax.legend(frameon=False, fontsize=7, **(dict(loc="center left", bbox_to_anchor=(0, .66)) if score == "cross"
                  else dict(loc="lower right") if b == "liquid" else dict(loc="center left", bbox_to_anchor=(0, .3))))
    ax = axes[3]
    for column in ("q3_avg", "q4_avg", "q6_avg"):
        name, color = SENSITIVITY_SERIES[column]
        for (a, b), style in ((("Ih", "liquid"), "-"), (("Ih", "Ic"), "--")):
            x, y = series("average", column, a, b)
            ax.plot(x, y, color=color, lw=2, ls=style, marker="o", ms=4,
                    label=f"{name}, {a} vs {b}")
    ax.axvline(3.5, color=GRID, lw=1, zorder=0)
    ax.set_xlabel("averaging radius (Å); q_lm at 3.5 Å")
    ax.set_title("(d) Averaging radius", loc="left")
    ax.legend(frameon=False, fontsize=6.5, ncol=1)
    ax = axes[4]
    for (kind, a, b), style in ((("thermal", "Ih", "liquid"), "-"), (("thermal", "Ih", "Ic"), "--"),
                                (("quenched", "Ih", "Ic"), ":")):
        x, y = series("entropy", "entropy_avg", a, b, kind=kind)
        ax.plot(x, y, color=INK2, lw=2, ls=style, marker="o", ms=4, label=f"{a} vs {b}, {kind}")
    ax.axvline(5., color=GRID, lw=1, zorder=0)
    ax.set_xlabel("entropy integration and averaging radius r_m (Å)")
    ax.set_title("(e) Averaged pair entropy s̄", loc="left")
    ax.legend(frameon=False, fontsize=7)
    ax = axes[5]
    labels = table[(table.kind == "thermal") & table.b.isna()]
    for method, name in METHOD_NAMES.items():
        if method == "ptm":
            continue
        for phase, style in (("Ih", "-"), ("Ic", "--")):
            part = labels[(labels.column == method) & (labels.a == phase)].sort_values("value")
            ax.plot(part.value, part.label_score, color=SERIES[method], lw=2, ls=style, marker="o", ms=4,
                    label=f"{name}, bulk {phase}")
    ax.axvline(3.5, color=GRID, lw=1, zorder=0)
    ax.set_xlabel("neighbor cutoff (Å)")
    ax.set_ylabel("fraction labeled with the correct polytype")
    ax.set_title("(f) Labelers, thermal bulk ice", loc="left")
    ax.legend(frameon=False, fontsize=7)
    for ax in axes:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in axes[:5]:
        ax.set_ylabel("held-out balanced accuracy" if ax is axes[2] else "oracle score")
    fig.tight_layout()
    save(fig, "F5_sensitivity")


def f6_cost():
    """Wall time per molecule against system size (cost.py)."""
    tables = [pd.read_csv(path) for path in sorted(RESULTS.glob("cost_jobs*.csv"))]
    table = pd.concat(tables, ignore_index=True)
    liquid = table[table.phase == "liquid"]
    largest = liquid.molecules.max()
    top = liquid[(liquid.molecules == largest) & (liquid.jobs == 1) & ~liquid.entry.isin(("snapshot", "all"))]
    order = top.sort_values("us_per_molecule").entry.tolist()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw=dict(width_ratios=(1.15, 1)))
    ax = axes[0]
    jobs = sorted(liquid.jobs.unique())
    colors = dict(zip(jobs, ("#2a78d6", "#eb6834")))
    for k, n in enumerate(jobs):
        part = liquid[(liquid.molecules == largest) & (liquid.jobs == n)].set_index("entry").loc[order]
        ax.barh(np.arange(len(order)) + (k - .5)*.4, part.us_per_molecule, height=.38, color=colors[n],
                label=f"{n} thread{'s' if n > 1 else ''}")
    ax.set_yticks(range(len(order)), order)
    ax.set_xscale("log")
    ax.set_xticks([2, 5, 10, 20, 50], ["2", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_xlabel(f"µs per molecule, liquid, {largest:,} molecules")
    ax.set_title("(a) Cost per registry entry", loc="left")
    ax.legend(frameon=False, loc="lower right")
    ax = axes[1]
    for (phase, n), part in table[table.entry == "all"].groupby(["phase", "jobs"]):
        part = part.sort_values("molecules")
        ax.plot(part.molecules, part.seconds, color=colors[n], ls="-" if phase == "liquid" else "--",
                lw=2, marker="o", ms=4, label=f"{phase}, {n} thread{'s' if n > 1 else ''}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([2000, 4000, 8000, 16000, 32000, 64000], ["2k", "4k", "8k", "16k", "32k", "64k"])
    ax.set_yticks([.3, 1, 3, 10], ["0.3", "1", "3", "10"])
    ax.minorticks_off()
    ax.set_xlabel("molecules")
    ax.set_ylabel("seconds per frame, all entries")
    ax.set_title("(b) Full registry against system size", loc="left")
    ax.legend(frameon=False)
    for ax in axes:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    save(fig, "F6_cost")


if __name__ == "__main__":
    f1_ideal()
    f1_ideal_counts()
    f1_thermal()
    f2_redundancy()
    f7_labels()
    s2_heldout()
    f7b_depth()
    s3_grouped()
    s4_model()
    f8_pressure()
    f5_sensitivity()
    f6_cost()
