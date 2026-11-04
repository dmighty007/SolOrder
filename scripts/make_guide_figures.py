"""Draw the figures of the beginner's guide from docs/guide/data into docs/guide/figures.

Every number in a figure comes from the same SolOrder calls shown in the guide.

  python scripts/make_guide_figures.py
"""

import warnings
from pathlib import Path

import matplotlib
import MDAnalysis as mda
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

from SolOrder import OrderParameter, eFingerprint  # noqa: E402

warnings.filterwarnings("ignore", module="MDAnalysis")

ROOT = Path(__file__).resolve().parents[1]
DATA, FIGURES = ROOT/"docs/guide/data", ROOT/"docs/guide/figures"
PHASES = {"ice_Ih": "ice Ih", "ice_Ic": "ice Ic", "liquid": "liquid"}
COLORS = {"ice_Ih": "#2a78d6", "ice_Ic": "#eb6834", "liquid": "#1baf7a"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight", "font.size": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": .6, "xtick.color": MUTED, "ytick.color": MUTED, "legend.frameon": False,
    "lines.linewidth": 2,
})


def load(name):
    return mda.Universe(str(DATA/f"{name}.gro"))


def distributions():
    """Four order parameters, one histogram panel each, ice Ih / Ic / liquid overlaid."""
    values = {}
    for name in PHASES:
        u = load(name)
        op = OrderParameter(u, cutoff=5.)
        values[name] = {
            "Tetrahedral order $q$": op.tetrahedral_order(),
            "Averaged Steinhardt $\\bar{q}_6$": op.steinhardt(6, averaged=True, cutoff=3.5),
            "Local structure index (Å$^2$)": op.local_structure_index(),
            "Local entropy $\\bar{s}$ ($k_B$)": eFingerprint(u, cutoff=5.).localEntropy(),
        }
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
    for ax, label in zip(axes.flat, values["liquid"]):
        data = np.concatenate([v[label] for v in values.values()])
        bins = np.linspace(*np.nanpercentile(data, (.5, 99.5)), 60)
        for name, v in values.items():
            ax.hist(v[label], bins=bins, density=True, histtype="step", color=COLORS[name], label=PHASES[name])
        ax.set_xlabel(label)
        ax.set_ylabel("probability density")
    axes[0, 0].legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(FIGURES/"distributions.png")
    plt.close(fig)
    return values


def scatter(values):
    """Two parameters together separate all three structures."""
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    x, y = "Tetrahedral order $q$", "Averaged Steinhardt $\\bar{q}_6$"
    for name, v in values.items():
        ax.scatter(v[x], v[y], s=6, alpha=.35, color=COLORS[name], edgecolors="none", label=PHASES[name])
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_xlim(.5, 1.)
    ax.legend(markerscale=3, loc="upper left")
    fig.savefig(FIGURES/"scatter.png")
    plt.close(fig)


def chill_fractions():
    """CHILL+ label fractions per structure as stacked horizontal bars."""
    names = ["hexagonal_ice", "cubic_ice", "interfacial_ice", "hydrate", "interfacial_hydrate", "other"]
    shown = ["hexagonal ice", "cubic ice", "interfacial ice", "hydrate", "interfacial hydrate", "other"]
    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#b5b4ae"]
    fractions = {}
    for name in PHASES:
        labels = OrderParameter(load(name), cutoff=3.5).chill_plus().names
        fractions[name] = [np.mean(labels == n) for n in names]
    fig, ax = plt.subplots(figsize=(7.2, 2.2))
    for row, name in enumerate(PHASES):
        left = 0.
        for value, color in zip(fractions[name], colors):
            ax.barh(row, value, left=left, color=color, edgecolor="white", linewidth=1.5, height=.6)
            if value > .06:
                ax.text(left + value/2, row, f"{100*value:.0f}%", ha="center", va="center", color="white",
                        fontsize=8, fontweight="bold")
            left += value
    ax.set_yticks(range(len(PHASES)), PHASES.values())
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of molecules")
    ax.grid(axis="y", visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors]
    ax.legend(handles, shown, ncol=6, loc="lower center", bbox_to_anchor=(.5, 1.02), fontsize=7.5,
              handlelength=1, columnspacing=1)
    fig.savefig(FIGURES/"chill_plus_fractions.png")
    plt.close(fig)


def largest_cluster(op, ice):
    """Mask of the largest connected group of ``ice`` molecules (neighbors within the op cutoff)."""
    graph = op.neighbor_list()
    source = np.repeat(np.arange(op.n), graph.counts)
    keep = ice[source] & ice[graph.indices]
    adjacency = coo_matrix((np.ones(keep.sum()), (source[keep], graph.indices[keep])), shape=(op.n, op.n))
    _, component = connected_components(adjacency, directed=False)
    sizes = np.bincount(component[ice])
    return ice & (component == np.argmax(sizes))


def seed_slice():
    """A 6 Å slab through the centre of the ice seed, coloured by CHILL+ label."""
    u = mda.Universe(str(DATA/"seed.gro.gz"))
    op = OrderParameter(u, cutoff=3.5)
    labels = op.chill_plus().labels
    nucleus = largest_cluster(op, np.isin(labels, (1, 2)))
    truth = np.load(DATA/"seed_truth.npz")
    centre, radius = truth["centre"], float(truth["R"])
    d = op.pos - centre
    d -= op.box[:3]*np.round(d/op.box[:3])
    slab = (np.abs(d[:, 2]) < 3.) & (np.abs(d[:, 0]) < 40) & (np.abs(d[:, 1]) < 40)
    groups = [("other", ~np.isin(labels, (1, 2, 3)), "#c9c8c2"),
              ("interfacial ice", labels == 3, "#eb6834"),
              (f"largest ice cluster ({nucleus.sum()} molecules)", nucleus, "#2a78d6")]
    fig, ax = plt.subplots(figsize=(4.8, 4.8))
    for label, mask, color in groups:
        m = slab & mask
        ax.scatter(d[m, 0], d[m, 1], s=14, color=color, edgecolors="white", linewidths=.4, label=label)
    ax.add_patch(plt.Circle((0, 0), radius, fill=False, ls="--", lw=1.2, color=INK))
    edge = radius/np.sqrt(2)
    ax.annotate(f"constructed seed\n{truth['seed'].sum()} molecules", (edge, edge), (24, 33), fontsize=7.5,
                color=INK, ha="left", va="bottom", arrowprops=dict(arrowstyle="-", color=INK, lw=.8),
                bbox=dict(boxstyle="round,pad=.25", fc="white", ec="none"))
    ax.set_aspect("equal")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.grid(False)
    ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.0), ncol=1, fontsize=8, markerscale=1.4)
    fig.savefig(FIGURES/"seed_slice.png")
    plt.close(fig)


def coordination_cutoff():
    """Mean coordination number against cutoff: choose the cutoff on the plateau."""
    cutoffs = np.arange(2.6, 4.81, .05)
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for name in PHASES:
        op = OrderParameter(load(name), cutoff=5.)
        ax.plot(cutoffs, [op.coordination_number(c).mean() for c in cutoffs], color=COLORS[name],
                label=PHASES[name])
    ax.axvline(3.5, color=MUTED, lw=1, ls="--")
    ax.text(3.53, .3, "3.5 Å", color=MUTED, fontsize=8)
    ax.set_xlabel("cutoff (Å)")
    ax.set_ylabel("mean coordination number")
    ax.set_ylim(0, 12)
    ax.legend(loc="upper left")
    fig.savefig(FIGURES/"coordination_cutoff.png")
    plt.close(fig)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    scatter(distributions())
    chill_fractions()
    seed_slice()
    coordination_cutoff()


if __name__ == "__main__":
    main()
