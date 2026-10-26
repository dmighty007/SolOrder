"""Create CSV summaries and standalone figures from sampled validation results."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


GROUPS = ("PathGennie_raw", "metadynamics", "MD_shooting", "mW_bruteforce", "PathGennie_quenched")
TOP_AUC = 15  # rows in the README association table


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=root/"docs/validation/nucleation")
    args = parser.parse_args()
    directory = args.directory
    manifest = json.loads((directory/"manifest.json").read_text())
    if "analyzed_frames" not in manifest:
        raise RuntimeError("Validation batch is incomplete; wait for the analysis process to finish")
    frames = pd.read_csv(directory/"frames.csv")
    if manifest["analyzed_frames"] != len(frames):
        raise RuntimeError("Manifest and frame table do not agree")
    checks = pd.read_csv(directory/"checks.csv")
    groups = pd.read_csv(directory/"phase_conditioned.csv")
    source_map = {s["name"]: s for s in manifest["sources"]}
    audits_path = directory/"independent_audits.csv"
    audits = pd.read_csv(audits_path) if audits_path.exists() else None
    summary = []
    for name, data in frames.groupby("source", sort=False):
        data = data.sort_values("frame")
        source = source_map[name]
        first, last = data.iloc[0], data.iloc[-1]
        summary.append(dict(source=name, family=source["family"], samples=len(data),
                            stored_frames=source["frames"], centers=int(first.waters), full_atom=source["full_atom"],
                            approximate_box=source["box_approximate"], recovered=source["kind"] == "zlib_recovery",
                            core_first=int(first.largest_Ih_Ic_core), core_last=int(last.largest_Ih_Ic_core),
                            core_max=int(data.largest_Ih_Ic_core.max()),
                            Q_first=first.Q_mean, Q_last=last.Q_mean,
                            entropy_first=first.entropy_local_mean, entropy_last=last.entropy_local_mean,
                            HB_first=first.get("HB_mean", np.nan), HB_last=last.get("HB_mean", np.nan)))
    pd.DataFrame(summary).to_csv(directory/"source_summary.csv", index=False)
    # Restrict the primary contrast ranking to readable, fixed-box raw PG paths.
    reliable = {s["name"] for s in manifest["sources"] if s["family"] == "PathGennie_raw" and
                s["kind"] == "hdf5" and not s["box_approximate"]}
    eligible = frames[frames.source.isin(reliable) & (frames.largest_Ih_Ic_core >= 20)].copy()
    ranking = []
    for column in frames.columns:
        if column.endswith("_core_auc"):
            values = eligible[column].dropna()
            if len(values):
                ranking.append(dict(parameter=column[:-9], frames=len(values), median_auc=float(values.median()),
                                    q25_auc=float(values.quantile(.25)), q75_auc=float(values.quantile(.75)),
                                    median_separation=float(np.median(np.abs(values-.5)))))
    ranking = sorted(ranking, key=lambda r: r["median_separation"], reverse=True)
    pd.DataFrame(ranking).to_csv(directory/"core_contrast.csv", index=False)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    for axis, family, title in zip(axes, ("PathGennie_raw", "metadynamics", "mW_bruteforce"),
                                   ("Raw PathGennie paths", "TIP4P/Ice metadynamics", "mW brute-force paths")):
        colors = {"PG_E": "tab:green", "PG_F": "tab:purple", "PG_G_Ic": "gray",
                  "PG_I_Ic": "tab:orange", "PG_J_Ih": "tab:red", "PG_K_Ih": "tab:brown",
                  "PG_M_Ih": "tab:blue", "PG_O_Ic": "tab:cyan"}
        for name, data in frames[frames.family == family].groupby("source", sort=False):
            source = source_map[name]
            style = "--" if source["box_approximate"] else ":" if source["kind"] == "zlib_recovery" else "-"
            axis.plot(data.progress, data.largest_Ih_Ic_core, style, marker=".", lw=1.4, color=colors.get(name), label=name)
        axis.set(title=title, xlabel="Stored-frame progress", ylabel="Largest Ih/Ic core (waters)", xlim=(0, 1))
        axis.legend(fontsize=6, loc="upper left", ncol=2 if family != "metadynamics" else 1)
        axis.grid(alpha=.2)
    fig.suptitle("Sampled nucleation configurations: strict CHILL+ Ih/Ic connected core\nDotted: recovered storage order; dashed O_Ic: reconstructed NPT boxes", fontsize=10)
    fig.savefig(directory/"core_evolution.png", dpi=180)
    fig.savefig(directory/"core_evolution.pdf")
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for axis, parameter, label in zip(axes.ravel(), ("Q", "entropy_local", "HB", "F4"),
                                     ("Orientational tetrahedral Q", "Local entropy (kB units)", "H-bonds per water", "F4")):
        for name, color in (("PG_I_Ic", "tab:orange"), ("PG_M_Ih", "tab:blue")):
            for group, style in (("largest_core", "-"), ("outside_core", "--")):
                data = groups[(groups.source == name) & (groups.parameter == parameter) & (groups.group == group)].sort_values("frame")
                # Tiny random cores are omitted from core curves.
                if group == "largest_core":
                    data = data[data.group_size >= 20]
                axis.plot(data.progress, data["mean"], linestyle=style, marker=".", color=color, label=f"{name}: {group}")
        axis.set(xlabel="Stored-frame progress", ylabel=label, xlim=(0, 1))
        axis.grid(alpha=.2)
        if parameter == "HB":
            axis.set_ylim(3.5, 4.1)
    axes[0, 0].legend(fontsize=7, handlelength=4)
    fig.suptitle("Raw TIP4P/Ice: nucleus-local OP contrast versus the rest of the box", fontsize=12)
    fig.savefig(directory/"core_op_contrast.png", dpi=180)
    fig.savefig(directory/"core_op_contrast.pdf")
    plt.close(fig)
    # The strongest associations, without gaps; core_contrast.csv has the rest.
    selected = ranking[:TOP_AUC]
    fig, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    y = np.arange(len(selected))
    median = np.array([r["median_auc"] for r in selected])
    lower = median-np.array([r["q25_auc"] for r in selected])
    upper = np.array([r["q75_auc"] for r in selected])-median
    axis.errorbar(median, y, xerr=[lower, upper], fmt="o", color="tab:blue", capsize=3)
    axis.axvline(.5, color="gray", linestyle="--", lw=1)
    axis.set(yticks=y, yticklabels=[r["parameter"] for r in selected], xlim=(0, 1),
             xlabel="Per-frame ROC AUC against largest-core membership (median and IQR)")
    axis.invert_yaxis()
    axis.set_title("Descriptive core association on raw E/F/I/M samples; core size >=20")
    fig.savefig(directory/"core_association.png", dpi=180)
    fig.savefig(directory/"core_association.pdf")
    plt.close(fig)
    failures = checks[(checks.label_mismatches > 0) | (checks.cn_mismatches > 0) |
                      (checks.hb_mismatches.fillna(0) > 0) | (checks.bound_violations > 0) |
                      (checks.inf_count > 0) | (checks.voronoi_relative_volume_error > 1e-5)]
    coverage = frames.groupby("family").agg(sources=("source", "nunique"), frames=("frame", "size"),
                                             centers=("waters", "sum"))
    checked_parallel = checks.parallel_max_abs_error.dropna()
    report = f'''# Validation on nucleation configurations

**{len(frames)} sampled frames from {frames.source.nunique()} trajectories / path files**
were analyzed from `{manifest['source_project']}`. This is a deterministic sampled
validation run, not an exhaustive calculation of every stored frame.
The analyzed configurations contain {int(frames.waters.sum()):,} oxygen-center instances.

| Dataset family | Sources | Sampled frames | Oxygen-center instances |
|---|---:|---:|---:|
'''
    for family in GROUPS:
        if family in coverage.index:
            r = coverage.loc[family]
            report += f"| {family} | {int(r.sources)} | {int(r.frames)} | {int(r.centers):,} |\n"
    report += f'''
## Numerical checks

- Frames with failed checks: **{len(failures)}**.
- Input coordinates and periodic boxes were verified after MDAnalysis precision
  conversion before calculation; coordinates were rechecked after OP calls.
- Independent project CHILL+ mismatches: **{int(checks.label_mismatches.sum())}** centers.
- Independent cKDTree coordination mismatches: **{int(checks.cn_mismatches.sum())}** centers.
- Independent donor-H-acceptor H-bond count mismatches: **{int(checks.hb_mismatches.fillna(0).sum())}** centers
  across {int(checks.hb_mismatches.notna().sum())} full-atom snapshots.
- Serial/parallel comparison: one complete snapshot per source ({len(checked_parallel)} snapshots);
  maximum absolute per-center difference **{checked_parallel.max():.3g}**.
- Maximum relative Voronoi volume-conservation error: **{checks.voronoi_relative_volume_error.max():.3g}**.
- Infinite OP outputs: **{int(checks.inf_count.sum())}**. NaN outputs (centers lacking the
  neighbors a parameter needs): **{int(checks.nan_count.sum())}**; per-parameter counts are in frames.csv.

The independent CHILL+ reference is the source project's `lib/order.py`, using
an explicit harmonic polynomial implementation and cKDTree neighbor search.
SolOrder uses SciPy spherical harmonics, MDAnalysis geometry, and Numba reductions.
The reference was evaluated on the same analysis coordinate precision and box.
This checks implementation agreement; it does not supply independent experimental
phase labels. Selected coordinate/box hashes, source file sizes and timestamps,
package versions and source-reference hash are in [manifest.json](manifest.json).

## Raw PathGennie paths

| Path | Samples / stored frames | Largest core: first -> last | Maximum sampled core | Mean Q: first -> last |
|---|---:|---:|---:|---:|
'''
    for row in summary:
        if row["family"] == "PathGennie_raw":
            markers = " (recovered storage order)" if row["recovered"] else " (approximate NPT boxes)" if row["approximate_box"] else ""
            report += f"| {row['source']}{markers} | {row['samples']} / {row['stored_frames']} | {row['core_first']} -> {row['core_last']} | {row['core_max']} | {row['Q_first']:.4f} -> {row['Q_last']:.4f} |\n"
    report += f'''
![Largest connected ice core](core_evolution.png)

Stored-frame progress is normalized independently per source. PathGennie accepted
paths and metadynamics are selected/biased ensembles; these curves are not an
unbiased nucleation rate or a common physical clock. Native MD/NPZ times are
preserved in the CSV where available. Shooting branches and all individual
sources are listed in [source_summary.csv](source_summary.csv).

## Local OP contrast

![Local OP means](core_op_contrast.png)

Core membership means the largest connected set of CHILL+ Ih/Ic centers, using
3.5-Angstrom O--O edges; interface labels are excluded from this strict core.
Outside-core averages include every remaining center. Tiny cores (<20 waters)
are omitted from the plotted core mean to avoid interpreting single random
ice-like molecules as a nucleus.

The following table summarizes signed ROC AUC against that core-membership mask
on readable fixed-box raw E/F/I/M path samples with cores of at least 20 waters.
Values above 0.5 indicate larger OP values inside the core; below 0.5 indicate
smaller values. Values are medians over frames, not pooled atom counts. This is
a descriptive association to a geometry-derived label, not held-out predictive
accuracy or independent phase-identification ground truth.

The label is not independent of every parameter. CHILL+ is built from q3
bond-orientational correlations over 3.5-Angstrom neighbors and labels only
four-coordinated centers, so high scores for averaged q3 (and, less directly,
other averaged bond-order parameters and crystalline connections) partly
restate the label, and coordination-based scores are constrained by it.
Parameter rankings against independent ground truth come from constructed
seeds instead (benchmarks/, seed test).

The table lists the {TOP_AUC} parameters with the largest median separation from
0.5, in order and without omissions.

| OP | Frames | Median AUC | Interquartile range |
|---|---:|---:|---:|
'''
    for row in selected:
        report += f"| {row['parameter']} | {row['frames']} | {row['median_auc']:.3f} | {row['q25_auc']:.3f}–{row['q75_auc']:.3f} |\n"
    report += '''
![Association with largest-core membership](core_association.png)

Full rankings are in [core_contrast.csv](core_contrast.csv). Near-constant model
moments, such as Q_T for rigid TIP4P/Ice, are not expected to distinguish lattice
order. Global OP averages can change little when the nucleus contains only a
small fraction of the simulation box.

## Coverage and limits

All 25 original calculations are evaluated on full-atom TIP4P/Ice inputs,
including local/global entropy variants, raw/normalized bond order, both CHILL
classifiers, Voronoi metrics, and the complete water-angle distribution.
Additional A/B/H parameters and geometry-derived TIP4P/Ice Q_T are also evaluated.
H-bonds, F4, water-H angles and Q_T are excluded from mW and oxygen-only movie
inputs because those inputs contain no hydrogens. The validation script does
not infer or synthesize hydrogen positions. Optional DScribe/PTM descriptors
were validated separately; this run focuses on structural OPs.

Periodic dimensions come from native DCD/NPZ metadata, the documented NVT box,
or the following explicitly approximate NPT source. Entropy spheres fit inside
half the smallest periodic height. Default OP search cutoff is 5 A; first-shell
coordination, CHILL+, F4 and bond-order comparisons use 3.5 A.

### Checkpoint and preprocessing caveats

- G/J/K checkpoint object headers are unreadable through HDF5. Complete, bounded
  zlib coordinate streams were recovered read-only, in the same physical-storage
  convention used by the source project. Logical HDF5 order and exact accepted
  path history are not guaranteed for those recovered inputs. Input files were
  not repaired or modified.
- O_Ic checkpoints contain no per-frame boxes. The source project's movie-info
  boxes were reconstructed as cycle-window medians; they are used and flagged
  as approximate, with no fallback to the NVT box. The results are not an exact
  NPT box/trajectory validation.
- Prepared movie files are FIRE-quenched, oxygen-only, and resampled. They are
  reported as a separate family. Some movies describe earlier checkpoints;
  they must not be interpreted as exactly paired thermal/quenched trajectories.
- Older duplicate checkpoints, seed-test/training arrays, and unrelated bulk
  controls are excluded. The mW T215/T220/T225 directories contain CV logs but no
  local coordinate trajectories, so no OP calculation is claimed for them.
- Finite sampling may miss transient maxima. No kinetic rate, critical nucleus
  size or universal OP threshold is estimated here.
'''
    for source in manifest["sources"]:
        missing = source.get("unrecovered_relative_to_movie", 0)
        if missing:
            report += f"- {source['name']}: recovered {source['frames']} complete coordinate chunks versus {source['source_movie_reported_frames']} frames in the saved movie-info table ({missing} fewer).\n"
    report += '''
## Reproduction and artifacts

```sh
python -m pip install -e ".[validation]"
MPLCONFIGDIR=/tmp/solorder-mpl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \\
  python scripts/validate_nucleation.py \\
  --data-root /home/dm/Dibyendu/Projects/ice_nucleation_pg
MPLCONFIGDIR=/tmp/solorder-mpl python scripts/audit_nucleation.py
MPLCONFIGDIR=/tmp/solorder-mpl python scripts/report_nucleation.py
```

- [Frame-level means, standard deviations, NaN counts and association scores](frames.csv)
- [Per-core/outside-core statistics](phase_conditioned.csv)
- [Independent checks per frame](checks.csv)
- [All source summaries](source_summary.csv)
- [Input provenance and selected-frame hashes](manifest.json)
- Standalone figures: [core evolution](core_evolution.pdf),
  [local OP contrast](core_op_contrast.pdf), [core association](core_association.pdf).
'''
    if audits is not None:
        maxima = audits.filter(regex="_max_error$").max()
        legendre = float(maxima.filter(like="Legendre").max())
        freud_columns = [c for c in maxima.index if c.startswith(("q3_", "q4_", "q6_", "w6_")) and "Legendre" not in c]
        section = f"\n## Independent numerical spot checks\n\n{len(audits)} endpoint snapshots were additionally audited against separate numerical\nreferences; **{int(audits.passed.sum())}/{len(audits)} passed**. Non-averaged q3/q4/q6\nagreed with the double-precision Legendre addition formula within {legendre:.3g}.\nDirect tetrahedral, radial and LSI formulas agreed within\n{max(maxima.get('Q_max_error', 0), maxima.get('TTO_max_error', 0), maxima.get('LSI_max_error', 0)):.3g}.\nIndependent MDAnalysis torsions/angles agreed with F4 and P(theta) within\n{max(maxima.get('F4_max_error', 0), maxima.get('Ptheta_max_error', 0)):.3g}.\n\nFreud 3.4 bond-order comparisons differed by up to {float(maxima[freud_columns].max()):.3g};\nthese backend comparisons use an explicit 1e-4 absolute tolerance. Direct\ngeometric formula checks use 1e-10, Legendre checks 1e-7, and angle/torsion\nchecks 1e-6. Backend differences are reported rather than treated as exact\nagreement. See [independent_audits.csv](independent_audits.csv).\n"
        report = report.replace("\n## Coverage and limits", section+"\n## Coverage and limits")
    report = report.replace("Full rankings are in", "Averaged q6 and averaged tetrahedral order show strong positive core association,\nwhile averaged entropy is lower in the core. H-bond counts are nearly four in\nboth populations and show little ranking separation here (median AUC about\n0.51). CHILL+ itself requires four-coordination, so this is a conditional result\non these selected configurations.\n\nFull rankings are in")
    unchanged = all(Path(source["path"]).stat().st_size == source["size_bytes"] and
                    Path(source["path"]).stat().st_mtime_ns == source["mtime_ns"]
                    for source in manifest["sources"])
    manifest["source_sizes_timestamps_unchanged"] = unchanged
    (directory/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    if not unchanged:
        raise RuntimeError("A source trajectory changed during validation")
    (directory/"README.md").write_text(report)
    print(f"Report written for {len(frames)} frames; {len(failures)} failed-check frames")


if __name__ == "__main__":
    main()
