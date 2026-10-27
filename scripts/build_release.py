"""Package the benchmark as an archival dataset: coordinates, code snapshot, metadata, splits, results.

Layout of release/solorder-ice-benchmark-v<VERSION>/:
  coordinates/   every structure and trajectory, in the layout of benchmarks/cache/,
                 so that the pipeline runs on it unchanged:
                   genice/      ideal GenIce2 structures (seed<N>/<phase>.gro)
                   seedtest/    run A (seeds, bulk) and run B (train_*.npz), the
                                TIP4P/Ice topology and the MD helpers they need
                   replicas/    runs C1-C5
                   conditions/  runs D (TIP4P/2005, TIP4P/Ice NPT control) and E
                                (high pressure); conditions/mw_inputs/ is mW
                   interface/   run F (seeds of 50-1600 molecules, slabs) and the
                                run G confirmation set (G_*, preparation p4)
                   applications/ per-water order-parameter tables of the two
                                unbinding-path applications (supplementary
                                material); the trajectories are not included
  code/          exact snapshot of the code that produced the results (package,
                 benchmark scripts, examples, tests, environment), with
                 code/VERSION.json (git commit, uncommitted files) and
                 code/uncommitted.diff
  TOPOLOGY.md    provenance of the force-field files
  LABELS.md      meaning of every reference label
  splits.json    the fixed fitting and evaluation splits used in the paper
  results/       every result table and figure
  MANIFEST.sha256
Per-molecule order-parameter tables of the ice benchmark are not included:
compute.py regenerates them from the coordinates (about 5 minutes on 10 cores).

  python scripts/build_release.py
"""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT/"benchmarks"
CACHE = BENCH/"cache"
SEEDTEST = Path(os.environ.get("SOLORDER_SEEDTEST", CACHE/"seedtest"))
VERSION = "1.0"
OUT = ROOT/"release"/f"solorder-ice-benchmark-v{VERSION}"

COORDINATES = {  # destination: (source, patterns)
    "genice": (CACHE/"genice", ("seed*/*.gro", "seed*/manifest.json")),
    "seedtest": (SEEDTEST, ("*.gro", "*.top", "*.npz", "*.py", "*.itp", "charmm36-mar2019.ff/*")),
    "replicas": (CACHE/"replicas", ("*.npz", "*.gro", "*.top", "*.json", "*.npy")),
    "conditions": (CACHE/"conditions", ("D_*.npz", "D_*.json", "E_*.npz", "E_*.json", "E_*.gro", "E_*.top",
                                        "mw_inputs/*.npz", "mw_inputs/*.json", "mw_inputs/*.gro",
                                        "mw_inputs/*.npy", "mw_inputs/*.py", "mw_inputs/*.sh")),
    "interface": (CACHE/"interface", ("F_seed_*", "F_slab_*", "F_liquid_*.npz", "G_seed_*", "G_prism_*",
                                      "G_long_*")),
    "applications": (CACHE/"applications", ("*_features.npz",)),
}
CODE = ("src/SolOrder/*.py", "benchmarks/*.py", "benchmarks/README.md", "scripts/*.py", "examples/*.py",
        "examples/README.md", "tests/*.py", "pyproject.toml", "environment.yml", "requirements.txt",
        "README.md", "MANIFEST.in", "LICENSE", ".gitignore")
DESCRIPTIONS = {
    "genice": "18 GenIce2 structures from 5 generator seeds (jitter is applied in compute.py)",
    "seedtest": "run A: TIP4P/Ice 220 K bulk crystals and seeds in liquid, 10 frames 1 ps apart; run B "
                "(train_*.npz): bulk Ih, Ic, Isd, liquid at 220 K, Ih, Ic at 250 K, 60 frames 5 ps apart",
    "replicas": "runs C1-C5: five independent preparations, TIP4P/Ice 220 K (Ih, Ic also 250 K)",
    "conditions": "run D: TIP4P/2005 (220, 200 K) and TIP4P/Ice NPT control on the C1-C3 structures, mW "
                  "(mw_inputs/); run E: ices II, III, V, VI, VII, VIII and pressure-matched liquids",
    "interface": "run F: Ih, Ic, Isd seeds of 50-1600 molecules and planar basal slabs, three preparations; "
                 "run G (G_*): seeds of 300 and 600 molecules, prism slabs and basal slabs followed to 500 ps, "
                 "a fourth preparation built after the calibration was frozen (confirmation_protocol.md)",
    "applications": "supplementary applications: water-only order parameters, geometric region labels and "
                    "split columns per water and frame for OAMe-G2 (oame_features.npz) and T4L L99A "
                    "(t4l_features.npz) unbinding paths; input of host_guest.py analyze and "
                    "confined_water.py analyze",
}

SPLITS = {
    "description": "Fixed splits used in the paper. A threshold or model is fitted on 'fit' and scored on 'test'.",
    "heldout_runs": {"fit": "seedtest frames_*.npz (run A)", "test": "seedtest train_*.npz (run B)"},
    "independent_preparations": {"fit": "four of replicas p1-p5 (C1-C5)", "test": "the fifth"},
    "frozen_thresholds": {"fit": "replicas p1-p3, thermal Ih versus liquid, 220 K",
                          "test": ["replicas p4, p5 liquid", "conditions D tip4p2005 liquids", "conditions mW liquids",
                                   "conditions E liquids p1, p2"]},
    "model_transfer": {"fit": "conditions D_tip4pice_*", "test": ["D_tip4p2005_*_220", "D_tip4p2005_*_200",
                       "mw_inputs/D_mW_*"], "report": "worst preparation"},
    "pressure": {"fit": "E p1 (p2)", "test": "E p2 (p1)", "report": "worse direction"},
    "ice_V": {"select": "E_V_p1 and E_liquid_p1_0.5GPa, frames 0-29 fit, 30-59 validate", "refit": "p1",
              "test": ["p2", "p3 (used only for the final test)"]},
    "nucleus_size": {"calibrate": "interface F_slab_* only (count dead layer, frozen before run G)",
                     "test": "every interface F_seed_*",
                     "confirm": "interface G_seed_* (built and simulated after the calibration was frozen)",
                     "facet_and_time": "interface G_prism_*, G_long_* (measured, never used to recalibrate)"},
    "interface_depth_models": {"fit": "F_seed_* of 100 and 200 molecules (and F_slab_*)",
                               "test": "F_seed_* of 50, 400, 800 and 1600 molecules"},
    "sensitivity": {"fit": "replicas p1 at the default setting", "test": "p2 and p3 at each varied setting"},
    "oame_g2": {"fold": "path number (column 'path', trajectory_<k>) modulo 2", "report": "worse direction"},
    "t4l_l99a": {"fit": "replicas r1, r2 (column 'replica')", "test": "r3, r4", "report": "worse direction"},
}

LABELS = """# Reference labels

All labels describe how a structure was constructed, not what a labeler says.

- `phase` (file names): the phase a run was started in. Crystals were checked
  for retention (oxygen displacement, dipole orientation, density); see the
  per-run .json files.
- `stack` (seedtest/truth.npz, replicas *.stack.npy, interface *.truth.npz):
  1 if the molecule's layer is hexagonally stacked (its neighboring layers sit
  in the same lateral position), 2 if cubic; 0 for liquid molecules.
- `seed` (seedtest/truth.npz, interface *.truth.npz): the molecule belongs to
  the constructed crystal; constructed molecules come first in each .gro.
- `interior` (seedtest, interface seeds): all four constructed lattice
  neighbors are also in the crystal (interface slabs: true for every crystal
  molecule).
- `depth` (interface): R - r for seeds, the distance below the outermost
  constructed oxygen; for slabs, the distance below the outermost oxygen plane
  of the nearer face (prism slabs G_prism_*: the {10-10} face). NaN for liquid.
- `time_ps` (interface *.npz of run G): time after release of each saved frame.
- Dead layer (benchmarks/nucleus_size.py): delta = Gamma/n, Gamma the number of
  crystal molecules a labeler misses per unit area of one slab face
  (sum of 1 - p_i/f over crystal molecules, per face and frame), f its recall
  deeper than 8 A.
- `R` (seedtest, interface): radius of the outermost constructed seed molecule.
  The equimolar radius used for nucleus sizes is (3N/4 pi n)^(1/3) with
  n = 3/(A^2 sqrt(3) B), A = B sqrt(8/3), B = 2.76 A (benchmarks/nucleus_size.py).
- Liquid reference in run A: molecules more than 8 A beyond R.
- Intact (computed by seed_audit.py and interface_depth.py, not stored):
  within 1 A of the constructed site after removing the crystal's drift, and
  all constructed bonds shorter than 3.5 A.
"""

TOPOLOGY = """# Force-field files

Every TIP4P/Ice and TIP4P/2005 system is built from the GROMACS topology in
`coordinates/seedtest/` (`*.top`, which includes
`charmm36-mar2019.ff/forcefield.itp` and `tip4p-sol.itp`) by
`coordinates/seedtest/system_ice.py` (OpenMM, PME, 0.9 nm cutoff, dispersion
correction, rigid water).

- `charmm36-mar2019.ff/` is the CHARMM36 (March 2019) port for GROMACS,
  **modified**: its `ffnonbonded.itp` entry for the water oxygen type `OWT4`
  carries the TIP4P/Ice Lennard-Jones parameters, sigma = 0.31668 nm and
  epsilon = 0.88211 kJ/mol (Abascal et al., J. Chem. Phys. 122, 234511, 2005).
  Do not replace it with an unmodified CHARMM36 folder. Compared with an
  unmodified copy of the same port (`forcefield_modification.diff`), the
  stock block "atomtypes for additional water models" (TIP4P OWT4 with sigma
  0.315365 nm and epsilon 0.648520 kJ/mol, HWT4, MWT4, SPC OW and HW) is
  replaced by TIP4P/Ice OWT4, HW and IW and three methane types (C4-C6, from
  an earlier methane-hydrate project, unused here), and `tip4p.itp` (not
  included by any topology here) also differs. A water-only system uses only
  OWT4, HW and IW, so no other entry can affect the results.
- `tip4p-sol.itp` defines the rigid four-site molecule with TIP4P/Ice charges
  (H +0.5897 e, M -1.1794 e), O-H 0.09572 nm, H-H 0.15139 nm and the virtual
  site M at weights a = b = 0.13458 (O-M 0.01577 nm).
- TIP4P/2005 runs (conditions/D_tip4p2005_*) rewrite charges, Lennard-Jones
  parameters and the M-site weights in code (benchmarks/md_conditions.py,
  set_model: H +0.5564 e, sigma 0.31589 nm, epsilon 0.774907 kJ/mol, O-M
  0.01546 nm).
- mW runs use no topology file (benchmarks/md_mw.py, Stillinger-Weber with
  the parameters of Molinero and Moore, J. Phys. Chem. B 113, 4008, 2009).

SHA-256 of the files above are in MANIFEST.sha256.
"""


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout


def code_snapshot():
    files = sorted({p for pattern in CODE for p in ROOT.glob(pattern) if p.is_file()})
    for f in files:
        copy(f, OUT/"code"/f.relative_to(ROOT))
    paths = [str(f.relative_to(ROOT)) for f in files]
    status = git("status", "--porcelain", "--", *paths).splitlines()
    (OUT/"code/uncommitted.diff").write_text(git("diff", "HEAD", "--", *paths))
    info = dict(commit=git("rev-parse", "HEAD").strip(), branch=git("rev-parse", "--abbrev-ref", "HEAD").strip(),
                uncommitted=status, files=len(files),
                note="code/ is the exact code used; uncommitted.diff is its difference from the commit "
                     "(untracked files, marked ??, are only in code/)")
    (OUT/"code/VERSION.json").write_text(json.dumps(info, indent=1))
    return info


def freeze(python):
    """`pip freeze` of an interpreter, or a note if it is not available."""
    try:
        result = subprocess.run([str(python), "-m", "pip", "freeze"], capture_output=True, text=True, timeout=120)
        return result.stdout if result.returncode == 0 else f"# pip freeze failed: {result.stderr.strip()}\n"
    except (OSError, subprocess.TimeoutExpired) as error:
        return f"# not recorded: {error}\n"


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    datasets = {}
    for name, (source, patterns) in COORDINATES.items():
        files = sorted({p for pattern in patterns for p in source.glob(pattern) if p.is_file()})
        for f in files:
            copy(f, OUT/"coordinates"/name/f.relative_to(source))
        datasets[name] = dict(description=DESCRIPTIONS[name], files=len(files))
    info = code_snapshot()
    shutil.copytree(BENCH/"results", OUT/"results", ignore=shutil.ignore_patterns("pilot_v1"))
    copy(BENCH/"cache/manifest.json", OUT/"metadata/input_hashes.json")
    (OUT/"splits.json").write_text(json.dumps(SPLITS, indent=1))
    (OUT/"LABELS.md").write_text(LABELS)
    (OUT/"TOPOLOGY.md").write_text(TOPOLOGY)
    copy(ROOT/"docs/forcefield_modification.diff", OUT/"forcefield_modification.diff")
    copy(BENCH/"confirmation_protocol.md", OUT/"confirmation_protocol.md")
    import SolOrder
    (OUT/"metadata/versions.json").write_text(json.dumps(dict(solorder=SolOrder.__version__, release=VERSION,
                                                              code=info, datasets=datasets), indent=1))
    (OUT/"metadata/environment_analysis.txt").write_text(f"# environment {Path(sys.executable).parents[1].name}\n" + freeze(sys.executable))
    md_python = os.environ.get("SOLORDER_MD_PYTHON")
    (OUT/"metadata/environment_simulation.txt").write_text(
        f"# environment {Path(md_python).parents[1].name}\n" + freeze(md_python) if md_python else
        "# set SOLORDER_MD_PYTHON to the OpenMM interpreter when building the release\n")
    (OUT/"README.md").write_text(readme(datasets, info))
    files = sorted(p for p in OUT.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    (OUT/"MANIFEST.sha256").write_text("".join(f"{sha256(p)}  {p.relative_to(OUT)}\n" for p in files))
    size = sum(p.stat().st_size for p in files)
    print(f"{OUT}: {len(files)} files, {size/1e9:.2f} GB; code {info['commit'][:10]} "
          f"with {len(info['uncommitted'])} uncommitted files")


def readme(datasets, info):
    lines = [f"# SolOrder ice order-parameter benchmark, v{VERSION}", "",
             "Coordinates, code, metadata, fixed evaluation splits and results for the benchmark in",
             "'Local order parameters for ice polymorphs: a unified comparison and practical guide'.",
             "Units: A, ps, K, GPa. Verify with `sha256sum -c MANIFEST.sha256`.", "",
             "## Contents", ""]
    lines += [f"- `coordinates/{k}/` ({v['files']} files): {v['description']}" for k, v in datasets.items()]
    lines += [f"- `code/`: the exact code used (git commit {info['commit'][:10]} plus `code/uncommitted.diff`)",
              "- `results/`: every result table and figure; `splits.json`, `LABELS.md`, `TOPOLOGY.md`",
              "- `confirmation_protocol.md`: the frozen nucleus-size calibration and pass criterion, written",
              "  before the run G confirmation set was simulated", "",
              "## Reproduce", "",
              "```sh",
              "conda env create -f code/environment.yml && conda activate solorder-dev",
              "pip install -e 'code[descriptors,benchmark]'",
              "mkdir -p code/benchmarks/cache && cp -r coordinates/* code/benchmarks/cache/",
              "mkdir -p code/benchmarks/results && cp -r results/. code/benchmarks/results/  # archived tables",
              "python code/examples/benchmark_quickstart.py        # one published score in seconds",
              "python code/benchmarks/figures.py                   # every figure from the archived tables",
              "cd code/benchmarks",
              "python compute.py && python metrics.py              # per-molecule tables, separability",
              "python seed_audit.py && python interface_depth.py compute && python interface_depth.py metrics",
              "python nucleus_size.py && python interface_depth.py compute --confirmation && python confirmation.py",
              "python false_nuclei.py && python ice_v.py",
              "python sensitivity.py compute && python sensitivity.py metrics",
              "python liquid_diagnostics.py && python redundancy_by_phase.py",
              "python host_guest.py analyze && python confined_water.py analyze  # from coordinates/applications",
              "OMP_NUM_THREADS=1 python cost.py --jobs 1 && python cost.py --jobs 8  # timings, machine-dependent",
              "python figures.py                                   # figures from the regenerated tables",
              "```", "",
              "Each step overwrites the archived table it regenerates in code/benchmarks/results/; tables",
              "that need OpenMM (quench_check.py quench, the md_*.py runs) stay as archived. The",
              "confined-water tables come from coordinates/applications/, which holds per-water order",
              "parameters and labels but not the unbinding trajectories (separate path-sampling",
              "campaigns); host_guest.py scan/compute and confined_water.py compute need those.", "",
              "Simulations (OpenMM with CUDA) are rerun with md_replicas.py, md_conditions.py, md_mw.py and",
              "md_interface.py; the seed test with coordinates/seedtest/build.py and md.py. They need a",
              "separate environment, e.g. `conda create -n solorder-md -c conda-forge python=3.11 openmm",
              "networkx scipy` followed by `pip install -e 'code[simulation]'`. Set SOLORDER_SEEDTEST to use",
              "a seed-test folder outside benchmarks/cache/.", "",
              "The exact package versions used are in metadata/environment_analysis.txt (analysis) and",
              "metadata/environment_simulation.txt (simulations).", "",
              "Data license: CC BY 4.0. Code license: MIT."]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
